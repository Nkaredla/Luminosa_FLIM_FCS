function out = Calc_mIRF_Global_GammaShifted_fast(head, tcspc, tau0, opts)
% Calc_mIRF_Global_GammaShifted_fast
% Faster global shifted-gamma IRF fit.
%
% Main speedups relative to the original function:
%   1) accepts a pre-collapsed 1D decay vector
%   2) uses fewer multistart evaluations by default
%   3) can collapse x/y/channel dimensions on GPU before fitting
%
% GPU only helps the decay reduction step. The simplex optimization itself
% is low-dimensional and remains CPU-side.
%
% Inputs
% ------
% head : PTU header struct
% tcspc: decay vector or reducible TCSPC array
% tau0 : initial lifetimes [ns]
% opts.useGPU : logical, default false
% opts.nCasc  : default 4
% opts.nSub   : default 6
%
% Output fields match Calc_mIRF_Global_GammaShifted.

if nargin < 3 || isempty(tau0)
    tau0 = [0.4 2.0];
end
if nargin < 4 || isempty(opts)
    opts = struct();
end
useGPU = get_opt(opts, 'useGPU', false);
nCasc  = get_opt(opts, 'nCasc', 4);
nSub   = get_opt(opts, 'nSub', 6);

tau0 = tau0(:);
Resolution = get_resolution_ns(head);
y = reduce_to_single_decay_fast(tcspc, useGPU);
y = double(y(:));
t = Resolution .* ((1:numel(y))' - 0.5);

if numel(y) ~= numel(t)
    error('Time axis length does not match reduced TCSPC decay length.');
end

[~, indMax] = max(y);
t0_0 = t(indMax);
alpha0 = 3.0;
beta0  = max(2*Resolution, 0.03);

p0 = [t0_0; alpha0; beta0; tau0];
pl = [t0_0 - 1.0; 0.5; max(Resolution/10, 0.001); max(0.05 * tau0, 0.01)];
pu = [t0_0 + 1.0; 20.0; 2.0; max(20 * tau0, 10)];

bestErr = Inf;
bestP = p0;

for casc = 1:nCasc
    for sub = 1:nSub
        rf = p0;
        rf(1) = p0(1) + (0.75 / casc) * (2*rand - 1);
        rf(2:end) = p0(2:end) .* 2.^(0.7 * (rand(size(p0(2:end))) - 0.5) / casc);
        rf = max(rf, pl);
        rf = min(rf, pu);

        if exist('Simplex', 'file') == 2
            pfit = Simplex('TCSPC_Fun_GammaShifted', rf, pl, pu, [], [], t, y, []);
        else
            pfit = fminsearch(@(p) boundedErrGamma(p, pl, pu, t, y), rf, optimset('Display','off','MaxIter',150));
            pfit = min(max(pfit(:), pl), pu);
        end
        err = TCSPC_Fun_GammaShifted(pfit, t, y, []);
        if err < bestErr
            bestErr = err;
            bestP = pfit(:);
        end
    end
end

[~, cfit, zz, zfit, IRF] = TCSPC_Fun_GammaShifted(bestP, t, y, []);
ndof = max(numel(y) - numel(bestP), 1);
chi = sum((y - zfit).^2 ./ max(abs(zfit), 10)) / ndof;

out = struct();
out.IRF    = IRF(:);
out.pIRF   = bestP(1:3);
out.tauFit = bestP(4:end);
out.A      = cfit;
out.zz     = zz;
out.fit    = zfit;
out.t      = t;
out.chi    = chi;
out.y      = y;
out.pAll   = bestP;
end

function err = boundedErrGamma(p, pl, pu, t, y)
p = min(max(p(:), pl), pu);
err = TCSPC_Fun_GammaShifted(p, t, y, []);
end

function y = reduce_to_single_decay_fast(tcspc, useGPU)
if isvector(tcspc)
    y = double(tcspc(:));
    return
end

sz = size(tcspc);
nd = ndims(tcspc);

if nd == 4
    % Expected layout: x, y, time, channel.
    if useGPU && gpuIsAvailableLocal()
        g = gpuArray(single(tcspc));
        y = gather(double(squeeze(sum(sum(sum(g, 1), 2), 4))));
    else
        y = squeeze(sum(sum(sum(double(tcspc), 1), 2), 4));
    end
elseif nd == 3
    % Expected layout: x, y, time.
    if useGPU && gpuIsAvailableLocal()
        g = gpuArray(single(tcspc));
        y = gather(double(squeeze(sum(sum(g, 1), 2))));
    else
        y = squeeze(sum(sum(double(tcspc), 1), 2));
    end
else
    x = squeeze(double(tcspc));
    if isvector(x)
        y = x(:);
    else
        % Fallback: assume time is the last dimension.
        tdim = numel(sz);
        dimsToSum = setdiff(1:tdim, tdim);
        y = double(tcspc);
        for k = sort(dimsToSum, 'descend')
            y = sum(y, k);
        end
        y = squeeze(y);
        if ~isvector(y)
            error('Could not reduce TCSPC input to a single decay vector safely.');
        end
    end
end

y = double(y(:));
end

function Resolution = get_resolution_ns(head)
if isfield(head, 'MeasDesc_Resolution') && ~isempty(head.MeasDesc_Resolution)
    Resolution = double(head.MeasDesc_Resolution) * 1e9;
elseif isfield(head, 'Resolution') && ~isempty(head.Resolution)
    Resolution = max(double(head.Resolution));
else
    error('Could not determine TCSPC resolution from head.');
end
end

function v = get_opt(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultVal;
end
end

function ok = gpuIsAvailableLocal()
ok = false;
try
    ok = gpuDeviceCount > 0;
catch
    ok = false;
end
end
