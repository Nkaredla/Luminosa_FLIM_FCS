function out = Calc_mIRF_Global_GammaShifted(head, tcspc, tau0)
% Calc_mIRF_Global_GammaShifted
%
% Fit one shared shifted-gamma IRF on the summed TCSPC decay.
%
% Inputs
% ------
% head : PTU header struct
% tcspc: decay vector, image cube, or any array reducible to one decay
% tau0 : initial lifetime guesses [ns], e.g. [0.35 1.8 4.2]
%
% Output fields
% -------------
% out.IRF      : fitted normalized IRF vector
% out.pIRF     : [t0 alpha beta]
% out.tauFit   : fitted lifetimes [ns]
% out.A        : fitted amplitudes (background first)
% out.zz       : basis matrix
% out.fit      : fitted decay
% out.t        : time axis [ns]
% out.chi      : Pearson-like chi^2 / dof
% out.y        : global decay used for fitting
% out.pAll     : full fitted parameter vector
%
% Requires
% --------
% Simplex.m
% TCSPC_Fun_GammaShifted.m
% IRF_GammaShifted.m
% Convol.m

if nargin < 3 || isempty(tau0)
    tau0 = [0.4 2.0];
end
tau0 = tau0(:);

Resolution = get_resolution_ns(head);
t = Resolution .* ((1:get_num_time_bins(tcspc))' - 0.5);

y = reduce_to_single_decay(tcspc);
y = double(y(:));

if numel(y) ~= numel(t)
    error('Time axis length does not match reduced TCSPC decay length.');
end

% Initial onset near histogram peak
[~, indMax] = max(y);
t0_0 = t(indMax);

% Conservative starting gamma parameters
alpha0 = 3.0;                      % moderate skew
beta0  = max(2*Resolution, 0.03);  % ns

p0 = [t0_0; alpha0; beta0; tau0];

% Bounds
pl = [t0_0 - 1.0; ...
      0.5; ...
      max(Resolution/10, 0.001); ...
      max(0.05 * tau0, 0.01)];

pu = [t0_0 + 1.0; ...
      20.0; ...
      2.0; ...
      max(20 * tau0, 10)];

bestErr = Inf;
bestP   = p0;

nCasc = 8;
nSub  = 12;

for casc = 1:nCasc
    for sub = 1:nSub
        rf = p0;

        % Additive randomization for onset
        rf(1) = p0(1) + (1.0 / casc) * (2*rand - 1);

        % Multiplicative randomization for positive params
        rf(2:end) = p0(2:end) .* 2.^(0.9 * (rand(size(p0(2:end))) - 0.5) / casc);

        rf = max(rf, pl);
        rf = min(rf, pu);

        pfit = Simplex('TCSPC_Fun_GammaShifted', rf, pl, pu, [], [], t, y, []);
        err  = TCSPC_Fun_GammaShifted(pfit, t, y, []);

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


function y = reduce_to_single_decay(tcspc)
% Reduce arbitrary TCSPC array to a single decay vector over time.

y = squeeze(double(tcspc));

if isvector(y)
    y = y(:);
    return
end

while ~isvector(y)
    y = sum(y, 1);
    y = squeeze(y);
end

y = y(:);
end


function nT = get_num_time_bins(tcspc)
x = squeeze(tcspc);
if isvector(x)
    nT = numel(x);
    return
end

sz = size(x);
nT = sz(end);
end


function Resolution = get_resolution_ns(head)
% Return TCSPC bin width in ns

if isfield(head, 'MeasDesc_Resolution') && ~isempty(head.MeasDesc_Resolution)
    Resolution = double(head.MeasDesc_Resolution) * 1e9;
elseif isfield(head, 'Resolution') && ~isempty(head.Resolution)
    Resolution = max(double(head.Resolution));
else
    error('Could not determine TCSPC resolution from head.');
end
end