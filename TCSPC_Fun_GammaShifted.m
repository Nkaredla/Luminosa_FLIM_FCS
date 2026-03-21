function [err, c, zz, z, IRF] = TCSPC_Fun_GammaShifted(p, t, y, para)
% TCSPC_Fun_GammaShifted
%
% Reconvolution fit model using a shifted-gamma IRF.
% Uses PIRLSnonneg / PIRLSnonneg_batch_gpu_matlab for the inner
% nonnegative coefficient fit instead of lsqnonneg.
%
% Model:
%   y(t) ~ c0*1 + sum_i ci * [ IRF_gamma (*) exp(-t/tau_i)/tau_i ]
%
% Parameters
% ----------
% If para is empty:
%   p = [t0; alpha; beta; tau1; tau2; ...]
%
% If para is provided:
%   para = [t0; alpha; beta] fixed IRF parameters
%   p    = [tau1; tau2; ...]
%
% Outputs
% -------
% err : objective value
% c   : nonnegative coefficients (background in row 1)
% zz  : basis matrix
% z   : fitted model
% IRF : normalized shifted-gamma IRF

p = p(:);
t = t(:);

if nargin > 3 && ~isempty(para)
    p = [para(:); p(:)];
end

nIRF = 3;
if numel(p) <= nIRF
    error('Parameter vector must contain [t0; alpha; beta; tau1; ...].');
end

if size(y,1) ~= numel(t)
    y = y(:);
end
if size(y,1) ~= numel(t)
    error('Length of t must match first dimension of y.');
end

tau = p(nIRF+1:end);
tau = tau(:);
tau(tau <= 0) = 1e-6;

IRF = IRF_GammaShifted(p(1:nIRF), t, []);

% Exponentials start at first TCSPC bin
t0 = t - t(1);
t0(t0 < 0) = 0;

nT   = numel(t);
nExp = numel(tau);

zz = zeros(nT, nExp + 1);
zz(:,1) = ones(nT,1);   % background column

for i = 1:nExp
    decay = exp(-t0 ./ tau(i)) ./ tau(i);
    tmp = Convol(IRF, decay);
    zz(:,1+i) = tmp(1:nT);
end

% Normalize columns for numerical conditioning
colsum = sum(zz,1);
colsum(colsum <= 0) = 1;
zz = zz ./ (ones(nT,1) * colsum);

% Force y to be 2D [Nt x Ny]
if isvector(y)
    y = y(:);
end

% Solve nonnegative coefficients with PIRLS
[c, z] = solve_nonneg_with_pirls(zz, y);

% Pearson-like objective
den = max(abs(z), 10);
err = sum(sum((y - z).^2 ./ den));
end


% ========================================================================
function [c, z] = solve_nonneg_with_pirls(M, Y)
% Solve min ||Y - M*C||^2 subject to C >= 0
% Preference order:
%   1) PIRLSnonneg_batch_gpu_matlab(M,Y)    for multiple RHS / GPU batch
%   2) PIRLSnonneg(M,Y)                     if your PIRLS accepts matrices
%   3) PIRLSnonneg(M,y_j) per column
%   4) fallback projected-gradient NNLS

[nT, nBasis] = size(M); %#ok<ASGLU>
nY = size(Y,2);

% Case 1: batch GPU PIRLS if available and useful
if nY > 1 && exist('PIRLSnonneg_batch_gpu_matlab', 'file') == 2
    try
        c = PIRLSnonneg_batch_gpu_matlab(M, Y);
        c = max(c, 0);
        z = M * c;
        return
    catch
        % fall through
    end
end

% Case 2: generic PIRLS that may accept matrix RHS
if exist('PIRLSnonneg', 'file') == 2
    try
        c = PIRLSnonneg(M, Y);
        c = max(c, 0);
        z = M * c;
        return
    catch
        % fall through to columnwise
    end

    try
        c = zeros(nBasis, nY);
        for j = 1:nY
            c(:,j) = PIRLSnonneg(M, Y(:,j));
        end
        c = max(c, 0);
        z = M * c;
        return
    catch
        % fall through
    end
end

% Case 3: local fallback
c = local_projected_gradient_nnls(M, Y);
c = max(c, 0);
z = M * c;
end


% ========================================================================
function C = local_projected_gradient_nnls(M, Y)
% Fallback NNLS if no PIRLS implementation is on the MATLAB path.
% This is NOT PIRLS; it is only a safety net.
%
% Solves: min ||M*C - Y||_F^2, C >= 0

[nT, nBasis] = size(M); %#ok<ASGLU>
nY = size(Y,2);

MtM = M.' * M;
MtY = M.' * Y;

% Step size from Lipschitz constant
L = norm(MtM, 2);
if ~isfinite(L) || L <= 0
    L = 1;
end
step = 1 / L;

C = max(MtY, 0);   % crude initialization
if isempty(C)
    C = zeros(nBasis, nY);
end

maxIter = 500;
tol = 1e-8;

for it = 1:maxIter
    Cprev = C;
    grad = MtM * C - MtY;
    C = C - step * grad;
    C(C < 0) = 0;

    denom = max(norm(Cprev(:)), 1);
    if norm(C(:) - Cprev(:)) / denom < tol
        break
    end
end
end