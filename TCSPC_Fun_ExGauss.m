function [err, c, zz, z, IRF] = TCSPC_Fun_ExGauss(p, t, y, para)
% TCSPC_Fun_ExGauss
%
% Reconvolution fit model using an ex-Gaussian IRF.
% Inner nonnegative amplitude solve uses PIRLSnonneg if available.
%
% Model:
%   y(t) ~ c0*1 + sum_i ci * [ IRF_exG (*) exp(-t/tau_i)/tau_i ]
%
% Parameters
% ----------
% If para is empty:
%   p = [t0; sigma; tauR; tau1; tau2; ...]
%
% If para is provided:
%   para = [t0; sigma; tauR] fixed IRF parameters
%   p    = [tau1; tau2; ...]
%
% Inputs
% ------
% t : time axis [ns]
% y : decay vector or matrix [Nt x Ny]
%
% Outputs
% -------
% err : objective value
% c   : nonnegative amplitudes (background in row 1)
% zz  : basis matrix
% z   : fitted model
% IRF : normalized ex-Gaussian IRF

p = p(:);
t = t(:);

if nargin > 3 && ~isempty(para)
    p = [para(:); p(:)];
end

nIRF = 3;
if numel(p) <= nIRF
    error('Parameter vector must contain [t0; sigma; tauR; tau1; ...].');
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

IRF = IRF_ExGauss(p(1:nIRF), t, []);

% Decay kernel starts at first TCSPC bin
tdec = t - t(1);
tdec(tdec < 0) = 0;

nT   = numel(t);
nExp = numel(tau);

zz = zeros(nT, nExp + 1);
zz(:,1) = ones(nT,1);   % background column

for i = 1:nExp
    decay = exp(-tdec ./ tau(i)) ./ tau(i);
    tmp = Convol(IRF, decay);
    zz(:,1+i) = tmp(1:nT);
end

% Normalize columns for conditioning
colsum = sum(zz,1);
colsum(colsum <= 0) = 1;
zz = zz ./ (ones(nT,1) * colsum);

% Solve amplitudes with PIRLS
[c, z] = solve_nonneg_with_pirls(zz, y);

% Pearson-like objective
den = max(abs(z), 10);
err = sum(sum((y - z).^2 ./ den));

end


function [c, z] = solve_nonneg_with_pirls(M, Y)
% Solve min ||Y - M*C||^2 with C >= 0
% Preference order:
%   1) PIRLSnonneg_batch_gpu_matlab(M,Y)
%   2) PIRLSnonneg(M,Y)
%   3) PIRLSnonneg(M,Y(:,j)) columnwise
%   4) fallback to lsqnonneg columnwise

nBasis = size(M,2);
nY = size(Y,2);

% 1) Batch GPU PIRLS if present
if nY > 1 && exist('PIRLSnonneg_batch_gpu_matlab', 'file') == 2
    try
        c = PIRLSnonneg_batch_gpu_matlab(M, Y);
        c = max(c, 0);
        z = M * c;
        return
    catch
    end
end

% 2) PIRLS that may accept matrix RHS
if exist('PIRLSnonneg', 'file') == 2
    try
        c = PIRLSnonneg(M, Y);
        c = max(c, 0);
        z = M * c;
        return
    catch
    end

    % 3) Columnwise PIRLS
    try
        c = zeros(nBasis, nY);
        for j = 1:nY
            c(:,j) = PIRLSnonneg(M, Y(:,j));
        end
        c = max(c, 0);
        z = M * c;
        return
    catch
    end
end

% 4) Fallback
c = zeros(nBasis, nY);
z = zeros(size(M,1), nY);
for j = 1:nY
    c(:,j) = lsqnonneg(M, Y(:,j));
    z(:,j) = M * c(:,j);
end
end