function [err, c, zz, z, IRF] = TCSPC_Fun_ExGauss(p, t, y, para)
% TCSPC_Fun_ExGauss
%
% Reconvolution fit model using an ex-Gaussian IRF.
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
%   p    = [tau1; tau2; ...]  or [other free params]
%
% Inputs
% ------
% t : time axis [ns]
% y : decay vector or matrix [Nt x Ny]
%
% Outputs
% -------
% err : objective value
% c   : NNLS amplitudes (background in row 1)
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

% Decay kernel starts at the first TCSPC bin
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

% Optional mild normalization of non-background columns
% Keeps scales manageable without changing span
colsum = sum(zz,1);
colsum(colsum <= 0) = 1;
zz = zz ./ (ones(nT,1) * colsum);

nY = size(y,2);
c  = zeros(size(zz,2), nY);
z  = zeros(nT, nY);

for j = 1:nY
    c(:,j) = lsqnonneg(zz, y(:,j));
    z(:,j) = zz * c(:,j);
end

% Pearson-like objective
den = max(abs(z), 10);
err = sum(sum((y - z).^2 ./ den));
end