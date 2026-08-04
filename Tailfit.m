function [tau, amp, offset, z, t, err, coeff] = Tailfit(y, dt, tau0, limits, mode, bild, maxSteps)
%TAILFIT Fast direct multi-exponential tail fit without IRF reconvolution.
%
% [tau, amp, offset, z, t, err] = Tailfit(y, dt, tau0, limits, mode, bild)
%
% y      decay tail data
% dt     TCSPC channel width in ns
% tau0   initial lifetimes in ns, default [0.4 2.0]
% limits optional [lower upper] lifetime limits in ns. Give either
%        2 scalars or an nTau-by-2 array.
% mode   'pirls', 'ls', or 'mle' amplitude solve, default 'pirls'
% bild   if true, plot tail data, fit, and residual

y = real(double(y(:)));
y(~isfinite(y)) = 0;
y = max(y, 0);
n = numel(y);
if n < 3
    error('Tailfit requires at least three tail samples.');
end

dt = double(dt);
if ~isfinite(dt) || dt <= 0
    error('Tailfit requires a positive dt in ns.');
end

if nargin < 3 || isempty(tau0)
    tau0 = [0.4 2.0];
end
tau0 = real(double(tau0(:)));
tau0 = tau0(isfinite(tau0) & tau0 > 0);
if isempty(tau0)
    tau0 = [0.4; 2.0];
end

if nargin < 4
    limits = [];
end
[lb, ub] = tailLifetimeBounds(tau0, limits, dt, n);

if nargin < 5 || isempty(mode)
    mode = 'pirls';
end
mode = lower(strrep(char(mode), '-', '_'));

if nargin < 6 || isempty(bild)
    bild = 0;
end
if nargin < 7 || isempty(maxSteps)
    maxSteps = max(40, 80*numel(tau0));
end

t = dt*(0:n-1)';
obj = @(p) tailfitObjective(p, t, y, mode);
if exist('Simplex', 'file') == 2
    tau = Simplex(obj, tau0, lb, ub, 1e-5, maxSteps);
else
    tau = fminsearch(@(p) obj(boundVector(p(:), lb, ub)), tau0, ...
        optimset('Display', 'off', 'MaxIter', maxSteps));
    tau = boundVector(tau(:), lb, ub);
end

[err, coeff, z] = tailfitObjective(tau, t, y, mode);
offset = coeff(1);
amp = coeff(2:end);

if bild
    figure;
    subplot(4,1,1:3);
    semilogy(t, y, 'ob', t, max(z, eps), 'r', 'LineWidth', 1.5);
    axis tight
    xlabel('time [ns]');
    ylabel('count');
    subplot(4,1,4);
    plot(t, (y-z)./sqrt(max(abs(z), eps)), t, 0*t, 'k:');
    axis tight
    xlabel('time [ns]');
    ylabel('weighted residual');
end
end

function [err, coeff, z] = tailfitObjective(tau, t, y, mode)
tau = real(double(tau(:)));
if isempty(tau) || any(~isfinite(tau)) || any(tau <= 0)
    err = Inf;
    coeff = [];
    z = [];
    return;
end

M = [ones(numel(t), 1) exp(-t ./ tau(:).')];
switch mode
    case {'pirls', 'pirlsnonneg', 'poisson_irls', 'mle', 'ml', 'maximum_likelihood'}
        coeff = PIRLSnonneg(M, y);
    case {'ls', 'least_squares'}
        coeff = lsqnonneg(M, y);
    otherwise
        error('Tailfit mode must be ''pirls'', ''mle'', or ''ls''.');
end
z = M*coeff;
err = sum((z-y).^2./max(abs(z), eps))/numel(y);
if ~isscalar(err) || ~isfinite(err) || abs(imag(err)) > 1e-10*max(1, abs(real(err)))
    err = Inf;
else
    err = real(err);
end
end

function [lb, ub] = tailLifetimeBounds(tau0, limits, dt, n)
nTau = numel(tau0);
if isempty(limits)
    lb = max(dt/4, 0.01) * ones(nTau, 1);
    ub = max([20*tau0(:), 10*ones(nTau, 1), dt*n*ones(nTau, 1)], [], 2);
    return;
end

limits = real(double(limits));
if isvector(limits) && numel(limits) == 2
    lb = limits(1) * ones(nTau, 1);
    ub = limits(2) * ones(nTau, 1);
elseif isvector(limits) && numel(limits) == 2*nTau
    limits = limits(:);
    lb = limits(1:nTau);
    ub = limits(nTau+1:end);
elseif isequal(size(limits), [nTau 2])
    lb = limits(:, 1);
    ub = limits(:, 2);
else
    error('Tailfit limits must be [lower upper], [lower... upper...], or nTau-by-2.');
end

lb(~isfinite(lb) | lb <= 0) = max(dt/4, 0.01);
badUpper = ~isfinite(ub) | ub <= lb;
if any(badUpper)
    ub(badUpper) = max([20*tau0(badUpper), dt*n*ones(sum(badUpper), 1)], [], 2);
end
end

function x = boundVector(x, lb, ub)
x = max(lb(:), min(ub(:), real(double(x(:)))));
end
