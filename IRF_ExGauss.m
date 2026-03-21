function z = IRF_ExGauss(p, t, pic)
% IRF_ExGauss
%
% Ex-Gaussian IRF model:
%   h(t) = Gaussian(t0, sigma) (*) one-sided exponential(tauR)
%
% Parameters
% ----------
% p(1) = t0    : temporal shift [ns]
% p(2) = sigma : Gaussian width [ns]
% p(3) = tauR  : exponential right-tail constant [ns]
%
% Input
% -----
% t   : time axis [ns]
% pic : optional plotting flag
%
% Output
% ------
% z   : normalized IRF column vector, sum(z)=1

t = t(:);
p = p(:);

t0    = p(1);
sigma = max(p(2), 1e-6);
tauR  = max(p(3), 1e-6);

x = t - t0;

% Ex-Gaussian PDF (Gaussian convolved with one-sided exponential)
if tauR < 1e-5
    h = exp(-0.5 * (x ./ sigma).^2);
else
    arg = (sigma.^2 ./ tauR - x) ./ (sqrt(2) * sigma);
    h = (1 ./ (2 * tauR)) .* exp((sigma.^2) ./ (2 * tauR.^2) - x ./ tauR) .* erfc(arg);
end

h(~isfinite(h)) = 0;
h(h < 0) = 0;

% Fallback if everything underflows
if sum(h) <= 0
    [~, ind] = min(abs(x));
    h = zeros(size(t));
    h(ind) = 1;
end

z = h(:) ./ sum(h);

if nargin > 2 && ~isempty(pic)
    if pic == 1
        plot(t, z, 'r', 'LineWidth', 1.5); drawnow
    else
        semilogy(t, max(z, eps), 'r', 'LineWidth', 1.5); drawnow
    end
end
end