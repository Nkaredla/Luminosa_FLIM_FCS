function z = IRF_GammaShifted(p, t, pic)
% IRF_GammaShifted
%
% Shifted gamma IRF model:
%   h(t) = H(t-t0) * ((t-t0)^(alpha-1) * exp(-(t-t0)/beta)) / (gamma(alpha)*beta^alpha)
%
% Parameters
% ----------
% p(1) = t0    : temporal shift [ns]
% p(2) = alpha : gamma shape (>0)
% p(3) = beta  : gamma scale [ns] (>0)
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
alpha = max(p(2), 1e-6);
beta  = max(p(3), 1e-6);

x = t - t0;

h = zeros(size(t));
ind = (x >= 0);

if any(ind)
    xx = x(ind);

    % Use log form for numerical stability
    logh = (alpha - 1) .* log(max(xx, 1e-15)) ...
         - xx ./ beta ...
         - gammaln(alpha) ...
         - alpha .* log(beta);

    h(ind) = exp(logh);
end

h(~isfinite(h)) = 0;
h(h < 0) = 0;

% Fallback if underflow wipes everything out
if sum(h) <= 0
    [~, idx] = min(abs(x));
    h = zeros(size(t));
    h(idx) = 1;
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