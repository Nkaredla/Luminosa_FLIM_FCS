function err = pirlsfit(param, irf, y, p)
%PIRLSFIT Fluofit objective using PIRLSnonneg for nonnegative amplitudes.
%   Used by Fluofit when fitMode is 'pirls'.

persistent nnlsOptions
if isempty(nnlsOptions)
    nnlsOptions = optimset('Display', 'off', 'TolX', 1e-8, 'MaxIter', 1000);
end

n = length(irf);
if ~isscalar(p) || ~isfinite(p) || abs(imag(p)) > 1e-10*max(1, abs(real(p))) || real(p) <= 0
    err = 1e300;
    return
end
if any(~isfinite(param(:))) || any(~isfinite(irf(:))) || any(~isfinite(y(:))) || ...
        any(abs(imag(param(:))) > 1e-10 * max(1, abs(real(param(:))))) || ...
        any(abs(imag(irf(:))) > 1e-10 * max(1, abs(real(irf(:))))) || ...
        any(abs(imag(y(:))) > 1e-10 * max(1, abs(real(y(:)))))
    err = 1e300;
    return
end

p = real(p);
param = real(param);
irf = real(irf);
y = real(y(:));

t = 1:n;
tp = (1:p)';
c = param(1);
tau = param(2:end);
tau = tau(:).';

if any(tau <= 0)
    err = 1e300;
    return
end

x = exp(-(tp-1)*(1./tau))*diag(1./(1-exp(-p./tau)));
irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + ...
    (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
M = Convol(irs, x);
M = [ones(size(M,1),1) M];

if any(~isfinite(M(:))) || any(abs(imag(M(:))) > 1e-10 * max(1, abs(real(M(:)))))
    err = 1e300;
    return
end
M = real(M);
if any(M(:) < -1e-12 * max(1, max(abs(M(:)))))
    err = 1e300;
    return
end
M(M < 0) = 0;

try
    if exist('PIRLSnonneg', 'file') == 2
        A = PIRLSnonneg(M, y, 10);
    else
        [A, ~, ~, exitflag] = lsqnonneg(M, y, nnlsOptions);
        if exitflag <= 0
            err = 1e300;
            return
        end
    end
catch
    err = 1e300;
    return
end

if any(~isfinite(A(:))) || any(abs(imag(A(:))) > 1e-10 * max(1, abs(real(A(:)))))
    err = 1e300;
    return
end
A = max(real(A), 0);

z = M * A;
if any(~isfinite(z(:))) || any(abs(imag(z(:))) > 1e-10 * max(1, abs(real(z(:)))))
    err = 1e300;
    return
end
z = real(z);

ind = y > 0;
if any(z(ind) <= 0) || any(z(~ind) < 0)
    err = 1e300;
    return
end

denom = max(1, n - numel(tau));
err = sum(y(ind).*log(y(ind)./z(ind)) - y(ind) + z(ind)) / denom;
end
