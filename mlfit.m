function err = mlfit(param, irf, y, p)
%	MLFIT(param, irf, y, p) returns the ML error between the data y
%	and the computed values.
%	MLFIT assumes a function of the form:
%
%	  y =  yoffset + A(1)*convol(irf,exp(-t/tau(1)/(1-exp(-p/tau(1)))) + ...
%
%	param(1) is the color shift value between irf and y.
%	param(2) is the irf offset.
%	param(3:...) are the decay times.
%	irf is the measured Instrumental Response Function.
%	y is the measured fluorescence decay curve.
%	p is the time between to laser excitations (in number of TCSPC channels).


% global bild Plothandle
bild = 0;
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
y = real(y);

t = 1:n;
tp = (1:p)';
c = param(1);
tau = param(2:length(param)); tau = tau(:)';

if any(tau <= 0)
    err = 1e300;
    return
end

x = exp(-(tp-1)*(1./tau))*diag(1./(1-exp(-p./tau)));
irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
z = Convol(irs, x);
z = [ones(size(z,1),1) z];

if any(~isfinite(z(:))) || any(abs(imag(z(:))) > 1e-10 * max(1, abs(real(z(:)))))
    err = 1e300;
    return
end
z = real(z);
if any(z(:) < -1e-12 * max(1, max(abs(z(:)))))
    err = 1e300;
    return
end
z(z < 0) = 0;

[A, ~, ~, exitflag] = lsqnonneg(z, y, nnlsOptions);
if exitflag <= 0
    err = 1e300;
    return
end
if any(~isfinite(A(:))) || any(abs(imag(A(:))) > 1e-10 * max(1, abs(real(A(:)))))
    err = 1e300;
    return
end
A = real(A);

z = z*A;
if any(~isfinite(z(:))) || any(abs(imag(z(:))) > 1e-10 * max(1, abs(real(z(:)))))
    err = 1e300;
    return
end
z = real(z);

if bild
    % 	set(Plothandle,'ydata',z)
    figure(Plothandle)
    hold all
    plot(t,log10(z));
    hold off
    drawnow
end
ind = y>0;
if any(z(ind) <= 0) || any(z(~ind) < 0)
    err = 1e300;
    return
end

denom = max(1, n-length(tau));
err = sum(y(ind).*log(y(ind)./z(ind))-y(ind)+z(ind))/denom;



