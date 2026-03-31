function [err, c, z, zz] = lsFCS_diff1_fit(p, t, xi, G, pixel_um)
% One-diffusion spatio-temporal correlation model for lsFCS / lsCrossRead.
%
% p(1) = tauD (s)
% p(2) = w0 (um), where w0 is the Gaussian standard deviation (sigma)
%
% t  : time axis (s), column or row
% xi : spatial lag (pixels), column or row
% G  : correlation map (numel(t) x numel(xi))
% pixel_um : pixel size (um)

tauD = p(1);
w0   = p(2); % sigma of Gaussian (not 1/e^2 radius)

if tauD <= 0 || w0 <= 0 || any(~isfinite(p))
    err = Inf;
    c = [];
    z = [];
    zz = [];
    return;
end

[xixi, tt] = meshgrid(xi * pixel_um, t);

den = 1 + tt ./ tauD;
tmp = exp(-(xixi.^2) ./ (2 * (w0.^2) .* den)) ./ den;

basis = tmp ./ max(tmp(:));
A = [basis(:), ones(numel(basis), 1)];

c = lsqnonneg(A, G(:));
z = c(1) * basis + c(2);
zz = cat(3, basis, ones(size(basis)));

resid = G(:) - z(:);
err = sum((resid.^2) ./ max(abs(z(:)), eps));
