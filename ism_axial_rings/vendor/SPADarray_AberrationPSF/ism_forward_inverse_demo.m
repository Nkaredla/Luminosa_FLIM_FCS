function ism_forward_inverse_demo(outDir)
% Physics-based forward/inverse ISM bead model with Enderlein/Bessel propagation.
% MATLAB companion to the Python script used for the manuscript figures.
%
% The forward model computes the raw detector-channel stack for the configured array.
% The inverse demo performs a local Gauss-Newton fit for a single active mode
% plus lateral bead shift, using the full detector-resolved stack.
%
% This file is intentionally written as a single self-contained script with
% local functions so that it is easy to adapt to existing ISM pipelines.

if nargin < 1 || isempty(outDir)
    outDir = 'output_matlab';
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

sim = defaultParams();
makeModeFigure(sim, fullfile(outDir, 'figure1_combined_matlab.png'));
makeRetrievalFigure(sim, fullfile(outDir, 'retrieval_validation_matlab.png'));

end

function sim = defaultParams()
sim.modeOrder = {'tilt_x','tilt_y','defocus','astig_x','astig_y','coma_x','coma_y','spherical'};
sim.lamExc = 0.488;
sim.lamEm = 0.520;
sim.lamRef = 0.520;
sim.NA = 1.2;
sim.nMedium = 1.33;
sim.fovXY = 1.8;
sim.nzRange = 0.8;
sim.nx = 27;
sim.nz = 5;
sim.beadRadius = 0.08;
sim.detPitch = 0.18;
sim.detSize = sim.detPitch;
sim.detectorLayout = 'honeycomb23';
sim.detectorPixelShape = 'hex';
sim.detectorHexRadius = sim.detSize / sqrt(3);
sim.Nr = 36;
sim.Nphi = 72;
sim.M = 5;
sim.x = linspace(-sim.fovXY/2, sim.fovXY/2, sim.nx);
sim.y = sim.x;
sim.z = linspace(-sim.nzRange/2, sim.nzRange/2, sim.nz);
sim.dx = abs(sim.x(2)-sim.x(1));
[sim.detXY, sim.detectorIndexGrid, sim.detectorLayoutInfo] = ...
    detectorLayout(sim.detectorLayout, sim.detPitch);
sim.nDet = size(sim.detXY, 1);
sim.detectorGridSize = size(sim.detectorIndexGrid);
sim.arrayN = sim.detectorGridSize(1);
sim.obj = beadObject3D(sim);
end

function obj = beadObject3D(sim)
[X, Y] = meshgrid(sim.x, sim.y);
obj = zeros(numel(sim.y), numel(sim.x), numel(sim.z));
for iz = 1:numel(sim.z)
    obj(:,:,iz) = double(X.^2 + Y.^2 + sim.z(iz).^2 <= sim.beadRadius^2);
end
obj = obj / sum(obj(:));
end

function Z = zernikeModes(rho, phi)
Z.tilt_x = 2*rho.*cos(phi);
Z.tilt_y = 2*rho.*sin(phi);
Z.astig_x = sqrt(6)*rho.^2.*cos(2*phi);
Z.astig_y = sqrt(6)*rho.^2.*sin(2*phi);
Z.defocus = sqrt(3)*(2*rho.^2 - 1);
Z.coma_x = sqrt(8)*(3*rho.^3 - 2*rho).*cos(phi);
Z.coma_y = sqrt(8)*(3*rho.^3 - 2*rho).*sin(phi);
Z.spherical = sqrt(5)*(6*rho.^4 - 6*rho.^2 + 1);
end

function coeffs = coeffStruct(sim, a)
coeffs = struct();
for k = 1:numel(sim.modeOrder)
    if abs(a(k)) > 1e-15
        coeffs.(sim.modeOrder{k}) = a(k);
    end
end
end

function [rho, ms, pm] = pupilFourierCoeffs(sim, coeffs, wavelength)
rho = linspace(0, 1, sim.Nr);
phi = linspace(0, 2*pi, sim.Nphi+1); phi(end) = [];
[R, P] = ndgrid(rho, phi);
W = zeros(size(R));
Z = zernikeModes(R, P);
fn = fieldnames(coeffs);
for k = 1:numel(fn)
    W = W + coeffs.(fn{k}) * Z.(fn{k});
end
smax = sim.NA / sim.nMedium;
cosT = sqrt(max(1 - (smax*R).^2, 0));
pupil = sqrt(cosT) .* exp(1i*2*pi*(sim.lamRef/wavelength).*W);
A = fft(pupil, [], 2) / sim.Nphi;
ms = -sim.M:sim.M;
pm = complex(zeros(sim.Nr, numel(ms)));
for ii = 1:numel(ms)
    idx = mod(ms(ii), sim.Nphi) + 1;
    pm(:,ii) = A(:,idx);
end
end

function psf = scalarPSFBessel(sim, coeffs, wavelength)
[rho, ms, pm] = pupilFourierCoeffs(sim, coeffs, wavelength);
[X, Y] = meshgrid(sim.x, sim.y);
r = sqrt(X.^2 + Y.^2);
psi = atan2(Y, X);
rv = r(:);
ur = unique(round(rv, 12));
[~, inv] = ismember(round(rv, 12), ur);
k0 = 2*pi / wavelength;
smax = sim.NA / sim.nMedium;
cosT = sqrt(max(1 - (smax*rho).^2, 1e-12));
weight = (smax^2) * rho ./ cosT;
arg = k0 * sim.NA * (ur * rho);
Jtab = cell(numel(ms),1);
for ii = 1:numel(ms)
    m = ms(ii);
    J = besselj(abs(m), arg);
    if m < 0
        J = ((-1)^(-m)) * J;
    end
    Jtab{ii} = J;
end
psf = zeros(numel(sim.y), numel(sim.x), numel(sim.z));
psiV = psi(:);
for iz = 1:numel(sim.z)
    axial = exp(1i * k0 * sim.nMedium * sim.z(iz) * cosT);
    radial = complex(zeros(numel(ur), numel(ms)));
    for ii = 1:numel(ms)
        basis = (pm(:,ii).').*axial.*weight;
        radial(:,ii) = trapz(rho, Jtab{ii} .* basis, 2);
    end
    U = complex(zeros(size(rv)));
    for ii = 1:numel(ms)
        m = ms(ii);
        U = U + (2*pi) * (1i^m) * exp(1i*m*psiV) .* radial(inv,ii);
    end
    psf(:,:,iz) = abs(reshape(U, size(r))).^2;
end
psf = psf / sum(psf(:));
end

function out = fourierShift2D(img, shiftXPx, shiftYPx)
[ny, nx] = size(img);
ux = ifftshift((-floor(nx/2):ceil(nx/2)-1)/nx);
uy = ifftshift((-floor(ny/2):ceil(ny/2)-1)/ny);
[Ux, Uy] = meshgrid(ux, uy);
phase = exp(-2i*pi*(shiftXPx*Ux + shiftYPx*Uy));
out = real(ifft2(fft2(img).*phase));
end

function out = boxBlurStack(stack, detPix)
if mod(detPix,2)==0, detPix = detPix+1; end
ker = ones(detPix, detPix);
ker = ker / sum(ker(:));
out = zeros(size(stack));
for iz = 1:size(stack,3)
    out(:,:,iz) = conv2(stack(:,:,iz), ker, 'same');
end
end

function stack = normalizedStack(sim, coeffs, x0, y0)
if nargin < 3, x0 = 0; end
if nargin < 4, y0 = 0; end
hExc = scalarPSFBessel(sim, coeffs, sim.lamExc);
hDet = scalarPSFBessel(sim, coeffs, sim.lamEm);
detPSF = detectorBlurStack(hDet, sim);
chans = zeros(sim.nx, sim.nx, size(sim.detXY,1));
for k = 1:size(sim.detXY,1)
    sx = sim.detXY(k,1); sy = sim.detXY(k,2);
    hEff = zeros(size(hExc));
    for iz = 1:size(hExc,3)
        shifted = max(fourierShift2D(detPSF(:,:,iz), sx/sim.dx, sy/sim.dx), 0);
        hEff(:,:,iz) = hExc(:,:,iz) .* shifted;
    end
    hEff = hEff / sum(hEff(:));
    img = zeros(sim.nx, sim.nx);
    for iz = 1:size(sim.obj,3)
        % IMPORTANT: use correlation, not direct convolution
        img = img + conv2(sim.obj(:,:,iz), rot90(hEff(:,:,iz),2), 'same');
    end
    chans(:,:,k) = img;
end
stack = chans / sum(chans(:));
if abs(x0) > 0 || abs(y0) > 0
    shifted = zeros(size(stack));
    for k = 1:size(stack,3)
        shifted(:,:,k) = max(fourierShift2D(stack(:,:,k), x0/sim.dx, y0/sim.dx), 0);
    end
    stack = shifted / sum(shifted(:));
end
end

function [dataN, beta] = normalizeDataStack(data)
border = [reshape(data(1,:,:), [], 1); reshape(data(end,:,:), [], 1); reshape(data(:,1,:), [], 1); reshape(data(:,end,:), [], 1)];
beta = max(median(border), 0);
d = max(data - beta, 0);
s = sum(d(:));
if s <= 0
    dataN = d;
else
    dataN = d / s;
end
end

function res = singleModeFit(sim, modeName, amp, x0, y0, photons, bg, maxIter)
if nargin < 4, x0 = 0.02; end
if nargin < 5, y0 = -0.015; end
if nargin < 6, photons = 2e5; end
if nargin < 7, bg = 0.1; end
if nargin < 8, maxIter = 12; end
midx = find(strcmp(sim.modeOrder, modeName));
trueA = zeros(1, numel(sim.modeOrder));
trueA(midx) = amp;
mTrue = normalizedStack(sim, coeffStruct(sim, trueA), x0, y0);
data = poissrnd(photons*mTrue + bg);
[dataN, beta] = normalizeDataStack(data);
p = zeros(1, numel(sim.modeOrder)+2);
fd = 0.01;
stepXY = sim.dx/4;
for it = 1:maxIter
    m0 = normalizedStack(sim, coeffStruct(sim, p(1:numel(sim.modeOrder))), p(end-1), p(end));
    r = dataN(:) - m0(:);
    pp = p; pm = p; pp(midx)=pp(midx)+fd; pm(midx)=pm(midx)-fd;
    mp = normalizedStack(sim, coeffStruct(sim, pp(1:numel(sim.modeOrder))), pp(end-1), pp(end));
    mm = normalizedStack(sim, coeffStruct(sim, pm(1:numel(sim.modeOrder))), pm(end-1), pm(end));
    J = zeros(numel(m0), 3);
    J(:,1) = (mp(:)-mm(:))/(2*fd);
    for q = 1:2
        pp = p; pm = p; pp(end-2+q)=pp(end-2+q)+stepXY; pm(end-2+q)=pm(end-2+q)-stepXY;
        mp = normalizedStack(sim, coeffStruct(sim, pp(1:numel(sim.modeOrder))), pp(end-1), pp(end));
        mm = normalizedStack(sim, coeffStruct(sim, pm(1:numel(sim.modeOrder))), pm(end-1), pm(end));
        J(:,1+q) = (mp(:)-mm(:))/(2*stepXY);
    end
    H = J.'*J + diag([1e-4 1e-5 1e-5]);
    g = J.'*r;
    if ~all(isfinite(H(:))) || ~all(isfinite(g(:)))
        error('Non-finite values encountered in Gauss-Newton system.');
    end
    delta = H \ g;
    p(midx) = p(midx) + delta(1);
    p(end-1) = p(end-1) + delta(2);
    p(end) = p(end) + delta(3);
end
res.modeName = modeName;
res.trueAmp = amp;
res.estAmp = p(midx);
res.trueXY = [x0 y0];
res.estXY = [p(end-1) p(end)];
res.beta = beta;
end

function makeModeFigure(sim, outFile)
scenarios = {
    'None', struct();
    'Defocus', struct('defocus', 0.18);
    'Astig x', struct('astig_x', 0.18);
    'Coma x', struct('coma_x', 0.18);
    'Spherical', struct('spherical', 0.18)};
figure('Color','w','Position',[100 100 1300 480]);
for j = 1:size(scenarios,1)
    stack = normalizedStack(sim, scenarios{j,2});
    central = stack(:,:,ceil(size(stack,3)/2));
    s = sum(stack,3);
    [~, idx] = max(s(:));
    [iy, ix] = ind2sub(size(s), idx);
    idxGrid = detectorIndexGridFromXY(sim.detXY, false, true, false);
    vals = squeeze(stack(iy, ix, :));
    micro = nan(size(idxGrid));
    validGrid = isfinite(idxGrid) & idxGrid >= 1 & idxGrid <= numel(vals);
    micro(validGrid) = vals(idxGrid(validGrid));
    subplot(2, size(scenarios,1), j);
    imagesc(sim.x, sim.y, central); axis image; title(scenarios{j,1}); xlabel('x (\mum)'); ylabel('y (\mum)'); colorbar;
    subplot(2, size(scenarios,1), j + size(scenarios,1));
    imagesc(micro); axis image; xlabel('detector x'); ylabel('detector y'); colorbar;
end
sgtitle('Detector-resolved signatures for representative aberrations (0.18 waves RMS)');
exportgraphics(gcf, outFile, 'Resolution', 180);
close(gcf);
end

function makeRetrievalFigure(sim, outFile)
modes = {'defocus','astig_x','coma_x','spherical'};
amp = 0.03;
ests = zeros(size(modes));
for k = 1:numel(modes)
    res = singleModeFit(sim, modes{k}, amp);
    if ismember(modes{k}, {'defocus','astig_x','spherical'})
        ests(k) = abs(res.estAmp);
    else
        ests(k) = res.estAmp;
    end
end
figure('Color','w','Position',[100 100 650 380]);
xs = 1:numel(modes);
bar(xs-0.17, repmat(amp, size(xs)), 0.34); hold on;
bar(xs+0.17, ests, 0.34);
set(gca, 'XTick', xs, 'XTickLabel', {'defocus','astig. x','coma x','spherical'});
ylabel('aberration amplitude (waves RMS)');
title(sprintf('Single-mode recovery from the raw %d-channel stack', size(sim.detXY,1)));
legend({'true','retrieved'}, 'Box', 'off', 'Location', 'northwest');
grid on;
exportgraphics(gcf, outFile, 'Resolution', 180);
close(gcf);
end
