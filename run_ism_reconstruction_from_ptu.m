function results = run_ism_reconstruction_from_ptu(ptuOut, params)
% RUN_ACO_ISM_FROM_PTU
%
% APR + ACO-ISM reconstruction compatible with PTU_MultiFrameScanRead output.
%
% Main features
% -------------
% 1) extracts detector image stack from ptuOut.tags or ptuOut.tag
% 2) estimates APR shifts by phase correlation
% 3) identifies the center detector automatically from the shift cloud
% 4) builds APR image (same size as raw images)
% 5) builds averaged autocorrelation from RAW detector images
% 6) masks the zero-shift spike
% 7) inverts the averaged autocorrelation with Schultz-Snyder,
%    initialized from APR
% 8) plots properly labeled OTF diagnostics (for reference only)
%
% INPUT
% -----
% ptuOut : struct from PTU_MultiFrameScanRead
%
% params fields (all optional)
% ------------------------------------------
% .imageSource          : 'tags' or 'tag'           (default 'tags')
% .frameIndices         : frames for imageSource='tag' (default all)
% .frameCombine         : 'sum' or 'mean'           (default 'sum')
% .dropEmptyChannels    : true/false                (default true)
%
% .centerDetectorIndex  : manually set center detector index (default auto)
% .tempReferenceIndex   : temporary reference for first-pass registration
%                         (default brightest channel)
%
% .smoothSigma          : Gaussian sigma before registration (default 1)
% .useWindow            : true/false Hann window for registration (default true)
% .normalizeImages      : true/false                (default true)
% .upsampleReg          : registration refinement factor (default 20)
%
% .nIter                : Schultz-Snyder iterations (default 1000)
% .checkEvery           : convergence check period   (default 25)
% .stopTol              : relative RMSE tolerance    (default 1e-7)
% .minIter              : minimum iterations         (default 50)
% .supportMask          : logical support mask, or [] (default [])
% .preserveFlux         : true/false normalize to APR flux (default true)
%
% .pixelSize            : object pixel size, e.g. nm (for OTF diagnostics)
% .lambda               : emission wavelength, same units as pixelSize
% .NA                   : objective NA
%
% .showPlots            : true/false                (default true)
%
% OUTPUT
% ------
% results struct with fields:
% .imgStack
% .channelIDs
% .tempReferenceIndex
% .centerDetectorIndex
% .shiftsToTempRef         [Nd x 2], shift-to-apply to align to temp ref
% .detectorPositions       [Nd x 2], relative detector geometry from shift cloud
% .shiftsToCenter          [Nd x 2], shift-to-apply to align each image to center
% .rawSum
% .aprImage
% .acoAverage
% .acoAverageMasked
% .acoImage
% .convHistory
% .otfWF
% .otfISMexp
% .otfISMideal

if nargin < 2 || isempty(params)
    params = struct();
end

params = setDefault(params, 'imageSource', 'tags');
params = setDefault(params, 'frameIndices', []);
params = setDefault(params, 'frameCombine', 'sum');
params = setDefault(params, 'dropEmptyChannels', true);

params = setDefault(params, 'centerDetectorIndex', NaN);
params = setDefault(params, 'tempReferenceIndex', NaN);

params = setDefault(params, 'smoothSigma', 1);
params = setDefault(params, 'useWindow', true);
params = setDefault(params, 'normalizeImages', true);
params = setDefault(params, 'upsampleReg', 20);

params = setDefault(params, 'nIter', 1000);
params = setDefault(params, 'checkEvery', 25);
params = setDefault(params, 'stopTol', 1e-7);
params = setDefault(params, 'minIter', 50);
params = setDefault(params, 'supportMask', []);
params = setDefault(params, 'preserveFlux', true);

params = setDefault(params, 'pixelSize', 20);
params = setDefault(params, 'lambda', 650);
params = setDefault(params, 'NA', 1.4);

params = setDefault(params, 'showPlots', true);

    function s = setDefault(s, fieldName, defaultValue)
        if ~isfield(s, fieldName) || isempty(s.(fieldName))
            s.(fieldName) = defaultValue;
        end
    end

% ------------------------------------------------------------
% 1) Extract detector image stack from PTU output
% ------------------------------------------------------------
[imgStack, channelIDs] = extract_img_stack_from_ptu(ptuOut, params);
[H, W, Nd] = size(imgStack);

% ------------------------------------------------------------
% 2) Raw detector sum
% ------------------------------------------------------------
rawSum = sum(imgStack, 3);

% ------------------------------------------------------------
% 3) First-pass registration against a temporary reference
%    Used only to recover detector geometry and identify center detector
% ------------------------------------------------------------
procStack = zeros(H, W, Nd);
for k = 1:Nd
    procStack(:,:,k) = preprocessForRegistration( ...
        imgStack(:,:,k), params.smoothSigma, params.useWindow, params.normalizeImages);
end

if isnan(params.tempReferenceIndex)
    tempRefIdx = chooseBrightestDetector(imgStack);
else
    tempRefIdx = params.tempReferenceIndex;
end

shiftsToTempRef = zeros(Nd, 2); % shift-to-apply to align channel k to temp ref
refProc = procStack(:,:,tempRefIdx);
for k = 1:Nd
    shiftsToTempRef(k,:) = phaseCorrShiftToApply(refProc, procStack(:,:,k), params.upsampleReg);
end
shiftsToTempRef(tempRefIdx,:) = [0 0];

% Relative detector positions up to translation/sign convention.
% Using the shift cloud directly is enough to find the center detector:
% the center is the point nearest the centroid.
detectorPositions = shiftsToTempRef;
centroid = mean(detectorPositions, 1);

if isnan(params.centerDetectorIndex)
    d2 = sum((detectorPositions - centroid).^2, 2);
    [~, centerIdx] = min(d2);
else
    centerIdx = params.centerDetectorIndex;
end

% ------------------------------------------------------------
% 4) Final APR shifts: register all channels to the center detector
% ------------------------------------------------------------
shiftsToCenter = zeros(Nd, 2); % shift-to-apply to align each image to center image
refCenterProc = procStack(:,:,centerIdx);
for k = 1:Nd
    shiftsToCenter(k,:) = phaseCorrShiftToApply(refCenterProc, procStack(:,:,k), params.upsampleReg);
end
shiftsToCenter(centerIdx,:) = [0 0];

% ------------------------------------------------------------
% 5) APR image (same size as raw)
%    This follows i_ISM(xs) = sum_d i(xs + mu(xd) | xd)
%    with mu estimated by phase correlation.
% ------------------------------------------------------------
aprImage = zeros(H, W);
for k = 1:Nd
    aprImage = aprImage + shiftImageFourier(double(imgStack(:,:,k)), ...
        shiftsToCenter(k,1), shiftsToCenter(k,2));
end
aprImage = max(aprImage, 0);

% ------------------------------------------------------------
% 6) Average autocorrelation of RAW detector images
%    ACO-ISM uses autocorrelation of each raw detector image,
%    averages them, then inverts the average autocorrelation.
% ------------------------------------------------------------
acoAverage = zeros(H, W);
for k = 1:Nd
    Ik = double(imgStack(:,:,k));
    Ik = Ik - min(Ik(:));   % as described in Ancora et al.
    Ik = max(Ik, 0);
    acoAverage = acoAverage + autocorr2_fft(Ik);
end
acoAverage = acoAverage / Nd;

% mask/replace zero-shift spike
acoAverageMasked = mask_zero_shift_pixel(acoAverage);

% ------------------------------------------------------------
% 7) Schultz-Snyder inversion initialized from APR
% ------------------------------------------------------------
[acoImage, convHistory] = schultz_snyder_from_apr( ...
    acoAverageMasked, aprImage, params);

% ------------------------------------------------------------
% 8) OTF diagnostics only (not used by ACO reconstruction)
%    Based on PNAS focal-plane approximation:
%    Uexp ~ D(q/2)^2, Uideal ~ D(q/2)
% ------------------------------------------------------------
otfWF = theoreticalIncoherentOTF(H, W, params.pixelSize, params.lambda, params.NA);
otfISMideal = rescaleOTFSupport(otfWF, 0.5);
otfISMexp = otfISMideal .^ 2;

otfWF = normalizeOTFdc(otfWF);
otfISMideal = normalizeOTFdc(otfISMideal);
otfISMexp = normalizeOTFdc(otfISMexp);

% ------------------------------------------------------------
% 9) Collect outputs
% ------------------------------------------------------------
results = struct();
results.imgStack = imgStack;
results.channelIDs = channelIDs;
results.tempReferenceIndex = tempRefIdx;
results.centerDetectorIndex = centerIdx;
results.shiftsToTempRef = shiftsToTempRef;
results.detectorPositions = detectorPositions;
results.shiftsToCenter = shiftsToCenter;
results.rawSum = rawSum;
results.aprImage = aprImage;
results.acoAverage = acoAverage;
results.acoAverageMasked = acoAverageMasked;
results.acoImage = acoImage;
results.convHistory = convHistory;
results.otfWF = otfWF;
results.otfISMexp = otfISMexp;
results.otfISMideal = otfISMideal;
results.paramsUsed = params;

if params.showPlots
    show_aco_ism_results(results);
end
end


% =========================================================================
% PTU extraction
% =========================================================================
function [imgStack, channelIDs] = extract_img_stack_from_ptu(ptuOut, params)

switch lower(params.imageSource)
    case 'tags'
        assert(isfield(ptuOut, 'tags'), 'ptuOut.tags not found.');
        imgStack = double(ptuOut.tags);

    case 'tag'
        assert(isfield(ptuOut, 'tag'), 'ptuOut.tag not found.');
        assert(ndims(ptuOut.tag) == 4, 'ptuOut.tag must be H x W x Nd x Nf.');

        nFrames = size(ptuOut.tag, 4);
        if isempty(params.frameIndices)
            frameIdx = 1:nFrames;
        else
            frameIdx = params.frameIndices(:)';
            assert(all(frameIdx >= 1 & frameIdx <= nFrames), 'frameIndices out of range.');
        end

        switch lower(params.frameCombine)
            case 'sum'
                imgStack = sum(double(ptuOut.tag(:,:,:,frameIdx)), 4);
            case 'mean'
                imgStack = mean(double(ptuOut.tag(:,:,:,frameIdx)), 4);
            otherwise
                error('frameCombine must be ''sum'' or ''mean''.');
        end

    otherwise
        error('imageSource must be ''tags'' or ''tag''.');
end

if isfield(ptuOut, 'dind')
    channelIDs = double(ptuOut.dind(:));
    if numel(channelIDs) ~= size(imgStack, 3)
        channelIDs = (1:size(imgStack,3)).';
    end
else
    channelIDs = (1:size(imgStack,3)).';
end

if params.dropEmptyChannels
    s = squeeze(sum(sum(imgStack, 1), 2));
    keep = s > 0;
    imgStack = imgStack(:,:,keep);
    channelIDs = channelIDs(keep);
end
end


% =========================================================================
% Registration / APR
% =========================================================================
function idx = chooseBrightestDetector(imgStack)
Nd = size(imgStack, 3);
s = zeros(Nd,1);
for k = 1:Nd
    s(k) = sum(imgStack(:,:,k), 'all');
end
[~, idx] = max(s);
end

function Iout = preprocessForRegistration(I, smoothSigma, useWindow, normalizeImages)
I = double(I);

if smoothSigma > 0
    I = imgaussfilt(I, smoothSigma);
end

I = I - median(I(:));

if normalizeImages
    sd = std(I(:));
    if sd > 0
        I = I / sd;
    end
end

if useWindow
    [H, W] = size(I);
    wy = hann(H, 'periodic');
    wx = hann(W, 'periodic');
    I = I .* (wy * wx.');
end

Iout = I;
end

function shiftToApply = phaseCorrShiftToApply(refImg, movImg, upsample)
% Returns the shift [dy dx] to APPLY to movImg in order to align it to refImg.
%
% Robustness trick:
%   use phase correlation peak for coarse estimate,
%   test both signs,
%   refine with fminsearch on the direct overlap objective.

F1 = fft2(refImg);
F2 = fft2(movImg);

CPS = F1 .* conj(F2);
CPS = CPS ./ max(abs(CPS), eps);

cc = real(ifft2(CPS));
cc = fftshift(cc);

[H, W] = size(cc);
[~, idx] = max(cc(:));
[py, px] = ind2sub([H, W], idx);

cy = floor(H/2) + 1;
cx = floor(W/2) + 1;

dy0 = py - cy;
dx0 = px - cx;

% subpixel peak interpolation
dySub = 0;
dxSub = 0;

if py > 1 && py < H
    c1 = cc(py-1, px);
    c2 = cc(py,   px);
    c3 = cc(py+1, px);
    den = c1 - 2*c2 + c3;
    if abs(den) > 1e-12
        dySub = 0.5 * (c1 - c3) / den;
    end
end

if px > 1 && px < W
    c1 = cc(py, px-1);
    c2 = cc(py, px);
    c3 = cc(py, px+1);
    den = c1 - 2*c2 + c3;
    if abs(den) > 1e-12
        dxSub = 0.5 * (c1 - c3) / den;
    end
end

p0 = [dy0 + dySub, dx0 + dxSub];

% choose sign that gives better overlap
sc1 = sum(refImg .* shiftImageFourier(movImg,  p0(1),  p0(2)), 'all');
sc2 = sum(refImg .* shiftImageFourier(movImg, -p0(1), -p0(2)), 'all');
if sc2 > sc1
    p0 = -p0;
end

if upsample > 1
    fun = @(p) -sum(refImg .* shiftImageFourier(movImg, p(1), p(2)), 'all');
    opts = optimset('Display', 'off', 'TolX', 1/upsample, 'TolFun', 1e-7);
    p = fminsearch(fun, p0, opts);
else
    p = p0;
end

shiftToApply = p(:).';
end


% =========================================================================
% ACO / Schultz-Snyder
% =========================================================================
function A = autocorr2_fft(I)
% 2D autocorrelation with zero-shift at (1,1)
F = fft2(I);
A = real(ifft2(abs(F).^2));
end

function Aout = mask_zero_shift_pixel(Ain)
% Replace the autocorrelation zero-shift spike.
% Paper uses biharmonic inpainting; here we use regionfill if available,
% otherwise a local 8-neighbor average.

Ashift = fftshift(Ain);
[H, W] = size(Ashift);
cy = floor(H/2) + 1;
cx = floor(W/2) + 1;

if exist('regionfill', 'file') == 2
    mask = false(H, W);
    mask(cy, cx) = true;
    Ashift = regionfill(Ashift, mask);
else
    y1 = max(1, cy-1); y2 = min(H, cy+1);
    x1 = max(1, cx-1); x2 = min(W, cx+1);
    patch = Ashift(y1:y2, x1:x2);
    m = true(size(patch));
    m(cy-y1+1, cx-x1+1) = false;
    vals = patch(m);
    Ashift(cy, cx) = mean(vals(:), 'omitnan');
end

Aout = ifftshift(Ashift);
end

function [x, history] = schultz_snyder_from_apr(Iaco, aprImage, params)
% Schultz-Snyder fixed-point iteration initialized from APR image.

x = max(double(aprImage), 0);
targetFlux = sum(x(:));

if targetFlux <= 0
    error('APR image has zero total intensity.');
end

% normalize initialization but keep same total intensity afterward
sx = sum(x(:));
if sx > 0
    x = x / sx * targetFlux;
end

if ~isempty(params.supportMask)
    supportMask = logical(params.supportMask);
    assert(isequal(size(supportMask), size(x)), 'supportMask size mismatch.');
else
    supportMask = [];
end

nIter = params.nIter;
checkEvery = params.checkEvery;
stopTol = params.stopTol;
minIter = params.minIter;

history.iter = [];
history.relRMSE = [];

for it = 1:nIter
    xx = autocorr2_fft(x);                 % x ⋆ x
    ratio = Iaco ./ max(xx, eps);

    % convolution: x * ratio
    termConv = real(ifft2(fft2(x) .* fft2(ratio)));

    % correlation: x ⋆ ratio
    termCorr = real(ifft2(conj(fft2(x)) .* fft2(ratio)));

    x = x .* max(termConv + termCorr, 0);
    x = max(x, 0);

    if ~isempty(supportMask)
        x(~supportMask) = 0;
    end

    if params.preserveFlux
        sx = sum(x(:));
        if sx > 0
            x = x / sx * targetFlux;
        end
    end

    if (mod(it, checkEvery) == 0) || (it == 1) || (it == nIter)
        xxNow = autocorr2_fft(x);
        relRMSE = sqrt(mean((xxNow(:) - Iaco(:)).^2)) / max(sqrt(mean(Iaco(:).^2)), eps);

        history.iter(end+1,1) = it; %#ok<AGROW>
        history.relRMSE(end+1,1) = relRMSE; %#ok<AGROW>

        if it >= minIter && relRMSE < stopTol
            break;
        end
    end
end
end


% =========================================================================
% Fourier / OTF diagnostics
% =========================================================================
function Ishift = shiftImageFourier(I, dy, dx)
[H, W] = size(I);
[KY, KX] = frequencyGrid(H, W);
ramp = exp(-1i * 2*pi * (KX * dx + KY * dy));
Ishift = real(ifft2(fft2(I) .* ramp));
end

function [KY, KX] = frequencyGrid(H, W)
fy = ifftshift((-floor(H/2):ceil(H/2)-1) / H);
fx = ifftshift((-floor(W/2):ceil(W/2)-1) / W);
[KX, KY] = meshgrid(fx, fy);
end

function otf = theoreticalIncoherentOTF(H, W, pixelSize, lambda, NA)
% 2D scalar incoherent circular-pupil OTF
[KYpx, KXpx] = frequencyGrid(H, W);
FX = KXpx / pixelSize;
FY = KYpx / pixelSize;
FR = sqrt(FX.^2 + FY.^2);

fc = 2 * NA / lambda;
rho = FR / fc;

otf = zeros(H, W);
m = rho <= 1;
r = rho(m);
otf(m) = (2/pi) * (acos(r) - r .* sqrt(1 - r.^2));
otf = normalizeOTFdc(otf);
end

function out = rescaleOTFSupport(otf, alpha)
% out(k) = otf(alpha * k)
[H, W] = size(otf);
[x, y] = meshgrid(1:W, 1:H);
cx = (W + 1) / 2;
cy = (H + 1) / 2;

xs = alpha * (x - cx) + cx;
ys = alpha * (y - cy) + cy;

otfShift = fftshift(otf);
outShift = interp2(otfShift, xs, ys, 'linear', 0);
out = ifftshift(outShift);
end

function otf = normalizeOTFdc(otf)
dc = otf(1,1);
if abs(dc) > 0
    otf = otf / dc;
end
end


% =========================================================================
% Plotting
% =========================================================================
function show_aco_ism_results(res)

imgStack = res.imgStack;
[H, W, ~] = size(imgStack);

figure('Color', 'w', 'Name', 'APR + ACO-ISM reconstruction');

subplot(2,3,1);
imagesc(res.rawSum); axis image off; colormap hot;
title('Raw detector sum');

subplot(2,3,2);
imagesc(res.aprImage); axis image off; colormap hot;
title('APR image');

subplot(2,3,3);
imagesc(res.acoImage); axis image off; colormap hot;
title('ACO-ISM image');

subplot(2,3,4);
plot(res.detectorPositions(:,2), -res.detectorPositions(:,1), 'o', 'LineWidth', 1.4);
hold on;
plot(res.detectorPositions(res.centerDetectorIndex,2), ...
    -res.detectorPositions(res.centerDetectorIndex,1), ...
    'rp', 'MarkerSize', 12, 'LineWidth', 2);
axis equal; grid on;
xlabel('\Delta x [px]');
ylabel('\Delta y [px]');
title('Recovered detector geometry');

subplot(2,3,5);
fx = (-floor(W/2):ceil(W/2)-1) / W;
crow = floor(H/2) + 1;
plot(fx, fftshift(abs(res.otfWF(crow,:))), 'LineWidth', 1.2); hold on;
plot(fx, fftshift(abs(res.otfISMexp(crow,:))), 'LineWidth', 1.2);
plot(fx, fftshift(abs(res.otfISMideal(crow,:))), '--', 'LineWidth', 1.2);
xlabel('Spatial frequency [cycles/pixel]');
ylabel('|OTF|');
legend('WF OTF', 'ISM experimental', 'ISM ideal', 'Location', 'best');
grid on;
title('OTF cross-sections');

subplot(2,3,6);
montage(normalizeMontage(imgStack), ...
    'Size', [ceil(sqrt(size(imgStack,3))) ceil(sqrt(size(imgStack,3)))]);
title('Detector images');

% convergence figure
if ~isempty(res.convHistory.iter)
    figure('Color', 'w', 'Name', 'ACO convergence');
    semilogy(res.convHistory.iter, res.convHistory.relRMSE, '-o', 'LineWidth', 1.3);
    xlabel('Iteration');
    ylabel('Relative RMSE of autocorrelation');
    grid on;
    title('Schultz-Snyder convergence');
end

% show masked autocorrelation target
figure('Color', 'w', 'Name', 'ACO target');
subplot(1,2,1);
imagesc(fftshift(res.acoAverage)); axis image off; colormap hot;
title('Average autocorrelation');

subplot(1,2,2);
imagesc(fftshift(res.acoAverageMasked)); axis image off; colormap hot;
title('Masked average autocorrelation');
end

function out = normalizeMontage(stack)
[H,W,N] = size(stack);
out = zeros(H,W,1,N);
for k = 1:N
    I = double(stack(:,:,k));
    I = I - min(I(:));
    mx = max(I(:));
    if mx > 0
        I = I / mx;
    end
    out(:,:,1,k) = I;
end
end

