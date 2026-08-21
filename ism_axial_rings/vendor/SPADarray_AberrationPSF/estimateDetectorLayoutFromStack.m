function [detXY, diagnostics] = estimateDetectorLayoutFromStack(stack, pixelSizeUm, opts)
%ESTIMATEDETECTORLAYOUTFROMSTACK Estimate detector positions from channel images.
%
%   [detXY, diagnostics] = estimateDetectorLayoutFromStack(stack, pixelSizeUm)
%
%   stack can be [ny nx nDetector] or [ny nx nDetector nPlane]. For 4-D data
%   only opts.planeIndex is used (default 1, the focal plane). Shifts are
%   never estimated independently from a defocused plane.
%
%   The estimate uses FFT phase correlation between each detector-channel
%   image and a reference channel. The measured shifts are ISM reassignment
%   shifts: shifts to apply to channel images, not physical detector offsets.
%   For matched excitation/detection PSFs the physical detector offset is
%   approximately -2 times the reassignment shift. The conversion is exposed
%   through opts.positionSign and opts.detectorScale.
%
%   Useful opts fields:
%       refChannel   : temporary detector reference [brightest]
%       centerDetectorIndex : final center detector [estimated from shift cloud]
%       positionSign : +1 or -1 multiplier for output coordinates [-1]
%       detectorScale: reassignment-to-detector scale [2]
%       scanAxisSigns : [x y] signs mapping increasing image column/row to
%                       physical positive scan x/y [[1 1]]
%       centerMode   : 'reference', 'signalWeighted', or 'mean' ['reference']
%       planeIndex   : plane used when stack is 4-D [1]
%       airyUnitUm   : Airy-disk diameter in sample space [NaN]
%       smoothSigma  : Gaussian sigma before registration [1]
%       useWindow    : apply a periodic Hann window [true]
%       normalizeImages : normalize each image by its standard deviation [true]
%       upsampleReg  : direct-overlap subpixel refinement factor [20]
%       logCompress  : use log1p intensity before correlation [false]

    if nargin < 2 || isempty(pixelSizeUm)
        pixelSizeUm = 1;
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    planeIndex = getOptionLocal(opts, 'planeIndex', 1);
    imgStack = collapseStack(stack, planeIndex);
    [ny, nx, nCh] = size(imgStack);

    logCompress = getOptionLocal(opts, 'logCompress', false);
    smoothSigma = getOptionLocal(opts, 'smoothSigma', 1);
    useWindow = getOptionLocal(opts, 'useWindow', true);
    normalizeImages = getOptionLocal(opts, 'normalizeImages', true);
    upsampleReg = getOptionLocal(opts, 'upsampleReg', 20);
    imgs = zeros(ny, nx, nCh);
    signal = zeros(nCh, 1);

    for k = 1:nCh
        img = double(imgStack(:,:,k));
        signalImage = img - borderMedian(img);
        signalImage(signalImage < 0) = 0;
        signal(k) = sum(signalImage(:));
        if logCompress
            img = log1p(img);
        end
        imgs(:,:,k) = preprocessForRegistration( ...
            img, smoothSigma, useWindow, normalizeImages);
    end

    tempRefChannel = getOptionLocal(opts, 'refChannel', []);
    if isempty(tempRefChannel)
        [~, tempRefChannel] = max(signal);
    end
    if ~isnumeric(tempRefChannel) || ~isscalar(tempRefChannel) || ...
            ~isfinite(tempRefChannel) || tempRefChannel ~= round(tempRefChannel) || ...
            tempRefChannel < 1 || tempRefChannel > nCh
        error('estimateDetectorLayoutFromStack:BadReference', ...
            'refChannel must be between 1 and %d.', nCh);
    end

    shiftsToTempRefYX = zeros(nCh, 2);
    tempPeakValue = zeros(nCh, 1);
    tempRef = imgs(:,:,tempRefChannel);
    for k = 1:nCh
        [shiftsToTempRefYX(k,:), tempPeakValue(k)] = phaseCorrShiftToApply( ...
            tempRef, imgs(:,:,k), upsampleReg);
    end
    shiftsToTempRefYX(tempRefChannel,:) = [0 0];

    centerDetectorIndex = getOptionLocal(opts, 'centerDetectorIndex', []);
    if isempty(centerDetectorIndex)
        shiftCloudCentroidYX = mean(shiftsToTempRefYX, 1);
        distanceSquared = sum((shiftsToTempRefYX - shiftCloudCentroidYX).^2, 2);
        [~, centerDetectorIndex] = min(distanceSquared);
    elseif ~isnumeric(centerDetectorIndex) || ~isscalar(centerDetectorIndex) || ...
            ~isfinite(centerDetectorIndex) || ...
            centerDetectorIndex ~= round(centerDetectorIndex) || ...
            centerDetectorIndex < 1 || centerDetectorIndex > nCh
        error('estimateDetectorLayoutFromStack:BadCenterDetector', ...
            'centerDetectorIndex must be between 1 and %d.', nCh);
    else
        shiftCloudCentroidYX = mean(shiftsToTempRefYX, 1);
    end

    shiftsToCenterYX = zeros(nCh, 2);
    peakValue = zeros(nCh, 1);
    centerRef = imgs(:,:,centerDetectorIndex);
    for k = 1:nCh
        [shiftsToCenterYX(k,:), peakValue(k)] = phaseCorrShiftToApply( ...
            centerRef, imgs(:,:,k), upsampleReg);
    end
    shiftsToCenterYX(centerDetectorIndex,:) = [0 0];
    shiftsArrayPx = [shiftsToCenterYX(:,2), shiftsToCenterYX(:,1)];
    scanAxisSigns = double(getOptionLocal(opts, 'scanAxisSigns', [1 1]));
    scanAxisSigns = scanAxisSigns(:).';
    if numel(scanAxisSigns) ~= 2 || ...
            any(~ismember(scanAxisSigns,[-1 1]))
        error('estimateDetectorLayoutFromStack:BadScanAxisSigns', ...
            'scanAxisSigns must be [signX signY] with values +1 or -1.');
    end
    shiftsPx = shiftsArrayPx .* scanAxisSigns;

    positionSign = getOptionLocal(opts, 'positionSign', -1);
    detectorScale = getOptionLocal(opts, 'detectorScale', 2);
    if ~isnumeric(positionSign) || ~isscalar(positionSign) || ...
            ~isfinite(positionSign) || ~ismember(positionSign, [-1 1])
        error('estimateDetectorLayoutFromStack:BadPositionSign', ...
            'positionSign must be +1 or -1.');
    end
    if ~isnumeric(detectorScale) || ~isscalar(detectorScale) || ...
            ~isfinite(detectorScale) || detectorScale <= 0
        error('estimateDetectorLayoutFromStack:BadDetectorScale', ...
            'detectorScale must be a finite positive scalar.');
    end
    detXY = positionSign * detectorScale * shiftsPx * pixelSizeUm;

    centerMode = lower(char(getOptionLocal(opts, 'centerMode', 'reference')));
    switch centerMode
        case 'reference'
            detXY = bsxfun(@minus, detXY, detXY(centerDetectorIndex,:));
        case 'signalweighted'
            w = signal(:);
            if sum(w) > 0
                detXY = bsxfun(@minus, detXY, sum(bsxfun(@times, detXY, w), 1) / sum(w));
            end
        case 'mean'
            detXY = bsxfun(@minus, detXY, mean(detXY, 1));
        otherwise
            error('estimateDetectorLayoutFromStack:BadCenterMode', ...
                'centerMode must be ''reference'', ''signalWeighted'', or ''mean''.');
    end

    airyUnitUm = getOptionLocal(opts, 'airyUnitUm', NaN);
    diagnostics.planeIndex = planeIndex;
    diagnostics.tempRefChannel = tempRefChannel;
    diagnostics.refChannel = centerDetectorIndex;
    diagnostics.centerDetectorIndex = centerDetectorIndex;
    diagnostics.shiftCloudCentroidYX = shiftCloudCentroidYX;
    diagnostics.shiftsToTempRefYX = shiftsToTempRefYX;
    diagnostics.tempPeakValue = tempPeakValue;
    diagnostics.shiftsToCenterYX = shiftsToCenterYX;
    diagnostics.scanAxisSigns = scanAxisSigns;
    diagnostics.shiftsArrayPx = shiftsArrayPx;
    diagnostics.pixelSizeUm = pixelSizeUm;
    diagnostics.pixelSizeNm = 1000 * pixelSizeUm;
    diagnostics.airyUnitUm = airyUnitUm;
    diagnostics.airyUnitDefinition = 'Airy disk diameter: 1 AU = 1.22*lambda_em/NA';
    diagnostics.positionSign = positionSign;
    diagnostics.detectorScale = detectorScale;
    diagnostics.smoothSigma = smoothSigma;
    diagnostics.useWindow = logical(useWindow);
    diagnostics.normalizeImages = logical(normalizeImages);
    diagnostics.upsampleReg = upsampleReg;
    diagnostics.logCompress = logical(logCompress);
    diagnostics.centerMode = centerMode;
    diagnostics.signal = signal;
    diagnostics.shiftsPx = shiftsPx;
    diagnostics.imagePeakDisplacementPx = -shiftsPx;
    diagnostics.shiftsUm = shiftsPx * pixelSizeUm;
    diagnostics.shiftsNm = 1000 * diagnostics.shiftsUm;
    diagnostics.detectorOffsetFromShiftUm = positionSign * detectorScale * diagnostics.shiftsUm;
    diagnostics.coordinateConvention = [ ...
        'shiftsPx is the physical-scan-frame image shift to apply; ' ...
        'image peak displacement is -shiftsPx; model detXY is ' ...
        'positionSign*detectorScale*shiftsPx'];
    diagnostics.detectorOffsetFromShiftNm = 1000 * diagnostics.detectorOffsetFromShiftUm;
    diagnostics.detXYUm = detXY;
    diagnostics.detXYNm = 1000 * detXY;
    if isfinitePositiveScalarLocal(airyUnitUm)
        diagnostics.shiftsAU = diagnostics.shiftsUm / airyUnitUm;
        diagnostics.shiftMagnitudeAU = hypot( ...
            diagnostics.shiftsAU(:,1), diagnostics.shiftsAU(:,2));
        diagnostics.detectorOffsetFromShiftAU = ...
            diagnostics.detectorOffsetFromShiftUm / airyUnitUm;
        diagnostics.detXYAU = detXY / airyUnitUm;
        diagnostics.detectorRadiusAU = hypot( ...
            diagnostics.detXYAU(:,1), diagnostics.detXYAU(:,2));
    else
        diagnostics.shiftsAU = nan(size(shiftsPx));
        diagnostics.shiftMagnitudeAU = nan(nCh, 1);
        diagnostics.detectorOffsetFromShiftAU = nan(size(shiftsPx));
        diagnostics.detXYAU = nan(size(detXY));
        diagnostics.detectorRadiusAU = nan(nCh, 1);
    end
    diagnostics.peakValue = peakValue;
end

function imgStack = collapseStack(stack, planeIndex)
    if ndims(stack) == 3
        imgStack = stack;
    elseif ndims(stack) == 4
        if ~isnumeric(planeIndex) || ~isscalar(planeIndex) || ...
                ~isfinite(planeIndex) || planeIndex < 1 || ...
                planeIndex > size(stack,4) || planeIndex ~= round(planeIndex)
            error('estimateDetectorLayoutFromStack:BadPlaneIndex', ...
                'planeIndex must be an integer between 1 and %d.', size(stack,4));
        end
        imgStack = stack(:,:,:,planeIndex);
    else
        error('estimateDetectorLayoutFromStack:BadStackShape', ...
            'stack must be [ny nx nDetector] or [ny nx nDetector nPlane].');
    end
end

function b = borderMedian(img)
    border = [img(1,:), img(end,:), img(:,1).', img(:,end).'];
    b = median(border(:));
end

function out = preprocessForRegistration(img, smoothSigma, useWindow, normalizeImages)
    out = double(img);
    if smoothSigma > 0
        out = gaussianSmooth(out, smoothSigma);
    end
    out = out - median(out(:));
    if normalizeImages
        sd = std(out(:));
        if sd > 0
            out = out / sd;
        end
    end
    if useWindow
        wy = periodicHann(size(out, 1));
        wx = periodicHann(size(out, 2));
        out = out .* (wy * wx.');
    end
end

function out = gaussianSmooth(img, sigma)
    if exist('imgaussfilt', 'file') == 2
        out = imgaussfilt(img, sigma);
        return;
    end
    radius = max(1, ceil(3*sigma));
    x = -radius:radius;
    kernel = exp(-0.5 * (x/sigma).^2);
    kernel = kernel / sum(kernel);
    out = conv2(conv2(img, kernel, 'same'), kernel.', 'same');
end

function w = periodicHann(n)
    if n <= 1
        w = ones(n, 1);
    else
        index = (0:n-1).';
        w = 0.5 - 0.5*cos(2*pi*index/n);
    end
end

function [shiftToApply, peak] = phaseCorrShiftToApply(refImg, movImg, upsample)
    F1 = fft2(refImg);
    F2 = fft2(movImg);
    CPS = F1 .* conj(F2);
    CPS = CPS ./ max(abs(CPS), eps);

    cc = fftshift(real(ifft2(CPS)));
    [height, width] = size(cc);
    [peak, idx] = max(cc(:));
    [py, px] = ind2sub([height, width], idx);
    cy = floor(height/2) + 1;
    cx = floor(width/2) + 1;

    dy0 = py - cy;
    dx0 = px - cx;
    dySub = 0;
    dxSub = 0;
    if py > 1 && py < height
        dySub = parabolicOffset(cc(py-1,px), cc(py,px), cc(py+1,px));
    end
    if px > 1 && px < width
        dxSub = parabolicOffset(cc(py,px-1), cc(py,px), cc(py,px+1));
    end

    p0 = [dy0 + dySub, dx0 + dxSub];
    scorePositive = overlapScore(refImg, movImg, p0);
    scoreNegative = overlapScore(refImg, movImg, -p0);
    if scoreNegative > scorePositive
        p0 = -p0;
    end
    if upsample > 1
        objective = @(p) -overlapScore(refImg, movImg, p);
        searchOpts = optimset('Display', 'off', 'TolX', 1/upsample, 'TolFun', 1e-7);
        p0 = fminsearch(objective, p0, searchOpts);
    end
    shiftToApply = p0(:).';
end

function score = overlapScore(ref, img, shift)
    shifted = shiftImageFourier(img, shift(1), shift(2));
    score = sum(ref(:) .* shifted(:));
end

function shifted = shiftImageFourier(img, dy, dx)
    [ny, nx] = size(img);
    fy = ifftshift((-floor(ny/2):ceil(ny/2)-1) / ny);
    fx = ifftshift((-floor(nx/2):ceil(nx/2)-1) / nx);
    [kx, ky] = meshgrid(fx, fy);
    ramp = exp(-1i * 2*pi * (kx*dx + ky*dy));
    shifted = real(ifft2(fft2(img) .* ramp));
end

function off = parabolicOffset(left, center, right)
    denom = left - 2*center + right;
    if abs(denom) < eps
        off = 0;
    else
        off = 0.5 * (left - right) / denom;
        off = max(min(off, 1), -1);
    end
end

function tf = isfinitePositiveScalarLocal(value)
    tf = isnumeric(value) && isscalar(value) && isfinite(value) && value > 0;
end

function v = getOptionLocal(opts, name, defaultValue)
    if isfield(opts, name) && ~isempty(opts.(name))
        v = opts.(name);
    else
        v = defaultValue;
    end
end
