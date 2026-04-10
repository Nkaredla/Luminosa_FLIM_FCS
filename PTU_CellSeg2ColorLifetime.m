function res = PTU_CellSeg2ColorLifetime(name, cfg)
% PTU_CellSeg2ColorLifetime
% Segment sparse cells from detector 1 / pulse 1 in a cnum=2 PTU scan,
% extract ROI TCSPC curves for detector 1 / pulse 1 and detector 2 / pulse 2,
% fit the detector 2 / pulse 2 curves for each ROI plus background, display
% an ROI-average lifetime image, and save ROI/background TCSPC outputs.
%
% Example
% -------
% cfg = struct();
% cfg.tau0 = [0.4 2.5];
% cfg.segDetectorIndex = 1;   % d1
% cfg.fitDetectorIndex = 2;   % d2
% cfg.segPulseIndex = 1;      % channel 1
% cfg.fitPulseIndex = 2;      % channel 2
% res = PTU_CellSeg2ColorLifetime('D:\data\RawImage.ptu', cfg);

    if nargin < 2
        cfg = struct();
    end
    cfg = setDefaultCfg(cfg);

    if ~(ischar(name) || isstring(name))
        error('Input "name" must be a PTU filename.');
    end
    name = char(name);
    if numel(name) < 4 || ~strcmpi(name(end-2:end), 'ptu')
        error('Input file must be a .ptu file.');
    end
    if cfg.cnum ~= 2
        error('This function is written for cnum = 2 measurements.');
    end

    logVerbose(cfg, 'PTU_CellSeg2ColorLifetime: %s', name);
    progressState = initProgressState(cfg, name);
    progressCleanup = onCleanup(@() closeProgressState(progressState, cfg)); %#ok<NASGU>

    head = PTU_Read_Head(name);
    if isempty(head)
        error('Could not read PTU header.');
    end
    pixelSizeUm = resolvePixelSizeUmFromHead(head);
    cfg = attachPhysicalSegmentationScale(cfg, pixelSizeUm);
    if isfinite(cfg.pixelSizeUm) && isfield(cfg, 'minPixelSizeUmToProcess') && ...
            isfinite(cfg.minPixelSizeUmToProcess) && cfg.minPixelSizeUmToProcess > 0 && ...
            cfg.pixelSizeUm < cfg.minPixelSizeUmToProcess
        error('Pixel size %.4f um is below the minimum allowed threshold %.4f um.', ...
            cfg.pixelSizeUm, cfg.minPixelSizeUmToProcess);
    end
    if isfinite(cfg.pixelSizeUm) && isfield(cfg, 'maxPixelSizeUmToProcess') && ...
            isfinite(cfg.maxPixelSizeUmToProcess) && cfg.maxPixelSizeUmToProcess > 0 && ...
            cfg.pixelSizeUm > cfg.maxPixelSizeUmToProcess
        error('Pixel size %.4f um is above the maximum allowed threshold %.4f um.', ...
            cfg.pixelSizeUm, cfg.maxPixelSizeUmToProcess);
    end
    updateProgressState(progressState, 'Header read', struct('head', head));
    logVerbose(cfg, 'Header read.');
    if isfinite(cfg.pixelSizeUm)
        logVerbose(cfg, 'Pixel size = %.4f um. Minimum kept cell diameter = %.2f um (%.2f px).', ...
            cfg.pixelSizeUm, cfg.minCellDiameterUm, cfg.minCellDiameterPixels);
    else
        logVerbose(cfg, 'Pixel size not found in header. Using pixel-only area threshold.');
    end

    rawDtNs = double(head.MeasDesc_Resolution) * 1e9;
    nNativeBins = ceil(double(head.MeasDesc_GlobalResolution) / max(double(head.MeasDesc_Resolution), eps)) + 1;

    readOpts = struct();
    readOpts.photonsPerChunk = cfg.photonsPerChunk;
    readOpts.computePerFrame = false;
    readOpts.storeTcspcPix = false;
    readOpts.storePhotonLists = true;
    readOpts.storeTimeCell = false;
    readOpts.showWaitbar = cfg.showWaitbar;
    readOpts.minLifetimeBin_ns = rawDtNs;
    readOpts.maxNgate = nNativeBins;

    logVerbose(cfg, 'Reading PTU photon lists...');
    out = PTU_FLIM_GPU(name, readOpts);
    if isempty(out.im_tcspc)
        error('Photon lists were not stored. PTU read failed to keep per-photon TCSPC data.');
    end
    updateProgressState(progressState, 'Photon lists loaded', struct('out', out));
    logVerbose(cfg, 'Photon list read complete. %d photons kept.', numel(out.im_tcspc));
    logVerbose(cfg, 'Photon vector lengths: tcspc=%d, chan=%d, col=%d, line=%d', ...
        numel(out.im_tcspc), numel(out.im_chan), numel(out.im_col), numel(out.im_line));

    detectorIDs = double(out.dind(:));
    if numel(detectorIDs) < max(cfg.segDetectorIndex, cfg.fitDetectorIndex)
        error('Requested detector indices exceed the number of active detectors in the PTU file.');
    end

    tcspcByDetector = buildDetectorHistogram(out);
    [gateStarts, gateStops, gateLen, gateInfo] = detectPulseGatesFromHistogram(tcspcByDetector, cfg.cnum, cfg);
    if numel(gateStarts) < cfg.cnum
        error('Could not detect %d pulse windows reliably.', cfg.cnum);
    end
    updateProgressState(progressState, 'Pulse gates detected', struct( ...
        'tcspcByDetector', tcspcByDetector, 'gateStarts', gateStarts, 'gateStops', gateStops));
    logVerbose(cfg, 'Pulse gates detected: [%s] bins, gate length %d bins.', num2str(gateStarts(:).'), gateLen);

    segDetectorID = detectorIDs(cfg.segDetectorIndex);
    fitDetectorID = detectorIDs(cfg.fitDetectorIndex);
    segGateStart = gateStarts(cfg.segPulseIndex);
    fitGateStart = gateStarts(cfg.fitPulseIndex);

    nx = double(out.head.ImgHdr_PixX);
    ny = double(out.head.ImgHdr_PixY);

    gatedCh1D1 = selectGatedPhotons(out, segDetectorID, segGateStart, gateLen, nx, ny);
    gatedCh2D2 = selectGatedPhotons(out, fitDetectorID, fitGateStart, gateLen, nx, ny);

    segImage = countsImageFromGated(gatedCh1D1, nx, ny);
    fitImage = countsImageFromGated(gatedCh2D2, nx, ny);
    updateProgressState(progressState, 'Gated images built', struct('segImage', segImage, 'fitImage', fitImage));
    logVerbose(cfg, 'Built gated images for d%d/pulse%d and d%d/pulse%d.', ...
        cfg.segDetectorIndex, cfg.segPulseIndex, cfg.fitDetectorIndex, cfg.fitPulseIndex);

    [labelImage, segMask, roiStats, segDebug] = segmentSparseCells(segImage, cfg);
    nCells = max(labelImage(:));
    updateProgressState(progressState, 'Cells segmented', struct('segImage', segImage, 'labelImage', labelImage, 'segMask', segMask));
    logVerbose(cfg, 'Segmentation found %d cells.', nCells);

    [roiCurvesCh1D1, bgCurveCh1D1, globalCurveCh1D1, roiPhotonCountsCh1D1, bgPhotonCountCh1D1] = curvesFromLabels(gatedCh1D1, labelImage, gateLen);
    [roiCurvesCh2D2, bgCurveCh2D2, globalCurveCh2D2, roiPhotonCountsCh2D2, bgPhotonCountCh2D2] = curvesFromLabels(gatedCh2D2, labelImage, gateLen);

    if sum(globalCurveCh2D2) < cfg.minPhotonsForIRF
        error('Not enough detector 2 / pulse 2 photons to estimate a global IRF.');
    end

    fitHead = struct('MeasDesc_Resolution', rawDtNs * 1e-9);
    irfOpts = struct('useGPU', logical(cfg.useGPU), 'nCasc', cfg.irfNCasc, 'nSub', cfg.irfNSub);
    logVerbose(cfg, 'Estimating global IRF from ch2/d2.');
    irfOut = Calc_mIRF_Global_GammaShifted_fast(fitHead, globalCurveCh2D2(:), cfg.tau0(:), irfOpts);
    updateProgressState(progressState, 'IRF estimated', struct('tNs', ((0:gateLen-1)' + 0.5) * rawDtNs, ...
        'globalCurveCh2D2', globalCurveCh2D2, 'irfOut', irfOut));

    tAxisNs = ((0:gateLen-1)' + 0.5) * rawDtNs;

    backgroundFit = fitCurveWithFixedIRF(bgCurveCh2D2, rawDtNs, irfOut.IRF(:), cfg.tau0, cfg);
    logVerbose(cfg, 'Background fit complete. tau_int = %.3f ns', backgroundFit.tauIntNs);

    roiFits = repmat(emptyFitResult(gateLen), nCells, 1);
    for k = 1:nCells
        logVerbose(cfg, 'Fitting ROI %d / %d', k, nCells);
        roiFits(k) = fitCurveWithFixedIRF(roiCurvesCh2D2(:,k), rawDtNs, irfOut.IRF(:), cfg.tau0, cfg);
        updateProgressState(progressState, sprintf('Fitting ROI %d/%d', k, nCells), struct( ...
            'tNs', tAxisNs, ...
            'currentCurve', roiCurvesCh2D2(:,k), ...
            'currentFit', roiFits(k).fitCounts, ...
            'labelImage', labelImage, ...
            'segImage', segImage, ...
            'fitImage', fitImage));
    end

    roiTauInt = nan(nCells, 1);
    for k = 1:nCells
        roiTauInt(k) = roiFits(k).tauIntNs;
    end

    componentAmpMaps = [];
    componentFracMaps = [];
    pixelFitInfo = struct();
    if cfg.usePixelwiseIntensityLifetimeMap
        if numel(cfg.tau0) ~= 2
            error('usePixelwiseIntensityLifetimeMap currently requires exactly two lifetime components in cfg.tau0.');
        end
        logVerbose(cfg, 'Building pixelwise component amplitudes inside ROIs.');
        cubeCh2D2 = tcspcCubeFromGated(gatedCh2D2, nx, ny, gateLen);
        [tauImage, componentAmpMaps, componentFracMaps, pixelFitInfo] = buildPixelwiseRoiLifetimeMap( ...
            cubeCh2D2, labelImage, roiFits, backgroundFit, rawDtNs, cfg, progressState, tAxisNs);
        updateProgressState(progressState, 'Pixelwise ROI lifetime map built', struct( ...
            'res', struct('images', struct('averageLifetimeImageNs', tauImage))));
    else
        tauImage = nan(nx, ny);
        validBg = fitImage > 0;
        if isfinite(backgroundFit.tauIntNs)
            tauImage(validBg) = backgroundFit.tauIntNs;
        end
        for k = 1:nCells
            if isfinite(roiTauInt(k))
                tauImage(labelImage == k) = roiTauInt(k);
            end
        end
        updateProgressState(progressState, 'ROI lifetime map built', struct( ...
            'res', struct('images', struct('averageLifetimeImageNs', tauImage))));
    end

    rois = repmat(struct( ...
        'index', [], ...
        'areaPixels', [], ...
        'centroidXY', [], ...
        'meanSegIntensity', [], ...
        'maxSegIntensity', [], ...
        'curveCh1D1', [], ...
        'curveCh2D2', [], ...
        'photonsCh1D1', [], ...
        'photonsCh2D2', [], ...
        'componentCoeffCh2D2', [], ...
        'componentFracCh2D2', [], ...
        'fitCh2D2', []), nCells, 1);
    for k = 1:nCells
        rois(k).index = k;
        rois(k).areaPixels = roiStats(k).Area;
        rois(k).centroidXY = roiStats(k).Centroid;
        rois(k).meanSegIntensity = roiStats(k).MeanIntensity;
        rois(k).maxSegIntensity = roiStats(k).MaxIntensity;
        rois(k).curveCh1D1 = roiCurvesCh1D1(:,k);
        rois(k).curveCh2D2 = roiCurvesCh2D2(:,k);
        rois(k).photonsCh1D1 = roiPhotonCountsCh1D1(k);
        rois(k).photonsCh2D2 = roiPhotonCountsCh2D2(k);
        rois(k).componentCoeffCh2D2 = roiFits(k).speciesCoeff;
        rois(k).componentFracCh2D2 = roiFits(k).speciesFrac;
        rois(k).fitCh2D2 = roiFits(k);
    end

    res = struct();
    res.name = name;
    res.config = cfg;
    res.head = head;
    res.pixelSizeUm = cfg.pixelSizeUm;
    res.rawResolutionNs = rawDtNs;
    res.detectorIDs = detectorIDs;
    res.segDetectorIndex = cfg.segDetectorIndex;
    res.fitDetectorIndex = cfg.fitDetectorIndex;
    res.segDetectorID = segDetectorID;
    res.fitDetectorID = fitDetectorID;
    res.gates = struct( ...
        'startsBin', gateStarts(:), ...
        'stopsBin', gateStops(:), ...
        'lengthBins', gateLen, ...
        'startsNs', (double(gateStarts(:)) - 1) * rawDtNs, ...
        'stopsNs', double(gateStops(:)) * rawDtNs, ...
        'info', gateInfo);
    res.images = struct( ...
        'segmentationSourceCh1D1', segImage, ...
        'analysisSourceCh2D2', fitImage, ...
        'segmentationMask', segMask, ...
        'labelImage', labelImage, ...
        'averageLifetimeImageNs', tauImage, ...
        'intensityAveragedLifetimeImageNs', tauImage, ...
        'componentAmplitudeMapsCh2D2', componentAmpMaps, ...
        'componentFractionMapsCh2D2', componentFracMaps, ...
        'segmentationDebug', segDebug);
    res.tcspc = struct( ...
        'tNs', tAxisNs, ...
        'globalCh1D1', globalCurveCh1D1, ...
        'globalCh2D2', globalCurveCh2D2, ...
        'backgroundCh1D1', bgCurveCh1D1, ...
        'backgroundCh2D2', bgCurveCh2D2, ...
        'roiCh1D1', roiCurvesCh1D1, ...
        'roiCh2D2', roiCurvesCh2D2);
    res.irf = irfOut;
    res.background = struct( ...
        'curveCh1D1', bgCurveCh1D1, ...
        'photonCountCh1D1', bgPhotonCountCh1D1, ...
        'curveCh2D2', bgCurveCh2D2, ...
        'photonCountCh2D2', bgPhotonCountCh2D2, ...
        'fitCh2D2', backgroundFit);
    res.pixelFit = pixelFitInfo;
    res.rois = rois;

    figHandles = gobjects(0);
    if cfg.showFigures || cfg.saveFigures
        visibleState = ternaryString(cfg.showFigures, 'on', 'off');
        figHandles = showSummaryFigures(res, visibleState);
    end

    if cfg.saveOutputs
        logVerbose(cfg, 'Saving outputs...');
        saved = saveResultOutputs(res, name, cfg, figHandles);
        res.savedOutputs = saved;
    else
        res.savedOutputs = struct();
    end

    updateProgressState(progressState, 'Done', struct('res', res));
    logVerbose(cfg, 'Done.');

    if ~isempty(figHandles) && cfg.saveFigures && cfg.closeFiguresAfterSave
        deleteValidFigures(figHandles);
    end
end

function cfg = setDefaultCfg(cfg)
    defaults.cnum = 2;
    defaults.segDetectorIndex = 1;
    defaults.fitDetectorIndex = 2;
    defaults.segPulseIndex = 1;
    defaults.fitPulseIndex = 2;
    defaults.photonsPerChunk = 5e6;
    defaults.showWaitbar = false;
    defaults.useGPU = false;

    defaults.segSmoothSigma = 1.2;
    defaults.segSmoothSigmaUm = 0.25;
    defaults.segTophatRadius = 8;
    defaults.segTophatRadiusUm = 1.0;
    defaults.segmentationMode = 'variance';
    defaults.segVarianceWindow = 3;
    defaults.segVarianceThresholdScale = 1.35;
    defaults.segVarianceBackgroundQuantile = 0.40;
    defaults.segRemoveBorderObjects = false;
    defaults.segBorderMarginPixels = 1;
    defaults.segFlatFilterSize = 9;
    defaults.segNormalizedClusterPeak = 250;
    defaults.segTophatRadiusPixelsMax = 12;
    defaults.segThresholdScale = 0.9;
    defaults.segThresholdMAD = 3.0;
    defaults.segUseAdaptiveThreshold = true;
    defaults.segAdaptiveSensitivity = 0.52;
    defaults.segAdaptiveNeighborhoodUm = 3.0;
    defaults.minCellArea = 12;
    defaults.maxCellArea = inf;
    defaults.minCellDiameterUm = 3.0;
    defaults.minPixelSizeUmToProcess = NaN;
    defaults.maxPixelSizeUmToProcess = NaN;
    defaults.segCloseRadius = 3;
    defaults.segCloseRadiusUm = 0.8;
    defaults.segCloseRadiusPixelsMax = 5;
    defaults.segPrefilterMinArea = 30;
    defaults.fillHoles = true;
    defaults.fillWholeCellArea = true;
    defaults.useWatershed = true;
    defaults.watershedH = 0.5;

    defaults.gateThresholdFrac = 0.15;
    defaults.gatePreBins = 100;
    defaults.minGateSeparationBins = 50;

    defaults.tau0 = [0.4 2.5];
    defaults.useGuiStyleMultiExpFit = false;
    defaults.usePixelwiseIntensityLifetimeMap = false;
    defaults.includeBackgroundInFit = true;
    defaults.pixelwiseIncludeBackground = true;
    defaults.pixelwiseMinPhotons = 1;
    defaults.minPhotonsForIRF = 500;
    defaults.minPhotonsPerFit = 100;
    defaults.optimizeTau = true;
    defaults.irfNCasc = 4;
    defaults.irfNSub = 6;
    defaults.fitMaxIter = 400;

    defaults.showFigures = true;
    defaults.saveOutputs = true;
    defaults.saveMat = true;
    defaults.saveCsv = true;
    defaults.saveFigures = true;
    defaults.closeFiguresAfterSave = true;
    defaults.outputDir = '';
    defaults.outputStem = '';
    defaults.verbose = false;
    defaults.plotProgress = false;
    defaults.closeProgressFigure = true;
    defaults.pixelSizeUm = NaN;
    defaults.minCellDiameterPixels = NaN;
    defaults.minCellAreaPixelsEffective = defaults.minCellArea;
    defaults.segSmoothSigmaPixelsEffective = defaults.segSmoothSigma;
    defaults.segTophatRadiusPixelsEffective = defaults.segTophatRadius;
    defaults.segCloseRadiusPixelsEffective = defaults.segCloseRadius;
    defaults.segAdaptiveNeighborhoodPixelsEffective = 31;

    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(cfg, names{k}) || isempty(cfg.(names{k}))
            cfg.(names{k}) = defaults.(names{k});
        end
    end
    cfg.tau0 = double(cfg.tau0(:));
end

function pixelSizeUm = resolvePixelSizeUmFromHead(head)
    pixelSizeUm = NaN;
    if ~isstruct(head) || ~isfield(head, 'ImgHdr_PixResol') || isempty(head.ImgHdr_PixResol)
        return;
    end
    pixelSizeUm = double(head.ImgHdr_PixResol);
    if ~isfinite(pixelSizeUm) || pixelSizeUm <= 0
        pixelSizeUm = NaN;
        return;
    end
    if pixelSizeUm < 1e-3
        pixelSizeUm = pixelSizeUm * 1e6;
    end
end

function cfg = attachPhysicalSegmentationScale(cfg, pixelSizeUm)
    cfg.pixelSizeUm = pixelSizeUm;
    cfg.minCellDiameterPixels = NaN;
    cfg.minCellAreaPixelsEffective = max(1, double(cfg.minCellArea));
    cfg.segSmoothSigmaPixelsEffective = double(cfg.segSmoothSigma);
    cfg.segTophatRadiusPixelsEffective = max(1, round(double(cfg.segTophatRadius)));
    cfg.segCloseRadiusPixelsEffective = max(1, round(double(cfg.segCloseRadius)));
    cfg.segAdaptiveNeighborhoodPixelsEffective = 31;

    if isfinite(pixelSizeUm) && pixelSizeUm > 0 && isfield(cfg, 'minCellDiameterUm') && ...
            ~isempty(cfg.minCellDiameterUm) && isfinite(cfg.minCellDiameterUm) && cfg.minCellDiameterUm > 0
        cfg.minCellDiameterPixels = double(cfg.minCellDiameterUm) / pixelSizeUm;
        cfg.minCellAreaPixelsEffective = max(cfg.minCellAreaPixelsEffective, ...
            pi * (0.5 * cfg.minCellDiameterPixels)^2);
    end

    if isfinite(pixelSizeUm) && pixelSizeUm > 0
        cfg.segSmoothSigmaPixelsEffective = resolvePhysicalScalePx(cfg.segSmoothSigma, cfg.segSmoothSigmaUm, pixelSizeUm, false);
        cfg.segTophatRadiusPixelsEffective = resolvePhysicalScalePx(cfg.segTophatRadius, cfg.segTophatRadiusUm, pixelSizeUm, true);
        cfg.segCloseRadiusPixelsEffective = resolvePhysicalScalePx(cfg.segCloseRadius, cfg.segCloseRadiusUm, pixelSizeUm, true);
        cfg.segAdaptiveNeighborhoodPixelsEffective = resolveAdaptiveNeighborhoodPx(cfg.segAdaptiveNeighborhoodUm, pixelSizeUm);
    end

    if isfield(cfg, 'segTophatRadiusPixelsMax') && isfinite(cfg.segTophatRadiusPixelsMax) && cfg.segTophatRadiusPixelsMax > 0
        cfg.segTophatRadiusPixelsEffective = min(cfg.segTophatRadiusPixelsEffective, double(cfg.segTophatRadiusPixelsMax));
    end
    if isfield(cfg, 'segCloseRadiusPixelsMax') && isfinite(cfg.segCloseRadiusPixelsMax) && cfg.segCloseRadiusPixelsMax > 0
        cfg.segCloseRadiusPixelsEffective = min(cfg.segCloseRadiusPixelsEffective, double(cfg.segCloseRadiusPixelsMax));
    end
end

function valPx = resolvePhysicalScalePx(defaultPx, scaleUm, pixelSizeUm, roundToInt)
    valPx = double(defaultPx);
    if nargin < 4
        roundToInt = false;
    end
    if isfinite(scaleUm) && scaleUm > 0 && isfinite(pixelSizeUm) && pixelSizeUm > 0
        valPx = scaleUm / pixelSizeUm;
    end
    if roundToInt
        valPx = max(1, round(valPx));
    else
        valPx = max(eps, valPx);
    end
end

function winPx = resolveAdaptiveNeighborhoodPx(neighborhoodUm, pixelSizeUm)
    winPx = 31;
    if isfinite(neighborhoodUm) && neighborhoodUm > 0 && isfinite(pixelSizeUm) && pixelSizeUm > 0
        winPx = max(9, round(neighborhoodUm / pixelSizeUm));
    end
    if mod(winPx, 2) == 0
        winPx = winPx + 1;
    end
end

function tcspcByDetector = buildDetectorHistogram(out)
    detectorIDs = double(out.dind(:));
    nDet = numel(detectorIDs);
    nBins = double(out.Ngate);

    lut = zeros(1, max(256, max(detectorIDs) + 1));
    lut(detectorIDs + 1) = 1:nDet;

    [t, chan] = getAlignedPhotonArrays(out, false);
    if isempty(t)
        tcspcByDetector = zeros(nBins, nDet);
        return;
    end
    det = zeros(size(chan));
    chanValid = chan >= 0 & chan <= (numel(lut) - 1);
    det(chanValid) = lut(chan(chanValid) + 1);
    valid = isfinite(t) & isfinite(det) & det >= 1 & det <= nDet & t >= 1 & t <= nBins;

    if any(valid)
        tcspcByDetector = accumarray([round(t(valid)), round(det(valid))], 1, [nBins, nDet], @sum, 0);
    else
        tcspcByDetector = zeros(nBins, nDet);
    end
end

function [gateStarts, gateStops, gateLen, info] = detectPulseGatesFromHistogram(tcspcByDetector, cnum, cfg)
    profile = mean(double(tcspcByDetector), 2);
    profile = max(profile(:), 0);

    if numel(profile) < 2 || all(profile <= 0)
        error('No usable TCSPC signal was found for pulse-gate detection.');
    end

    smoothBins = max(3, min(21, 2 * floor(numel(profile) / 200) + 1));
    kernel = ones(smoothBins, 1) / smoothBins;
    profileSmooth = conv(profile, kernel, 'same');

    baseline = meanOfLowestFraction(profileSmooth, 0.2);
    peakVal = max(profileSmooth);
    thr = baseline + cfg.gateThresholdFrac * max(peakVal - baseline, 0);
    mask = profileSmooth > thr;

    [startsRaw, stopsRaw, scoresRaw] = contiguousRuns(mask, profileSmooth);

    if numel(startsRaw) < cnum
        [startsRaw, stopsRaw, scoresRaw] = fallbackPeakWindows(profileSmooth, cnum, cfg.minGateSeparationBins);
    end
    if numel(startsRaw) < cnum
        error('Could not identify %d pulse windows from the global TCSPC profile.', cnum);
    end

    [~, orderByScore] = sort(scoresRaw(:), 'descend');
    keep = sort(orderByScore(1:cnum));
    startsRaw = startsRaw(keep);
    stopsRaw = stopsRaw(keep);
    scoresRaw = scoresRaw(keep);

    [startsRaw, orderByStart] = sort(startsRaw(:));
    stopsRaw = stopsRaw(orderByStart);
    scoresRaw = scoresRaw(orderByStart);

    gateStarts = max(1, startsRaw - cfg.gatePreBins);
    nextStart = [gateStarts(2:end) - 1; numel(profileSmooth)];
    gateLen = min(nextStart - gateStarts + 1);
    gateLen = max(1, gateLen);
    gateStops = gateStarts + gateLen - 1;

    info = struct();
    info.profile = profile;
    info.profileSmooth = profileSmooth;
    info.threshold = thr;
    info.rawStarts = startsRaw;
    info.rawStops = stopsRaw;
    info.rawScores = scoresRaw;
end

function m = meanOfLowestFraction(x, frac)
    x = sort(double(x(:)));
    n = max(1, round(frac * numel(x)));
    m = mean(x(1:n));
end

function [starts, stops, scores] = contiguousRuns(mask, scoreTrace)
    mask = logical(mask(:));
    d = diff([false; mask; false]);
    starts = find(d == 1);
    stops = find(d == -1) - 1;
    scores = zeros(numel(starts), 1);
    for k = 1:numel(starts)
        scores(k) = sum(scoreTrace(starts(k):stops(k)));
    end
end

function [starts, stops, scores] = fallbackPeakWindows(profile, cnum, minSep)
    y = double(profile(:));
    n = numel(y);
    starts = zeros(0,1);
    stops = zeros(0,1);
    scores = zeros(0,1);

    isPeak = false(size(y));
    if n >= 3
        isPeak(2:end-1) = y(2:end-1) >= y(1:end-2) & y(2:end-1) >= y(3:end);
    end
    if n >= 1
        isPeak(1) = y(1) >= y(min(2,n));
        isPeak(end) = y(end) >= y(max(1,n-1));
    end

    peakIdx = find(isPeak & y > 0);
    if isempty(peakIdx)
        [~, peakIdx] = max(y);
    end

    [~, ord] = sort(y(peakIdx), 'descend');
    peakIdx = peakIdx(ord);

    picked = zeros(0,1);
    for k = 1:numel(peakIdx)
        if isempty(picked) || all(abs(peakIdx(k) - picked) >= minSep)
            picked(end+1,1) = peakIdx(k); %#ok<AGROW>
        end
        if numel(picked) >= cnum
            break;
        end
    end
    if numel(picked) < cnum
        return;
    end

    picked = sort(picked);
    halfWidth = max(5, floor(minSep / 2));
    starts = max(1, picked - halfWidth);
    stops = min(n, picked + halfWidth);
    scores = y(picked);
end

function gated = selectGatedPhotons(out, detectorID, gateStart, gateLen, nx, ny)
    [t, det, x, y] = getAlignedPhotonArrays(out, true);

    gateStop = gateStart + gateLen - 1;
    keep = det == double(detectorID) & ...
        t >= double(gateStart) & t <= double(gateStop) & ...
        x >= 1 & x <= double(nx) & ...
        y >= 1 & y <= double(ny);

    gated = struct();
    gated.x = x(keep);
    gated.y = y(keep);
    gated.t = t(keep) - double(gateStart) + 1;
    gated.pix = gated.x + (gated.y - 1) * double(nx);
end

function [t, chan, x, y, nAligned] = getAlignedPhotonArrays(out, needXY)
    if nargin < 2
        needXY = true;
    end

    t = double(out.im_tcspc(:));
    chan = double(out.im_chan(:));
    if needXY
        x = double(out.im_col(:));
        y = double(out.im_line(:));
        nAligned = min([numel(t), numel(chan), numel(x), numel(y)]);
        x = x(1:nAligned);
        y = y(1:nAligned);
    else
        x = [];
        y = [];
        nAligned = min([numel(t), numel(chan)]);
    end

    if nAligned <= 0
        t = zeros(0,1);
        chan = zeros(0,1);
        x = zeros(0,1);
        y = zeros(0,1);
        nAligned = 0;
        return;
    end

    t = t(1:nAligned);
    chan = chan(1:nAligned);
end

function img = countsImageFromGated(gated, nx, ny)
    if isempty(gated.pix)
        img = zeros(nx, ny);
        return;
    end
    cnt = accumarray(round(gated.pix(:)), 1, [nx * ny, 1], @sum, 0);
    img = reshape(cnt, [nx, ny]);
end

function [labelImage, segMask, roiStats, segDebug] = segmentSparseCells(segImage, cfg)
    im = double(segImage);
    im(isnan(im)) = 0;
    im = max(im, 0);

    if ~any(im(:) > 0)
        segMask = false(size(im));
        labelImage = zeros(size(im));
        roiStats = repmat(struct('Area', [], 'Centroid', [], 'MeanIntensity', [], 'MaxIntensity', []), 0, 1);
        segDebug = struct( ...
            'sourceImage', im, ...
            'varianceImage', zeros(size(im)), ...
            'varianceBackground', NaN, ...
            'varianceThreshold', NaN, ...
            'seedMask', false(size(im)), ...
            'filledMask', false(size(im)), ...
            'blurredMaskedImage', zeros(size(im)), ...
            'normalizedClusterImage', zeros(size(im)), ...
            'clusterPeakBlur', zeros(0, 1));
        return;
    end

    [seedMask, varianceDebug] = buildVarianceSegmentationMask(im, cfg);
    minSeedAreaPx = max(1, round(double(cfg.segPrefilterMinArea)));
    seedMask = keepConnectedComponentsMinAreaLocal(seedMask, minSeedAreaPx);

    if cfg.fillHoles && exist('imfill', 'file') == 2
        filledMask = imfill(seedMask, 'holes');
    else
        filledMask = seedMask;
    end

    if ~cfg.fillWholeCellArea
        segMask = seedMask;
    else
        segMask = filledMask;
    end

    segMask = keepConnectedComponentsMinAreaLocal(segMask, minSeedAreaPx);

    labelImage = bwlabel(segMask);
    segMask = applyObjectSizeLimits(segMask, labelImage, cfg);
    labelImage = bwlabel(segMask);

    [blurredMaskedImage, normalizedClusterImage, clusterPeakBlur] = ...
        buildNormalizedClusterImageLocal(im, labelImage, cfg);

    if max(labelImage(:)) > 0
        roiStats = regionprops(labelImage, segImage, 'Area', 'Centroid', 'MeanIntensity', 'MaxIntensity');
    else
        roiStats = repmat(struct('Area', [], 'Centroid', [], 'MeanIntensity', [], 'MaxIntensity', []), 0, 1);
    end

    segDebug = struct( ...
        'sourceImage', im, ...
        'varianceImage', varianceDebug.varianceImage, ...
        'varianceBackground', varianceDebug.backgroundVariance, ...
        'varianceThreshold', varianceDebug.threshold, ...
        'seedMask', seedMask, ...
        'filledMask', filledMask, ...
        'blurredMaskedImage', blurredMaskedImage, ...
        'normalizedClusterImage', normalizedClusterImage, ...
        'clusterPeakBlur', clusterPeakBlur);
end

function segMask = buildAdaptiveBrightMask(imSmooth, cfg)
    segMask = false(size(imSmooth));
    imSmooth = double(imSmooth);
    if ~any(imSmooth(:) > 0)
        return;
    end

    imNorm = imSmooth - min(imSmooth(:));
    maxVal = max(imNorm(:));
    if maxVal <= 0
        return;
    end
    imNorm = imNorm ./ maxVal;

    if exist('adaptthresh', 'file') == 2 && exist('imbinarize', 'file') == 2
        try
            winPx = cfg.segAdaptiveNeighborhoodPixelsEffective;
            T = adaptthresh(imNorm, cfg.segAdaptiveSensitivity, 'ForegroundPolarity', 'bright', ...
                'NeighborhoodSize', [winPx winPx], 'Statistic', 'gaussian');
            segMask = imbinarize(imNorm, T);
            return;
        catch
        end
    end

    if exist('imbinarize', 'file') == 2
        try
            segMask = imbinarize(imNorm, localGraythresh(imNorm) * max(0.5, cfg.segAdaptiveSensitivity));
        catch
            segMask = imNorm >= (mean(imNorm(:)) + 0.5 * std(imNorm(:)));
        end
    else
        segMask = imNorm >= (mean(imNorm(:)) + 0.5 * std(imNorm(:)));
    end
end

function [segMask, debug] = buildVarianceSegmentationMask(im, cfg)
    segMask = false(size(im));
    debug = struct('varianceImage', zeros(size(im)), 'backgroundVariance', NaN, 'threshold', NaN);

    if ~any(im(:) > 0)
        return;
    end

    winPx = 3;
    if isfield(cfg, 'segVarianceWindow') && isfinite(cfg.segVarianceWindow) && cfg.segVarianceWindow > 0
        winPx = max(1, round(double(cfg.segVarianceWindow)));
    end
    if mod(winPx, 2) == 0
        winPx = winPx + 1;
    end

    if exist('stdfilt', 'file') == 2
        varianceImage = stdfilt(double(im), ones(winPx)).^2;
    else
        varianceImage = localVarianceImageFallback(double(im), winPx);
    end
    varianceImage(~isfinite(varianceImage)) = 0;
    debug.varianceImage = varianceImage;

    vals = sort(varianceImage(:));
    if isempty(vals)
        return;
    end
    q = 0.60;
    if isfield(cfg, 'segVarianceBackgroundQuantile') && isfinite(cfg.segVarianceBackgroundQuantile)
        q = min(max(double(cfg.segVarianceBackgroundQuantile), 0.01), 0.99);
    end

    idx = min(numel(vals), max(1, ceil(numel(vals) * q)));
    backgroundVariance = vals(idx);
    thresholdScale = 2;
    if isfield(cfg, 'segVarianceThresholdScale') && isfinite(cfg.segVarianceThresholdScale)
        thresholdScale = max(double(cfg.segVarianceThresholdScale), 0);
    end

    thr = thresholdScale * backgroundVariance;

    segMask = varianceImage > thr;
    debug.backgroundVariance = backgroundVariance;
    debug.threshold = thr;
end

function segMask = keepConnectedComponentsMinAreaLocal(segMask, minAreaPx)
    segMask = logical(segMask);
    if ~any(segMask(:))
        return;
    end

    cc = bwconncomp(segMask);
    if cc.NumObjects <= 0
        segMask(:) = false;
        return;
    end

    keepMask = false(size(segMask));
    for k = 1:cc.NumObjects
        if numel(cc.PixelIdxList{k}) >= minAreaPx
            keepMask(cc.PixelIdxList{k}) = true;
        end
    end
    segMask = keepMask;
end

function [blurredMaskedImage, normalizedClusterImage, clusterPeakBlur] = buildNormalizedClusterImageLocal(im, labelImage, cfg)
    blurredMaskedImage = zeros(size(im));
    normalizedClusterImage = zeros(size(im));
    nObj = max(labelImage(:));
    clusterPeakBlur = zeros(nObj, 1);
    if nObj <= 0
        return;
    end

    flatSize = 9;
    if isfield(cfg, 'segFlatFilterSize') && isfinite(cfg.segFlatFilterSize) && cfg.segFlatFilterSize > 0
        flatSize = max(1, round(double(cfg.segFlatFilterSize)));
    end
    flatFilter = ones(flatSize, flatSize) / max(flatSize * flatSize, 1);
    if exist('imfilter', 'file') == 2
        blurredMaskedImage = imfilter(double(im) .* double(labelImage > 0), flatFilter, 'conv');
    else
        blurredMaskedImage = conv2(double(im) .* double(labelImage > 0), flatFilter, 'same');
    end

    targetPeak = 250;
    if isfield(cfg, 'segNormalizedClusterPeak') && isfinite(cfg.segNormalizedClusterPeak) && cfg.segNormalizedClusterPeak > 0
        targetPeak = double(cfg.segNormalizedClusterPeak);
    end

    for k = 1:nObj
        objMask = (labelImage == k);
        peakVal = max(blurredMaskedImage(objMask));
        if ~isfinite(peakVal) || peakVal <= 0
            peakVal = 1;
        end
        clusterPeakBlur(k) = peakVal;
        normalizedClusterImage(objMask) = blurredMaskedImage(objMask) .* (targetPeak / peakVal);
    end
end

function mode = resolveSegmentationModeLocal(cfg)
    mode = 'variance';
    if isfield(cfg, 'segmentationMode') && ~isempty(cfg.segmentationMode)
        mode = lower(strtrim(char(cfg.segmentationMode)));
    end
    if ~ismember(mode, {'variance', 'threshold', 'hybrid'})
        mode = 'variance';
    end
end

function varImage = localVarianceImageFallback(im, winPx)
    kernel = ones(winPx, winPx) / max(winPx * winPx, 1);
    localMean = conv2(im, kernel, 'same');
    localMeanSq = conv2(im.^2, kernel, 'same');
    varImage = max(localMeanSq - localMean.^2, 0);
end

function segMask = applyObjectSizeLimits(segMask, labelImage, cfg)
    segMask = logical(segMask);
    if nargin < 2 || isempty(labelImage)
        labelImage = bwlabel(segMask);
    end
    nObj = max(labelImage(:));
    if nObj <= 0
        return;
    end

    stats = regionprops(labelImage, 'Area', 'EquivDiameter');
    remove = false(nObj, 1);

    if isfinite(cfg.maxCellArea)
        areas = [stats.Area];
        remove = remove | (areas(:) > double(cfg.maxCellArea));
    end

    if isfield(cfg, 'minCellDiameterPixels') && isfinite(cfg.minCellDiameterPixels) && cfg.minCellDiameterPixels > 0
        eqd = [stats.EquivDiameter];
        remove = remove | (eqd(:) < double(cfg.minCellDiameterPixels));
    else
        areas = [stats.Area];
        remove = remove | (areas(:) < double(max(1, round(cfg.minCellAreaPixelsEffective))));
    end

    if any(remove)
        removeIdx = find(remove);
        for k = 1:numel(removeIdx)
            segMask(labelImage == removeIdx(k)) = false;
        end
    end
end

function segMask = solidifySegmentationMask(segMask, cfg)
    segMask = logical(segMask);
    if ~any(segMask(:))
        return;
    end

    if isfield(cfg, 'segCloseRadiusPixelsEffective') && cfg.segCloseRadiusPixelsEffective > 0 && exist('imclose', 'file') == 2
        try
            se = strel('disk', round(cfg.segCloseRadiusPixelsEffective), 0);
            segMask = imclose(segMask, se);
        catch
        end
    end

    if isfield(cfg, 'fillHoles') && cfg.fillHoles && exist('imfill', 'file') == 2
        segMask = imfill(segMask, 'holes');
    end
end

function segMask = clearBorderObjectsLocal(segMask, cfg)
    segMask = logical(segMask);
    if ~any(segMask(:))
        return;
    end
    if ~isfield(cfg, 'segRemoveBorderObjects') || ~cfg.segRemoveBorderObjects
        return;
    end

    margin = 1;
    if isfield(cfg, 'segBorderMarginPixels') && isfinite(cfg.segBorderMarginPixels) && cfg.segBorderMarginPixels >= 0
        margin = max(1, round(double(cfg.segBorderMarginPixels)));
    end

    if exist('imclearborder', 'file') == 2 && margin <= 1
        segMask = imclearborder(segMask);
        return;
    end

    labelImage = bwlabel(segMask);
    nObj = max(labelImage(:));
    if nObj <= 0
        return;
    end

    edgeMask = false(size(segMask));
    margin = min([margin, size(segMask, 1), size(segMask, 2)]);
    edgeMask(1:margin, :) = true;
    edgeMask(end-margin+1:end, :) = true;
    edgeMask(:, 1:margin) = true;
    edgeMask(:, end-margin+1:end) = true;

    removeIds = unique(labelImage(edgeMask));
    removeIds(removeIds <= 0) = [];
    for k = 1:numel(removeIds)
        segMask(labelImage == removeIds(k)) = false;
    end
end

function segMask = fillLabeledObjects(labelImage, cfg)
    segMask = false(size(labelImage));
    nObj = max(labelImage(:));
    if nObj <= 0
        return;
    end

    for k = 1:nObj
        objMask = (labelImage == k);
        objMask = solidifySegmentationMask(objMask, cfg);
        if exist('imfill', 'file') == 2
            objMask = imfill(objMask, 'holes');
        end
        segMask = segMask | objMask;
    end
end

function level = localGraythresh(im)
    im = double(im);
    im = im - min(im(:));
    maxVal = max(im(:));
    if maxVal <= 0
        level = 0;
        return;
    end
    im = im / maxVal;
    if exist('graythresh', 'file') == 2
        level = graythresh(im);
    else
        level = 0.5;
    end
end

function v = localMAD(x)
    x = double(x(:));
    med = median(x);
    v = median(abs(x - med));
end

function [roiCurves, bgCurve, globalCurve, roiPhotonCounts, bgPhotonCount] = curvesFromLabels(gated, labelImage, gateLen)
    nRoi = max(labelImage(:));
    roiCurves = zeros(gateLen, nRoi);
    bgCurve = zeros(gateLen, 1);
    globalCurve = zeros(gateLen, 1);
    roiPhotonCounts = zeros(nRoi, 1);
    bgPhotonCount = 0;

    if isempty(gated.t)
        return;
    end

    t = round(gated.t(:));
    pix = round(gated.pix(:));
    roiId = double(labelImage(pix));

    globalCurve = accumarray(t, 1, [gateLen, 1], @sum, 0);

    inRoi = roiId > 0;
    if any(inRoi)
        roiCurves = accumarray([t(inRoi), roiId(inRoi)], 1, [gateLen, nRoi], @sum, 0);
        roiPhotonCounts = sum(roiCurves, 1).';
    end

    bgCurve = accumarray(t(~inRoi), 1, [gateLen, 1], @sum, 0);
    bgPhotonCount = sum(bgCurve);
end

function fitRes = emptyFitResult(nBins)
    fitRes = struct( ...
        'ok', false, ...
        'tauFitNs', [], ...
        'tauAvgNs', NaN, ...
        'tauAmpNs', NaN, ...
        'tauIntNs', NaN, ...
        'coeff', [], ...
        'backgroundCoeff', NaN, ...
        'speciesCoeff', [], ...
        'speciesFrac', [], ...
        'includeBG', true, ...
        'modelMatrix', [], ...
        'fitCounts', zeros(nBins, 1), ...
        'chi2red', NaN, ...
        'nPhotons', 0, ...
        'message', '');
end

function fitRes = fitCurveWithFixedIRF(counts, dtNs, irf, tau0, cfg)
    counts = double(counts(:));
    nBins = numel(counts);
    fitRes = emptyFitResult(nBins);
    fitRes.nPhotons = sum(counts);
    fitRes.includeBG = logical(cfg.includeBackgroundInFit);

    if fitRes.nPhotons < cfg.minPhotonsPerFit
        fitRes.message = sprintf('Too few photons for fitting (%d < %d).', round(fitRes.nPhotons), round(cfg.minPhotonsPerFit));
        return;
    end

    tau0 = double(tau0(:));
    irf = max(double(irf(:)), 0);
    irf = irf ./ max(sum(irf), eps);

    if cfg.useGuiStyleMultiExpFit
        [tauSorted, coeffSorted, bestFit, stats] = fitDecayIRFModelMultiStartLocal( ...
            counts, dtNs, tau0, irf, fitRes.includeBG, cfg.optimizeTau, cfg.fitMaxIter);
    else
        [tauSorted, ~] = sortLifetimesAndAmpsLocal(tau0(:), []);
        [~, coeffSorted, bestFit] = roiTcspcErrRawIRFLocal(tauSorted(:), counts, dtNs, irf, fitRes.includeBG);
        stats = calcDecayFitStatsLocal(counts, bestFit, numel(tauSorted), fitRes.includeBG);
    end

    speciesCoeff = coeffSorted(1 + double(fitRes.includeBG):end);
    speciesCoeff = max(speciesCoeff(:), 0);
    speciesFrac = speciesCoeff ./ max(sum(speciesCoeff), eps);
    tauAmpNs = weightedTauAverage(speciesCoeff, tauSorted(:), 'amplitude');
    tauIntNs = weightedTauAverage(speciesCoeff, tauSorted(:), 'intensity');
    modelMatrix = buildGuiDecayModelFromIRFLocal(irf, dtNs, tauSorted(:), fitRes.includeBG);

    fitRes.ok = true;
    fitRes.tauFitNs = tauSorted(:).';
    fitRes.tauAvgNs = tauIntNs;
    fitRes.tauAmpNs = tauAmpNs;
    fitRes.tauIntNs = tauIntNs;
    fitRes.coeff = coeffSorted(:).';
    if fitRes.includeBG && ~isempty(coeffSorted)
        fitRes.backgroundCoeff = coeffSorted(1);
    else
        fitRes.backgroundCoeff = 0;
    end
    fitRes.speciesCoeff = speciesCoeff(:).';
    fitRes.speciesFrac = speciesFrac(:).';
    fitRes.modelMatrix = modelMatrix;
    fitRes.fitCounts = bestFit(:);
    fitRes.chi2red = stats.chi2red;
    fitRes.message = '';
end

function [tauFit, coeff, fitCountsRawFull, stats] = fitDecayIRFModelMultiStartLocal(countsRawFull, dtNs, tau0, irf, includeBG, optimizeTau, maxIter)
    tau0 = double(tau0(:));
    countsRawFull = double(countsRawFull(:));
    irf = max(double(irf(:)), 0);
    irf = irf ./ max(sum(irf), eps);
    seedMat = buildTauSeedMatrixLocal(tau0);
    best.err = inf;
    best.tauFit = tau0(:).';
    best.coeff = zeros(numel(tau0) + double(includeBG), 1);
    best.fitCountsRawFull = zeros(size(countsRawFull));
    best.stats = struct('chi2red', NaN);
    for iseed = 1:size(seedMat,2)
        tauSeed = seedMat(:,iseed);
        if optimizeTau
            p0 = tauSeed;
            xmin = max(0.03, p0 / 10);
            xmax = max(p0 * 10, p0 + 0.05);
            tol = 1e-5;
            steps = max(maxIter, 180 * numel(p0));
            try
                if exist('Simplex', 'file') == 2
                    pfit = Simplex(@roiTcspcErrRawIRFLocal, p0, xmin, xmax, tol, steps, countsRawFull, dtNs, irf, includeBG);
                else
                    pfit = fminsearch(@(p) boundedRoiTcspcErrRawIRFLocal(p, xmin, xmax, countsRawFull, dtNs, irf, includeBG), ...
                        p0, optimset('Display', 'off', 'MaxIter', steps));
                    pfit = min(max(pfit(:), xmin), xmax);
                end
            catch
                pfit = p0;
            end
        else
            pfit = tauSeed;
        end
        [err, coeffCand, fitCand] = roiTcspcErrRawIRFLocal(pfit, countsRawFull, dtNs, irf, includeBG);
        [tauSorted, ~] = sortLifetimesAndAmpsLocal(pfit(:), []);
        [err, coeffCand, fitCand] = roiTcspcErrRawIRFLocal(tauSorted(:), countsRawFull, dtNs, irf, includeBG);
        statsCand = calcDecayFitStatsLocal(countsRawFull, fitCand, numel(tauSorted), includeBG);
        if err < best.err
            best.err = err;
            best.tauFit = tauSorted(:).';
            best.coeff = coeffCand(:);
            best.fitCountsRawFull = fitCand(:);
            best.stats = statsCand;
        end
    end
    [tauFit, coeff] = sortLifetimesAndAmpsLocal(best.tauFit(:), best.coeff(:));
    [~, coeff, fitCountsRawFull] = roiTcspcErrRawIRFLocal(tauFit(:), countsRawFull, dtNs, irf, includeBG);
    stats = calcDecayFitStatsLocal(countsRawFull, fitCountsRawFull, numel(tauFit), includeBG);
end

function err = boundedRoiTcspcErrRawIRFLocal(tau, tauLo, tauHi, counts, dtNs, irf, includeBG)
    tau = min(max(double(tau(:)), tauLo(:)), tauHi(:));
    err = roiTcspcErrRawIRFLocal(tau, counts, dtNs, irf, includeBG);
end

function [err, coeff, fitCounts] = roiTcspcErrRawIRFLocal(tau, counts, dtNs, irf, includeBG)
    tau = max(double(tau(:)), 1e-6);
    counts = max(double(counts(:)), 0);
    M = buildGuiDecayModelFromIRFLocal(irf, dtNs, tau, includeBG);
    coeff = lsqnonneg(M, counts);
    fitCounts = M * coeff;
    err = sum((counts - fitCounts).^2 ./ max(fitCounts, 1));
end

function M = buildGuiDecayModelFromIRFLocal(irf, dtNs, tauVec, includeBG)
    irf = max(double(irf(:)), 0);
    irf = irf ./ max(sum(irf), eps);
    tauVec = max(double(tauVec(:)), eps);
    nBins = numel(irf);
    nExp = numel(tauVec);
    t0 = ((0:nBins-1)') * dtNs;
    t1 = ((1:nBins)') * dtNs;
    M = zeros(nBins, nExp + double(includeBG));
    if includeBG
        M(:,1) = 1;
    end
    for k = 1:nExp
        tk = tauVec(k);
        decayBin = tk * (exp(-t0 ./ tk) - exp(-t1 ./ tk));
        convSig = conv(irf, decayBin, 'full');
        M(:, double(includeBG) + k) = convSig(1:nBins);
    end
end

function stats = calcDecayFitStatsLocal(counts, fitCounts, nExp, includeBG)
    counts = max(double(counts(:)), 0);
    fitCounts = max(double(fitCounts(:)), eps);
    n = numel(counts);
    k = nExp + double(includeBG);
    chi2 = sum((counts - fitCounts).^2 ./ max(fitCounts, 1));
    ndf = max(n - k, 1);
    chi2red = chi2 / ndf;
    bic = chi2 + k * log(max(n,2));
    aicc = chi2 + 2*k + (2*k*(k+1)) / max(n-k-1, 1);
    stats = struct('chi2', chi2, 'chi2red', chi2red, 'bic', bic, 'aicc', aicc, ...
        'n', n, 'k', k, 'includeBG', logical(includeBG));
end

function tauOut = weightedTauAverage(coeff, tau, modeName)
    coeff = max(double(coeff(:)), 0);
    tau = max(double(tau(:)), eps);
    switch lower(modeName)
        case 'intensity'
            denom = sum(coeff .* tau);
            if denom > 0
                tauOut = sum(coeff .* tau.^2) / denom;
            else
                tauOut = NaN;
            end
        otherwise
            denom = sum(coeff);
            if denom > 0
                tauOut = sum(coeff .* tau) / denom;
            else
                tauOut = NaN;
            end
    end
end

function seedMat = buildTauSeedMatrixLocal(tau0)
    tau0 = double(tau0(:));
    scales = [1.0 0.75 1.35];
    seedMat = zeros(numel(tau0), numel(scales));
    for k = 1:numel(scales)
        seedMat(:,k) = max(0.03, tau0 * scales(k));
    end
end

function [tauSorted, coeffSorted] = sortLifetimesAndAmpsLocal(tau, coeff)
    tau = double(tau(:));
    [tauSorted, ord] = sort(tau, 'ascend');
    coeffSorted = coeff;
    if ~isempty(coeff)
        coeff = double(coeff(:));
        if numel(coeff) == numel(tau) + 1
            coeffSorted = [coeff(1); coeff(1 + ord)];
        else
            coeffSorted = coeff;
        end
    end
end

function cube = tcspcCubeFromGated(gated, nx, ny, gateLen)
    cube = zeros(nx, ny, gateLen);
    if isempty(gated) || ~isfield(gated, 'pix') || isempty(gated.pix) || ~isfield(gated, 't') || isempty(gated.t)
        return;
    end

    pix = round(double(gated.pix(:)));
    t = round(double(gated.t(:)));
    valid = pix >= 1 & pix <= (nx * ny) & t >= 1 & t <= gateLen;
    if ~any(valid)
        return;
    end

    counts = accumarray([pix(valid), t(valid)], 1, [nx * ny, gateLen], @sum, 0);
    cube = reshape(counts, [nx, ny, gateLen]);
end

function [tauImage, componentAmpMaps, componentFracMaps, pixelFitInfo] = buildPixelwiseRoiLifetimeMap(cube, labelImage, roiFits, backgroundFit, dtNs, cfg, progressState, tAxisNs)
    %#ok<INUSD>
    [nx, ny, gateLen] = size(cube);
    nPix = nx * ny;
    nRoi = max(labelImage(:));
    nComp = numel(cfg.tau0);

    cube2d = reshape(double(cube), [nPix, gateLen]);
    fitImageFlat = sum(cube2d, 2);
    labelFlat = double(labelImage(:));

    tauFlat = nan(nPix, 1);
    ampFlat = nan(nPix, nComp);
    fracFlat = nan(nPix, nComp);

    bgMask = (labelFlat == 0) & (fitImageFlat > 0);
    if isfinite(backgroundFit.tauIntNs)
        tauFlat(bgMask) = backgroundFit.tauIntNs;
    end

    roiInfoTemplate = struct( ...
        'index', [], ...
        'tauFitNs', [], ...
        'tauAmpNs', NaN, ...
        'tauIntNs', NaN, ...
        'speciesCoeff', [], ...
        'speciesFrac', [], ...
        'nPixelsTotal', 0, ...
        'nPixelsEligible', 0, ...
        'nPixelsFit', 0, ...
        'nPixelsFallback', 0, ...
        'pixelPhotonCounts', []);
    pixelFitInfo = struct();
    pixelFitInfo.nComponents = nComp;
    pixelFitInfo.pixelwiseMinPhotons = cfg.pixelwiseMinPhotons;
    pixelFitInfo.pixelwiseIncludeBackground = logical(cfg.pixelwiseIncludeBackground);
    pixelFitInfo.backgroundTauIntNs = backgroundFit.tauIntNs;
    pixelFitInfo.roi = repmat(roiInfoTemplate, nRoi, 1);

    progressStride = max(1, ceil(max(nRoi, 1) / 10));
    for k = 1:nRoi
        roiMask = (labelFlat == k);
        pixIdx = find(roiMask);
        nPixRoi = numel(pixIdx);

        pixelFitInfo.roi(k).index = k;
        pixelFitInfo.roi(k).nPixelsTotal = nPixRoi;

        if nPixRoi == 0
            continue;
        end

        tauFallback = roiFits(k).tauIntNs;
        if ~isfinite(tauFallback)
            tauFallback = backgroundFit.tauIntNs;
        end
        if isfinite(tauFallback)
            tauFlat(pixIdx) = tauFallback;
        end

        pixelFitInfo.roi(k).tauFitNs = roiFits(k).tauFitNs;
        pixelFitInfo.roi(k).tauAmpNs = roiFits(k).tauAmpNs;
        pixelFitInfo.roi(k).tauIntNs = roiFits(k).tauIntNs;
        pixelFitInfo.roi(k).speciesCoeff = roiFits(k).speciesCoeff;
        pixelFitInfo.roi(k).speciesFrac = roiFits(k).speciesFrac;

        if ~(roiFits(k).ok && numel(roiFits(k).tauFitNs) == nComp && ~isempty(roiFits(k).modelMatrix))
            pixelFitInfo.roi(k).nPixelsFallback = nPixRoi;
            continue;
        end

        roiCounts = cube2d(pixIdx, :);
        photonTotals = sum(roiCounts, 2);
        eligible = photonTotals >= cfg.pixelwiseMinPhotons;
        pixelFitInfo.roi(k).pixelPhotonCounts = photonTotals;
        pixelFitInfo.roi(k).nPixelsEligible = sum(eligible);
        pixelFitInfo.roi(k).nPixelsFallback = nPixRoi - pixelFitInfo.roi(k).nPixelsEligible;

        if ~any(eligible)
            continue;
        end

        M = roiFits(k).modelMatrix;
        if roiFits(k).includeBG && ~cfg.pixelwiseIncludeBackground
            M = M(:, 2:end);
        end
        tauComp = double(roiFits(k).tauFitNs(:));

        eligibleIdx = find(eligible);
        for p = 1:numel(eligibleIdx)
            rowIdx = eligibleIdx(p);
            y = double(roiCounts(rowIdx, :).');
            coeffPix = lsqnonneg(M, y);

            if roiFits(k).includeBG && cfg.pixelwiseIncludeBackground
                speciesCoeffPix = coeffPix(2:end);
            else
                speciesCoeffPix = coeffPix(:);
            end
            speciesCoeffPix = max(double(speciesCoeffPix(:)), 0);

            if numel(speciesCoeffPix) ~= nComp
                pixelFitInfo.roi(k).nPixelsFallback = pixelFitInfo.roi(k).nPixelsFallback + 1;
                continue;
            end

            tauPix = weightedTauAverage(speciesCoeffPix, tauComp, 'intensity');
            if isfinite(tauPix)
                tauFlat(pixIdx(rowIdx)) = tauPix;
                ampFlat(pixIdx(rowIdx), :) = speciesCoeffPix(:).';
                fracFlat(pixIdx(rowIdx), :) = (speciesCoeffPix ./ max(sum(speciesCoeffPix), eps)).';
                pixelFitInfo.roi(k).nPixelsFit = pixelFitInfo.roi(k).nPixelsFit + 1;
            else
                pixelFitInfo.roi(k).nPixelsFallback = pixelFitInfo.roi(k).nPixelsFallback + 1;
            end
        end

        if mod(k, progressStride) == 0 || k == nRoi
            updateProgressState(progressState, sprintf('Pixelwise ROI %d/%d', k, nRoi), struct( ...
                'currentCurve', sum(roiCounts, 1).', ...
                'currentFit', roiFits(k).fitCounts, ...
                'tNs', tAxisNs, ...
                'res', struct('images', struct('averageLifetimeImageNs', reshape(tauFlat, [nx, ny])))));
        end
    end

    tauImage = reshape(tauFlat, [nx, ny]);
    componentAmpMaps = reshape(ampFlat, [nx, ny, nComp]);
    componentFracMaps = reshape(fracFlat, [nx, ny, nComp]);
end

function figHandles = showSummaryFigures(res, visibleState)
    if nargin < 2 || isempty(visibleState)
        visibleState = 'on';
    end
    segImage = res.images.segmentationSourceCh1D1;
    fitImage = res.images.analysisSourceCh2D2;
    labelImage = res.images.labelImage;
    tauImage = res.images.averageLifetimeImageNs;
    roiStats = res.rois;
    bgTau = res.background.fitCh2D2.tauIntNs;
    tNs = res.tcspc.tNs;

    f1 = figure('Name', '2-color cell segmentation and ROI lifetime map', 'Color', 'w', 'Visible', visibleState);
    tl = tiledlayout(f1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tl, 1);
    imagesc(ax1, segImage');
    axis(ax1, 'image');
    colormap(ax1, 'gray');
    colorbar(ax1);
    title(ax1, sprintf('Segmentation source: detector %d (ID %d), pulse %d', ...
        res.segDetectorIndex, res.segDetectorID, res.config.segPulseIndex));
    hold(ax1, 'on');
    contour(ax1, labelImage' > 0, [0.5 0.5], 'y', 'LineWidth', 0.8);
    hold(ax1, 'off');

    ax2 = nexttile(tl, 2);
    imagesc(ax2, labelImage' > 0);
    axis(ax2, 'image');
    colormap(ax2, gray);
    colorbar(ax2);
    title(ax2, sprintf('Segmented cells (n = %d)', max(labelImage(:))));
    hold(ax2, 'on');
    for k = 1:numel(roiStats)
        c = roiStats(k).centroidXY;
        text(ax2, c(2), c(1), sprintf('%d', k), 'Color', 'y', ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    hold(ax2, 'off');

    ax3 = nexttile(tl, 3);
    imagesc(ax3, fitImage');
    axis(ax3, 'image');
    colormap(ax3, 'gray');
    colorbar(ax3);
    title(ax3, sprintf('Analysis source: detector %d (ID %d), pulse %d', ...
        res.fitDetectorIndex, res.fitDetectorID, res.config.fitPulseIndex));

    ax4 = nexttile(tl, 4);
    imagesc(ax4, tauImage', 'AlphaData', ~isnan(tauImage'));
    axis(ax4, 'image');
    colormap(ax4, parula);
    cb = colorbar(ax4);
    ylabel(cb, 'Intensity-averaged lifetime (ns)');
    title(ax4, 'Pixel/ROI intensity-averaged lifetime image');
    hold(ax4, 'on');
    for k = 1:numel(roiStats)
        c = roiStats(k).centroidXY;
        if isfinite(roiStats(k).fitCh2D2.tauIntNs)
            lbl = sprintf('%.2f', roiStats(k).fitCh2D2.tauIntNs);
        else
            lbl = 'NaN';
        end
        text(ax4, c(2), c(1), lbl, 'Color', 'w', ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    if isfinite(bgTau)
        text(ax4, 5, 8, sprintf('BG %.2f ns', bgTau), 'Color', 'w', ...
            'FontWeight', 'bold', 'BackgroundColor', 'k', 'Margin', 2, ...
            'VerticalAlignment', 'top');
    end
    hold(ax4, 'off');

    nCurves = numel(roiStats) + 1;
    nCols = min(4, max(1, ceil(sqrt(nCurves))));
    nRows = ceil(nCurves / nCols);
    f2 = figure('Name', 'Detector 2 / pulse 2 TCSPC fits', 'Color', 'w', 'Visible', visibleState);
    tl2 = tiledlayout(f2, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile(tl2);
    semilogy(tNs, max(res.background.curveCh2D2, 1e-3), 'k.', 'MarkerSize', 8);
    hold on;
    semilogy(tNs, max(res.background.fitCh2D2.fitCounts, 1e-3), 'r-', 'LineWidth', 1.1);
    hold off;
    grid on;
    xlabel('Delay time (ns)');
    ylabel('Counts');
    title(sprintf('Background, <\\tau>_{int} = %.2f ns', res.background.fitCh2D2.tauIntNs));

    for k = 1:numel(roiStats)
        nexttile(tl2);
        semilogy(tNs, max(roiStats(k).curveCh2D2, 1e-3), 'k.', 'MarkerSize', 8);
        hold on;
        semilogy(tNs, max(roiStats(k).fitCh2D2.fitCounts, 1e-3), 'r-', 'LineWidth', 1.1);
        hold off;
        grid on;
        xlabel('Delay time (ns)');
        ylabel('Counts');
        title(sprintf('Cell %d, <\\tau>_{int} = %.2f ns', k, roiStats(k).fitCh2D2.tauIntNs));
    end

    figHandles = [f1; f2];
end

function saved = saveResultOutputs(res, name, cfg, figHandles)
    [outputDir, outputStem] = resolveOutputSpec(name, cfg);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    saved = struct();
    saved.outputDir = outputDir;

    if cfg.saveMat
        saved.matPath = fullfile(outputDir, [outputStem '_cellseg2color.mat']);
        resToSave = res; %#ok<NASGU>
        save(saved.matPath, 'resToSave', '-v7.3');
    else
        saved.matPath = '';
    end

    if cfg.saveCsv
        tcspcCh1 = buildTcspcTable(res.tcspc.tNs, res.background.curveCh1D1, res.tcspc.roiCh1D1, 'ch1d1');
        tcspcCh2 = buildTcspcTable(res.tcspc.tNs, res.background.curveCh2D2, res.tcspc.roiCh2D2, 'ch2d2');
        roiSummary = buildRoiSummaryTable(res.rois, numel(res.config.tau0));

        saved.tcspcCh1Csv = fullfile(outputDir, [outputStem '_tcspc_ch1d1.csv']);
        saved.tcspcCh2Csv = fullfile(outputDir, [outputStem '_tcspc_ch2d2.csv']);
        saved.roiSummaryCsv = fullfile(outputDir, [outputStem '_roi_summary.csv']);
        localWriteTable(tcspcCh1, saved.tcspcCh1Csv);
        localWriteTable(tcspcCh2, saved.tcspcCh2Csv);
        localWriteTable(roiSummary, saved.roiSummaryCsv);
    else
        saved.tcspcCh1Csv = '';
        saved.tcspcCh2Csv = '';
        saved.roiSummaryCsv = '';
    end

    if cfg.saveFigures && ~isempty(figHandles)
        saved.summaryFigurePng = fullfile(outputDir, [outputStem '_summary.png']);
        saved.fitFigurePng = fullfile(outputDir, [outputStem '_fits.png']);
        saveFigurePng(figHandles(1), saved.summaryFigurePng);
        if numel(figHandles) >= 2
            saveFigurePng(figHandles(2), saved.fitFigurePng);
        else
            saved.fitFigurePng = '';
        end
    else
        saved.summaryFigurePng = '';
        saved.fitFigurePng = '';
    end
end

function [outputDir, outputStem] = resolveOutputSpec(name, cfg)
    [parentDir, stem] = fileparts(name);
    outputStem = stem;
    if isfield(cfg, 'outputStem') && ~isempty(cfg.outputStem)
        outputStem = char(cfg.outputStem);
    end

    if isfield(cfg, 'outputDir') && ~isempty(cfg.outputDir)
        outputDir = char(cfg.outputDir);
    else
        outputDir = fullfile(parentDir, [stem '_CellSeg2Color_results']);
    end
end

function tbl = buildTcspcTable(tNs, bgCurve, roiCurves, curveTag)
    tbl = table(double(tNs(:)), double(bgCurve(:)), 'VariableNames', {'time_ns', ['background_' curveTag]});
    nRoi = size(roiCurves, 2);
    for k = 1:nRoi
        tbl.(['roi_' sprintf('%03d', k) '_' curveTag]) = double(roiCurves(:,k));
    end
end

function tbl = buildRoiSummaryTable(rois, nTau)
    nRoi = numel(rois);
    roiIndex = zeros(nRoi, 1);
    areaPixels = zeros(nRoi, 1);
    centroidX = zeros(nRoi, 1);
    centroidY = zeros(nRoi, 1);
    photonsCh1D1 = zeros(nRoi, 1);
    photonsCh2D2 = zeros(nRoi, 1);
    tauAvgNs = nan(nRoi, 1);
    tauAmpNs = nan(nRoi, 1);
    tauIntNs = nan(nRoi, 1);
    chi2red = nan(nRoi, 1);
    tauFit = nan(nRoi, nTau);
    componentCoeff = nan(nRoi, nTau);
    componentFrac = nan(nRoi, nTau);

    for k = 1:nRoi
        roiIndex(k) = rois(k).index;
        areaPixels(k) = rois(k).areaPixels;
        centroidX(k) = rois(k).centroidXY(1);
        centroidY(k) = rois(k).centroidXY(2);
        photonsCh1D1(k) = rois(k).photonsCh1D1;
        photonsCh2D2(k) = rois(k).photonsCh2D2;
        tauAvgNs(k) = rois(k).fitCh2D2.tauAvgNs;
        tauAmpNs(k) = rois(k).fitCh2D2.tauAmpNs;
        tauIntNs(k) = rois(k).fitCh2D2.tauIntNs;
        chi2red(k) = rois(k).fitCh2D2.chi2red;
        tauHere = rois(k).fitCh2D2.tauFitNs;
        coeffHere = rois(k).componentCoeffCh2D2;
        fracHere = rois(k).componentFracCh2D2;
        nHere = min(nTau, numel(tauHere));
        if nHere > 0
            tauFit(k, 1:nHere) = tauHere(1:nHere);
        end
        nCoeff = min(nTau, numel(coeffHere));
        if nCoeff > 0
            componentCoeff(k, 1:nCoeff) = coeffHere(1:nCoeff);
        end
        nFrac = min(nTau, numel(fracHere));
        if nFrac > 0
            componentFrac(k, 1:nFrac) = fracHere(1:nFrac);
        end
    end

    tbl = table(roiIndex, areaPixels, centroidX, centroidY, photonsCh1D1, photonsCh2D2, tauAvgNs, tauAmpNs, tauIntNs, chi2red);
    for j = 1:nTau
        tbl.(['tauFit_' num2str(j) '_ns']) = tauFit(:,j);
        tbl.(['componentCoeff_' num2str(j)]) = componentCoeff(:,j);
        tbl.(['componentFrac_' num2str(j)]) = componentFrac(:,j);
    end
end

function localWriteTable(tbl, outPath)
    if exist('writetable', 'file') == 2
        writetable(tbl, outPath);
    else
        writecell([tbl.Properties.VariableNames; table2cell(tbl)], outPath);
    end
end

function saveFigurePng(figHandle, outPath)
    if isempty(figHandle) || ~ishandle(figHandle)
        return;
    end
    try
        exportgraphics(figHandle, outPath, 'Resolution', 200);
    catch
        print(figHandle, outPath, '-dpng', '-r200');
    end
end

function out = ternaryString(cond, trueVal, falseVal)
    if cond
        out = trueVal;
    else
        out = falseVal;
    end
end

function deleteValidFigures(figHandles)
    for k = 1:numel(figHandles)
        if ~isempty(figHandles(k)) && ishandle(figHandles(k))
            delete(figHandles(k));
        end
    end
end

function logVerbose(cfg, varargin)
    if isfield(cfg, 'verbose') && cfg.verbose
        fprintf('[CellSeg2Color] ');
        fprintf(varargin{:});
        fprintf('\n');
    end
end

function state = initProgressState(cfg, name)
    state = struct('enabled', false, 'fig', [], 'ax1', [], 'ax2', [], 'ax3', [], 'ax4', [], 'titleText', '');
    if ~(isfield(cfg, 'plotProgress') && cfg.plotProgress)
        return;
    end

    state.enabled = true;
    state.fig = figure('Name', sprintf('CellSeg progress: %s', name), 'Color', 'w', 'Visible', 'on');
    tl = tiledlayout(state.fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    state.ax1 = nexttile(tl, 1);
    state.ax2 = nexttile(tl, 2);
    state.ax3 = nexttile(tl, 3);
    state.ax4 = nexttile(tl, 4);
    state.titleText = title(tl, 'Starting...');
    drawnow;
end

function updateProgressState(state, stageName, data)
    if nargin < 3
        data = struct();
    end
    if ~isstruct(state) || ~isfield(state, 'enabled') || ~state.enabled
        return;
    end
    if isempty(state.fig) || ~ishandle(state.fig)
        return;
    end

    if ~isempty(state.titleText) && isgraphics(state.titleText)
        state.titleText.String = stageName;
    end

    if isfield(data, 'tcspcByDetector')
        cla(state.ax1);
        semilogy(state.ax1, max(double(data.tcspcByDetector), 1e-3), 'LineWidth', 1.0);
        xlabel(state.ax1, 'TCSPC bin');
        ylabel(state.ax1, 'Counts');
        title(state.ax1, 'Detector histograms');
        grid(state.ax1, 'on');
        if isfield(data, 'gateStarts') && isfield(data, 'gateStops')
            hold(state.ax1, 'on');
            for k = 1:numel(data.gateStarts)
                xline(state.ax1, data.gateStarts(k), 'g--');
                xline(state.ax1, data.gateStops(k), 'r--');
            end
            hold(state.ax1, 'off');
        end
    elseif isfield(data, 'segImage')
        cla(state.ax1);
        imagesc(state.ax1, data.segImage');
        axis(state.ax1, 'image');
        colormap(state.ax1, 'gray');
        title(state.ax1, 'Segmentation source');
        colorbar(state.ax1);
    end

    if isfield(data, 'labelImage')
        cla(state.ax2);
        imagesc(state.ax2, (data.labelImage > 0)');
        axis(state.ax2, 'image');
        title(state.ax2, 'Segmentation mask');
        colorbar(state.ax2);
    elseif isfield(data, 'fitImage')
        cla(state.ax2);
        imagesc(state.ax2, data.fitImage');
        axis(state.ax2, 'image');
        colormap(state.ax2, 'gray');
        title(state.ax2, 'Analysis source');
        colorbar(state.ax2);
    end

    if isfield(data, 'globalCurveCh2D2') && isfield(data, 'tNs')
        cla(state.ax3);
        semilogy(state.ax3, data.tNs, max(double(data.globalCurveCh2D2(:)), 1e-3), 'k-', 'LineWidth', 1.0);
        hold(state.ax3, 'on');
        if isfield(data, 'irfOut') && isstruct(data.irfOut) && isfield(data.irfOut, 'IRF') && ~isempty(data.irfOut.IRF)
            irf = double(data.irfOut.IRF(:));
            irf = irf ./ max(irf, eps) * max(double(data.globalCurveCh2D2(:)));
            semilogy(state.ax3, data.tNs, max(irf, 1e-3), 'r--', 'LineWidth', 1.0);
        end
        hold(state.ax3, 'off');
        xlabel(state.ax3, 'Delay time (ns)');
        ylabel(state.ax3, 'Counts');
        title(state.ax3, 'Global ch2/d2');
        grid(state.ax3, 'on');
    elseif isfield(data, 'currentCurve') && isfield(data, 'tNs')
        cla(state.ax3);
        semilogy(state.ax3, data.tNs, max(double(data.currentCurve(:)), 1e-3), 'k.', 'MarkerSize', 8);
        hold(state.ax3, 'on');
        if isfield(data, 'currentFit') && ~isempty(data.currentFit)
            semilogy(state.ax3, data.tNs, max(double(data.currentFit(:)), 1e-3), 'r-', 'LineWidth', 1.0);
        end
        hold(state.ax3, 'off');
        xlabel(state.ax3, 'Delay time (ns)');
        ylabel(state.ax3, 'Counts');
        title(state.ax3, 'Current ROI fit');
        grid(state.ax3, 'on');
    end

    if isfield(data, 'fitImage')
        cla(state.ax4);
        imagesc(state.ax4, data.fitImage');
        axis(state.ax4, 'image');
        colormap(state.ax4, 'gray');
        title(state.ax4, 'ch2/d2 image');
        colorbar(state.ax4);
    elseif isfield(data, 'res') && isstruct(data.res)
        cla(state.ax4);
        imagesc(state.ax4, data.res.images.averageLifetimeImageNs', 'AlphaData', ~isnan(data.res.images.averageLifetimeImageNs'));
        axis(state.ax4, 'image');
        colormap(state.ax4, parula);
        title(state.ax4, 'Lifetime image');
        colorbar(state.ax4);
    end

    drawnow limitrate;
end

function closeProgressState(state, cfg)
    if ~isstruct(state) || ~isfield(state, 'enabled') || ~state.enabled
        return;
    end
    if ~isfield(cfg, 'closeProgressFigure') || ~cfg.closeProgressFigure
        return;
    end
    if ~isempty(state.fig) && ishandle(state.fig)
        delete(state.fig);
    end
end
