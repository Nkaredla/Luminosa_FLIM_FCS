function result = estimateFocalPlaneCenterPointISMPoisson(stackInput, varargin)
%ESTIMATEFOCALPLANECENTERPOINTISMPOISSON Fit focal center microimage by Poisson MLE.
%
%   result = estimateFocalPlaneCenterPointISMPoisson(stackInput)
%
%   This is a single-plane experimental estimator. It loads a detector-resolved
%   ISM z stack, selects the focal plane, fits the bead center in the
%   detector-summed focal image, extracts the 23 raw detector counts at that
%   scan pixel, and estimates aberration coefficients by Poisson likelihood:
%
%   ITERATIVE CENTERING
%   The detector-summed focal spot is the aberrated PSF, so asymmetric
%   aberrations (notably coma) shift its peak away from the emitter's true
%   on-axis position, while the microimage model assumes an on-axis emitter.
%   With 'iterativeCenter' true (default) the center and the aberration fit
%   are solved by coordinate descent: fit aberrations at the current center,
%   then relocate the (sub-pixel) sampling center to the position whose data
%   microimage best matches the on-axis model microimage (minimum Poisson
%   deviance), and repeat. Do not combine this with 'fitXY' or with tilt in
%   'fitModes' (both re-introduce the lateral-offset/tilt degeneracy). Set
%   'iterativeCenter', false for the legacy single-pass Gaussian-center
%   behavior. For sub-pixel recentering use 'centerSampleMode', 'subpixel'.
%
%       mu_k = photonScale * model_k(coeffs) + backgroundScale * dark_k
%
%   No defocused/diversity plane is used. Even-mode signs are therefore not
%   resolved by axial diversity and the result should be treated as a focal
%   center-microimage diagnostic.
%
%   A diagnostic figure is shown (and saved) with the observed, fitted, and
%   Poisson-residual detector micro-images drawn on the honeycomb hex layout,
%   the observed-vs-expected detector counts, and the fitted Zernike
%   aberration coefficients with +/-1 sigma uncertainty error bars (the
%   Fisher/Cramer-Rao standard deviations in result.fit.coeffStdWaves /
%   coeffStdNm). Set 'showFigures', false to suppress the on-screen figure
%   and 'writeOutputs', false to skip the saved PNG/CSVs.
%
%   For raw Luminosa scan folders, the function first looks for existing batch
%   alignment outputs. If none are found and autoPreprocessRawStack is true, it
%   runs batchAnalyzeLuminosaPSFs and uses the generated alignment CSV.
%
%   POLARIZATION
%   Excitation is modeled as circularly polarized by default
%   ('excitationPolarizationMode' = 'circular', a rotationally symmetric
%   focus). Set it to 'x'/'y' for a linearly polarized excitation laser. A
%   polarization mismatch (e.g. modeling linear when the rig is circular)
%   injects a fixed lateral anisotropy that the fit absorbs as spurious
%   astigmatism, so match this to the experiment.

    if nargin < 1 || isempty(stackInput)
        stackInput = 'D:\Luminosa\Aberration_ISM\PDA23_centered_bead_point_10s_20260702-124209';
    end

    opts = parseOptions(varargin{:});
    addRequiredPaths(opts);

    [data, inputInfo] = loadInputData(stackInput, opts);
    data = ensureDataFields(data, opts);

    [focusIdx, focusSelection] = chooseFocusPlane(data, opts);
    [focusImage, rawFocus, backgroundCounts] = focalPlaneForCentering(data, focusIdx);
    [includeCh, channelInfo] = selectedChannels( ...
        data.channelIDs(:), backgroundCounts, opts);

    [sim, detectorDiagnostics] = configureFocalSimulation(data, focusIdx, opts);
    p0 = initialParameterVector(sim, opts);

    % Initial center from the (aberration-biased) detector-summed image, then
    % refine it by coordinate descent so the sampling position is consistent
    % with the on-axis aberration model (see fitWithIterativeCenter).
    [initialCenterXY, ~, centerFitInfo] = estimateCenterXY(focusImage, opts);
    [fit, centerXY, obs, centerIteration] = fitWithIterativeCenter( ...
        rawFocus, backgroundCounts, includeCh, data.channelIDs(:), ...
        initialCenterXY, sim, opts, p0);
    fit = attachFitDiagnostics(fit, obs, sim, opts);

    centerCounts = obs.counts;
    centerInfo = buildFocalCenterInfo(focusImage, centerXY, initialCenterXY, ...
        centerFitInfo, data, opts, centerIteration);

    result = struct();
    result.method = 'focal-plane center-microimage Poisson MLE';
    result.input = inputInfo;
    result.data = data;
    result.focusSelection = focusSelection;
    result.focusIndex = focusIdx;
    result.focusZUm = data.stageZUm(focusIdx);
    result.centerCounts = centerCounts(:);
    result.backgroundCounts = backgroundCounts(:);
    result.centerInfo = centerInfo;
    result.centerIteration = centerIteration;
    result.channelInfo = channelInfo;
    result.includedChannels = find(includeCh(:)).';
    result.includedChannelIDs = data.channelIDs(includeCh);
    result.excludedChannelIDs = data.channelIDs(~includeCh);
    result.fit = fit;
    result.detectorLayoutDiagnostics = detectorDiagnostics;
    result.validity = struct( ...
        'correctionReady', false, ...
        'recommendedUse', 'diagnostic focal-plane center-microimage fit only', ...
        'reason', ['Only the focal plane is used, so axial-diversity sign ' ...
        'checks and full-stack residual validation are not available.']);
    result.options = opts;
    result.outputDir = resolveOutputDir(opts, focusSelection, inputInfo);

    if opts.writeOutputs
        writeOutputs(result);
    end
    if opts.verbose
        printSummary(result);
    end
    if opts.showFigures
        result.summaryFigure = buildSummaryFigure(result, 'on');
    end
end

function opts = parseOptions(varargin)
    p = inputParser;
    p.FunctionName = 'estimateFocalPlaneCenterPointISMPoisson';

    addParameter(p, 'stageZUm', []);
    addParameter(p, 'xyPixelSizeUm', []);
    addParameter(p, 'zStepUm', 0.05);
    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'channelOrder', []);
    addParameter(p, 'ptuReaderFolder', defaultPtuReaderFolder());
    addParameter(p, 'ptuPhotonsPerChunk', 1e6);
    addParameter(p, 'boundaryWidthPx', 3);
    addParameter(p, 'darkFile', 'D:\Luminosa\Data\ISMdark_counts.ptu');
    addParameter(p, 'darkScale', 1);
    addParameter(p, 'darkPerPixel', []);
    addParameter(p, 'darkMeasurementMode', 'tttr');
    addParameter(p, 'backgroundMode', 'auto');
    addParameter(p, 'flatField', []);
    addParameter(p, 'flatFieldVariable', '');
    addParameter(p, 'cacheFile', '');
    addParameter(p, 'reuseCache', true);

    addParameter(p, 'autoPreprocessRawStack', true);
    addParameter(p, 'batchResultsRoot', '');
    addParameter(p, 'batchInputSource', 'auto');
    addParameter(p, 'batchOptions', struct());

    addParameter(p, 'focusIndex', []);
    addParameter(p, 'focusZUm', []);
    addParameter(p, 'focusMetric', 'signal');
    addParameter(p, 'centerXY', []);
    addParameter(p, 'centerMode', 'gaussian');
    addParameter(p, 'centerThresholdFraction', 0.20);
    addParameter(p, 'centerSampleMode', 'nearest');
    addParameter(p, 'iterativeCenter', true);
    addParameter(p, 'maxCenterIterations', 5);
    addParameter(p, 'centerSearchRadiusPx', 2);
    addParameter(p, 'centerSearchStepPx', 0.25);
    addParameter(p, 'centerRefineStepPx', 0.05);
    addParameter(p, 'centerTolerancePx', 0.05);

    addParameter(p, 'sim', []);
    addParameter(p, 'objectiveNA', 1.2);
    addParameter(p, 'objectiveMagnification', 60);
    addParameter(p, 'emissionWavelengthUm', []);
    addParameter(p, 'excitationWavelengthUm', []);
    addParameter(p, 'mediumRefractiveIndex', []);
    addParameter(p, 'sampleGeometry', 'airOnGlass');
    addParameter(p, 'interfaceRadialWeightMode', 'sampleSolidAngle');
    addParameter(p, 'diffractionModel', 'vectorial');
    addParameter(p, 'vectorialPolarizationMode', 'xyAverage');
    addParameter(p, 'excitationPolarizationMode', 'circular');
    addParameter(p, 'collectionPolarizationMode', 'xyAverage');
    addParameter(p, 'airInterfaceStageMedium', 'immersion');
    addParameter(p, 'immersionRefractiveIndex', 1.33);
    addParameter(p, 'glassRefractiveIndex', 1.518);
    addParameter(p, 'sampleRefractiveIndex', 1.0003);
    addParameter(p, 'designGlassRefractiveIndex', 1.518);
    addParameter(p, 'coverslipThicknessUm', 190);
    addParameter(p, 'designCoverslipThicknessUm', 190);
    addParameter(p, 'beadBottomHeightUm', 0);
    addParameter(p, 'airBeadAxialSamples', 3);
    addParameter(p, 'beadRadiusUm', 0.08);
    addParameter(p, 'beadSubsamples', [3 3 3]);
    addParameter(p, 'modelDzUm', 0.025);
    addParameter(p, 'modelZPaddingUm', 0.50);
    addParameter(p, 'modelSampleXY', [0 0]);

    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detectorPixelShape', 'hex');
    addParameter(p, 'detectorSubsamples', []);
    addParameter(p, 'estimateDetectorLayout', false);
    addParameter(p, 'detectorXYUm', []);
    addParameter(p, 'detectorLayoutPositionSign', -1);
    addParameter(p, 'scanAxisSigns', [1 1]);
    addParameter(p, 'detectorLayoutScale', 2);
    addParameter(p, 'detectorLayoutCenterMode', 'reference');
    addParameter(p, 'detectorCenterIndex', []);
    addParameter(p, 'detectorShiftSmoothSigma', 1);
    addParameter(p, 'detectorShiftUseWindow', true);
    addParameter(p, 'detectorShiftNormalizeImages', true);
    addParameter(p, 'detectorShiftUpsample', 20);
    addParameter(p, 'detectorPitchSampleUm', []);
    addParameter(p, 'detectorHardwarePitchUm', []);
    addParameter(p, 'detectorTotalMagnification', []);
    addParameter(p, 'detectorQE', []);

    addParameter(p, 'fitModes', {'astig_x','astig_y','coma_x','coma_y','spherical'});
    addParameter(p, 'initialCoeffs', struct());
    addParameter(p, 'fitXY', false);
    addParameter(p, 'initialXY', [0 0]);
    addParameter(p, 'fitDetectorPitchScale', false);
    addParameter(p, 'initialDetectorPitchScale', 1);
    addParameter(p, 'fitBackgroundScale', true);
    addParameter(p, 'backgroundScale', 1);
    addParameter(p, 'backgroundScaleBounds', [0 5]);
    addParameter(p, 'maxIter', 8);
    addParameter(p, 'jacobianScheme', 'forward');
    addParameter(p, 'fdCoeff', 0.01);
    addParameter(p, 'fdXY', []);
    addParameter(p, 'fdDetectorPitchScale', 0.02);
    addParameter(p, 'fdBackgroundScale', 0.05);
    addParameter(p, 'regCoeff', 1e-5);
    addParameter(p, 'regXY', 1e-4);
    addParameter(p, 'regDetectorPitchScale', 1e-4);
    addParameter(p, 'regBackgroundScale', 1e-4);
    addParameter(p, 'maxCoeffStep', 0.04);
    addParameter(p, 'maxXYStep', 0.03);
    addParameter(p, 'maxDetectorPitchScaleStep', 0.05);
    addParameter(p, 'maxBackgroundScaleStep', 0.25);
    addParameter(p, 'coefficientBoundsWaves', [-0.30 0.30]);
    addParameter(p, 'detectorPitchScaleBounds', [0.5 2.0]);
    addParameter(p, 'tolStep', 1e-5);
    addParameter(p, 'minExpectedCount', 1e-9);
    addParameter(p, 'excludeChannelIDs', []);
    addParameter(p, 'excludeHotDarkPixels', false);
    addParameter(p, 'hotPixelFactor', 8);

    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'showFigures', true);
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;

    opts.channelIDs = double(opts.channelIDs(:)).';
    opts.channelOrder = double(opts.channelOrder(:)).';
    opts.ptuReaderFolder = char(opts.ptuReaderFolder);
    opts.darkFile = char(opts.darkFile);
    opts.darkMeasurementMode = lower(char(opts.darkMeasurementMode));
    opts.backgroundMode = lower(char(opts.backgroundMode));
    opts.cacheFile = char(opts.cacheFile);
    opts.flatFieldVariable = char(opts.flatFieldVariable);
    opts.batchResultsRoot = char(opts.batchResultsRoot);
    opts.batchInputSource = lower(char(opts.batchInputSource));
    opts.focusMetric = lower(char(opts.focusMetric));
    opts.centerMode = lower(char(opts.centerMode));
    opts.centerSampleMode = lower(char(opts.centerSampleMode));
    opts.sampleGeometry = char(opts.sampleGeometry);
    opts.interfaceRadialWeightMode = char(opts.interfaceRadialWeightMode);
    opts.diffractionModel = char(opts.diffractionModel);
    opts.vectorialPolarizationMode = char(opts.vectorialPolarizationMode);
    opts.excitationPolarizationMode = char(opts.excitationPolarizationMode);
    opts.collectionPolarizationMode = char(opts.collectionPolarizationMode);
    opts.airInterfaceStageMedium = char(opts.airInterfaceStageMedium);
    opts.detectorLayout = char(opts.detectorLayout);
    opts.detectorPixelShape = char(opts.detectorPixelShape);
    opts.detectorLayoutCenterMode = lower(char(opts.detectorLayoutCenterMode));
    opts.jacobianScheme = lower(char(opts.jacobianScheme));
    opts.outputDir = char(opts.outputDir);

    if ischar(opts.fitModes) || isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    elseif ~iscell(opts.fitModes)
        error('estimateFocalPlaneCenterPointISMPoisson:BadFitModes', ...
            'fitModes must be a character vector, string array, or cell array.');
    end
    if ~ismember(lower(opts.sampleGeometry), {'aironglass','homogeneous'})
        error('estimateFocalPlaneCenterPointISMPoisson:BadSampleGeometry', ...
            'sampleGeometry must be airOnGlass or homogeneous.');
    end
    if ~ismember(opts.centerSampleMode, {'nearest','round','pixel','subpixel','linear','bilinear','interp'})
        error('estimateFocalPlaneCenterPointISMPoisson:BadCenterSampleMode', ...
            'centerSampleMode must be nearest or subpixel.');
    end
    opts.detectorQE = validateDetectorQE(opts.detectorQE, opts.channelIDs);
    opts.backgroundScaleBounds = validateTwoElementBounds( ...
        opts.backgroundScaleBounds, 'backgroundScaleBounds');
    opts.coefficientBoundsWaves = validateTwoElementBounds( ...
        opts.coefficientBoundsWaves, 'coefficientBoundsWaves');
    opts.detectorPitchScaleBounds = validateTwoElementBounds( ...
        opts.detectorPitchScaleBounds, 'detectorPitchScaleBounds');
    if ~isstruct(opts.batchOptions) || ~isscalar(opts.batchOptions)
        error('estimateFocalPlaneCenterPointISMPoisson:BadBatchOptions', ...
            'batchOptions must be a scalar struct.');
    end
end

function folder = defaultPtuReaderFolder()
    thisDir = fileparts(mfilename('fullpath'));
    folder = fullfile(fileparts(fileparts(thisDir)), 'Luminosa_FLIM_FCS');
end

function addRequiredPaths(opts)
    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    parentDir = fileparts(thisDir);
    if exist(parentDir, 'dir') == 7
        addpath(parentDir);
    end
    if ~isempty(opts.ptuReaderFolder) && exist(opts.ptuReaderFolder, 'dir') == 7
        addpath(opts.ptuReaderFolder);
    end
end

function [data, inputInfo] = loadInputData(stackInput, opts)
    inputInfo = struct('originalInput', stackInput, 'resolvedInput', '', ...
        'preprocessedRawFolder', false);
    if isstruct(stackInput)
        data = stackInput;
        inputInfo.resolvedInput = 'preloaded data struct';
        return;
    end

    [resolvedInput, preprocessedRawFolder] = resolveStackInput(stackInput, opts);
    inputInfo.resolvedInput = resolvedInput;
    inputInfo.preprocessedRawFolder = preprocessedRawFolder;
    data = loadFullStackISMData(resolvedInput, ...
        'stageZUm', opts.stageZUm, ...
        'xyPixelSizeUm', opts.xyPixelSizeUm, ...
        'channelIDs', opts.channelIDs, ...
        'channelOrder', opts.channelOrder, ...
        'ptuReaderFolder', opts.ptuReaderFolder, ...
        'ptuPhotonsPerChunk', opts.ptuPhotonsPerChunk, ...
        'boundaryWidthPx', opts.boundaryWidthPx, ...
        'darkFile', opts.darkFile, ...
        'darkScale', opts.darkScale, ...
        'darkPerPixel', opts.darkPerPixel, ...
        'darkMeasurementMode', opts.darkMeasurementMode, ...
        'backgroundMode', opts.backgroundMode, ...
        'flatField', opts.flatField, ...
        'flatFieldVariable', opts.flatFieldVariable, ...
        'cacheFile', opts.cacheFile, ...
        'reuseCache', opts.reuseCache, ...
        'verbose', opts.verbose);
end

function [resolvedInput, preprocessedRawFolder] = resolveStackInput(stackInput, opts)
    preprocessedRawFolder = false;
    if isnumeric(stackInput)
        resolvedInput = stackInput;
        return;
    end
    inputPath = char(stackInput);
    if exist(inputPath, 'dir') ~= 7
        resolvedInput = inputPath;
        return;
    end

    csvFile = findBatchAlignmentCsvForScanFolder(inputPath, opts);
    if ~isempty(csvFile) && exist(csvFile, 'file') == 2
        resolvedInput = csvFile;
        return;
    end

    if opts.autoPreprocessRawStack && hasRawStackFiles(inputPath)
        preprocessRawScanFolder(inputPath, opts);
        preprocessedRawFolder = true;
        csvFile = findBatchAlignmentCsvForScanFolder(inputPath, opts);
        if ~isempty(csvFile) && exist(csvFile, 'file') == 2
            resolvedInput = csvFile;
            return;
        end
    end

    resolvedInput = inputPath;
end

function data = ensureDataFields(data, opts)
    if ~isfield(data, 'rawCounts')
        error('estimateFocalPlaneCenterPointISMPoisson:BadDataStruct', ...
            'Input data struct must contain rawCounts.');
    end
    data.rawCounts = double(data.rawCounts);
    [ny, nx, nCh, nPlane] = size(data.rawCounts);

    if ~isfield(data, 'channelIDs') || isempty(data.channelIDs)
        data.channelIDs = opts.channelIDs(1:nCh);
    end
    data.channelIDs = double(data.channelIDs(:)).';

    if ~isfield(data, 'stageZUm') || isempty(data.stageZUm)
        if numel(opts.stageZUm) == nPlane
            data.stageZUm = double(opts.stageZUm(:)).';
        else
            data.stageZUm = ((1:nPlane) - (nPlane+1)/2) * opts.zStepUm;
        end
    end

    if ~isfield(data, 'xyPixelSizeUm') || ~isfinitePositiveScalar(data.xyPixelSizeUm)
        if isfinitePositiveScalar(opts.xyPixelSizeUm)
            data.xyPixelSizeUm = double(opts.xyPixelSizeUm);
        else
            data.xyPixelSizeUm = 1;
        end
    end
    if ~isfield(data, 'xUm') || numel(data.xUm) ~= nx
        data.xUm = centeredAxis(nx, data.xyPixelSizeUm);
    end
    if ~isfield(data, 'yUm') || numel(data.yUm) ~= ny
        data.yUm = centeredAxis(ny, data.xyPixelSizeUm);
    end
    if ~isfield(data, 'head') || ~isstruct(data.head)
        data.head = struct();
    end
    if ~isfield(data, 'backgroundPerPixel') || isempty(data.backgroundPerPixel)
        data.backgroundPerPixel = zeros(1, 1, nCh, nPlane);
        data.backgroundSource = 'zero background';
        data.backgroundIndependent = false;
        data.backgroundExposureCalibrated = false;
    end
    if ~isfield(data, 'backgroundSource')
        data.backgroundSource = 'unknown';
    end
    if ~isfield(data, 'flatFieldGain') || numel(data.flatFieldGain) ~= nCh
        data.flatFieldGain = ones(1, nCh);
    end
    if ~isfield(data, 'flatFieldCalibrated')
        data.flatFieldCalibrated = false;
    end
end

function axisUm = centeredAxis(n, pixelSizeUm)
    axisUm = ((1:n) - (n+1)/2) * pixelSizeUm;
end

function [focusIdx, selection] = chooseFocusPlane(data, opts)
    signal = focalSignalTrace(data, opts.focusMetric);
    nPlane = numel(signal);
    if ~isempty(opts.focusIndex)
        focusIdx = round(opts.focusIndex);
        if focusIdx < 1 || focusIdx > nPlane
            error('estimateFocalPlaneCenterPointISMPoisson:BadFocusIndex', ...
                'focusIndex must be between 1 and %d.', nPlane);
        end
    elseif ~isempty(opts.focusZUm)
        [~, focusIdx] = min(abs(data.stageZUm - double(opts.focusZUm)));
    else
        [~, focusIdx] = max(signal);
    end

    selection = struct();
    selection.focusIndex = focusIdx;
    selection.focusZUm = data.stageZUm(focusIdx);
    selection.stageZUm = data.stageZUm(:).';
    selection.signalTrace = signal(:).';
    selection.focusMetric = opts.focusMetric;
end

function signal = focalSignalTrace(data, metric)
    [ny, nx, nCh, nPlane] = size(data.rawCounts);
    signal = zeros(1, nPlane);
    for iz = 1:nPlane
        bg = reshape(data.backgroundPerPixel(:,:,:,iz), 1, 1, nCh);
        plane = max(data.rawCounts(:,:,:,iz) - repmat(bg, ny, nx, 1), 0);
        switch metric
            case {'signal','sum','total'}
                signal(iz) = sum(plane(:));
            case {'peak','max'}
                summed = sum(plane, 3);
                signal(iz) = max(summed(:));
            otherwise
                error('estimateFocalPlaneCenterPointISMPoisson:BadFocusMetric', ...
                    'focusMetric must be signal or peak.');
        end
    end
end

function [focusImage, rawFocus, background] = focalPlaneForCentering(data, focusIdx)
%FOCALPLANEFORCENTERING Focal detector-summed image, raw detector plane, and background.
    [ny, nx, nCh, ~] = size(data.rawCounts);
    rawFocus = double(data.rawCounts(:,:,:,focusIdx));
    bg = reshape(data.backgroundPerPixel(:,:,:,focusIdx), 1, 1, nCh);
    signalFocus = max(rawFocus - repmat(bg, ny, nx, 1), 0);
    focusImage = sum(signalFocus, 3);
    background = reshape(bg, nCh, 1);
end

function [fit, centerXY, obs, iterInfo] = fitWithIterativeCenter( ...
        rawFocus, background, include, channelIDs, centerXY0, sim, opts, p0)
%FITWITHITERATIVECENTER Coordinate descent between the aberration fit and the
%   data sampling center. Each pass fits aberrations at the current center,
%   then relocates the center to the (sub-pixel) sampling position whose data
%   microimage best matches the on-axis model microimage (minimum Poisson
%   deviance). Both steps reduce the same objective, so the loop is stable and
%   the converged center is the emitter's on-axis position rather than the
%   aberration-shifted summed-image peak.

    centerXY = double(centerXY0(:)).';
    if opts.iterativeCenter
        maxIt = max(1, round(opts.maxCenterIterations));
    else
        maxIt = 1;
    end

    history = nan(maxIt, 4);   % [x_px, y_px, deviance, move_px]
    fit = struct();
    obs = struct();
    it = 0;
    for it = 1:maxIt
        counts = sampleDetectorPlaneAtXY(rawFocus, centerXY, opts.centerSampleMode);
        obs = struct('counts', counts(:), 'background', background(:), ...
            'include', include(:), 'channelIDs', channelIDs(:), ...
            'minExpectedCount', opts.minExpectedCount);
        fit = fitFocalPoisson(obs, sim, opts, p0);
        history(it, 1:3) = [centerXY, fit.deviance];

        if ~opts.iterativeCenter || it == maxIt
            history(it, 4) = 0;
            break;
        end

        newCenter = recenterByModelMatch(rawFocus, fit, obs, opts, centerXY);
        move = norm(newCenter - centerXY);
        history(it, 4) = move;
        if move < opts.centerTolerancePx
            break;
        end
        centerXY = newCenter;
    end

    converged = (~opts.iterativeCenter) || (it < maxIt);
    iterInfo = struct( ...
        'enabled', logical(opts.iterativeCenter), ...
        'nIterations', it, ...
        'converged', converged, ...
        'initialCenterXY', double(centerXY0(:)).', ...
        'finalCenterXY', centerXY, ...
        'totalShiftPx', norm(centerXY - double(centerXY0(:)).'), ...
        'history', history(1:it, :));
end

function newCenter = recenterByModelMatch(rawFocus, fit, obs, opts, centerNow)
%RECENTERBYMODELMATCH Find the data sampling center best matching the fixed
%   on-axis model microimage (coarse grid then local refinement).
    include = obs.include(:);
    m = fit.modelProbability(:);
    if numel(m) ~= numel(include)
        newCenter = centerNow;   % model/observation mismatch; leave center as is
        return;
    end
    m = m(include);
    b = max(fit.backgroundScale * obs.background(include), 0);

    coarse = searchCenter(rawFocus, centerNow, opts.centerSearchRadiusPx, ...
        opts.centerSearchStepPx, m, b, include, opts);
    newCenter = searchCenter(rawFocus, coarse, opts.centerSearchStepPx, ...
        opts.centerRefineStepPx, m, b, include, opts);
end

function best = searchCenter(rawFocus, c0, radius, step, m, b, include, opts)
    if ~(step > 0), step = max(radius, eps); end
    offsets = -radius:step:radius;
    if isempty(offsets), offsets = 0; end
    best = c0;
    bestDeviance = inf;
    bestDist = inf;
    for dy = offsets
        for dx = offsets
            cand = c0 + [dx, dy];
            counts = sampleDetectorPlaneAtXY(rawFocus, cand, opts.centerSampleMode);
            y = counts(include);
            scale = profileSignalScale(y, m, b);
            mu = max(scale * m + b, opts.minExpectedCount);
            deviance = sum(poissonDevianceResidual(y, mu).^2);
            dist = hypot(dx, dy);
            % Prefer lower deviance; break near-ties (e.g. same rounded pixel
            % in 'nearest' mode) toward the smallest move to avoid drift.
            if deviance < bestDeviance - 1e-9 || ...
                    (abs(deviance - bestDeviance) <= 1e-9 && dist < bestDist)
                bestDeviance = deviance;
                bestDist = dist;
                best = cand;
            end
        end
    end
end

function info = buildFocalCenterInfo(focusImage, centerXY, initialCenterXY, ...
        centerFitInfo, data, opts, iterInfo)
    ix = min(max(round(centerXY(1)), 1), size(focusImage, 2));
    iy = min(max(round(centerXY(2)), 1), size(focusImage, 1));
    info = struct();
    info.focusImage = focusImage;
    info.centerXY = centerXY(:).';
    info.initialCenterXY = double(initialCenterXY(:)).';
    info.centerPixelYX = [iy, ix];
    info.centerOffsetFromRoundedXY = centerXY(:).' - [ix, iy];
    info.centerMode = opts.centerMode;
    info.centerSampleMode = opts.centerSampleMode;
    info.centerFit = centerFitInfo;
    info.iterativeCenter = iterInfo;
    info.backgroundSource = data.backgroundSource;
end

function [centerXY, centerYX, fitInfo] = estimateCenterXY(img, opts)
    img = double(img);
    fitInfo = struct('method', 'explicit centerXY option');
    if ~isempty(opts.centerXY)
        centerXY = double(opts.centerXY(:)).';
        if numel(centerXY) ~= 2
            error('estimateFocalPlaneCenterPointISMPoisson:BadCenterXY', ...
                'centerXY must be [x y].');
        end
    else
        switch opts.centerMode
            case 'peak'
                [~, idx] = max(img(:));
                [cy, cx] = ind2sub(size(img), idx);
                centerXY = [cx, cy];
                fitInfo.method = 'summed focal image peak';
            case {'centroid','weighted'}
                [centerXY, fitInfo] = centerOfMassXY(img, opts.centerThresholdFraction);
            case {'gaussian','gaussiancom','gaussian_com','comgaussian'}
                [centerXY, fitInfo] = gaussianCenterFromCom(img, opts.centerThresholdFraction);
            otherwise
                error('estimateFocalPlaneCenterPointISMPoisson:BadCenterMode', ...
                    'centerMode must be gaussian, centroid, or peak.');
        end
    end
    ix = min(max(round(centerXY(1)), 1), size(img, 2));
    iy = min(max(round(centerXY(2)), 1), size(img, 1));
    centerYX = [iy, ix];
end

function [centerXY, info] = centerOfMassXY(img, thresholdFraction)
    positive = max(double(img) - min(img(:)), 0);
    threshold = thresholdFraction * max(positive(:));
    weights = positive;
    weights(weights < threshold) = 0;
    if sum(weights(:)) <= 0
        [~, idx] = max(img(:));
        [cy, cx] = ind2sub(size(img), idx);
        centerXY = [cx, cy];
        info = struct('method', 'summed focal image peak fallback');
    else
        [yy, xx] = ndgrid(1:size(img,1), 1:size(img,2));
        mass = sum(weights(:));
        centerXY = [sum(xx(:).*weights(:)), sum(yy(:).*weights(:))] / mass;
        info = struct('method', 'summed focal image center of mass', ...
            'mass', mass, 'threshold', threshold);
    end
end

function [centerXY, info] = gaussianCenterFromCom(img, thresholdFraction)
    [seedXY, comInfo] = centerOfMassXY(img, thresholdFraction);
    img = double(img);
    [ny, nx] = size(img);
    cx0 = min(max(seedXY(1), 1), nx);
    cy0 = min(max(seedXY(2), 1), ny);
    radius = max(4, ceil(min([ny nx]) / 6));
    xlo = max(1, floor(cx0 - radius));
    xhi = min(nx, ceil(cx0 + radius));
    ylo = max(1, floor(cy0 - radius));
    yhi = min(ny, ceil(cy0 + radius));
    patch = img(ylo:yhi, xlo:xhi);
    [yy, xx] = ndgrid(ylo:yhi, xlo:xhi);

    bg0 = median(patch(:));
    amp0 = max(patch(:)) - bg0;
    if ~isfinite(amp0) || amp0 <= 0
        centerXY = seedXY;
        info = comInfo;
        info.method = 'summed focal image center of mass fallback';
        return;
    end

    sigma0 = max(1.5, radius / 2);
    p0 = [bg0, amp0, cx0, cy0, log(sigma0), log(sigma0)];
    objective = @(p) gaussianSse(p, xx, yy, patch);
    options = optimset('Display', 'off', 'MaxIter', 200, 'MaxFunEvals', 800);
    p = fminsearch(objective, p0, options);
    centerXY = [min(max(p(3), 1), nx), min(max(p(4), 1), ny)];
    info = struct('method', 'summed focal image Gaussian fit from COM seed', ...
        'centerOfMassXY', seedXY, 'fitSigmaXY', exp(p(5:6)), ...
        'fitAmplitude', p(2), 'fitBackground', p(1), ...
        'fitWindowYX', [ylo yhi xlo xhi]);
end

function sse = gaussianSse(p, xx, yy, patch)
    sx = max(exp(p(5)), 0.5);
    sy = max(exp(p(6)), 0.5);
    model = p(1) + p(2) * exp(-0.5 * (((xx - p(3)) / sx).^2 + ...
        ((yy - p(4)) / sy).^2));
    residual = model - patch;
    sse = sum(residual(:).^2);
    if ~isfinite(sse)
        sse = realmax;
    end
end

function values = sampleDetectorPlaneAtXY(plane, centerXY, mode)
    nCh = size(plane, 3);
    if size(plane, 1) == 1 && size(plane, 2) == 1
        values = reshape(double(plane(1, 1, :)), nCh, 1);
        return;
    end

    x = min(max(double(centerXY(1)), 1), size(plane, 2));
    y = min(max(double(centerXY(2)), 1), size(plane, 1));
    ix = min(max(round(x), 1), size(plane, 2));
    iy = min(max(round(y), 1), size(plane, 1));

    switch lower(char(mode))
        case {'nearest','round','pixel'}
            values = reshape(double(plane(iy, ix, :)), nCh, 1);
        case {'subpixel','linear','bilinear','interp'}
            values = zeros(nCh, 1);
            for ch = 1:nCh
                img = double(plane(:,:,ch));
                v = interp2(img, x, y, 'linear', NaN);
                if ~isfinite(v)
                    v = img(iy, ix);
                end
                values(ch) = v;
            end
        otherwise
            error('estimateFocalPlaneCenterPointISMPoisson:BadCenterSampleMode', ...
                'centerSampleMode must be nearest or subpixel.');
    end
end

function [includeCh, info] = selectedChannels(channelIDs, backgroundCounts, opts)
    includeCh = true(numel(channelIDs), 1);
    if ~isempty(opts.excludeChannelIDs)
        includeCh = includeCh & ~ismember(channelIDs(:), opts.excludeChannelIDs(:));
    end
    if opts.excludeHotDarkPixels
        bg = double(backgroundCounts(:));
        medBg = median(bg(isfinite(bg)));
        if isempty(medBg) || ~isfinite(medBg)
            medBg = 0;
        end
        hot = bg > opts.hotPixelFactor * max(medBg, eps);
        includeCh = includeCh & ~hot;
    end
    if nnz(includeCh) < max(6, ceil(numel(includeCh)/3))
        warning('estimateFocalPlaneCenterPointISMPoisson:FewChannels', ...
            'Too few channels remained after exclusions; using all channels.');
        includeCh(:) = true;
    end
    info = struct();
    info.includeMask = includeCh;
    info.excludeChannelIDs = opts.excludeChannelIDs;
    info.excludeHotDarkPixels = opts.excludeHotDarkPixels;
    info.hotPixelFactor = opts.hotPixelFactor;
end

function [sim, diagnostics] = configureFocalSimulation(data, focusIdx, opts)
    diagnostics = [];
    if ~isempty(opts.sim)
        sim = opts.sim;
    else
        sim = defaultParams();
        sim.detectorLayout = opts.detectorLayout;
        sim.detectorPixelShape = opts.detectorPixelShape;
        [sim.detXY, sim.detectorIndexGrid, sim.detectorLayoutInfo] = ...
            detectorLayout(sim.detectorLayout, sim.detPitch);
        sim.nDet = size(sim.detXY, 1);
        sim.detectorGridSize = size(sim.detectorIndexGrid);
        sim.arrayN = sim.detectorGridSize(1);
        if strcmpi(sim.detectorPixelShape, 'hex')
            sim.detectorHexRadius = sim.detSize / sqrt(3);
        end
    end

    sim = applyOpticsOptions(sim, data.head, opts);
    sim.sampleGeometry = opts.sampleGeometry;
    sim.interfaceRadialWeightMode = opts.interfaceRadialWeightMode;
    sim.airInterfaceStageMedium = opts.airInterfaceStageMedium;
    sim.diffractionModel = opts.diffractionModel;
    sim.vectorialPolarizationMode = opts.vectorialPolarizationMode;
    sim.excitationPolarizationMode = opts.excitationPolarizationMode;
    sim.collectionPolarizationMode = opts.collectionPolarizationMode;
    sim.includesVectorialPolarization = usesVectorialPSF(sim);
    sim.nImmersion = opts.immersionRefractiveIndex;
    sim.nGlass = opts.glassRefractiveIndex;
    sim.nSample = opts.sampleRefractiveIndex;
    sim.nDesignGlass = opts.designGlassRefractiveIndex;
    sim.coverslipThicknessUm = opts.coverslipThicknessUm;
    sim.designCoverslipThicknessUm = opts.designCoverslipThicknessUm;
    sim.beadBottomHeightUm = opts.beadBottomHeightUm;
    sim.airBeadAxialSamples = opts.airBeadAxialSamples;
    sim.nMedium = sim.nImmersion;
    sim.beadRadius = opts.beadRadiusUm;
    sim.beadSubsamples = opts.beadSubsamples;
    sim.x = data.xUm(:).';
    sim.y = data.yUm(:).';
    sim.nx = numel(sim.x);
    sim.ny = numel(sim.y);
    sim.dx = data.xyPixelSizeUm;
    sim.xyPixelSizeUm = data.xyPixelSizeUm;
    sim.fovX = max(sim.x) - min(sim.x);
    sim.fovY = max(sim.y) - min(sim.y);
    sim.fovXY = max(sim.fovX, sim.fovY);
    if ~isempty(opts.detectorSubsamples)
        sim.detectorSubsamples = opts.detectorSubsamples;
    end
    validateScalarOptics(sim);

    airyUnitUm = 1.22 * sim.lamEm / effectivePropagatingNA(sim);
    if ~isempty(opts.detectorXYUm)
        detXY = validateDetectorXY(opts.detectorXYUm, size(data.rawCounts, 3));
        diagnostics = detectorDiagnostics(detXY, data.xyPixelSizeUm, airyUnitUm, ...
            'explicit detectorXYUm option');
    elseif opts.estimateDetectorLayout
        layoutOpts = struct();
        layoutOpts.positionSign = opts.detectorLayoutPositionSign;
        layoutOpts.scanAxisSigns = opts.scanAxisSigns;
        layoutOpts.detectorScale = opts.detectorLayoutScale;
        layoutOpts.centerMode = opts.detectorLayoutCenterMode;
        layoutOpts.centerDetectorIndex = opts.detectorCenterIndex;
        layoutOpts.smoothSigma = opts.detectorShiftSmoothSigma;
        layoutOpts.useWindow = opts.detectorShiftUseWindow;
        layoutOpts.normalizeImages = opts.detectorShiftNormalizeImages;
        layoutOpts.upsampleReg = opts.detectorShiftUpsample;
        layoutOpts.airyUnitUm = airyUnitUm;
        stackForLayout = backgroundSubtractedPlane(data, focusIdx);
        [detXY, diagnostics] = estimateDetectorLayoutFromStack( ...
            stackForLayout, data.xyPixelSizeUm, layoutOpts);
        diagnostics.source = 'background-subtracted focal-plane ISM phase correlation';
    else
        detXY = sim.detXY;
        diagnostics = detectorDiagnostics(detXY, data.xyPixelSizeUm, airyUnitUm, ...
            sprintf('fixed regular %s detector layout', sim.detectorLayout));
    end

    calibratedPitch = resolveDetectorPitchSampleUm(opts);
    recoveredPitch = medianNearestDistance(detXY);
    if isfinitePositiveScalar(calibratedPitch) && ...
            isfinitePositiveScalar(recoveredPitch)
        detXY = detXY * (calibratedPitch / recoveredPitch);
        diagnostics.detXYBeforePitchCalibrationUm = diagnostics.detXYUm;
        diagnostics.detXYUm = detXY;
        diagnostics.detXYNm = 1000 * detXY;
        diagnostics.pitchCalibrationScale = calibratedPitch / recoveredPitch;
        diagnostics.calibratedPitchUm = calibratedPitch;
    end
    diagnostics.detectorPitchSource = detectorPitchSource(opts);
    sim = applyDetectorXY(sim, detXY);
    sim.detectorPhysicalPitchCalibrated = isfinitePositiveScalar(calibratedPitch) || ...
        ~isempty(opts.detectorXYUm);
    sim.detectorCalibratedPitchUm = calibratedPitch;

    sim = prepareFocalZGrid(sim, opts);
    if ~strcmpi(sim.sampleGeometry, 'airOnGlass')
        sim.obj = beadObject3D(sim);
    end
end

function sim = applyOpticsOptions(sim, head, opts)
    if isfinitePositiveScalar(opts.objectiveNA)
        sim.NA = double(opts.objectiveNA);
    end
    if isfinitePositiveScalar(opts.objectiveMagnification)
        sim.objectiveMagnification = double(opts.objectiveMagnification);
    end
    if isfinitePositiveScalar(opts.excitationWavelengthUm)
        sim.lamExc = double(opts.excitationWavelengthUm);
    end
    if isfinitePositiveScalar(opts.emissionWavelengthUm)
        sim.lamEm = double(opts.emissionWavelengthUm);
        sim.lamRef = double(opts.emissionWavelengthUm);
    end
    if isfinitePositiveScalar(opts.mediumRefractiveIndex)
        sim.nMedium = double(opts.mediumRefractiveIndex);
        sim.nImmersion = double(opts.mediumRefractiveIndex);
    end
    if ~isstruct(head)
        return;
    end
end

function plane = backgroundSubtractedPlane(data, focusIdx)
    [ny, nx, nCh, ~] = size(data.rawCounts);
    bg = reshape(data.backgroundPerPixel(:,:,:,focusIdx), 1, 1, nCh);
    plane = max(double(data.rawCounts(:,:,:,focusIdx)) - repmat(bg, ny, nx, 1), 0);
end

function sim = prepareFocalZGrid(sim, opts)
    pad = max(opts.modelZPaddingUm, sim.beadRadius);
    zMin = -pad - sim.beadRadius;
    zMax = pad + sim.beadRadius;
    nZ = ceil((zMax - zMin) / opts.modelDzUm) + 1;
    if mod(nZ, 2) == 0
        nZ = nZ + 1;
    end
    sim.z = linspace(zMin, zMax, nZ);
    sim.nz = numel(sim.z);
    sim.nzRange = zMax - zMin;
end

function sim = applyDetectorXY(sim, detXY)
    oldPitch = numericField(sim, 'detPitch', NaN);
    oldSize = numericField(sim, 'detSize', NaN);
    pitch = medianNearestDistance(detXY);
    sim.detXY = detXY;
    sim.nDet = size(detXY, 1);
    if isfinitePositiveScalar(pitch)
        fillRatio = 1;
        if isfinitePositiveScalar(oldPitch) && isfinitePositiveScalar(oldSize)
            fillRatio = oldSize / oldPitch;
        end
        sim.detPitch = pitch;
        sim.detSize = fillRatio * pitch;
        if isfield(sim, 'detectorPixelShape') && strcmpi(sim.detectorPixelShape, 'hex')
            sim.detectorHexRadius = sim.detSize / sqrt(3);
        end
    end
end

function p = initialParameterVector(sim, opts)
    p = zeros(1, numel(opts.fitModes));
    coeffs0 = opts.initialCoeffs;
    if isnumeric(coeffs0)
        coeffs0 = coeffStruct(sim, coeffs0);
    end
    for k = 1:numel(opts.fitModes)
        if isstruct(coeffs0) && isfield(coeffs0, opts.fitModes{k})
            p(k) = coeffs0.(opts.fitModes{k});
        end
    end
    if opts.fitXY
        p = [p double(opts.initialXY(:).')];
    end
    if opts.fitDetectorPitchScale
        p = [p double(opts.initialDetectorPitchScale)];
    end
    if opts.fitBackgroundScale
        p = [p double(opts.backgroundScale)];
    end
    p = clampParameters(p, opts);
    if ~isfield(sim, 'modeOrder')
        error('estimateFocalPlaneCenterPointISMPoisson:BadSimulation', ...
            'Simulation structure must contain modeOrder.');
    end
end

function fit = fitFocalPoisson(obs, sim, opts, p0)
    p = clampParameters(p0, opts);
    nParam = numel(p);
    steps = finiteDifferenceSteps(nParam, sim, opts);
    reg = regularizationVector(nParam, opts);
    maxStep = maxStepVector(nParam, opts);
    history = nan(opts.maxIter, 4);
    converged = false;
    terminationReason = 'maximum iterations reached';
    nForwardEvaluations = 0;

    [r, deviance, scale, mu, modelProb, backgroundScale] = ...
        evaluateParameters(p, obs, sim, opts);
    nForwardEvaluations = nForwardEvaluations + 1;

    for it = 1:opts.maxIter
        J = zeros(numel(r), nParam);
        for q = 1:nParam
            pp = p;
            pp(q) = pp(q) + steps(q);
            pp = clampParameters(pp, opts);
            dq = pp(q) - p(q);
            if abs(dq) < eps
                continue;
            end
            rp = evaluateParameters(pp, obs, sim, opts);
            nForwardEvaluations = nForwardEvaluations + 1;
            if strcmp(opts.jacobianScheme, 'central')
                pm = p;
                pm(q) = pm(q) - steps(q);
                pm = clampParameters(pm, opts);
                dm = p(q) - pm(q);
                if abs(dm) >= eps
                    rm = evaluateParameters(pm, obs, sim, opts);
                    nForwardEvaluations = nForwardEvaluations + 1;
                    J(:,q) = (rp - rm) / (dq + dm);
                else
                    J(:,q) = (rp - r) / dq;
                end
            else
                J(:,q) = (rp - r) / dq;
            end
        end

        H = J.'*J + diag(reg);
        g = -J.'*r;
        if rcond(H) < 1e-12
            delta = pinv(H) * g;
        else
            delta = H \ g;
        end
        delta = max(min(delta(:).', maxStep), -maxStep);
        if norm(delta) < opts.tolStep
            converged = true;
            terminationReason = 'Gauss-Newton step below tolerance';
            history = history(1:max(it-1, 0), :);
            break;
        end

        accepted = false;
        lineScale = 1;
        for ls = 1:10
            trial = clampParameters(p + lineScale * delta, opts);
            [rt, dt, st, mut, modelt, bgt] = evaluateParameters(trial, obs, sim, opts);
            nForwardEvaluations = nForwardEvaluations + 1;
            if dt < deviance
                p = trial;
                r = rt;
                deviance = dt;
                scale = st;
                mu = mut;
                modelProb = modelt;
                backgroundScale = bgt;
                accepted = true;
                break;
            end
            lineScale = lineScale / 2;
        end

        history(it,:) = [deviance, norm(lineScale * delta), scale, accepted];
        if opts.verbose
            fprintf(['[focal-center Poisson] iter %02d: deviance %.5g, ' ...
                'step %.3g, photons %.4g, bgScale %.3f\n'], ...
                it, deviance, norm(lineScale * delta), scale, backgroundScale);
        end
        if ~accepted
            terminationReason = 'line search failed';
            history = history(1:it, :);
            break;
        end
        if norm(lineScale * delta) < opts.tolStep
            converged = true;
            terminationReason = 'accepted step below tolerance';
            history = history(1:it, :);
            break;
        end
        if it == opts.maxIter
            history = history(1:it, :);
        end
    end

    [coeffs, xy, detectorPitchScale, backgroundScale] = unpackParams(sim, opts, p);
    simFit = applyDetectorPitchScaleToSim(sim, detectorPitchScale);
    fit = struct();
    fit.objective = 'poisson';
    fit.paramVector = p;
    fit.paramNames = buildParamNames(opts);
    fit.fitModes = opts.fitModes;
    fit.estCoeffs = coeffs;
    fit.estCoeffVector = coeffStructToVector(sim, coeffs);
    fit.estCoeffNm = fit.estCoeffVector * sim.lamRef * 1000;
    fit.estXY = xy;
    fit.estDetectorPitchScale = detectorPitchScale;
    fit.estDetectorPitchUm = numericField(simFit, 'detPitch', NaN);
    fit.backgroundScale = backgroundScale;
    fit.photonScale = scale;
    fit.expectedCounts = mu(:);
    fit.modelProbability = modelProb(:);
    fit.devianceResidual = r(:);
    fit.deviance = deviance;
    fit.reducedDeviance = deviance / max(nnz(obs.include) - nParam, 1);
    fit.residualNorm = sqrt(deviance);
    fit.history = history;
    fit.converged = converged;
    fit.terminationReason = terminationReason;
    fit.nForwardEvaluations = nForwardEvaluations;
    fit.sim = simFit;
    fit.nominalSim = sim;
end

function [residual, deviance, scale, mu, modelProb, backgroundScale] = ...
        evaluateParameters(p, obs, sim, opts)
    modelProb = modelFocalProbability(sim, opts, p);
    [~, ~, ~, backgroundScale] = unpackParams(sim, opts, p);
    idx = obs.include(:);
    y = double(obs.counts(idx));
    m = double(modelProb(idx));
    b = max(backgroundScale * double(obs.background(idx)), 0);
    scale = profileSignalScale(y, m, b);
    mu = max(scale * m + b, obs.minExpectedCount);
    residual = poissonDevianceResidual(y, mu);
    deviance = sum(residual.^2);
end

function modelProb = modelFocalProbability(sim, opts, p)
    [coeffs, xy, detectorPitchScale] = unpackParams(sim, opts, p);
    simEval = applyDetectorPitchScaleToSim(sim, detectorPitchScale);
    if strcmpi(simEval.sampleGeometry, 'airOnGlass')
        stack = normalizedStackAirInterfaceZPlanes(simEval, coeffs, 0, xy(1), xy(2), 0);
    else
        stack = normalizedStackExplicitDetectorZPlanes(simEval, coeffs, 0, xy(1), xy(2), 0);
    end
    values = sampleModelStackAtXY(stack, simEval, opts.modelSampleXY, ...
        opts.centerSampleMode);
    if ~isempty(opts.detectorQE)
        values = values .* reshape(opts.detectorQE, [], 1);
    end
    values = max(values(:), 0);
    total = sum(values);
    if total <= 0 || ~isfinite(total)
        modelProb = ones(simEval.nDet, 1) / simEval.nDet;
    else
        modelProb = values / total;
    end
end

function values = sampleModelStackAtXY(stack, sim, sampleXY, mode)
    x = double(sampleXY(1));
    y = double(sampleXY(2));
    x = min(max(x, min(sim.x(:))), max(sim.x(:)));
    y = min(max(y, min(sim.y(:))), max(sim.y(:)));
    switch lower(char(mode))
        case {'nearest','round','pixel'}
            [~, ix] = min(abs(sim.x - x));
            [~, iy] = min(abs(sim.y - y));
            values = reshape(double(stack(iy, ix, :, 1)), sim.nDet, 1);
        case {'subpixel','linear','bilinear','interp'}
            values = zeros(sim.nDet, 1);
            for ch = 1:sim.nDet
                img = double(stack(:,:,ch,1));
                values(ch) = interp2(sim.x, sim.y, img, x, y, 'linear', 0);
            end
        otherwise
            error('estimateFocalPlaneCenterPointISMPoisson:BadCenterSampleMode', ...
                'centerSampleMode must be nearest or subpixel.');
    end
end

function scale = profileSignalScale(y, m, b)
    m = max(double(m(:)), 0);
    y = double(y(:));
    b = max(double(b(:)), 0);
    if sum(m) <= 0
        scale = 0;
        return;
    end
    signal = max(sum(y - b), 0);
    guess = max(signal / max(sum(m), eps), 1);
    derivative = @(a) sum(m .* (1 - y ./ max(a * m + b, realmin)));
    lo = 0;
    hi = guess;
    while derivative(hi) < 0 && hi < 1e18
        hi = 2 * hi;
    end
    for k = 1:60
        mid = (lo + hi) / 2;
        if derivative(mid) > 0
            hi = mid;
        else
            lo = mid;
        end
    end
    scale = (lo + hi) / 2;
end

function r = poissonDevianceResidual(y, mu)
    y = double(y(:));
    mu = max(double(mu(:)), realmin);
    term = mu - y;
    positive = y > 0;
    term(positive) = y(positive) .* log(y(positive) ./ mu(positive)) - ...
        (y(positive) - mu(positive));
    term = max(term, 0);
    r = sign(y - mu) .* sqrt(2 * term);
end

function fit = attachFitDiagnostics(fit, obs, sim, opts)
    p = fit.paramVector;
    r0 = fit.devianceResidual(:);
    nParam = numel(p);
    steps = finiteDifferenceSteps(nParam, sim, opts);
    J = zeros(numel(r0), nParam);
    for q = 1:nParam
        pp = clampParameters(p + unitStep(q, nParam) * steps(q), opts);
        dq = pp(q) - p(q);
        if abs(dq) < eps
            continue;
        end
        rp = evaluateParameters(pp, obs, sim, opts);
        J(:,q) = (rp - r0) / dq;
    end
    s = svd(J, 'econ');
    tol = max(size(J)) * eps(max([s(:); 1]));
    rankJ = sum(s > tol);
    fisher = J.' * J;
    covar = pinv(fisher);
    stdParam = sqrt(max(diag(covar), 0));
    fit.jacobian = J;
    fit.fisher = fisher;
    fit.paramStd = stdParam(:).';
    fit.coeffStdWaves = stdParam(1:numel(opts.fitModes)).';
    fit.coeffStdNm = fit.coeffStdWaves * sim.lamRef * 1000;
    fit.identifiability = struct('nObservations', nnz(obs.include), ...
        'nParameters', nParam, 'rank', rankJ, ...
        'isFullRank', rankJ == nParam, 'singularValues', s, ...
        'tolerance', tol, 'conditionNumber', conditionFromSingularValues(s, rankJ));
end

function e = unitStep(q, n)
    e = zeros(1, n);
    e(q) = 1;
end

function value = conditionFromSingularValues(s, rankJ)
    if rankJ <= 0
        value = Inf;
    else
        value = s(1) / s(rankJ);
    end
end

function [coeffs, xy, detectorPitchScale, backgroundScale] = unpackParams(sim, opts, p)
    fullVec = zeros(1, numel(sim.modeOrder));
    for k = 1:numel(opts.fitModes)
        idx = find(strcmp(sim.modeOrder, opts.fitModes{k}), 1, 'first');
        if isempty(idx)
            error('estimateFocalPlaneCenterPointISMPoisson:UnknownMode', ...
                'Unknown fit mode "%s".', opts.fitModes{k});
        end
        fullVec(idx) = p(k);
    end
    coeffs = coeffStruct(sim, fullVec);
    next = numel(opts.fitModes) + 1;
    xy = [0 0];
    if opts.fitXY
        xy = [p(next) p(next+1)];
        next = next + 2;
    end
    detectorPitchScale = 1;
    if opts.fitDetectorPitchScale
        detectorPitchScale = p(next);
        next = next + 1;
    end
    backgroundScale = opts.backgroundScale;
    if opts.fitBackgroundScale
        backgroundScale = p(next);
    end
end

function p = clampParameters(p, opts)
    nMode = numel(opts.fitModes);
    b = opts.coefficientBoundsWaves;
    p(1:nMode) = min(max(p(1:nMode), b(1)), b(2));
    next = nMode + 1;
    if opts.fitXY
        p(next:next+1) = min(max(p(next:next+1), -0.30), 0.30);
        next = next + 2;
    end
    if opts.fitDetectorPitchScale
        b = opts.detectorPitchScaleBounds;
        p(next) = min(max(p(next), b(1)), b(2));
        next = next + 1;
    end
    if opts.fitBackgroundScale
        b = opts.backgroundScaleBounds;
        p(next) = min(max(p(next), b(1)), b(2));
    end
end

function step = finiteDifferenceSteps(nParam, sim, opts)
    step = opts.fdCoeff * ones(1, numel(opts.fitModes));
    if opts.fitXY
        fdXY = opts.fdXY;
        if isempty(fdXY) || ~isfinitePositiveScalar(fdXY)
            fdXY = max(sim.dx / 4, 0.005);
        end
        step = [step fdXY fdXY];
    end
    if opts.fitDetectorPitchScale
        step = [step opts.fdDetectorPitchScale];
    end
    if opts.fitBackgroundScale
        step = [step opts.fdBackgroundScale];
    end
    step = step(1:nParam);
    step(step <= 0 | ~isfinite(step)) = opts.fdCoeff;
end

function reg = regularizationVector(nParam, opts)
    reg = opts.regCoeff * ones(1, numel(opts.fitModes));
    if opts.fitXY
        reg = [reg opts.regXY opts.regXY];
    end
    if opts.fitDetectorPitchScale
        reg = [reg opts.regDetectorPitchScale];
    end
    if opts.fitBackgroundScale
        reg = [reg opts.regBackgroundScale];
    end
    reg = reg(1:nParam);
end

function values = maxStepVector(nParam, opts)
    values = opts.maxCoeffStep * ones(1, numel(opts.fitModes));
    if opts.fitXY
        values = [values opts.maxXYStep opts.maxXYStep];
    end
    if opts.fitDetectorPitchScale
        values = [values opts.maxDetectorPitchScaleStep];
    end
    if opts.fitBackgroundScale
        values = [values opts.maxBackgroundScaleStep];
    end
    values = values(1:nParam);
end

function names = buildParamNames(opts)
    names = opts.fitModes;
    if opts.fitXY
        names = [names {'x_shift','y_shift'}];
    end
    if opts.fitDetectorPitchScale
        names = [names {'detector_pitch_scale'}];
    end
    if opts.fitBackgroundScale
        names = [names {'background_scale'}];
    end
end

function sim = applyDetectorPitchScaleToSim(sim, scale)
    if isempty(scale) || ~isfinitePositiveScalar(scale) || abs(scale - 1) < eps
        return;
    end
    if isfield(sim, 'detXY') && ~isempty(sim.detXY)
        sim.detXY = double(sim.detXY) * double(scale);
    end
    if isfield(sim, 'detPitch') && isfinitePositiveScalar(sim.detPitch)
        sim.detPitch = double(sim.detPitch) * double(scale);
    end
    if isfield(sim, 'detSize') && isfinitePositiveScalar(sim.detSize)
        sim.detSize = double(sim.detSize) * double(scale);
    elseif isfield(sim, 'detPitch') && isfinitePositiveScalar(sim.detPitch)
        sim.detSize = sim.detPitch;
    end
    if isfield(sim, 'detectorPixelShape') && strcmpi(sim.detectorPixelShape, 'hex') && ...
            isfield(sim, 'detSize') && isfinitePositiveScalar(sim.detSize)
        sim.detectorHexRadius = sim.detSize / sqrt(3);
    end
end

function vec = coeffStructToVector(sim, coeffs)
    vec = zeros(1, numel(sim.modeOrder));
    for k = 1:numel(sim.modeOrder)
        if isfield(coeffs, sim.modeOrder{k})
            vec(k) = coeffs.(sim.modeOrder{k});
        end
    end
end

function outDir = resolveOutputDir(opts, focusSelection, inputInfo)
    if ~isempty(opts.outputDir)
        outDir = opts.outputDir;
        return;
    end
    thisDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(thisDir);
    label = localFileStem(inputInfo.resolvedInput);
    if isempty(label)
        label = 'focal_stack';
    end
    label = sprintf('%s_focus%03d', label, focusSelection.focusIndex);
    outDir = fullfile(rootDir, 'output_matlab', ...
        'focal_center_point_ism_poisson', sanitizeFileName(label));
end

function writeOutputs(result)
    outDir = result.outputDir;
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    save(fullfile(outDir, 'focal_center_point_ism_poisson_fit.mat'), ...
        'result', '-v7.3');
    writetable(coefficientTable(result), ...
        fullfile(outDir, 'focal_center_point_poisson_coefficients.csv'));
    writetable(microimageTable(result), ...
        fullfile(outDir, 'focal_center_point_poisson_counts.csv'));
    writetable(focusSelectionTable(result), ...
        fullfile(outDir, 'focal_plane_selection.csv'));
    if ~isempty(result.detectorLayoutDiagnostics)
        writetable(detectorLayoutTable(result), ...
            fullfile(outDir, 'focal_center_detector_geometry.csv'));
    end
    writeSummaryFigure(result, ...
        fullfile(outDir, 'focal_center_point_poisson_summary.png'));
end

function T = coefficientTable(result)
    mode = result.fit.sim.modeOrder(:);
    coeffWavesRMS = result.fit.estCoeffVector(:);
    coeffNmRMS = coeffWavesRMS * result.fit.sim.lamRef * 1000;
    isFitted = ismember(mode, result.fit.fitModes(:));
    T = table(mode, coeffWavesRMS, coeffNmRMS, isFitted);
end

function T = microimageTable(result)
    nCh = numel(result.centerCounts);
    detectorIndex = (1:nCh).';
    channelID = result.data.channelIDs(:);
    observedCounts = result.centerCounts(:);
    backgroundCounts = result.backgroundCounts(:);
    expectedCounts = nan(nCh, 1);
    modelProbability = nan(nCh, 1);
    residual = nan(nCh, 1);
    includeMask = result.channelInfo.includeMask(:);
    expectedCounts(includeMask) = result.fit.expectedCounts(:);
    modelProbability(includeMask) = result.fit.modelProbability(includeMask);
    residual(includeMask) = result.fit.devianceResidual(:);
    T = table(detectorIndex, channelID, includeMask, observedCounts, ...
        backgroundCounts, modelProbability, expectedCounts, residual);
end

function T = focusSelectionTable(result)
    zUm = result.focusSelection.stageZUm(:);
    signal = result.focusSelection.signalTrace(:);
    isFocus = false(numel(zUm), 1);
    isFocus(result.focusIndex) = true;
    T = table((1:numel(zUm)).', zUm, signal, isFocus, ...
        'VariableNames', {'planeIndex','zUm','signal','isFocus'});
end

function T = detectorLayoutTable(result)
    d = result.detectorLayoutDiagnostics;
    nCh = size(d.detXYUm, 1);
    detectorIndex = (1:nCh).';
    channelID = result.data.channelIDs(:);
    detXUm = d.detXYUm(:,1);
    detYUm = d.detXYUm(:,2);
    detXNm = 1000 * detXUm;
    detYNm = 1000 * detYUm;
    detectorRadiusUm = hypot(detXUm, detYUm);
    source = repmat({d.source}, nCh, 1);
    T = table(detectorIndex, channelID, detXUm, detYUm, detXNm, detYNm, ...
        detectorRadiusUm, source);
end

function printSummary(result)
    fprintf('\n[estimateFocalPlaneCenterPointISMPoisson]\n');
    fprintf('  focus plane: %d at z = %.4f um\n', ...
        result.focusIndex, result.focusZUm);
    fprintf('  bead center: x = %.3f px, y = %.3f px (%s)\n', ...
        result.centerInfo.centerXY(1), result.centerInfo.centerXY(2), ...
        result.centerInfo.centerFit.method);
    ic = result.centerIteration;
    if isfield(ic, 'enabled') && ic.enabled
        if ic.converged, status = 'converged'; else, status = 'max iters'; end
        fprintf(['  iterative recenter: %d iters (%s), moved %.3f px from ' ...
            'Gaussian seed (%.3f, %.3f)\n'], ic.nIterations, status, ...
            ic.totalShiftPx, ic.initialCenterXY(1), ic.initialCenterXY(2));
    end
    fprintf('  background: %s; fitted bg scale %.3f\n', ...
        result.centerInfo.backgroundSource, result.fit.backgroundScale);
    fprintf('  Poisson deviance: %.4g (reduced %.3g), photons %.4g\n', ...
        result.fit.deviance, result.fit.reducedDeviance, result.fit.photonScale);
    fprintf('  identifiability rank: %d/%d, condition %.3g\n', ...
        result.fit.identifiability.rank, result.fit.identifiability.nParameters, ...
        result.fit.identifiability.conditionNumber);
    stdWaves = result.fit.coeffStdWaves(:);
    stdNm = result.fit.coeffStdNm(:);
    for k = 1:numel(result.options.fitModes)
        mode = result.options.fitModes{k};
        value = result.fit.estCoeffs.(mode);
        sw = uncertaintyAt(stdWaves, k);
        sn = uncertaintyAt(stdNm, k);
        fprintf('  %-12s %+8.4f +/- %6.4f waves RMS (%+7.1f +/- %5.1f nm)\n', ...
            mode, value, sw, value * result.fit.sim.lamRef * 1000, sn);
    end
    fprintf('  output: %s\n', result.outputDir);
end

function value = uncertaintyAt(vec, k)
    if numel(vec) >= k && isfinite(vec(k))
        value = vec(k);
    else
        value = NaN;
    end
end

function writeSummaryFigure(result, fileName)
    try
        fig = buildSummaryFigure(result, 'off');
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, fileName, '-dpng', '-r180');
        close(fig);
    catch err
        warning('estimateFocalPlaneCenterPointISMPoisson:FigureFailed', ...
            'Could not write summary figure: %s', err.message);
    end
end

function fig = buildSummaryFigure(result, visibleState)
%BUILDSUMMARYFIGURE Assemble the focal-center diagnostic figure.
%
%   Panels: focal summed image, observed / fitted / residual honeycomb hex
%   detector maps, the observed-vs-expected detector counts, and the fitted
%   aberration coefficients with +/-1 sigma uncertainty error bars.

    if nargin < 2 || isempty(visibleState)
        visibleState = 'on';
    end

    nCh = numel(result.centerCounts);
    includeMask = result.channelInfo.includeMask(:);
    expectedCounts = nan(nCh, 1);
    devianceResidual = nan(nCh, 1);
    expectedCounts(includeMask) = result.fit.expectedCounts(:);
    devianceResidual(includeMask) = result.fit.devianceResidual(:);
    fittedBackground = result.fit.backgroundScale * result.backgroundCounts(:);
    observedSignal = max(result.centerCounts(:) - fittedBackground, 0);
    expectedSignal = expectedCounts - fittedBackground;
    expectedSignal(~isfinite(expectedSignal)) = NaN;

    fig = figure('Visible', visibleState, 'Color', 'w', ...
        'Position', [100 100 1500 850]);

    subplot(2, 3, 1);
    imagesc(result.centerInfo.focusImage);
    axis image;
    colormap(gca, 'parula');
    colorbar;
    hold on;
    if isfield(result.centerInfo, 'initialCenterXY')
        plot(result.centerInfo.initialCenterXY(1), ...
            result.centerInfo.initialCenterXY(2), ...
            'g+', 'MarkerSize', 13, 'LineWidth', 1.5);
    end
    plot(result.centerInfo.centerXY(1), result.centerInfo.centerXY(2), ...
        'rx', 'MarkerSize', 14, 'LineWidth', 2);
    plot(result.centerInfo.centerXY(1), result.centerInfo.centerXY(2), ...
        'ro', 'MarkerSize', 14, 'LineWidth', 1.5);
    if isfield(result.centerInfo, 'initialCenterXY')
        legend({'Gaussian seed', 'iterative center'}, 'Location', 'best', ...
            'TextColor', 'w', 'Color', [0 0 0 0.3], 'EdgeColor', 'none');
    end
    title(sprintf('Focal summed image, z = %.4f um', result.focusZUm), ...
        'Interpreter', 'none');
    xlabel('x pixel');
    ylabel('y pixel');

    subplot(2, 3, 2);
    plotDetectorMap(result, observedSignal, 'Observed signal counts');

    subplot(2, 3, 3);
    plotDetectorMap(result, expectedSignal, 'Fitted signal counts');

    subplot(2, 3, 4);
    plotDetectorMap(result, devianceResidual, 'Poisson deviance residual');

    subplot(2, 3, 5);
    x = 1:nCh;
    bar(x - 0.18, result.centerCounts(:), 0.36, ...
        'FaceColor', [0.20 0.44 0.74], 'EdgeColor', 'none');
    hold on;
    bar(x + 0.18, expectedCounts(:), 0.36, ...
        'FaceColor', [0.85 0.33 0.10], 'EdgeColor', 'none');
    plot(x, fittedBackground, 'k.-', 'LineWidth', 1.2, 'MarkerSize', 10);
    grid on;
    xlim([0.5 nCh + 0.5]);
    xlabel('detector index');
    ylabel('counts');
    legend({'observed raw', 'Poisson expected', 'fitted background'}, ...
        'Location', 'best');
    title(sprintf('Poisson deviance %.4g, reduced %.3g, photons %.4g', ...
        result.fit.deviance, result.fit.reducedDeviance, ...
        result.fit.photonScale), 'Interpreter', 'none');

    subplot(2, 3, 6);
    plotFittedCoefficients(gca, result);
end

function plotFittedCoefficients(ax, result)
%PLOTFITTEDCOEFFICIENTS Bar chart of fitted Zernike coefficients with +/-1 sigma error bars.
    modes = result.fit.fitModes(:);
    nM = numel(modes);
    valuesNm = zeros(nM, 1);
    for k = 1:nM
        valuesNm(k) = result.fit.estCoeffs.(modes{k}) * result.fit.sim.lamRef * 1000;
    end
    if isfield(result.fit, 'coeffStdNm') && numel(result.fit.coeffStdNm) == nM
        errNm = double(result.fit.coeffStdNm(:));
    else
        errNm = nan(nM, 1);
    end

    x = (1:nM).';
    bar(ax, x, valuesNm, 0.6, 'FaceColor', [0.35 0.55 0.75], 'EdgeColor', 'none');
    hold(ax, 'on');
    errorbar(ax, x, valuesNm, errNm, 'k', 'LineStyle', 'none', ...
        'LineWidth', 1.1, 'CapSize', 8);
    yline(ax, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    hold(ax, 'off');
    set(ax, 'XTick', x, 'XTickLabel', prettyModeLabels(modes));
    try, xtickangle(ax, 30); catch, end
    xlim(ax, [0.5 nM + 0.5]);
    ylabel(ax, 'coefficient [nm RMS]');
    grid(ax, 'on');
    ident = result.fit.identifiability;
    title(ax, sprintf('Fitted aberrations \\pm1\\sigma (rank %d/%d, cond %.2g)', ...
        ident.rank, ident.nParameters, ident.conditionNumber), ...
        'Interpreter', 'tex');
end

function labels = prettyModeLabels(modes)
    labels = cell(numel(modes), 1);
    for k = 1:numel(modes)
        switch lower(modes{k})
            case 'astig_x',   labels{k} = 'astig 0/90';
            case 'astig_y',   labels{k} = 'astig 45';
            case 'coma_x',    labels{k} = 'coma x';
            case 'coma_y',    labels{k} = 'coma y';
            case 'trefoil_x', labels{k} = 'trefoil x';
            case 'trefoil_y', labels{k} = 'trefoil y';
            case 'spherical', labels{k} = 'spherical';
            case 'defocus',   labels{k} = 'defocus';
            otherwise,        labels{k} = strrep(modes{k}, '_', ' ');
        end
    end
end

function plotDetectorMap(result, values, titleText)
    values = double(values(:));
    if isfield(result.detectorLayoutDiagnostics, 'detXYUm') && ...
            size(result.detectorLayoutDiagnostics.detXYUm, 1) == numel(values)
        xy = result.detectorLayoutDiagnostics.detXYUm;
        ax = gca;
        plotDetectorHexMap(xy, values, 'Parent', ax, ...
            'CellScale', 1.01, 'EdgeColor', [0.25 0.25 0.25], ...
            'LineWidth', 0.8);
        hold(ax, 'on');
        for k = 1:numel(values)
            text(ax, xy(k,1), xy(k,2), sprintf('%d', result.data.channelIDs(k)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 8, 'Color', 'w', 'FontWeight', 'bold');
        end
        hold(ax, 'off');
    else
        bar(values);
        xlabel('detector index');
        ylabel('value');
    end
    colorbar;
    grid on;
    title(titleText, 'Interpreter', 'none');
end

function tf = hasRawStackFiles(scanFolder)
    tf = ~isempty(sortedSeriesFiles(scanFolder)) || ...
        exist(fullfile(scanFolder, 'IntensityImage.pqdat'), 'file') == 2 || ...
        ~isempty(dir(fullfile(scanFolder, 'IntensityImage*.pqdat')));
end

function preprocessRawScanFolder(scanFolder, opts)
    batchOpts = opts.batchOptions;
    batchOpts.inputSource = resolveBatchInputSource(scanFolder, opts);
    batchOpts.zStepUm = opts.zStepUm;
    batchOpts.verbose = opts.verbose;
    if ~isempty(opts.ptuReaderFolder)
        batchOpts.ptuReaderFolder = opts.ptuReaderFolder;
    end
    if isfinitePositiveScalar(opts.xyPixelSizeUm)
        batchOpts.xyPixelSizeUm = opts.xyPixelSizeUm;
    end
    if isempty(opts.batchResultsRoot)
        batchAnalyzeLuminosaPSFs(scanFolder, [], batchOpts);
    else
        batchAnalyzeLuminosaPSFs(scanFolder, opts.batchResultsRoot, batchOpts);
    end
end

function source = resolveBatchInputSource(scanFolder, opts)
    source = lower(char(opts.batchInputSource));
    if strcmp(source, 'auto')
        if ~isempty(sortedSeriesFiles(scanFolder))
            source = 'ptu';
        elseif exist(fullfile(scanFolder, 'IntensityImage.pqdat'), 'file') == 2 || ...
                ~isempty(dir(fullfile(scanFolder, 'IntensityImage*.pqdat')))
            source = 'pqdat';
        else
            source = 'ptu';
        end
    end
    if ~ismember(source, {'ptu','pqdat'})
        error('estimateFocalPlaneCenterPointISMPoisson:BadBatchInputSource', ...
            'batchInputSource must be auto, ptu, or pqdat.');
    end
end

function files = sortedSeriesFiles(folderPath)
    files = dir(fullfile(folderPath, 'Series*.ptu'));
    if isempty(files)
        return;
    end
    idx = zeros(numel(files), 1);
    for k = 1:numel(files)
        tok = regexp(files(k).name, 'Series_(\d+)\.ptu$', 'tokens', 'once');
        if isempty(tok)
            idx(k) = inf;
        else
            idx(k) = str2double(tok{1});
        end
    end
    [~, order] = sort(idx);
    files = files(order);
end

function csvFile = findBatchAlignmentCsvForScanFolder(scanFolder, opts)
    candidates = batchAlignmentCandidates(scanFolder, opts);
    csvFile = '';
    newestTime = -inf;
    for k = 1:numel(candidates)
        if exist(candidates{k}, 'file') ~= 2
            continue;
        end
        info = dir(candidates{k});
        if ~isempty(info) && info.datenum > newestTime
            newestTime = info.datenum;
            csvFile = candidates{k};
        end
    end
end

function candidates = batchAlignmentCandidates(scanFolder, opts)
    scanFolder = stripTrailingFilesep(scanFolder);
    [dataRoot, scanName] = fileparts(scanFolder);
    [dataParent, datasetName] = fileparts(dataRoot);
    stem = scanOutputStem(scanName);

    parentOutputRoot = fullfile(dataParent, 'PSF_batch_outputs', ...
        sanitizeFileName(datasetName));
    directOutputRoot = fullfile(dataRoot, 'PSF_batch_outputs', ...
        sanitizeFileName(scanName));
    pqdatOutputRoot = fullfile(scanFolder, 'results_psf3d');
    roots = {scanFolder, parentOutputRoot, directOutputRoot, pqdatOutputRoot};
    if ~isempty(opts.batchResultsRoot)
        roots = [{opts.batchResultsRoot}, roots];
    end

    candidates = {};
    for r = 1:numel(roots)
        root = roots{r};
        candidates{end+1} = fullfile(root, 'frame_alignment.csv'); %#ok<AGROW>
        candidates{end+1} = fullfile(root, 'xz_yz_plots', ...
            sprintf('%s_frame_alignment.csv', stem)); %#ok<AGROW>
        candidates{end+1} = fullfile(root, sanitizeFileName(scanName), ...
            'frame_alignment.csv'); %#ok<AGROW>
        hits = dir(fullfile(root, '*', 'frame_alignment.csv'));
        scanKey = lower(sanitizeFileName(scanName));
        stemKey = lower(stem);
        for k = 1:numel(hits)
            if ~hits(k).isdir
                [~, folderName] = fileparts(hits(k).folder);
                folderKey = lower(folderName);
                if contains(folderKey, scanKey) || contains(folderKey, stemKey)
                    candidates{end+1} = fullfile(hits(k).folder, hits(k).name); %#ok<AGROW>
                end
            end
        end
    end
end

function detXY = validateDetectorXY(detXY, nCh)
    detXY = double(detXY);
    if ~ismatrix(detXY) || size(detXY, 2) ~= 2 || size(detXY, 1) ~= nCh || ...
            any(~isfinite(detXY(:)))
        error('estimateFocalPlaneCenterPointISMPoisson:BadDetectorXY', ...
            'detectorXYUm must be a finite [%d x 2] array.', nCh);
    end
end

function diagnostics = detectorDiagnostics(detXY, pixelSizeUm, airyUnitUm, source)
    diagnostics = struct();
    diagnostics.source = source;
    diagnostics.pixelSizeUm = pixelSizeUm;
    diagnostics.pixelSizeNm = 1000 * pixelSizeUm;
    diagnostics.airyUnitUm = airyUnitUm;
    diagnostics.detXYUm = detXY;
    diagnostics.detXYNm = 1000 * detXY;
    diagnostics.detXYAU = detXY / airyUnitUm;
    diagnostics.detectorRadiusAU = hypot(diagnostics.detXYAU(:,1), ...
        diagnostics.detXYAU(:,2));
end

function value = effectivePropagatingNA(sim)
    value = sim.NA;
    if isfield(sim, 'sampleGeometry') && strcmpi(sim.sampleGeometry, 'airOnGlass')
        value = min(sim.NA, sim.nSample);
    end
end

function pitch = resolveDetectorPitchSampleUm(opts)
    pitch = NaN;
    if isfinitePositiveScalar(opts.detectorPitchSampleUm)
        pitch = double(opts.detectorPitchSampleUm);
        return;
    end
    if isfinitePositiveScalar(opts.detectorHardwarePitchUm) && ...
            isfinitePositiveScalar(opts.detectorTotalMagnification)
        pitch = double(opts.detectorHardwarePitchUm) / ...
            double(opts.detectorTotalMagnification);
    end
end

function source = detectorPitchSource(opts)
    if isfinitePositiveScalar(opts.detectorPitchSampleUm)
        source = 'detectorPitchSampleUm';
    elseif isfinitePositiveScalar(opts.detectorHardwarePitchUm) && ...
            isfinitePositiveScalar(opts.detectorTotalMagnification)
        source = 'detectorHardwarePitchUm / detectorTotalMagnification';
    else
        source = 'model default';
    end
end

function pitch = medianNearestDistance(xy)
    n = size(xy, 1);
    nearest = nan(n, 1);
    for k = 1:n
        d = hypot(xy(:,1) - xy(k,1), xy(:,2) - xy(k,2));
        d(k) = inf;
        nearest(k) = min(d);
    end
    pitch = median(nearest(isfinite(nearest)));
end

function qe = validateDetectorQE(qe, channelIDs)
    if isempty(qe)
        return;
    end
    qe = double(qe(:));
    if numel(qe) ~= numel(channelIDs) || any(~isfinite(qe)) || any(qe < 0)
        error('estimateFocalPlaneCenterPointISMPoisson:BadDetectorQE', ...
            'detectorQE must contain one nonnegative finite value per channel.');
    end
end

function bounds = validateTwoElementBounds(bounds, name)
    bounds = double(bounds(:)).';
    if numel(bounds) ~= 2 || any(~isfinite(bounds)) || bounds(2) <= bounds(1)
        error('estimateFocalPlaneCenterPointISMPoisson:BadBounds', ...
            '%s must be a finite increasing two-element vector.', name);
    end
end

function value = numericField(s, name, fallback)
    value = fallback;
    if isstruct(s) && isfield(s, name)
        candidate = s.(name);
        if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
            value = double(candidate);
        end
    end
end

function tf = isfinitePositiveScalar(v)
    tf = isnumeric(v) && isscalar(v) && isfinite(v) && v > 0;
end

function stem = localFileStem(pathName)
    if isempty(pathName) || isnumeric(pathName)
        stem = 'numeric_stack';
        return;
    end
    pathName = stripTrailingFilesep(char(pathName));
    [~, stem] = fileparts(pathName);
    if isempty(stem)
        stem = 'stack';
    end
end

function s = stripTrailingFilesep(s)
    while ~isempty(s) && (s(end) == filesep || s(end) == '/' || s(end) == '\')
        s(end) = [];
    end
end

function stem = sanitizeFileName(name)
    stem = regexprep(char(name), '[^A-Za-z0-9._-]+', '_');
    stem = regexprep(stem, '_+', '_');
    stem = regexprep(stem, '^_+|_+$', '');
    if isempty(stem)
        stem = 'stack';
    end
end

function stem = scanOutputStem(label)
    stem = regexprep(char(label), '[^A-Za-z0-9]+', '_');
    stem = regexprep(stem, '_+', '_');
    stem = regexprep(stem, '^_+|_+$', '');
    if isempty(stem)
        stem = 'scan';
    end
    if ~isletter(stem(1))
        stem = ['x', stem];
    end
end
