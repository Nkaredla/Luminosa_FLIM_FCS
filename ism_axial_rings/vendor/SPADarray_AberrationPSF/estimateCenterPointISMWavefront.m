function result = estimateCenterPointISMWavefront(focusInput, defocusInput, varargin)
%ESTIMATECENTERPOINTISMWAVEFRONT Fit aberrations from ISM center micro-images.
%
%   result = estimateCenterPointISMWavefront(focusInput, defocusInput)
%
%   This estimator uses only the detector intensity distribution sampled at
%   the bead centre fitted in the detector-summed focal image: 23 detector
%   values at focus plus 23 detector values at an axial diversity plane.
%
%   Inputs can be:
%       - focus and defocus detector images [y x 23]
%       - a single two-plane stack [y x 23 x 2] or [y x 2 x 23]
%       - center vectors [23 x 1] and [23 x 1]
%       - PTU/MAT files containing detector-resolved data
%
%   The default plane spacing is [0 0.5] um, i.e. focus plus 500 nm.
%
%   Example:
%       res = estimateCenterPointISMWavefront(focusPtu, defocusPtu, ...
%           'planeZ', [0 0.5]);
%
%       res = estimateCenterPointISMWavefront(twoPlaneRaw4, [], ...
%           'centerNormalization', 'perPlane');

    if nargin < 1 || isempty(focusInput)
        error('estimateCenterPointISMWavefront:MissingInput', ...
            'Provide focusInput, or a single [y x 23 x 2] two-plane input.');
    end
    if nargin < 2
        defocusInput = [];
    end

    opts = parseOptions(varargin{:});
    addRequiredPaths(opts);

    [raw4, inputMeta] = loadTwoPlaneInput(focusInput, defocusInput, opts);
    raw4 = standardizeTwoPlaneNumeric(raw4, opts.channelIDs);

    [corrected4, correction] = applyIntensityCorrections(raw4, inputMeta.channelIDs, opts);
    [centerValues, centerInfo] = extractCenterDetectorValues(corrected4, opts);
    [dataN, normInfo] = normalizeCenterValues(centerValues, opts.centerNormalization);

    [sim, detectorLayoutDiagnostics] = configureCenterPointSim(raw4, inputMeta, opts);
    [fit, focusAnchored] = fitCenterPointModel(dataN, sim, opts.planeZ, opts);
    diagOpts = opts;
    diagOpts.planeWeights = fit.planeWeights;
    sufficiency = centerPointIdentifiability(sim, opts.fitModes, opts.planeZ, ...
        fit.paramVector, diagOpts);

    result = struct();
    result.rawData = raw4;
    result.correctedData = corrected4;
    result.planeZ = opts.planeZ(:).';
    result.centerValues = centerValues;
    result.normalizedCenterValues = dataN;
    result.centerInfo = centerInfo;
    result.normalization = normInfo;
    result.correction = correction;
    result.inputMeta = inputMeta;
    result.detectorLayoutDiagnostics = detectorLayoutDiagnostics;
    result.detectorQE = opts.detectorQE;
    result.fit = fit;
    result.fitStrategy = opts.fitStrategy;
    result.focusAnchored = focusAnchored;
    result.sufficiency = sufficiency;
    result.validity = struct( ...
        'correctionReady', false, ...
        'recommendedUse', 'diagnostic or full-stack initialization only', ...
        'reason', ['The normalized center-point fit discards spatial and ' ...
        'absolute axial photon information and has no bootstrap/model-error ' ...
        'acceptance gate. Use estimateFullStackISMWavefront for correction.']);
    result.options = opts;
    result.outputDir = resolveOutputDir(opts);

    if opts.writeOutputs
        writeCenterPointOutputs(result);
    end

    if opts.verbose
        printCenterPointSummary(result);
    end
end

function opts = parseOptions(varargin)
    p = inputParser;
    p.FunctionName = 'estimateCenterPointISMWavefront';

    addParameter(p, 'planeZ', [0 0.5]);
    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'channelOrder', []);
    addParameter(p, 'ptuReaderFolder', defaultPtuReaderFolder());
    addParameter(p, 'ptuPhotonsPerChunk', 1e6);
    addParameter(p, 'frameIndices', []);
    addParameter(p, 'focusFrameIndices', []);
    addParameter(p, 'defocusFrameIndices', []);
    addParameter(p, 'frameCombine', 'sum');
    addParameter(p, 'matVariable', '');

    addParameter(p, 'darkFile', 'D:\Luminosa\Data\ISMdark_counts.ptu');
    addParameter(p, 'darkMode', 'auto');
    addParameter(p, 'darkScale', 1);
    addParameter(p, 'darkPerChannel', []);
    addParameter(p, 'subtractBoundary', true);
    addParameter(p, 'boundaryWidthPx', 3);

    addParameter(p, 'centerXY', []);
    addParameter(p, 'centerMode', 'gaussian');
    addParameter(p, 'centerThresholdFraction', 0.20);
    addParameter(p, 'centerSampleMode', 'subpixel');
    addParameter(p, 'centerNormalization', 'perPlane');
    addParameter(p, 'modelSampleXY', [0 0]);

    addParameter(p, 'objectiveNA', 1.2);
    addParameter(p, 'objectiveMagnification', 60);
    addParameter(p, 'emissionWavelengthUm', []);
    addParameter(p, 'excitationWavelengthUm', []);
    addParameter(p, 'mediumRefractiveIndex', []);
    addParameter(p, 'sampleGeometry', 'airOnGlass');
    addParameter(p, 'interfaceRadialWeightMode', 'sampleSolidAngle');
    addParameter(p, 'diffractionModel', 'vectorial');
    addParameter(p, 'vectorialPolarizationMode', 'xyAverage');
    addParameter(p, 'excitationPolarizationMode', 'x');
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
    addParameter(p, 'sim', []);
    addParameter(p, 'xyPixelSizeUm', []);
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
    addParameter(p, 'airyUnitUm', []);
    addParameter(p, 'detectorPitchSampleUm', []);
    addParameter(p, 'detectorHardwarePitchUm', []);
    addParameter(p, 'detectorTotalMagnification', []);
    addParameter(p, 'detectorQE', []);

    addParameter(p, 'fitModes', {'astig_x','astig_y','coma_x','coma_y','spherical'});
    addParameter(p, 'fitStrategy', 'joint');
    addParameter(p, 'planeWeights', []);
    addParameter(p, 'focusWeight', 1);
    addParameter(p, 'diversityWeight', 0.20);
    addParameter(p, 'focusSeedAmplitude', 0.05);
    addParameter(p, 'focusMaxStarts', 4);
    addParameter(p, 'focusCoarseMaxIter', 2);
    addParameter(p, 'focusRefineMaxIter', 2);
    addParameter(p, 'focusCoarseDetectorSubsamples', 3);
    addParameter(p, 'signSelectionModes', ...
        {'astig_x','astig_y','spherical'});
    addParameter(p, 'fitXY', false);
    addParameter(p, 'fitZ', false);
    addParameter(p, 'fitDetectorPitchScale', false);
    addParameter(p, 'initialCoeffs', struct());
    addParameter(p, 'initialXY', [0 0]);
    addParameter(p, 'initialZOffset', 0);
    addParameter(p, 'initialDetectorPitchScale', 1);
    addParameter(p, 'maxIter', 6);
    addParameter(p, 'jacobianScheme', 'forward');
    addParameter(p, 'fdCoeff', 0.01);
    addParameter(p, 'fdXY', []);
    addParameter(p, 'fdZ', 0.02);
    addParameter(p, 'fdDetectorPitchScale', 0.02);
    addParameter(p, 'regCoeff', 1e-5);
    addParameter(p, 'regXY', 1e-4);
    addParameter(p, 'regZ', 1e-4);
    addParameter(p, 'regDetectorPitchScale', 1e-4);
    addParameter(p, 'maxCoeffStep', 0.04);
    addParameter(p, 'maxXYStep', 0.03);
    addParameter(p, 'maxZStep', 0.03);
    addParameter(p, 'maxDetectorPitchScaleStep', 0.05);
    addParameter(p, 'detectorPitchScaleBounds', [0.5 2.0]);
    addParameter(p, 'tolStep', 1e-5);
    addParameter(p, 'beadSubsamples', [3 3 3]);
    addParameter(p, 'modelDz', 0.025);
    addParameter(p, 'modelZPadding', 0.50);

    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;

    opts.planeZ = double(opts.planeZ(:)).';
    if numel(opts.planeZ) ~= 2
        error('estimateCenterPointISMWavefront:PlaneZ', ...
            'This center-point estimator requires exactly two planeZ values.');
    end
    opts.frameCombine = lower(char(opts.frameCombine));
    opts.darkMode = lower(char(opts.darkMode));
    opts.centerMode = lower(char(opts.centerMode));
    opts.centerSampleMode = lower(char(opts.centerSampleMode));
    opts.centerNormalization = lower(char(opts.centerNormalization));
    opts.sampleGeometry = char(opts.sampleGeometry);
    opts.interfaceRadialWeightMode = char(opts.interfaceRadialWeightMode);
    opts.diffractionModel = char(opts.diffractionModel);
    opts.vectorialPolarizationMode = char(opts.vectorialPolarizationMode);
    opts.excitationPolarizationMode = char(opts.excitationPolarizationMode);
    opts.collectionPolarizationMode = char(opts.collectionPolarizationMode);
    opts.airInterfaceStageMedium = char(opts.airInterfaceStageMedium);
    if ~ismember(lower(opts.sampleGeometry),{'aironglass','homogeneous'})
        error('estimateCenterPointISMWavefront:BadSampleGeometry', ...
            'sampleGeometry must be airOnGlass or homogeneous.');
    end
    opts.detectorLayoutCenterMode = lower(char(opts.detectorLayoutCenterMode));
    opts.detectorLayoutScale = validatePositiveScalar( ...
        opts.detectorLayoutScale, 'detectorLayoutScale');
    opts.detectorQE = validateDetectorQE(opts.detectorQE, opts.channelIDs);
    if ~isnumeric(opts.detectorLayoutPositionSign) || ...
            ~isscalar(opts.detectorLayoutPositionSign) || ...
            ~isfinite(opts.detectorLayoutPositionSign) || ...
            ~ismember(opts.detectorLayoutPositionSign, [-1 1])
        error('estimateCenterPointISMWavefront:BadDetectorLayoutPositionSign', ...
            'detectorLayoutPositionSign must be +1 or -1.');
    end
    opts.fitStrategy = lower(char(opts.fitStrategy));
    opts.planeWeights = validatePlaneWeights(opts.planeWeights);
    opts.focusWeight = validateNonnegativeScalar(opts.focusWeight, 'focusWeight');
    opts.diversityWeight = validateNonnegativeScalar(opts.diversityWeight, 'diversityWeight');
    if opts.focusWeight <= 0 || opts.diversityWeight <= 0
        error('estimateCenterPointISMWavefront:BadFocusAnchoredWeights', ...
            'focusWeight and diversityWeight must both be positive.');
    end
    opts.focusSeedAmplitude = validateNonnegativeScalar( ...
        opts.focusSeedAmplitude, 'focusSeedAmplitude');
    opts.focusMaxStarts = validatePositiveInteger(opts.focusMaxStarts, ...
        'focusMaxStarts');
    opts.focusCoarseMaxIter = validateNonnegativeInteger( ...
        opts.focusCoarseMaxIter, 'focusCoarseMaxIter');
    opts.focusRefineMaxIter = validateNonnegativeInteger( ...
        opts.focusRefineMaxIter, 'focusRefineMaxIter');
    opts.focusCoarseDetectorSubsamples = validatePositiveInteger( ...
        opts.focusCoarseDetectorSubsamples, 'focusCoarseDetectorSubsamples');
    opts.jacobianScheme = lower(char(opts.jacobianScheme));
    if ~ismember(opts.jacobianScheme, {'forward','central'})
        error('estimateCenterPointISMWavefront:BadJacobianScheme', ...
            'jacobianScheme must be ''forward'' or ''central''.');
    end
    if ischar(opts.signSelectionModes) || isstring(opts.signSelectionModes)
        opts.signSelectionModes = cellstr(opts.signSelectionModes);
    end
    opts.signSelectionModes = opts.signSelectionModes(:).';
    opts.fitDetectorPitchScale = logical(opts.fitDetectorPitchScale);
    opts.initialDetectorPitchScale = validatePositiveScalar( ...
        opts.initialDetectorPitchScale, 'initialDetectorPitchScale');
    opts.fdDetectorPitchScale = validatePositiveScalar( ...
        opts.fdDetectorPitchScale, 'fdDetectorPitchScale');
    opts.regDetectorPitchScale = validateNonnegativeScalar( ...
        opts.regDetectorPitchScale, 'regDetectorPitchScale');
    opts.maxDetectorPitchScaleStep = validatePositiveScalar( ...
        opts.maxDetectorPitchScaleStep, 'maxDetectorPitchScaleStep');
    opts.detectorPitchScaleBounds = validateDetectorPitchScaleBounds( ...
        opts.detectorPitchScaleBounds);
    opts.modelSampleXY = double(opts.modelSampleXY(:)).';
    if numel(opts.modelSampleXY) ~= 2 || any(~isfinite(opts.modelSampleXY))
        error('estimateCenterPointISMWavefront:BadModelSampleXY', ...
            'modelSampleXY must be finite [x y] coordinates in micrometers.');
    end
    if ischar(opts.fitModes) || isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    end
    opts.fitModes = opts.fitModes(:).';
end

function folder = defaultPtuReaderFolder()
    thisDir = fileparts(mfilename('fullpath'));
    luminosaRoot = fileparts(fileparts(thisDir));
    folder = fullfile(luminosaRoot, 'Luminosa_FLIM_FCS');
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

function [raw4, meta] = loadTwoPlaneInput(focusInput, defocusInput, opts)
    meta = emptyInputMeta();

    if isempty(defocusInput)
        if isnumeric(focusInput)
            raw4 = focusInput;
            meta.source = 'numeric two-plane or center-vector input';
            meta.channelIDs = defaultChannelIDsFromData(raw4, opts.channelIDs);
            return;
        end

        keepFrames = ~isempty(opts.frameIndices);
        [stack, oneMeta] = loadDetectorInput(focusInput, opts, opts.frameIndices, keepFrames);
        raw4 = stack;
        meta = oneMeta;
        meta.source = 'single two-plane input';
        return;
    end

    if isnumeric(focusInput) && isnumeric(defocusInput)
        focusStack = standardizeSinglePlaneOrVector(focusInput, opts.channelIDs);
        defocusStack = standardizeSinglePlaneOrVector(defocusInput, opts.channelIDs);
        assertSameCenterInputSize(focusStack, defocusStack);
        raw4 = cat(4, focusStack, defocusStack);
        meta.source = 'numeric focus and defocus inputs';
        meta.channelIDs = defaultChannelIDsFromData(raw4, opts.channelIDs);
        return;
    end

    [focusStack, focusMeta] = loadDetectorInput(focusInput, opts, opts.focusFrameIndices, false);
    [defocusStack, defocusMeta] = loadDetectorInput(defocusInput, opts, opts.defocusFrameIndices, false);
    focusStack = standardizeSinglePlaneOrVector(focusStack, opts.channelIDs);
    defocusStack = standardizeSinglePlaneOrVector(defocusStack, opts.channelIDs);
    assertSameCenterInputSize(focusStack, defocusStack);

    raw4 = cat(4, focusStack, defocusStack);
    meta = focusMeta;
    meta.source = 'focus and defocus inputs';
    meta.focus = focusMeta;
    meta.defocus = defocusMeta;
    if ~isempty(focusMeta.channelIDs) && ~isempty(defocusMeta.channelIDs) && ...
            ~isequal(focusMeta.channelIDs(:), defocusMeta.channelIDs(:))
        warning('estimateCenterPointISMWavefront:ChannelIDMismatch', ...
            'Focus and defocus channel IDs differ; using the focus order.');
    end
end

function meta = emptyInputMeta()
    meta = struct();
    meta.source = '';
    meta.file = '';
    meta.files = {};
    meta.channelIDs = [];
    meta.xyPixelSizeUm = NaN;
    meta.head = struct();
end

function ids = defaultChannelIDsFromData(data, requested)
    nCh = 23;
    dims = size(data);
    if isvector(data) && numel(data) == numel(requested)
        nCh = numel(requested);
    elseif ndims(data) >= 3
        if any(dims == numel(requested))
            nCh = numel(requested);
        elseif any(dims == 23)
            nCh = 23;
        end
    end
    ids = requested(:);
    if numel(ids) ~= nCh
        ids = (1:nCh).';
    end
end

function [stack, meta] = loadDetectorInput(inputValue, opts, frameIndices, keepFrames)
    if isnumeric(inputValue)
        stack = inputValue;
        meta = emptyInputMeta();
        meta.source = 'numeric input';
        meta.channelIDs = defaultChannelIDsFromData(stack, opts.channelIDs);
        return;
    end

    fileName = resolveInputFile(char(inputValue));
    [~,~,ext] = fileparts(fileName);
    switch lower(ext)
        case '.ptu'
            [stack, meta] = readPtuDetectorStack(fileName, opts, frameIndices, keepFrames);
        case '.mat'
            [stack, meta] = readMatDetectorStack(fileName, opts);
        otherwise
            error('estimateCenterPointISMWavefront:BadInputFile', ...
                'Unsupported input file "%s". Use PTU, MAT, or numeric data.', fileName);
    end
end

function fileName = resolveInputFile(inputPath)
    if exist(inputPath, 'dir') == 7
        candidates = { ...
            fullfile(inputPath, 'RawImage.ptu'), ...
            fullfile(inputPath, 'RawData.ptu')};
        for k = 1:numel(candidates)
            if exist(candidates{k}, 'file') == 2
                fileName = candidates{k};
                return;
            end
        end

        files = dir(fullfile(inputPath, '*.ptu'));
        if isempty(files)
            files = dir(fullfile(inputPath, '*.mat'));
        end
        if isempty(files)
            error('estimateCenterPointISMWavefront:NoInputFile', ...
                'No PTU or MAT file found in %s.', inputPath);
        end
        [~, newest] = max([files.datenum]);
        fileName = fullfile(files(newest).folder, files(newest).name);
        return;
    end

    if exist(inputPath, 'file') ~= 2
        error('estimateCenterPointISMWavefront:MissingInputFile', ...
            'Input file was not found: %s', inputPath);
    end
    fileName = inputPath;
end

function [stack, meta] = readMatDetectorStack(fileName, opts)
    S = load(fileName);
    varName = char(opts.matVariable);
    if ~isempty(varName)
        if ~isfield(S, varName)
            error('estimateCenterPointISMWavefront:MissingMatVariable', ...
                'Variable "%s" was not found in %s.', varName, fileName);
        end
        stack = S.(varName);
    else
        stack = chooseMatDetectorVariable(S);
    end

    meta = emptyInputMeta();
    meta.source = 'MAT file';
    meta.file = fileName;
    meta.channelIDs = defaultChannelIDsFromData(stack, opts.channelIDs);
    if isfield(S, 'channelIDs')
        meta.channelIDs = double(S.channelIDs(:));
    elseif isfield(S, 'dind')
        meta.channelIDs = double(S.dind(:));
    end
    meta.xyPixelSizeUm = xyPixelSizeFromLoadedStruct(S);
end

function xyPixelSizeUm = xyPixelSizeFromLoadedStruct(S)
    xyPixelSizeUm = NaN;
    candidates = {};
    if isfield(S, 'xyPixelSizeUm')
        candidates{end+1} = S.xyPixelSizeUm; %#ok<AGROW>
    end
    if isfield(S, 'opts') && isstruct(S.opts) && isfield(S.opts, 'xyPixelSizeUm')
        candidates{end+1} = S.opts.xyPixelSizeUm; %#ok<AGROW>
    end
    if isfield(S, 'scale') && isstruct(S.scale) && isfield(S.scale, 'xyPixelSizeUm')
        candidates{end+1} = S.scale.xyPixelSizeUm; %#ok<AGROW>
    end
    if isfield(S, 'scans') && isstruct(S.scans) && isfield(S.scans, 'xyPixelSizeUm')
        candidates{end+1} = [S.scans.xyPixelSizeUm]; %#ok<AGROW>
    end
    for k = 1:numel(candidates)
        values = double(candidates{k}(:));
        values = values(isfinite(values) & values > 0);
        if ~isempty(values)
            xyPixelSizeUm = values(1);
            return;
        end
    end
end

function stack = chooseMatDetectorVariable(S)
    preferred = {'raw4','rawData','detectorVolume','detectorStack', ...
        'channelVolume','channelVolumes','spadVolume','spadStack', ...
        'centerValues','detectorValues'};
    for k = 1:numel(preferred)
        name = preferred{k};
        if isfield(S, name) && isnumeric(S.(name)) && ndims(S.(name)) >= 2
            value = S.(name);
            if isDetectorLike(value)
                stack = value;
                return;
            end
        end
    end

    names = fieldnames(S);
    for k = 1:numel(names)
        value = S.(names{k});
        if isnumeric(value) && isDetectorLike(value)
            stack = value;
            return;
        end
    end

    error('estimateCenterPointISMWavefront:NoMatDetectorData', ...
        'No detector-resolved variable was found in the MAT file.');
end

function tf = isDetectorLike(value)
    tf = false;
    if ~isnumeric(value)
        return;
    end
    dims = size(value);
    if isvector(value) && numel(value) == 23
        tf = true;
    elseif ndims(value) == 2 && any(dims == 23) && any(dims == 2)
        tf = true;
    elseif ndims(value) >= 3 && any(dims == 23)
        tf = true;
    end
end

function [stack, meta] = readPtuDetectorStack(fileName, opts, frameIndices, keepFrames)
    hasFastReader = exist('PTU_MultiFrameScanReadFast', 'file') == 2;
    hasSlowReader = exist('PTU_MultiFrameScanRead', 'file') == 2;
    if ~hasFastReader && ~hasSlowReader
        error('estimateCenterPointISMWavefront:MissingPtuReader', ...
            'PTU_MultiFrameScanReadFast/PTU_MultiFrameScanRead are not on the MATLAB path.');
    end

    waitbarCleanup = suppressPtuWaitbars(); %#ok<NASGU>
    ptuOut = [];
    fastMessage = 'fast reader not available';
    if hasFastReader
        try
            ptuOut = PTU_MultiFrameScanReadFast(fileName, opts.ptuPhotonsPerChunk, false, false);
        catch fastErr
            fastMessage = fastErr.message;
        end
    end
    if isempty(ptuOut)
        if ~hasSlowReader
            error('estimateCenterPointISMWavefront:PtuReadFailed', ...
                'Could not read %s as a detector scan. Fast: %s Slow: slow reader not available', ...
                fileName, fastMessage);
        end
        try
            ptuOut = PTU_MultiFrameScanRead(fileName, opts.ptuPhotonsPerChunk);
        catch slowErr
            error('estimateCenterPointISMWavefront:PtuReadFailed', ...
                'Could not read %s as a detector scan. Fast: %s Slow: %s', ...
                fileName, fastMessage, slowErr.message);
        end
    end

    [stack, channelIDs] = ptuOutputToDetectorStack(ptuOut, opts, frameIndices, keepFrames);

    meta = emptyInputMeta();
    meta.source = 'PTU detector scan';
    meta.file = fileName;
    meta.files = {fileName};
    meta.channelIDs = channelIDs;
    if isfield(ptuOut, 'head')
        meta.head = ptuOut.head;
        meta.xyPixelSizeUm = numericField(ptuOut.head, 'ImgHdr_PixResol', NaN);
    end
end

function cleanup = suppressPtuWaitbars()
    oldVisible = get(groot, 'DefaultFigureVisible');
    closePtuWaitbars();
    try
        set(groot, 'DefaultFigureVisible', 'off');
    catch
    end
    cleanup = onCleanup(@() restorePtuWaitbars(oldVisible));
end

function restorePtuWaitbars(oldVisible)
    try
        set(groot, 'DefaultFigureVisible', oldVisible);
    catch
    end
    closePtuWaitbars();
end

function closePtuWaitbars()
    try
        h = findall(groot, 'Type', 'figure', 'Tag', 'TMWWaitbar');
        if ~isempty(h)
            delete(h(ishghandle(h)));
        end
    catch
    end
end

function [stackYX, channelIDs] = ptuOutputToDetectorStack(ptuOut, opts, frameIndices, keepFrames)
    if isfield(ptuOut, 'dind')
        channelIDs = double(ptuOut.dind(:));
    else
        channelIDs = [];
    end

    if keepFrames
        if ~isfield(ptuOut, 'tag') || isempty(ptuOut.tag)
            error('estimateCenterPointISMWavefront:NoPtuFrames', ...
                'The PTU output does not contain per-frame tag data.');
        end
        frameIndices = validateFrameIndices(frameIndices, size(ptuOut.tag, 4));
        stack = double(ptuOut.tag(:,:,:,frameIndices));
        stackYX = permute(stack, [2 1 3 4]);
    else
        if isempty(frameIndices)
            if isfield(ptuOut, 'tags') && ~isempty(ptuOut.tags)
                stack = double(ptuOut.tags);
            elseif isfield(ptuOut, 'tag') && ~isempty(ptuOut.tag)
                stack = sum(double(ptuOut.tag), 4);
            else
                error('estimateCenterPointISMWavefront:NoPtuImages', ...
                    'No detector image stack was found in the PTU output.');
            end
        else
            if isfield(ptuOut, 'tag') && ~isempty(ptuOut.tag)
                frameIndices = validateFrameIndices(frameIndices, size(ptuOut.tag, 4));
                switch opts.frameCombine
                    case 'sum'
                        stack = sum(double(ptuOut.tag(:,:,:,frameIndices)), 4);
                    case 'mean'
                        stack = mean(double(ptuOut.tag(:,:,:,frameIndices)), 4);
                    otherwise
                        error('estimateCenterPointISMWavefront:BadFrameCombine', ...
                            'frameCombine must be ''sum'' or ''mean''.');
                end
            elseif isfield(ptuOut, 'tags') && ~isempty(ptuOut.tags)
                frameIndices = validateFrameIndices(frameIndices, 1);
                stack = double(ptuOut.tags);
            else
                error('estimateCenterPointISMWavefront:NoPtuFrames', ...
                    'Frame selection requires ptuOut.tag or single-frame ptuOut.tags.');
            end
        end
        stackYX = permute(stack, [2 1 3]);
    end

    if isempty(channelIDs)
        channelIDs = (1:size(stackYX, 3)).';
    end

    [stackYX, channelIDs] = selectAndOrderChannels(stackYX, channelIDs, opts);
end

function frameIndices = validateFrameIndices(frameIndices, nFrame)
    if isempty(frameIndices)
        frameIndices = 1:nFrame;
    end
    frameIndices = round(frameIndices(:)).';
    if any(frameIndices < 1) || any(frameIndices > nFrame)
        error('estimateCenterPointISMWavefront:BadFrameIndices', ...
            'frameIndices must be within 1:%d.', nFrame);
    end
end

function [stack, channelIDs] = selectAndOrderChannels(stack, channelIDs, opts)
    requested = opts.channelIDs(:);
    if ~isempty(requested)
        [present, loc] = ismember(double(requested), double(channelIDs(:)));
        if all(present)
            stack = stack(:,:,loc,:);
            channelIDs = channelIDs(loc);
        elseif size(stack, 3) == numel(requested)
            channelIDs = requested;
        else
            error('estimateCenterPointISMWavefront:MissingChannels', ...
                'Could not find all requested detector channel IDs in PTU output.');
        end
    end

    order = opts.channelOrder;
    if ~isempty(order)
        order = double(order(:));
        if numel(order) ~= size(stack, 3)
            error('estimateCenterPointISMWavefront:BadChannelOrder', ...
                'channelOrder must contain one entry per selected detector channel.');
        end

        if all(ismember(order, channelIDs))
            [~, loc] = ismember(order, channelIDs);
        elseif all(order >= 1 & order <= size(stack, 3))
            loc = order;
        else
            error('estimateCenterPointISMWavefront:BadChannelOrder', ...
                'channelOrder must be either selected channel IDs or 1-based indices.');
        end
        stack = stack(:,:,loc,:);
        channelIDs = channelIDs(loc);
    end
end

function raw4 = standardizeTwoPlaneNumeric(raw, channelIDs)
    raw = double(raw);
    nCh = numel(channelIDs);

    if isvector(raw) && numel(raw) == 2*nCh
        raw4 = reshape(raw(:), 1, 1, nCh, 2);
        return;
    end

    if ndims(raw) == 2 && any(size(raw) == nCh) && any(size(raw) == 2)
        if size(raw, 1) == nCh
            raw4 = reshape(raw, 1, 1, nCh, 2);
        else
            raw4 = reshape(raw.', 1, 1, nCh, 2);
        end
        return;
    end

    if ndims(raw) == 3 && size(raw, 3) == nCh
        error('estimateCenterPointISMWavefront:NeedTwoPlanes', ...
            'A single input produced one detector plane. Provide defocusInput or frameIndices.');
    end

    if ndims(raw) ~= 4
        error('estimateCenterPointISMWavefront:BadTwoPlaneData', ...
            'Two-plane input must be [y x 23 x 2], [y x 2 x 23], or [23 x 2].');
    end

    dims = size(raw);
    if dims(3) == nCh && dims(4) == 2
        raw4 = raw;
    elseif dims(3) == 2 && dims(4) == nCh
        raw4 = permute(raw, [1 2 4 3]);
    else
        error('estimateCenterPointISMWavefront:BadTwoPlaneData', ...
            'Could not identify detector and two-plane dimensions.');
    end
end

function stack = standardizeSinglePlaneOrVector(stack, channelIDs)
    stack = double(stack);
    nCh = numel(channelIDs);

    if isvector(stack) && numel(stack) == nCh
        stack = reshape(stack(:), 1, 1, nCh);
        return;
    end

    if ndims(stack) == 2 && any(size(stack) == nCh) && any(size(stack) == 1)
        stack = reshape(stack(:), 1, 1, nCh);
        return;
    end

    if ndims(stack) == 4
        if size(stack, 4) == 1
            stack = stack(:,:,:,1);
        elseif size(stack, 3) == 1 && size(stack, 4) == nCh
            stack = permute(stack(:,:,1,:), [1 2 4 3]);
            stack = stack(:,:,:,1);
        else
            error('estimateCenterPointISMWavefront:BadSinglePlaneData', ...
                'A focus or defocus input must contain one detector plane.');
        end
    end

    if ndims(stack) ~= 3 || size(stack, 3) ~= nCh
        error('estimateCenterPointISMWavefront:BadSinglePlaneData', ...
            'A focus or defocus detector input must be [y x 23] or [23 x 1].');
    end
end

function assertSameCenterInputSize(a, b)
    if ~isequal(size(a), size(b))
        error('estimateCenterPointISMWavefront:PlaneSizeMismatch', ...
            'Focus and defocus detector inputs must have the same size.');
    end
end

function [corrected, correction] = applyIntensityCorrections(raw4, channelIDs, opts)
    corrected = double(raw4);
    correction = struct();
    correction.darkMethod = 'none';
    correction.darkImage = [];
    correction.darkPerChannel = zeros(size(raw4, 3), 1);
    correction.boundaryBackground = [];

    if strcmpi(opts.darkMode, 'none')
        % keep raw data
    elseif ~isempty(opts.darkPerChannel)
        correction.darkMethod = 'darkPerChannel option';
        correction.darkPerChannel = double(opts.darkPerChannel(:));
        corrected = subtractChannelVector(corrected, correction.darkPerChannel, opts.darkScale);
    elseif any(strcmpi(opts.darkMode, {'auto','file'})) && ~isempty(opts.darkFile)
        try
            [darkImage, darkVector] = readDarkEstimate(opts.darkFile, size(raw4), channelIDs, opts);
            correction.darkImage = darkImage;
            correction.darkPerChannel = darkVector(:);
            if ~isempty(darkImage)
                correction.darkMethod = 'dark PTU image';
                corrected = max(corrected - opts.darkScale * darkImage, 0);
            elseif ~isempty(darkVector)
                correction.darkMethod = 'dark PTU per-channel vector';
                corrected = subtractChannelVector(corrected, darkVector, opts.darkScale);
            end
        catch err
            if strcmpi(opts.darkMode, 'file')
                rethrow(err);
            end
            warning('estimateCenterPointISMWavefront:DarkFileFallback', ...
                'Could not use dark file; falling back to boundary background. %s', err.message);
        end
    end

    if opts.subtractBoundary && size(corrected, 1) > 1 && size(corrected, 2) > 1
        bg = boundaryBackground(corrected, opts.boundaryWidthPx);
        correction.boundaryBackground = bg;
        corrected = max(corrected - bg, 0);
    end
end

function corrected = subtractChannelVector(data, values, scale)
    values = double(values(:));
    if numel(values) ~= size(data, 3)
        error('estimateCenterPointISMWavefront:BadDarkVector', ...
            'darkPerChannel must contain one value per selected detector channel.');
    end
    corrected = max(data - scale * reshape(values, 1, 1, [], 1), 0);
end

function [darkImage, darkVector] = readDarkEstimate(darkFile, rawSize, channelIDs, opts)
    darkImage = [];
    darkVector = [];
    if exist(darkFile, 'file') ~= 2
        error('estimateCenterPointISMWavefront:MissingDarkFile', ...
            'Dark-count file was not found: %s', darkFile);
    end

    darkOpts = opts;
    darkOpts.darkMode = 'none';
    [darkStack, darkMeta] = readPtuDetectorStack(darkFile, darkOpts, [], false);
    if ~isempty(darkMeta.channelIDs)
        darkStack = alignStackChannels(darkStack, darkMeta.channelIDs, channelIDs);
    end

    if isequal(size(darkStack, 1), rawSize(1)) && ...
            isequal(size(darkStack, 2), rawSize(2)) && ...
            isequal(size(darkStack, 3), rawSize(3))
        darkImage = reshape(darkStack, size(darkStack,1), size(darkStack,2), size(darkStack,3), 1);
        darkImage = repmat(darkImage, 1, 1, 1, rawSize(4));
    else
        darkVector = squeeze(median(median(max(darkStack, 0), 1), 2));
        darkVector = double(darkVector(:));
    end
end

function aligned = alignStackChannels(stack, sourceIDs, targetIDs)
    if isempty(sourceIDs) || isempty(targetIDs) || isequal(sourceIDs(:), targetIDs(:))
        aligned = stack;
        return;
    end
    if numel(sourceIDs) == numel(targetIDs) && isequal(targetIDs(:), (1:numel(targetIDs)).')
        aligned = stack;
        return;
    end
    [present, loc] = ismember(double(targetIDs(:)), double(sourceIDs(:)));
    if ~all(present)
        error('estimateCenterPointISMWavefront:DarkChannelMismatch', ...
            'The dark PTU does not contain all detector channels used by the data.');
    end
    aligned = stack(:,:,loc,:);
end

function bg = boundaryBackground(data, widthPx)
    widthPx = max(1, min(round(widthPx), floor(min(size(data,1), size(data,2))/2)));
    bg = zeros(1, 1, size(data, 3), size(data, 4));
    for ip = 1:size(data, 4)
        for ch = 1:size(data, 3)
            img = data(:,:,ch,ip);
            border = [reshape(img(1:widthPx,:), [], 1); ...
                      reshape(img(end-widthPx+1:end,:), [], 1); ...
                      reshape(img(:,1:widthPx), [], 1); ...
                      reshape(img(:,end-widthPx+1:end), [], 1)];
            border = border(isfinite(border));
            if isempty(border)
                bg(1,1,ch,ip) = 0;
            else
                bg(1,1,ch,ip) = median(border);
            end
        end
    end
end

function [values, info] = extractCenterDetectorValues(raw4, opts)
    focusImage = sum(max(raw4(:,:,:,1), 0), 3);
    selectionMethod = 'summed focus image';
    [centerXY, centerYX, centerFit] = estimateCenterXY(focusImage, opts);
    actualSampleMode = opts.centerSampleMode;
    if ~isempty(opts.centerXY)
        selectionMethod = 'explicit centerXY option';
    elseif isfield(centerFit, 'method')
        selectionMethod = centerFit.method;
    end
    values = sampleDetectorStackAtXY(raw4, centerXY, actualSampleMode);

    info = struct();
    info.focusImage = focusImage;
    info.centerXY = centerXY;
    info.centerPixelYX = centerYX;
    info.centerOffsetFromRoundedXY = centerXY - [centerYX(2), centerYX(1)];
    info.centerMode = opts.centerMode;
    info.centerSampleMode = actualSampleMode;
    info.selectionMethod = selectionMethod;
    info.centerFit = centerFit;
end

function [centerXY, centerYX, fitInfo] = estimateCenterXY(img, opts)
    img = double(img);
    fitInfo = struct('method', 'explicit centerXY option');
    if ~isempty(opts.centerXY)
        centerXY = double(opts.centerXY(:)).';
        if numel(centerXY) ~= 2
            error('estimateCenterPointISMWavefront:BadCenterXY', ...
                'centerXY must be [x y].');
        end
    else
        switch opts.centerMode
            case 'peak'
                [~, idx] = max(img(:));
                [cy, cx] = ind2sub(size(img), idx);
                centerXY = [cx, cy];
                fitInfo.method = 'summed focus image peak';
            case {'centroid','weighted'}
                [centerXY, fitInfo] = centerOfMassXY(img, opts.centerThresholdFraction);
            case {'gaussian','gaussiancom','gaussian_com','comgaussian'}
                [centerXY, fitInfo] = gaussianCenterFromCom(img, opts.centerThresholdFraction);
            otherwise
                error('estimateCenterPointISMWavefront:BadCenterMode', ...
                    'centerMode must be ''gaussian'', ''centroid'', or ''peak''.');
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
        info = struct('method', 'summed focus image peak fallback');
    else
        [yy, xx] = ndgrid(1:size(img,1), 1:size(img,2));
        mass = sum(weights(:));
        centerXY = [sum(xx(:).*weights(:)), sum(yy(:).*weights(:))] / mass;
        info = struct('method', 'summed focus image center of mass', ...
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
        info.method = 'summed focus image center of mass fallback';
        return;
    end

    sigma0 = max(1.5, radius / 2);
    p0 = [bg0, amp0, cx0, cy0, log(sigma0), log(sigma0)];
    objective = @(p) gaussianSse(p, xx, yy, patch);
    options = optimset('Display', 'off', 'MaxIter', 200, 'MaxFunEvals', 800);
    p = fminsearch(objective, p0, options);
    centerXY = [min(max(p(3), 1), nx), min(max(p(4), 1), ny)];
    info = struct('method', 'summed focus image Gaussian fit from COM seed', ...
        'centerOfMassXY', seedXY, 'fitSigmaXY', exp(p(5:6)), ...
        'fitAmplitude', p(2), 'fitBackground', p(1), ...
        'fitWindowYX', [ylo yhi xlo xhi]);
end

function sse = gaussianSse(p, xx, yy, patch)
    sx = max(exp(p(5)), 0.5);
    sy = max(exp(p(6)), 0.5);
    model = p(1) + p(2) * exp(-0.5 * (((xx - p(3)) / sx).^2 + ((yy - p(4)) / sy).^2));
    residual = model - patch;
    sse = sum(residual(:).^2);
    if ~isfinite(sse)
        sse = realmax;
    end
end

function values = sampleDetectorStackAtXY(raw4, centerXY, mode)
    nCh = size(raw4, 3);
    nPlane = size(raw4, 4);

    if size(raw4, 1) == 1 && size(raw4, 2) == 1
        values = squeeze(raw4(1, 1, :, :));
        values = reshape(double(values), nCh, nPlane);
        return;
    end

    x = min(max(double(centerXY(1)), 1), size(raw4, 2));
    y = min(max(double(centerXY(2)), 1), size(raw4, 1));
    ix = min(max(round(x), 1), size(raw4, 2));
    iy = min(max(round(y), 1), size(raw4, 1));

    switch lower(char(mode))
        case {'nearest','round','pixel'}
            values = squeeze(raw4(iy, ix, :, :));
            values = reshape(double(values), nCh, nPlane);
        case {'subpixel','linear','bilinear','interp'}
            values = zeros(nCh, nPlane);
            for ip = 1:nPlane
                for ch = 1:nCh
                    img = double(raw4(:,:,ch,ip));
                    v = interp2(img, x, y, 'linear', NaN);
                    if ~isfinite(v)
                        v = img(iy, ix);
                    end
                    values(ch, ip) = v;
                end
            end
        otherwise
            error('estimateCenterPointISMWavefront:BadCenterSampleMode', ...
                'centerSampleMode must be ''subpixel'' or ''nearest''.');
    end
end

function [dataN, info] = normalizeCenterValues(values, mode)
    values = max(double(values), 0);
    dataN = values;
    info = struct();
    info.mode = mode;
    info.scale = [];

    switch lower(mode)
        case {'perplane','plane','eachplane'}
            scale = sum(values, 1);
            scale(scale <= 0 | ~isfinite(scale)) = 1;
            dataN = values ./ reshape(scale, 1, []);
            info.scale = scale;
        case {'global','all'}
            scale = sum(values(:));
            if scale <= 0 || ~isfinite(scale)
                scale = 1;
            end
            dataN = values / scale;
            info.scale = scale;
        case {'none','raw'}
            info.scale = ones(1, size(values, 2));
        otherwise
            error('estimateCenterPointISMWavefront:BadNormalization', ...
                'centerNormalization must be ''perPlane'', ''global'', or ''none''.');
    end
end

function [sim, diagnostics] = configureCenterPointSim(raw4, inputMeta, opts)
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
    if ~isempty(opts.detectorSubsamples)
        sim.detectorSubsamples = opts.detectorSubsamples;
    end

    sim = applyOpticsOptions(sim, inputMeta, opts);
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
    validateScalarOptics(sim);
    sim.beadSubsamples = opts.beadSubsamples;
    airyUnitUm = opts.airyUnitUm;
    if isempty(airyUnitUm) || ~isfinitePositiveScalar(airyUnitUm)
        airyUnitUm = 1.22 * sim.lamEm / effectivePropagatingNALocal(sim);
    end
    sim.airyUnitUm = airyUnitUm;
    sim.airyUnitDefinition = ...
        'Airy diameter using propagating sample-side NA';

    xyPixelSizeUm = opts.xyPixelSizeUm;
    if isempty(xyPixelSizeUm) || ~isfinitePositiveScalar(xyPixelSizeUm)
        xyPixelSizeUm = inputMeta.xyPixelSizeUm;
    end
    if size(raw4, 1) > 1 && size(raw4, 2) > 1 && isfinitePositiveScalar(xyPixelSizeUm)
        ny = size(raw4, 1);
        nx = size(raw4, 2);
        sim.fovX = xyPixelSizeUm * (nx - 1);
        sim.fovY = xyPixelSizeUm * (ny - 1);
        sim.fovXY = max(sim.fovX, sim.fovY);
        sim.xyPixelSizeUm = xyPixelSizeUm;
        sim.nx = nx;
        sim.ny = ny;
        sim.x = linspace(-sim.fovX/2, sim.fovX/2, nx);
        sim.y = linspace(-sim.fovY/2, sim.fovY/2, ny);
        if nx > 1
            sim.dx = abs(sim.x(2) - sim.x(1));
        end
    end

    if ~isempty(opts.detectorXYUm)
        detXY = validateDetectorXY(opts.detectorXYUm, size(raw4, 3));
        sim = applyDetectorXY(sim, detXY);
        diagnostics = explicitDetectorLayoutDiagnostics( ...
            detXY, xyPixelSizeUm, airyUnitUm);
    elseif opts.estimateDetectorLayout && size(raw4, 1) > 1 && ...
            size(raw4, 2) > 1 && isfinitePositiveScalar(xyPixelSizeUm)
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
        layoutOpts.planeIndex = 1;
        layoutOpts.airyUnitUm = airyUnitUm;
        [detXY, diagnostics] = estimateDetectorLayoutFromStack( ...
            raw4(:,:,:,1), xyPixelSizeUm, layoutOpts);
        sim = applyDetectorXY(sim, detXY);
        diagnostics.source = 'raw focus-plane ISM phase correlation';
    else
        diagnostics = fixedDetectorLayoutDiagnostics( ...
            sim.detXY, xyPixelSizeUm, airyUnitUm, sim.detectorLayout);
    end
    if isstruct(diagnostics)
        diagnostics.detectorPitchSource = detectorPitchSourceLocal(opts);
    end

    calibratedPitch = resolveDetectorPitchSampleUmLocal(opts);
    if isempty(opts.detectorXYUm) && isfinite(calibratedPitch) && ...
            isfield(sim,'detXY')
        recoveredPitch = medianNearestDistance(sim.detXY);
        if isfinitePositiveScalar(recoveredPitch)
            detXY = sim.detXY*(calibratedPitch/recoveredPitch);
            sim = applyDetectorXY(sim,detXY);
            if isstruct(diagnostics)
                diagnostics.detXYBeforePitchCalibrationUm = diagnostics.detXYUm;
                diagnostics.detXYUm = detXY;
                diagnostics.detXYNm = 1000*detXY;
                diagnostics.pitchCalibrationScale = ...
                    calibratedPitch/recoveredPitch;
                diagnostics.calibratedPitchUm = calibratedPitch;
                diagnostics.detectorPitchSource = detectorPitchSourceLocal(opts);
            end
        end
    end
    sim.detectorPhysicalPitchCalibrated = isfinite(calibratedPitch) || ...
        ~isempty(opts.detectorXYUm);
    sim.detectorCalibratedPitchUm = calibratedPitch;

    sim = prepareCenterPointZGrid(sim, opts.planeZ, opts);
    if ~strcmpi(sim.sampleGeometry,'airOnGlass')
        sim.obj = beadObject3D(sim);
    end
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

function sim = applyDetectorPitchScaleToSim(sim, scale)
    if isempty(scale) || ~isnumeric(scale) || ~isscalar(scale) || ...
            ~isfinite(scale) || scale <= 0 || abs(scale - 1) < eps
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

function detXY = validateDetectorXY(detXY, nCh)
    detXY = double(detXY);
    if ~ismatrix(detXY) || size(detXY, 2) ~= 2 || size(detXY, 1) ~= nCh || ...
            any(~isfinite(detXY(:)))
        error('estimateCenterPointISMWavefront:BadDetectorXY', ...
            'detectorXYUm must be a finite [%d x 2] array in sample-equivalent um.', nCh);
    end
end

function diagnostics = explicitDetectorLayoutDiagnostics(detXY, pixelSizeUm, airyUnitUm)
    diagnostics = struct();
    diagnostics.source = 'explicit detectorXYUm option';
    diagnostics.pixelSizeUm = pixelSizeUm;
    diagnostics.pixelSizeNm = 1000 * pixelSizeUm;
    diagnostics.airyUnitUm = airyUnitUm;
    diagnostics.airyUnitDefinition = 'Airy disk diameter: 1 AU = 1.22*lambda_em/NA';
    diagnostics.positionSign = NaN;
    diagnostics.detectorScale = NaN;
    diagnostics.centerMode = 'explicit';
    [~, diagnostics.centerDetectorIndex] = min(sum(detXY.^2, 2));
    diagnostics.signal = nan(size(detXY,1), 1);
    diagnostics.peakValue = nan(size(detXY,1), 1);
    diagnostics.shiftsPx = nan(size(detXY));
    diagnostics.shiftsUm = nan(size(detXY));
    diagnostics.shiftsNm = nan(size(detXY));
    diagnostics.detectorOffsetFromShiftUm = detXY;
    diagnostics.detectorOffsetFromShiftNm = 1000 * detXY;
    diagnostics.detXYUm = detXY;
    diagnostics.detXYNm = 1000 * detXY;
    diagnostics.shiftsAU = nan(size(detXY));
    diagnostics.shiftMagnitudeAU = nan(size(detXY,1), 1);
    diagnostics.detectorOffsetFromShiftAU = detXY / airyUnitUm;
    diagnostics.detXYAU = detXY / airyUnitUm;
    diagnostics.detectorRadiusAU = hypot( ...
        diagnostics.detXYAU(:,1), diagnostics.detXYAU(:,2));
end

function diagnostics = fixedDetectorLayoutDiagnostics(detXY, pixelSizeUm, airyUnitUm, layoutName)
    diagnostics = explicitDetectorLayoutDiagnostics(detXY, pixelSizeUm, airyUnitUm);
    diagnostics.source = sprintf('fixed regular %s detector layout', char(layoutName));
    diagnostics.centerMode = 'fixed';
end

function sim = applyOpticsOptions(sim, inputMeta, opts)
    sim.opticsSource = 'defaultParams';
    if isfield(inputMeta, 'head') && isstruct(inputMeta.head)
        head = inputMeta.head;
        na = firstHeaderNumericByRegex(head, ...
            {'(^|_)NA($|_)', 'Objective.*NA', 'Obj.*NA', 'Numerical.*Aperture'}, 'na');
        magnification = firstHeaderNumericByRegex(head, ...
            {'Magnification', 'Objective.*Mag', 'Obj.*Mag', '(^|_)Mag($|_)'}, 'magnification');
        lamExc = firstHeaderNumericByRegex(head, ...
            {'Exc.*Wave', 'Excitation.*Wave', 'Laser.*Wave', 'Lambda.*Exc', 'Wavelength.*Exc'}, 'wavelength');
        lamEm = firstHeaderNumericByRegex(head, ...
            {'Em.*Wave', 'Emission.*Wave', 'Detection.*Wave', 'Lambda.*Em', 'Wavelength.*Em'}, 'wavelength');

        if isfinitePositiveScalar(na)
            sim.NA = na;
            sim.opticsSource = 'PTU header';
        end
        if isfinitePositiveScalar(magnification)
            sim.objectiveMagnification = magnification;
            sim.opticsSource = 'PTU header';
        end
        if isfinitePositiveScalar(lamExc)
            sim.lamExc = lamExc;
            sim.opticsSource = 'PTU header';
        end
        if isfinitePositiveScalar(lamEm)
            sim.lamEm = lamEm;
            sim.lamRef = lamEm;
            sim.opticsSource = 'PTU header';
        end
    end

    if isfinitePositiveScalar(opts.objectiveNA)
        sim.NA = double(opts.objectiveNA);
        sim.opticsSource = 'explicit options';
    end
    if isfinitePositiveScalar(opts.objectiveMagnification)
        sim.objectiveMagnification = double(opts.objectiveMagnification);
        sim.opticsSource = 'explicit options';
    end
    if isfinitePositiveScalar(opts.excitationWavelengthUm)
        sim.lamExc = double(opts.excitationWavelengthUm);
        sim.opticsSource = 'explicit options';
    end
    if isfinitePositiveScalar(opts.emissionWavelengthUm)
        sim.lamEm = double(opts.emissionWavelengthUm);
        sim.lamRef = double(opts.emissionWavelengthUm);
        sim.opticsSource = 'explicit options';
    end
    if isfinitePositiveScalar(opts.mediumRefractiveIndex)
        sim.nMedium = double(opts.mediumRefractiveIndex);
        sim.opticsSource = 'explicit options';
    end
end

function sim = prepareCenterPointZGrid(sim, planeZ, opts)
    dzTarget = opts.modelDz;
    zPadding = opts.modelZPadding;
    zMin = min([planeZ(:); 0]) - zPadding;
    zMax = max([planeZ(:); 0]) + zPadding;
    nZ = ceil((zMax - zMin) / dzTarget) + 1;
    if mod(nZ, 2) == 0
        nZ = nZ + 1;
    end
    sim.z = linspace(zMin, zMax, nZ);
    sim.nz = numel(sim.z);
    sim.nzRange = zMax - zMin;
end

function [fit, strategy] = fitCenterPointModel(dataN, sim, planeZ, opts)
    strategy = struct();
    switch opts.fitStrategy
        case {'joint','twoplane','two-plane'}
            fit = fitCenterPointSingleStart(dataN, sim, planeZ, opts);
            strategy.name = 'joint';
        case {'focusanchored','focus-anchored','focus'}
            [fit, strategy] = fitFocusAnchoredModel(dataN, sim, planeZ, opts);
        otherwise
            error('estimateCenterPointISMWavefront:BadFitStrategy', ...
                'fitStrategy must be ''joint'' or ''focusAnchored''.');
    end
end

function [fit, strategy] = fitFocusAnchoredModel(dataN, sim, planeZ, opts)
    stageClock = tic;
    focusOpts = opts;
    focusOpts.fitStrategy = 'joint';
    focusOpts.planeWeights = [1 0];
    focusOpts.maxIter = min(opts.maxIter, opts.focusCoarseMaxIter);
    focusOpts.verbose = false;

    starts = focusInitialVectors(sim, focusOpts);
    if size(starts, 1) > opts.focusMaxStarts
        starts = starts(1:opts.focusMaxStarts, :);
    end
    nStarts = size(starts, 1);
    startIndex = (1:nStarts).';
    focusResidualNorm = nan(nStarts, 1);
    focusWeightedResidualNorm = nan(nStarts, 1);
    coefficientNorm = nan(nStarts, 1);
    convergedIterations = nan(nStarts, 1);
    forwardEvaluations = nan(nStarts, 1);
    startFits = cell(nStarts, 1);
    focusSim = sim;
    focusSim.detectorSubsamples = min(detectorSubsampleCount(sim), ...
        opts.focusCoarseDetectorSubsamples);
    if opts.verbose
        fprintf(['[estimateCenterPointISMWavefront] fast focus initialization: ' ...
            '%d starts, %d iterations, %d detector subsamples.\n'], ...
            nStarts, focusOpts.maxIter, focusSim.detectorSubsamples);
    end

    for k = 1:nStarts
        oneOpts = optionsWithInitialVector(focusOpts, focusSim, starts(k,:));
        oneFit = fitCenterPointSingleStart(dataN, focusSim, planeZ, oneOpts);
        startFits{k} = oneFit;
        focusResidualNorm(k) = norm(oneFit.residual(:,1));
        focusWeightedResidualNorm(k) = oneFit.weightedResidualNorm;
        coefficientNorm(k) = norm(oneFit.paramVector(1:numel(opts.fitModes)));
        convergedIterations(k) = size(oneFit.history, 1);
        forwardEvaluations(k) = oneFit.nForwardEvaluations;
        if opts.verbose
            fprintf('  focus start %d/%d: residual %.4e\n', ...
                k, nStarts, focusResidualNorm(k));
        end
    end

    score = focusWeightedResidualNorm + 1e-10 * coefficientNorm;
    [~, bestStart] = min(score);
    coarseFocusFit = startFits{bestStart};
    refineOpts = optionsWithInitialVector(focusOpts, sim, ...
        coarseFocusFit.paramVector);
    refineOpts.maxIter = min(opts.maxIter, opts.focusRefineMaxIter);
    if opts.verbose
        fprintf('[estimateCenterPointISMWavefront] refining best focus start at full detector sampling.\n');
    end
    focusFit = fitCenterPointSingleStart(dataN, sim, planeZ, refineOpts);
    focusDiagOpts = focusOpts;
    focusDiagOpts.planeWeights = [1 0];
    focusSufficiency = centerPointIdentifiability(sim, opts.fitModes, planeZ, ...
        focusFit.paramVector, focusDiagOpts);

    signModes = intersect(opts.signSelectionModes, opts.fitModes, 'stable');
    signCandidates = buildSignCandidates(focusFit.paramVector, opts.fitModes, signModes);
    nSign = size(signCandidates, 1);
    diversityResidualNorm = nan(nSign, 1);
    totalResidualNorm = nan(nSign, 1);
    signForwardEvaluations = nan(nSign, 1);
    signLabel = cell(nSign, 1);

    diversityOpts = opts;
    diversityOpts.fitStrategy = 'joint';
    diversityOpts.planeWeights = [0 1];
    diversityOpts.maxIter = 0;
    diversityOpts.verbose = false;
    if opts.verbose
        fprintf('[estimateCenterPointISMWavefront] testing %d even-mode sign combinations at %.0f nm.\n', ...
            nSign, 1000*abs(planeZ(2)-planeZ(1)));
    end
    for k = 1:nSign
        oneOpts = optionsWithInitialVector(diversityOpts, sim, signCandidates(k,:));
        oneFit = fitCenterPointSingleStart(dataN, sim, planeZ, oneOpts);
        diversityResidualNorm(k) = norm(oneFit.residual(:,2));
        totalResidualNorm(k) = oneFit.residualNorm;
        signForwardEvaluations(k) = oneFit.nForwardEvaluations;
        signLabel{k} = signCandidateLabel(signCandidates(k,:), opts.fitModes, signModes);
    end
    [~, bestSign] = min(diversityResidualNorm);
    selectedP = signCandidates(bestSign,:);
    sortedSignResidual = sort(diversityResidualNorm);
    if numel(sortedSignResidual) >= 2
        signResidualMargin = sortedSignResidual(2) - sortedSignResidual(1);
        signResidualRatio = sortedSignResidual(2) / max(sortedSignResidual(1), eps);
    else
        signResidualMargin = NaN;
        signResidualRatio = NaN;
    end

    finalOpts = optionsWithInitialVector(opts, sim, selectedP);
    finalOpts.fitStrategy = 'joint';
    finalOpts.planeWeights = [opts.focusWeight opts.diversityWeight];
    if opts.verbose
        fprintf('[estimateCenterPointISMWavefront] final focus-dominant refinement.\n');
    end
    fit = fitCenterPointSingleStart(dataN, sim, planeZ, finalOpts);
    fit.strategy = 'focusAnchored';
    fit.focusOnlyFit = focusFit;
    fit.selectedSignIndex = bestSign;

    strategy = struct();
    strategy.name = 'focusAnchored';
    strategy.focusPlaneWeights = [1 0];
    strategy.finalPlaneWeights = fit.planeWeights;
    strategy.nCoarseStarts = nStarts;
    strategy.coarseMaxIter = focusOpts.maxIter;
    strategy.refineMaxIter = refineOpts.maxIter;
    strategy.coarseDetectorSubsamples = focusSim.detectorSubsamples;
    strategy.nForwardEvaluations = sum(forwardEvaluations) + ...
        focusFit.nForwardEvaluations + sum(signForwardEvaluations) + ...
        fit.nForwardEvaluations;
    strategy.elapsedSeconds = toc(stageClock);
    strategy.focusFit = focusFit;
    strategy.focusSufficiency = focusSufficiency;
    strategy.selectedSignIndex = bestSign;
    strategy.selectedSignLabel = signLabel{bestSign};
    strategy.signResidualMargin = signResidualMargin;
    strategy.signResidualRatio = signResidualRatio;
    strategy.focusStartTable = table(startIndex, focusResidualNorm, ...
        focusWeightedResidualNorm, coefficientNorm, convergedIterations, ...
        forwardEvaluations);
    isSelected = false(nSign, 1);
    isSelected(bestSign) = true;
    strategy.signSelectionTable = table((1:nSign).', signLabel, isSelected, ...
        diversityResidualNorm, totalResidualNorm, ...
        'VariableNames', {'candidateIndex','signLabel','isSelected', ...
        'diversityResidualNorm','totalResidualNorm'});
end

function fit = fitCenterPointSingleStart(dataN, sim, planeZ, opts)
    if isempty(opts.fdXY)
        opts.fdXY = sim.dx / 4;
    end
    planeWeights = resolvedPlaneWeights(opts.planeWeights);

    paramNames = buildParamNames(opts.fitModes, opts);
    p = initialParameterVector(sim, opts);
    step = finiteDifferenceSteps(numel(opts.fitModes), opts, sim);
    reg = regularizationVector(numel(opts.fitModes), opts);
    maxStep = maxUpdateVector(numel(opts.fitModes), opts);
    nForwardEvaluations = 0;

    history = zeros(opts.maxIter, 2);
    for it = 1:opts.maxIter
        m0 = modelCenterValues(sim, opts.fitModes, planeZ, p, opts);
        nForwardEvaluations = nForwardEvaluations + 1;
        r = weightedResidualVector(dataN, m0, planeWeights);
        J = zeros(numel(m0), numel(p));

        for q = 1:numel(p)
            pp = p;
            pp(q) = pp(q) + step(q);
            mp = modelCenterValues(sim, opts.fitModes, planeZ, pp, opts);
            nForwardEvaluations = nForwardEvaluations + 1;
            if strcmp(opts.jacobianScheme, 'central')
                pm = p;
                pm(q) = pm(q) - step(q);
                mm = modelCenterValues(sim, opts.fitModes, planeZ, pm, opts);
                nForwardEvaluations = nForwardEvaluations + 1;
                J(:,q) = weightedResidualVector(mp, mm, planeWeights) / ...
                    (2*step(q));
            else
                J(:,q) = weightedResidualVector(mp, m0, planeWeights) / step(q);
            end
        end

        H = J.'*J + diag(reg);
        g = J.'*r;
        if any(~isfinite(H(:))) || any(~isfinite(g(:)))
            error('estimateCenterPointISMWavefront:NonFiniteSystem', ...
                'Non-finite values encountered in center-point Gauss-Newton system.');
        end

        if rcond(H) < 1e-12
            delta = pinv(H) * g;
        else
            delta = H \ g;
        end

        delta = max(min(delta(:).', maxStep), -maxStep);
        p = clampCenterPointParams(p + delta, opts);

        history(it,1) = norm(r);
        history(it,2) = norm(delta);
        if opts.verbose
            fprintf('[estimateCenterPointISMWavefront] iter %02d  residual %.4e  step %.4e\n', ...
                it, history(it,1), history(it,2));
        end

        if norm(delta) < opts.tolStep
            history = history(1:it,:);
            break;
        end
    end

    modelN = modelCenterValues(sim, opts.fitModes, planeZ, p, opts);
    nForwardEvaluations = nForwardEvaluations + 1;
    [estCoeffs, estXY, estZOffset, estDetectorPitchScale] = ...
        unpackParams(sim, opts.fitModes, p, opts);
    fitSim = applyDetectorPitchScaleToSim(sim, estDetectorPitchScale);

    fit = struct();
    fit.dataN = dataN;
    fit.modelN = modelN;
    fit.residual = dataN - modelN;
    fit.residualNorm = norm(fit.residual(:));
    fit.weightedResidualNorm = norm(weightedResidualVector(dataN, modelN, planeWeights));
    fit.planeResidualNorm = sqrt(sum(fit.residual.^2, 1));
    fit.planeWeights = planeWeights;
    fit.fitModes = opts.fitModes;
    fit.paramNames = paramNames;
    fit.paramVector = p;
    fit.estCoeffs = estCoeffs;
    fit.estCoeffVector = coeffStructToVector(sim, estCoeffs);
    fit.estXY = estXY;
    fit.estZOffset = estZOffset;
    fit.estDetectorPitchScale = estDetectorPitchScale;
    fit.estDetectorPitchUm = numericField(fitSim, 'detPitch', NaN);
    fit.estDetectorTotalMagnification = detectorTotalMagnificationFromPitch( ...
        fit.estDetectorPitchUm, opts);
    fit.history = history;
    fit.nForwardEvaluations = nForwardEvaluations;
    fit.sim = fitSim;
    fit.nominalSim = sim;
    fit.nominalDetectorPitchUm = numericField(sim, 'detPitch', NaN);
    fit.phase = zernikePhaseMap(sim, estCoeffs, sim.lamEm);
end

function starts = focusInitialVectors(sim, opts)
    base = initialParameterVector(sim, opts);
    nModes = numel(opts.fitModes);
    amp = opts.focusSeedAmplitude;
    if amp <= 0 || nModes == 0
        starts = base;
        return;
    end

    starts = zeros(0, numel(base));
    patterns = [ones(1,nModes); -ones(1,nModes)];
    if nModes > 1
        alternating = ones(1,nModes);
        alternating(2:2:end) = -1;
        patterns = [patterns; alternating; -alternating]; %#ok<AGROW>
    end

    signIdx = find(ismember(opts.fitModes, opts.signSelectionModes));
    if ~isempty(signIdx)
        for k = 1:numel(signIdx)
            pattern = zeros(1,nModes);
            pattern(signIdx(k)) = 1;
            patterns = [patterns; pattern]; %#ok<AGROW>
            patterns = [patterns; -pattern]; %#ok<AGROW>
        end
    end

    patterns = unique(patterns, 'rows', 'stable');
    for k = 1:size(patterns,1)
        one = base;
        one(1:nModes) = one(1:nModes) + amp * patterns(k,:);
        starts(end+1,:) = one; %#ok<AGROW>
    end
end

function candidates = buildSignCandidates(p, fitModes, signModes)
    if isempty(signModes)
        candidates = p;
        return;
    end
    idx = zeros(1, numel(signModes));
    for k = 1:numel(signModes)
        idx(k) = find(strcmp(fitModes, signModes{k}), 1, 'first');
    end
    combos = binarySignMatrix(numel(idx));
    candidates = repmat(p, size(combos,1), 1);
    magnitude = abs(p(idx));
    for k = 1:size(combos,1)
        candidates(k,idx) = combos(k,:) .* magnitude;
    end
end

function signs = binarySignMatrix(n)
    if n <= 0
        signs = zeros(1,0);
        return;
    end
    nRows = 2^n;
    signs = ones(nRows, n);
    for row = 0:nRows-1
        for col = 1:n
            signs(row+1,col) = 2*bitget(row, col) - 1;
        end
    end
end

function label = signCandidateLabel(p, fitModes, signModes)
    parts = cell(1, numel(signModes));
    for k = 1:numel(signModes)
        idx = find(strcmp(fitModes, signModes{k}), 1, 'first');
        if p(idx) < 0
            signChar = '-';
        else
            signChar = '+';
        end
        parts{k} = sprintf('%s=%s', signModes{k}, signChar);
    end
    if isempty(parts)
        label = 'no sign-selected modes';
    else
        label = strjoin(parts, ', ');
    end
end

function opts = optionsWithInitialVector(opts, sim, p)
    nModes = numel(opts.fitModes);
    coeffs = struct();
    for k = 1:nModes
        coeffs.(opts.fitModes{k}) = p(k);
    end
    opts.initialCoeffs = coeffs;
    next = nModes + 1;
    if opts.fitXY
        opts.initialXY = p(next:next+1);
        next = next + 2;
    end
    if opts.fitZ
        opts.initialZOffset = p(next);
        next = next + 1;
    end
    if opts.fitDetectorPitchScale
        opts.initialDetectorPitchScale = p(next);
    end
    if ~isfield(sim, 'modeOrder')
        error('estimateCenterPointISMWavefront:MissingModeOrder', ...
            'Simulation structure must contain modeOrder.');
    end
end

function r = weightedResidualVector(a, b, planeWeights)
    residual = double(a) - double(b);
    residual = residual .* reshape(sqrt(planeWeights), 1, []);
    r = residual(:);
end

function weights = resolvedPlaneWeights(weights)
    if isempty(weights)
        weights = [1 1];
    end
    weights = double(weights(:)).';
end

function n = detectorSubsampleCount(sim)
    n = 5;
    if isstruct(sim) && isfield(sim, 'detectorSubsamples') && ...
            isnumeric(sim.detectorSubsamples) && ...
            isscalar(sim.detectorSubsamples) && ...
            isfinite(sim.detectorSubsamples) && sim.detectorSubsamples >= 1
        n = round(double(sim.detectorSubsamples));
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
    if opts.fitZ
        p = [p double(opts.initialZOffset)];
    end
    if opts.fitDetectorPitchScale
        p = [p double(opts.initialDetectorPitchScale)];
    end
    p = clampCenterPointParams(p, opts);
end

function p = clampCenterPointParams(p, opts)
    if opts.fitDetectorPitchScale
        idx = numel(opts.fitModes) + 1;
        if opts.fitXY
            idx = idx + 2;
        end
        if opts.fitZ
            idx = idx + 1;
        end
        b = opts.detectorPitchScaleBounds;
        p(idx) = max(min(p(idx), b(2)), b(1));
    end
end

function names = buildParamNames(fitModes, opts)
    names = fitModes;
    if opts.fitXY
        names = [names {'x_shift','y_shift'}];
    end
    if opts.fitZ
        names = [names {'z_offset'}];
    end
    if opts.fitDetectorPitchScale
        names = [names {'detector_pitch_scale'}];
    end
end

function step = finiteDifferenceSteps(nModes, opts, sim)
    step = opts.fdCoeff * ones(1, nModes);
    if opts.fitXY
        fdXY = opts.fdXY;
        if isempty(fdXY) || ~isnumeric(fdXY) || ~isscalar(fdXY) || ...
                ~isfinite(fdXY) || fdXY <= 0
            fdXY = sim.dx / 4;
        end
        step = [step fdXY fdXY];
    end
    if opts.fitZ
        fdZ = opts.fdZ;
        if isempty(fdZ) || ~isnumeric(fdZ) || ~isscalar(fdZ) || ...
                ~isfinite(fdZ) || fdZ <= 0
            fdZ = max(sim.dx / 4, 0.01);
        end
        step = [step fdZ];
    end
    if opts.fitDetectorPitchScale
        step = [step opts.fdDetectorPitchScale];
    end
    step(step <= 0) = sim.dx / 4;
end

function reg = regularizationVector(nModes, opts)
    reg = opts.regCoeff * ones(1, nModes);
    if opts.fitXY
        reg = [reg opts.regXY opts.regXY];
    end
    if opts.fitZ
        reg = [reg opts.regZ];
    end
    if opts.fitDetectorPitchScale
        reg = [reg opts.regDetectorPitchScale];
    end
end

function maxStep = maxUpdateVector(nModes, opts)
    maxStep = opts.maxCoeffStep * ones(1, nModes);
    if opts.fitXY
        maxStep = [maxStep opts.maxXYStep opts.maxXYStep];
    end
    if opts.fitZ
        maxStep = [maxStep opts.maxZStep];
    end
    if opts.fitDetectorPitchScale
        maxStep = [maxStep opts.maxDetectorPitchScaleStep];
    end
end

function values = modelCenterValues(sim, fitModes, planeZ, p, opts)
    [coeffs, xy, zOffset, detectorPitchScale] = ...
        unpackParams(sim, fitModes, p, opts);
    simEval = applyDetectorPitchScaleToSim(sim, detectorPitchScale);
    activePlane = true(size(planeZ));
    planeWeights = resolvedPlaneWeights(opts.planeWeights);
    if isPerPlaneNormalization(opts.centerNormalization) && ...
            numel(planeWeights) == numel(planeZ) && any(planeWeights == 0)
        activePlane = planeWeights > 0;
    end

    if strcmpi(simEval.sampleGeometry,'airOnGlass')
        stack = normalizedStackAirInterfaceZPlanes( ...
            simEval,coeffs,planeZ(activePlane),xy(1),xy(2),zOffset);
    else
        stack = normalizedStackExplicitDetectorZPlanes( ...
            simEval,coeffs,planeZ(activePlane),xy(1),xy(2),zOffset);
    end

    activeValues = sampleModelStackAtXY(stack, simEval, opts.modelSampleXY, ...
        opts.centerSampleMode, nnz(activePlane));
    if ~isempty(opts.detectorQE)
        activeValues = activeValues .* reshape(opts.detectorQE, [], 1);
    end
    values = zeros(simEval.nDet, numel(planeZ));
    values(:,activePlane) = activeValues;
    values = normalizeCenterValues(values, opts.centerNormalization);
end

function tf = isPerPlaneNormalization(mode)
    tf = ismember(lower(char(mode)), {'perplane','plane','eachplane'});
end

function values = sampleModelStackAtXY(stack, sim, sampleXY, mode, nPlane)
    if nargin < 5
        nPlane = size(stack, 4);
    end

    x = double(sampleXY(1));
    y = double(sampleXY(2));
    x = min(max(x, min(sim.x(:))), max(sim.x(:)));
    y = min(max(y, min(sim.y(:))), max(sim.y(:)));

    switch lower(char(mode))
        case {'nearest','round','pixel'}
            [~, ix] = min(abs(sim.x - x));
            [~, iy] = min(abs(sim.y - y));
            values = squeeze(stack(iy, ix, :, :));
            values = reshape(double(values), sim.nDet, nPlane);
        case {'subpixel','linear','bilinear','interp'}
            values = zeros(sim.nDet, nPlane);
            for ip = 1:nPlane
                for ch = 1:sim.nDet
                    img = double(stack(:,:,ch,ip));
                    values(ch, ip) = interp2(sim.x, sim.y, img, x, y, 'linear', 0);
                end
            end
        otherwise
            error('estimateCenterPointISMWavefront:BadCenterSampleMode', ...
                'centerSampleMode must be ''subpixel'' or ''nearest''.');
    end
end

function [coeffs, xy, zOffset, detectorPitchScale] = unpackParams(sim, fitModes, p, opts)
    fullVec = zeros(1, numel(sim.modeOrder));
    for k = 1:numel(fitModes)
        idx = find(strcmp(sim.modeOrder, fitModes{k}), 1, 'first');
        if isempty(idx)
            error('estimateCenterPointISMWavefront:UnknownMode', ...
                'Unknown fit mode "%s".', fitModes{k});
        end
        fullVec(idx) = p(k);
    end
    coeffs = coeffStruct(sim, fullVec);
    next = numel(fitModes) + 1;
    xy = [0 0];
    if opts.fitXY
        xy = [p(next) p(next+1)];
        next = next + 2;
    end
    zOffset = 0;
    if opts.fitZ
        zOffset = p(next);
        next = next + 1;
    end
    detectorPitchScale = 1;
    if opts.fitDetectorPitchScale
        detectorPitchScale = p(next);
    end
end

function diag = centerPointIdentifiability(sim, fitModes, planeZ, p, opts)
    step = finiteDifferenceSteps(numel(fitModes), opts, sim);
    m0 = modelCenterValues(sim, fitModes, planeZ, p, opts);
    planeWeights = resolvedPlaneWeights(opts.planeWeights);
    J = zeros(numel(m0), numel(p));
    for q = 1:numel(p)
        pp = p;
        pp(q) = pp(q) + step(q);
        mp = modelCenterValues(sim, fitModes, planeZ, pp, opts);
        if strcmp(opts.jacobianScheme, 'central')
            pm = p;
            pm(q) = pm(q) - step(q);
            mm = modelCenterValues(sim, fitModes, planeZ, pm, opts);
            J(:,q) = weightedResidualVector(mp, mm, planeWeights) / ...
                (2*step(q));
        else
            J(:,q) = weightedResidualVector(mp, m0, planeWeights) / step(q);
        end
    end
    s = svd(J, 'econ');
    tol = max(size(J)) * eps(max(s));
    rankJ = sum(s > tol);

    diag = struct();
    if nnz(planeWeights > 0) == 1
        diag.label = 'center point, one weighted plane';
    else
        diag.label = 'center point, two weighted planes';
    end
    diag.nObservations = size(m0,1) * nnz(planeWeights > 0);
    diag.nParameters = numel(p);
    diag.planeWeights = planeWeights;
    diag.parameterNames = buildParamNames(fitModes, opts);
    diag.rank = rankJ;
    diag.isFullRank = rankJ == numel(p);
    diag.singularValues = s;
    diag.tolerance = tol;
    if rankJ > 0
        diag.conditionNumber = s(1) / s(rankJ);
        diag.minSingularValue = s(rankJ);
    else
        diag.conditionNumber = Inf;
        diag.minSingularValue = 0;
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

function outDir = resolveOutputDir(opts)
    if ~isempty(opts.outputDir)
        outDir = char(opts.outputDir);
        return;
    end
    thisDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(thisDir);
    outDir = fullfile(rootDir, 'output_matlab', 'center_point_ism_wavefront');
end

function writeCenterPointOutputs(result)
    outDir = result.outputDir;
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    save(fullfile(outDir, 'center_point_ism_wavefront_fit.mat'), 'result', '-v7.3');
    writetable(coefficientTable(result.fit), fullfile(outDir, 'center_point_wavefront_coefficients.csv'));
    writetable(centerIntensityTable(result), fullfile(outDir, 'center_point_detector_intensities.csv'));
    writetable(sufficiencyTable(result.sufficiency), fullfile(outDir, 'center_point_sufficiency_diagnostic.csv'));
    if ~isempty(result.detectorLayoutDiagnostics)
        writetable(detectorLayoutTable(result), ...
            fullfile(outDir, 'center_point_detector_geometry.csv'));
        d = result.detectorLayoutDiagnostics;
        if isfield(d, 'shiftsAU') && any(isfinite(d.shiftsAU(:)))
            writeDetectorShiftFigure(result, ...
                fullfile(outDir, 'center_point_detector_shift_maps_AU.png'));
        end
    end
    if isstruct(result.focusAnchored) && ...
            isfield(result.focusAnchored, 'focusStartTable')
        writetable(result.focusAnchored.focusStartTable, ...
            fullfile(outDir, 'focus_only_multistart_diagnostic.csv'));
    end
    if isstruct(result.focusAnchored) && ...
            isfield(result.focusAnchored, 'focusSufficiency')
        writetable(sufficiencyTable(result.focusAnchored.focusSufficiency), ...
            fullfile(outDir, 'focus_only_sufficiency_diagnostic.csv'));
    end
    if isstruct(result.focusAnchored) && ...
            isfield(result.focusAnchored, 'signSelectionTable')
        writetable(result.focusAnchored.signSelectionTable, ...
            fullfile(outDir, 'diversity_plane_sign_selection.csv'));
    end
    writeSummaryFigure(result, fullfile(outDir, 'center_point_ism_wavefront_fit_summary.png'));
end

function T = coefficientTable(fit)
    mode = fit.sim.modeOrder(:);
    coeffWavesRMS = fit.estCoeffVector(:);
    coeffNmRMS = coeffWavesRMS * fit.sim.lamRef * 1000;
    T = table(mode, coeffWavesRMS, coeffNmRMS);
end

function T = centerIntensityTable(result)
    nCh = size(result.centerValues, 1);
    detectorIndex = repmat((1:nCh).', 2, 1);
    planeIndex = [ones(nCh,1); 2*ones(nCh,1)];
    zUm = [repmat(result.planeZ(1), nCh, 1); repmat(result.planeZ(2), nCh, 1)];
    rawIntensity = [result.centerValues(:,1); result.centerValues(:,2)];
    normalizedIntensity = [result.normalizedCenterValues(:,1); result.normalizedCenterValues(:,2)];
    fittedIntensity = [result.fit.modelN(:,1); result.fit.modelN(:,2)];
    residual = normalizedIntensity - fittedIntensity;
    T = table(planeIndex, zUm, detectorIndex, rawIntensity, ...
        normalizedIntensity, fittedIntensity, residual);
end

function T = detectorLayoutTable(result)
    d = result.detectorLayoutDiagnostics;
    nCh = size(d.detXYUm, 1);
    detectorIndex = (1:nCh).';
    channelID = detectorIndex;
    if isfield(result.inputMeta, 'channelIDs') && numel(result.inputMeta.channelIDs) == nCh
        channelID = result.inputMeta.channelIDs(:);
    end
    signal = columnOrNaN(d, 'signal', nCh);
    phaseCorrPeak = columnOrNaN(d, 'peakValue', nCh);
    isCenterDetector = false(nCh, 1);
    if isfield(d, 'centerDetectorIndex') && isfinite(d.centerDetectorIndex)
        isCenterDetector(d.centerDetectorIndex) = true;
    end
    shiftXPixel = matrixColumnOrNaN(d, 'shiftsPx', 1, nCh);
    shiftYPixel = matrixColumnOrNaN(d, 'shiftsPx', 2, nCh);
    shiftXAU = matrixColumnOrNaN(d, 'shiftsAU', 1, nCh);
    shiftYAU = matrixColumnOrNaN(d, 'shiftsAU', 2, nCh);
    shiftMagnitudeAU = columnOrNaN(d, 'shiftMagnitudeAU', nCh);
    detectorXAU = matrixColumnOrNaN(d, 'detXYAU', 1, nCh);
    detectorYAU = matrixColumnOrNaN(d, 'detXYAU', 2, nCh);
    detectorRadiusAU = columnOrNaN(d, 'detectorRadiusAU', nCh);
    scanPixelSizeUm = repmat(d.pixelSizeUm, nCh, 1);
    scanPixelSizeNm = 1000 * scanPixelSizeUm;
    airyUnitUm = repmat(d.airyUnitUm, nCh, 1);
    detectorScale = repmat(d.detectorScale, nCh, 1);
    positionSign = repmat(d.positionSign, nCh, 1);
    detectorPitchUm = repmat(medianNearestDistance(d.detXYUm), nCh, 1);
    detectorPitchNm = 1000 * detectorPitchUm;
    T = table(detectorIndex, channelID, isCenterDetector, signal, phaseCorrPeak, ...
        shiftXPixel, shiftYPixel, shiftXAU, shiftYAU, shiftMagnitudeAU, ...
        detectorXAU, detectorYAU, detectorRadiusAU, airyUnitUm, ...
        scanPixelSizeUm, scanPixelSizeNm, detectorScale, positionSign, ...
        detectorPitchUm, detectorPitchNm);
end

function v = columnOrNaN(s, fieldName, n)
    v = nan(n, 1);
    if isfield(s, fieldName) && numel(s.(fieldName)) == n
        raw = s.(fieldName);
        v = double(raw(:));
    end
end

function v = matrixColumnOrNaN(s, fieldName, col, n)
    v = nan(n, 1);
    if isfield(s, fieldName) && size(s.(fieldName),1) == n && ...
            size(s.(fieldName),2) >= col
        raw = s.(fieldName);
        v = double(raw(:,col));
    end
end

function pitch = medianNearestDistance(xy)
    xy = double(xy);
    n = size(xy, 1);
    nearest = inf(n, 1);
    for i = 1:n
        delta = xy - xy(i,:);
        distance = hypot(delta(:,1), delta(:,2));
        distance(i) = inf;
        nearest(i) = min(distance);
    end
    nearest = nearest(isfinite(nearest) & nearest > 0);
    if isempty(nearest)
        pitch = NaN;
    else
        pitch = median(nearest);
    end
end

function T = sufficiencyTable(diag)
    label = {diag.label};
    nObservations = diag.nObservations;
    nParameters = diag.nParameters;
    rank = diag.rank;
    isFullRank = diag.isFullRank;
    conditionNumber = diag.conditionNumber;
    minSingularValue = diag.minSingularValue;
    focusWeight = diag.planeWeights(1);
    diversityWeight = diag.planeWeights(2);
    T = table(label, nObservations, nParameters, rank, isFullRank, ...
        conditionNumber, minSingularValue, focusWeight, diversityWeight);
end

function writeSummaryFigure(result, outFile)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [60 60 1250 760]);
    tl = tiledlayout(fig, 3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    planeNames = {'focus','defocus'};
    for ip = 1:2
        ax = nexttile(tl, ip);
        plotDetectorValues(ax, result.fit.sim, result.normalizedCenterValues(:,ip), false);
        title(ax, sprintf('Measured %s, z=%.3f um', planeNames{ip}, result.planeZ(ip)));
        colorbar(ax);

        ax = nexttile(tl, ip + 3);
        plotDetectorValues(ax, result.fit.sim, result.fit.modelN(:,ip), false);
        title(ax, sprintf('Fitted %s', planeNames{ip}));
        colorbar(ax);

        ax = nexttile(tl, ip + 6);
        plotDetectorValues(ax, result.fit.sim, result.fit.residual(:,ip), true);
        title(ax, sprintf('Residual %s', planeNames{ip}));
        colorbar(ax);
    end

    ax = nexttile(tl, 3);
    imagesc(ax, result.centerInfo.focusImage);
    axis(ax, 'image');
    colormap(ax, 'hot');
    hold(ax, 'on');
    plot(ax, result.centerInfo.centerXY(1), result.centerInfo.centerXY(2), ...
        'co', 'MarkerSize', 8, 'LineWidth', 1.2);
    hold(ax, 'off');
    title(ax, 'Focal bead image and selected scan pixel');
    xlabel(ax, 'x pixel');
    ylabel(ax, 'y pixel');
    colorbar(ax);

    ax = nexttile(tl, 6);
    bar(ax, result.fit.estCoeffVector);
    set(ax, 'XTick', 1:numel(result.fit.sim.modeOrder), ...
        'XTickLabel', result.fit.sim.modeOrder, 'XTickLabelRotation', 45);
    ylabel(ax, 'waves RMS');
    title(ax, 'Recovered coefficients');
    grid(ax, 'on');

    ax = nexttile(tl, 9);
    plot(ax, result.fit.history(:,1), '-o', 'LineWidth', 1.1);
    xlabel(ax, 'iteration');
    ylabel(ax, 'residual norm');
    title(ax, sprintf('Rank %d/%d, cond %.2g', result.sufficiency.rank, ...
        result.sufficiency.nParameters, result.sufficiency.conditionNumber));
    grid(ax, 'on');

    exportFigure(fig, outFile);
end

function writeDetectorShiftFigure(result, outFile)
    d = result.detectorLayoutDiagnostics;
    nCh = size(d.shiftsAU, 1);
    displayXY = regularDetectorDisplayXY(result.fit.sim, nCh);
    channelIDs = (1:nCh).';
    if isfield(result.inputMeta, 'channelIDs') && numel(result.inputMeta.channelIDs) == nCh
        channelIDs = result.inputMeta.channelIDs(:);
    end

    shiftX = d.shiftsAU(:,1);
    shiftY = d.shiftsAU(:,2);
    shiftMagnitude = d.shiftMagnitudeAU;
    signedLimit = max(abs([shiftX(:); shiftY(:)]), [], 'omitnan');
    magnitudeLimit = max(shiftMagnitude, [], 'omitnan');

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1350 450]);
    tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    values = {shiftX, shiftY, shiftMagnitude};
    labels = {'Shift to center, x [AU]', 'Shift to center, y [AU]', ...
        'Shift to center, magnitude [AU]'};
    for k = 1:3
        ax = nexttile(tl);
        plotDetectorHexMap(displayXY, values{k}, 'Parent', ax);
        if k < 3 && isfinite(signedLimit) && signedLimit > 0
            caxis(ax, [-signedLimit signedLimit]);
            colormap(ax, redBlueMapLocal(256));
        else
            if isfinite(magnitudeLimit) && magnitudeLimit > 0
                caxis(ax, [0 magnitudeLimit]);
            end
            colormap(ax, parula);
        end
        overlayDetectorLabels(ax, displayXY, values{k}, channelIDs);
        colorbar(ax);
        title(ax, labels{k});
    end
    title(tl, sprintf(['Focus-plane ISM shifts reused for both planes; ' ...
        '1 AU = %.4f um = 1.22 lambda_{em}/NA'], d.airyUnitUm), ...
        'FontWeight', 'bold');
    exportFigure(fig, outFile);
end

function overlayDetectorLabels(ax, displayXY, values, channelIDs)
    hold(ax, 'on');
    for k = 1:numel(values)
        text(ax, displayXY(k,1), displayXY(k,2), ...
            sprintf('%g\n%.3f', channelIDs(k), values(k)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 7, 'Color', [0 0 0], 'FontWeight', 'bold');
    end
    hold(ax, 'off');
end

function plotDetectorValues(ax, sim, values, signedMap)
    displayXY = regularDetectorDisplayXY(sim, numel(values));
    if exist('plotDetectorHexMap', 'file') == 2
        plotDetectorHexMap(displayXY, values, 'Parent', ax);
    else
        scatter(ax, displayXY(:,1), displayXY(:,2), 300, values, 'filled');
        axis(ax, 'equal');
        axis(ax, 'off');
    end
    if signedMap
        lim = max(abs(values(:)));
        if isfinite(lim) && lim > 0
            caxis(ax, [-lim lim]);
        end
        colormap(ax, 'parula');
    else
        colormap(ax, 'parula');
    end
end

function displayXY = regularDetectorDisplayXY(sim, nCh)
    layoutName = 'honeycomb23';
    if isstruct(sim) && isfield(sim, 'detectorLayout') && ~isempty(sim.detectorLayout)
        layoutName = sim.detectorLayout;
    end
    try
        [displayXY, ~] = detectorLayout(layoutName, 1);
    catch
        displayXY = [];
    end
    if size(displayXY,1) ~= nCh
        nCol = ceil(sqrt(nCh));
        nRow = ceil(nCh/nCol);
        displayXY = zeros(nCh,2);
        for k = 1:nCh
            row = floor((k-1)/nCol);
            col = mod(k-1,nCol);
            displayXY(k,:) = [col-(nCol-1)/2, (nRow-1)/2-row];
        end
    end
end

function cmap = redBlueMapLocal(n)
    if nargin < 1
        n = 256;
    end
    x = linspace(-1, 1, n).';
    cmap = [min(1, 1+x), 1-abs(x), min(1, 1-x)];
    cmap = max(0, min(1, cmap));
end

function exportFigure(fig, outFile)
    try
        exportgraphics(fig, outFile, 'Resolution', 180);
    catch
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, outFile, '-dpng', '-r180');
    end
    close(fig);
end

function printCenterPointSummary(result)
    fit = result.fit;
    fprintf('\nCenter-point ISM wavefront fit\n');
    fprintf('  correction readiness: NO (%s)\n',result.validity.recommendedUse);
    fprintf('  center sample: x=%.3f px, y=%.3f px (%s; nearest row %d, col %d)\n', ...
        result.centerInfo.centerXY(1), result.centerInfo.centerXY(2), ...
        result.centerInfo.centerSampleMode, ...
        result.centerInfo.centerPixelYX(1), result.centerInfo.centerPixelYX(2));
    fprintf('    center selection: %s\n', result.centerInfo.selectionMethod);
    fprintf('  model sample: x=%.4g um, y=%.4g um\n', ...
        result.options.modelSampleXY(1), result.options.modelSampleXY(2));
    fprintf('  planeZ: [%.4g %.4g] um\n', result.planeZ(1), result.planeZ(2));
    fprintf('  fit strategy: %s, plane weights [%.3g %.3g]\n', ...
        result.fitStrategy, fit.planeWeights(1), fit.planeWeights(2));
    if isstruct(result.focusAnchored) && ...
            isfield(result.focusAnchored, 'selectedSignLabel')
        fprintf('  selected symmetric-mode signs: %s\n', ...
            result.focusAnchored.selectedSignLabel);
        fprintf('  sign-selection residual ratio (2nd/best): %.4g\n', ...
            result.focusAnchored.signResidualRatio);
        fprintf('  staged forward-model evaluations: %d\n', ...
            result.focusAnchored.nForwardEvaluations);
        fprintf('  staged fit elapsed time: %.1f s\n', ...
            result.focusAnchored.elapsedSeconds);
    end
    fprintf('  normalization: %s\n', result.normalization.mode);
    fprintf('  scan grid: %d y x %d x pixels', ...
        size(result.rawData,1), size(result.rawData,2));
    if isfield(fit.sim, 'xyPixelSizeUm') && isfinitePositiveScalar(fit.sim.xyPixelSizeUm)
        fprintf(', %.2f nm/pixel (FOV %.3f x %.3f um)\n', ...
            1000*fit.sim.xyPixelSizeUm, ...
            max(fit.sim.x)-min(fit.sim.x), max(fit.sim.y)-min(fit.sim.y));
    else
        fprintf(', physical pixel size unavailable\n');
    end
    if ~isempty(result.detectorLayoutDiagnostics)
        d = result.detectorLayoutDiagnostics;
        pitchNm = 1000 * medianNearestDistance(d.detXYUm);
        fprintf('  detector geometry: %s\n', d.source);
        if isfinite(d.positionSign) && isfinite(d.detectorScale)
            fprintf('  scan pixel: %.2f nm; reassignment conversion: sign %+.0f, scale %.3g\n', ...
                d.pixelSizeNm, d.positionSign, d.detectorScale);
        else
            fprintf('  scan pixel: %.2f nm; detector positions are not inferred from scan-image shifts\n', ...
                d.pixelSizeNm);
        end
        fprintf('  Airy unit: %.4f um (diameter convention)\n', d.airyUnitUm);
        fprintf('  median detector pitch: %.3f AU\n', ...
            pitchNm / (1000*d.airyUnitUm));
        if isfield(d, 'detectorPitchSource')
            fprintf('  detector pitch source: %s\n', d.detectorPitchSource);
        end
        if isfield(d, 'shiftsAU') && any(isfinite(d.shiftsAU(:)))
            fprintf('  shifts estimated once from focus plane; center detector index %d\n', ...
                d.centerDetectorIndex);
        else
            fprintf('  regular detector center index %d\n', d.centerDetectorIndex);
        end
    else
        fprintf('  detector geometry: configured/ideal sim.detXY\n');
    end
    if result.options.fitXY
        fprintf('  fitted bead offset: x=%+.1f nm, y=%+.1f nm\n', ...
            1000*fit.estXY(1), 1000*fit.estXY(2));
    else
        fprintf('  fitted bead offset: disabled; center taken from detector-summed scan image\n');
    end
    if isfield(fit, 'estDetectorPitchScale')
        fprintf('  detector pitch scale: %.5g; pitch %.4g um', ...
            fit.estDetectorPitchScale, fit.estDetectorPitchUm);
        if isfinite(fit.estDetectorTotalMagnification)
            fprintf('; total magnification %.4g x', fit.estDetectorTotalMagnification);
        end
        fprintf('\n');
    end
    if isfield(fit.sim, 'objectiveMagnification')
        fprintf('  objective: %.3g x, NA %.3g\n', fit.sim.objectiveMagnification, fit.sim.NA);
    else
        fprintf('  objective NA: %.3g\n', fit.sim.NA);
    end
    fprintf('  residual norm: %.4e\n', fit.residualNorm);
    fprintf('  local rank: %d/%d, condition %.3g\n', ...
        result.sufficiency.rank, result.sufficiency.nParameters, ...
        result.sufficiency.conditionNumber);
    for k = 1:numel(fit.sim.modeOrder)
        fprintf('  %-10s %+8.4f waves (%+7.1f nm)\n', ...
            fit.sim.modeOrder{k}, fit.estCoeffVector(k), ...
            fit.estCoeffVector(k) * fit.sim.lamRef * 1000);
    end
    fprintf('  outputs: %s\n\n', result.outputDir);
end

function value = firstHeaderNumericByRegex(head, patterns, kind)
    value = NaN;
    if ~isstruct(head)
        return;
    end
    names = fieldnames(head);
    for p = 1:numel(patterns)
        for k = 1:numel(names)
            if isempty(regexpi(names{k}, patterns{p}, 'once'))
                continue;
            end
            v = headerNumericValue(head.(names{k}));
            if isPlausibleHeaderValue(v, kind)
                value = normalizeHeaderValue(v, kind);
                return;
            end
        end
    end
end

function v = headerNumericValue(raw)
    v = NaN;
    if isnumeric(raw) && ~isempty(raw)
        v = double(raw(1));
    elseif isstruct(raw)
        candidates = {'value','Value','data','Data','TagValue','tagValue'};
        for k = 1:numel(candidates)
            name = candidates{k};
            if isfield(raw, name) && isnumeric(raw.(name)) && ~isempty(raw.(name))
                v = double(raw.(name)(1));
                return;
            end
        end
    end
end

function tf = isPlausibleHeaderValue(v, kind)
    tf = isfinitePositiveScalar(v);
    if ~tf
        return;
    end
    switch lower(kind)
        case 'na'
            tf = v >= 0.1 && v <= 2.0;
        case 'wavelength'
            vu = normalizeHeaderValue(v, kind);
            tf = vu >= 0.3 && vu <= 1.1;
        case 'magnification'
            tf = v >= 1 && v <= 200;
    end
end

function v = normalizeHeaderValue(v, kind)
    if strcmpi(kind, 'wavelength') && v > 10
        v = v / 1000;
    end
end

function value = effectivePropagatingNALocal(sim)
    value = sim.NA;
    if isfield(sim,'sampleGeometry') && strcmpi(sim.sampleGeometry,'airOnGlass')
        value = min(sim.NA,sim.nSample);
    end
end

function pitch = resolveDetectorPitchSampleUmLocal(opts)
    pitch = NaN;
    if isnumeric(opts.detectorPitchSampleUm) && ...
            isscalar(opts.detectorPitchSampleUm) && ...
            isfinite(opts.detectorPitchSampleUm) && ...
            opts.detectorPitchSampleUm > 0
        pitch = double(opts.detectorPitchSampleUm);
        return;
    end
    physical = opts.detectorHardwarePitchUm;
    magnification = opts.detectorTotalMagnification;
    if isnumeric(physical) && isscalar(physical) && isfinite(physical) && ...
            physical > 0 && isnumeric(magnification) && ...
            isscalar(magnification) && isfinite(magnification) && ...
            magnification > 0
        pitch = double(physical)/double(magnification);
    end
end

function source = detectorPitchSourceLocal(opts)
    if isnumeric(opts.detectorPitchSampleUm) && ...
            isscalar(opts.detectorPitchSampleUm) && ...
            isfinite(opts.detectorPitchSampleUm) && ...
            opts.detectorPitchSampleUm > 0
        source = 'detectorPitchSampleUm';
    elseif isnumeric(opts.detectorHardwarePitchUm) && ...
            isscalar(opts.detectorHardwarePitchUm) && ...
            isfinite(opts.detectorHardwarePitchUm) && ...
            opts.detectorHardwarePitchUm > 0 && ...
            isnumeric(opts.detectorTotalMagnification) && ...
            isscalar(opts.detectorTotalMagnification) && ...
            isfinite(opts.detectorTotalMagnification) && ...
            opts.detectorTotalMagnification > 0
        source = 'detectorHardwarePitchUm / detectorTotalMagnification';
    else
        source = 'defaultParams.detPitch';
    end
end

function magnification = detectorTotalMagnificationFromPitch(pitchUm, opts)
    magnification = NaN;
    physical = opts.detectorHardwarePitchUm;
    if isnumeric(physical) && isscalar(physical) && isfinite(physical) && ...
            physical > 0 && isnumeric(pitchUm) && isscalar(pitchUm) && ...
            isfinite(pitchUm) && pitchUm > 0
        magnification = double(physical) / double(pitchUm);
    end
end

function weights = validatePlaneWeights(weights)
    if isempty(weights)
        return;
    end
    weights = double(weights(:)).';
    if numel(weights) ~= 2 || any(~isfinite(weights)) || ...
            any(weights < 0) || all(weights == 0)
        error('estimateCenterPointISMWavefront:BadPlaneWeights', ...
            'planeWeights must contain two finite nonnegative values, not both zero.');
    end
end

function qe = validateDetectorQE(qe, channelIDs)
    if isempty(qe)
        return;
    end
    qe = double(qe(:));
    nCh = numel(channelIDs);
    if numel(qe) ~= nCh || any(~isfinite(qe)) || any(qe <= 0)
        error('estimateCenterPointISMWavefront:BadDetectorQE', ...
            'detectorQE must contain %d finite positive relative gains.', nCh);
    end
    qe = qe / mean(qe);
end

function value = validateNonnegativeScalar(value, name)
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value < 0
        error('estimateCenterPointISMWavefront:BadScalarOption', ...
            '%s must be a finite nonnegative scalar.', name);
    end
    value = double(value);
end

function value = validatePositiveScalar(value, name)
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
        error('estimateCenterPointISMWavefront:BadPositiveScalar', ...
            '%s must be a finite positive scalar.', name);
    end
    value = double(value);
end

function bounds = validateDetectorPitchScaleBounds(bounds)
    if ~isnumeric(bounds) || numel(bounds) ~= 2 || ...
            any(~isfinite(bounds(:))) || any(bounds(:) <= 0)
        error('estimateCenterPointISMWavefront:BadDetectorPitchScaleBounds', ...
            'detectorPitchScaleBounds must be two finite positive values.');
    end
    bounds = sort(double(bounds(:))).';
    if bounds(1) == bounds(2)
        error('estimateCenterPointISMWavefront:BadDetectorPitchScaleBounds', ...
            'detectorPitchScaleBounds must span a nonzero range.');
    end
end

function value = validatePositiveInteger(value, name)
    value = validateNonnegativeInteger(value, name);
    if value < 1
        error('estimateCenterPointISMWavefront:BadIntegerOption', ...
            '%s must be at least 1.', name);
    end
end

function value = validateNonnegativeInteger(value, name)
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
            value < 0 || abs(value - round(value)) > 1e-12
        error('estimateCenterPointISMWavefront:BadIntegerOption', ...
            '%s must be a nonnegative integer.', name);
    end
    value = round(double(value));
end

function tf = isfinitePositiveScalar(v)
    tf = isnumeric(v) && isscalar(v) && isfinite(v) && v > 0;
end

function v = numericField(s, name, defaultValue)
    v = defaultValue;
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name))
        v = double(s.(name)(1));
    end
end
