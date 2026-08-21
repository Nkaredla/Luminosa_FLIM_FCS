function result = estimateTwoPlaneISMWavefront(focusInput, defocusInput, varargin)
%ESTIMATETWOPLANEISMWAVEFRONT Estimate low-order wavefront aberrations from two ISM planes.
%
%   result = estimateTwoPlaneISMWavefront(focusInput, defocusInput)
%
%   focusInput and defocusInput can be detector-resolved numeric stacks
%   [y x 23], PTU files, folders containing a PTU file, or MAT files
%   containing detector-resolved data. The function builds a two-plane raw
%   stack [y x 23 2], subtracts a dark-count PTU, fits the low-order Zernike
%   coefficients with phaseRetrieval3DBead, and reports a local numerical
%   rank diagnostic for the two-plane inverse problem.
%
%   A single PTU/MAT/numeric [y x 23 2] input is also accepted:
%
%       result = estimateTwoPlaneISMWavefront(twoPlaneInput, []);
%
%   For a single multi-frame PTU, pass two frame indices:
%
%       result = estimateTwoPlaneISMWavefront('scan.ptu', [], ...
%           'frameIndices', [1 2]);
%
%   The default axial diversity is [0 0.3] um, i.e. focus plus a 300 nm
%   plane. If the stage metadata or experiment truly uses 300 um, pass
%   'planeZ',[0 300], but that produces a very large axial model grid.
%
%   Common options:
%       'planeZ'        : recorded plane z positions in um, default [0 0.3]
%       'darkFile'      : dark-count PTU, default D:\Luminosa\Data\ISMdark_counts.ptu
%       'darkScale'     : multiplier applied to the dark image/counts, default 1
%       'channelIDs'    : PTU routing IDs for PDA-23, default 9:31
%       'channelOrder'  : optional explicit order by channel ID or index
%       'xyPixelSizeUm' : scan pixel size in um, otherwise PTU header/default
%       'fitModes'      : Zernike modes to fit
%       'outputDir'     : output directory; default output_matlab/two_plane_ism_wavefront
%
%   The 23 detector images in two planes are locally sufficient for the
%   configured fit when result.sufficiency.twoPlane.isFullRank is true.

    if nargin < 1 || isempty(focusInput)
        error('estimateTwoPlaneISMWavefront:MissingInput', ...
            'Provide focusInput, or a single [y x 23 2] two-plane input.');
    end
    if nargin < 2
        defocusInput = [];
    end

    opts = parseOptions(varargin{:});
    addRequiredPaths(opts);

    [raw4, inputMeta] = loadTwoPlaneInput(focusInput, defocusInput, opts);
    assertTwoPlaneDetectorStack(raw4);

    planeZ = opts.planeZ(:).';
    if numel(planeZ) ~= 2
        error('estimateTwoPlaneISMWavefront:PlaneZ', ...
            'This two-plane estimator requires exactly two planeZ values.');
    end
    warnIfLargePlaneSeparation(planeZ);

    simForDark = configureDetectorOnlySim(opts);
    [dark, corrected4] = applyDarkAndHotPixelCorrection(raw4, inputMeta.channelIDs, simForDark, opts);

    [raw4, squareInfo] = makeRaw4Square(raw4, opts.squareDataMode);
    [corrected4, correctedSquareInfo] = makeRaw4Square(corrected4, opts.squareDataMode);
    squareInfo.corrected = correctedSquareInfo;

    sim = configureFitSim(corrected4, inputMeta, opts);
    if opts.estimateDetectorLayout
        [sim, detectorLayoutDiagnostics] = overrideDetectorLayoutFromData(raw4, sim, inputMeta, opts);
    else
        detectorLayoutDiagnostics = [];
    end

    fitOpts = buildPhaseRetrievalOptions(sim, corrected4, opts);
    fit = phaseRetrieval3DBead(true, corrected4, planeZ, fitOpts);

    sufficiency = struct();
    if opts.diagnoseSufficiency
        sufficiency.twoPlane = localIdentifiabilityDiagnostic( ...
            fit.sim, opts.fitModes, planeZ, fit.paramVector, opts, 'two planes');
        sufficiency.focusPlaneOnly = localIdentifiabilityDiagnostic( ...
            fit.sim, opts.fitModes, planeZ(1), fit.paramVector, opts, 'focus plane only');
    end

    result = struct();
    result.rawData = raw4;
    result.correctedData = corrected4;
    result.planeZ = planeZ;
    result.dark = dark;
    result.inputMeta = inputMeta;
    result.squareInfo = squareInfo;
    result.fit = fit;
    result.sufficiency = sufficiency;
    result.detectorLayoutDiagnostics = detectorLayoutDiagnostics;
    result.options = opts;
    result.outputDir = resolveOutputDir(opts);

    if opts.writeOutputs
        writeTwoPlaneOutputs(result);
    end

    if opts.verbose
        printTwoPlaneSummary(result);
    end
end

function opts = parseOptions(varargin)
    p = inputParser;
    p.FunctionName = 'estimateTwoPlaneISMWavefront';

    addParameter(p, 'planeZ', [0 0.3]);
    addParameter(p, 'darkFile', 'D:\Luminosa\Data\ISMdark_counts.ptu');
    addParameter(p, 'darkScale', 1);
    addParameter(p, 'darkMode', 'auto');
    addParameter(p, 'darkFrameIndices', []);
    addParameter(p, 'darkPerChannel', []);
    addParameter(p, 'useBoundaryDarkFallback', true);

    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'channelOrder', []);
    addParameter(p, 'ptuReaderFolder', defaultPtuReaderFolder());
    addParameter(p, 'ptuPhotonsPerChunk', 1e6);
    addParameter(p, 'frameIndices', []);
    addParameter(p, 'focusFrameIndices', []);
    addParameter(p, 'defocusFrameIndices', []);
    addParameter(p, 'frameCombine', 'sum');

    addParameter(p, 'matVariable', '');
    addParameter(p, 'squareDataMode', 'crop');
    addParameter(p, 'xyPixelSizeUm', []);
    addParameter(p, 'objectiveNA', []);
    addParameter(p, 'objectiveMagnification', []);
    addParameter(p, 'emissionWavelengthUm', []);
    addParameter(p, 'excitationWavelengthUm', []);
    addParameter(p, 'mediumRefractiveIndex', []);
    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detectorPixelShape', 'hex');
    addParameter(p, 'detectorSubsamples', []);
    addParameter(p, 'estimateDetectorLayout', false);
    addParameter(p, 'detectorLayoutPositionSign', -1);
    addParameter(p, 'detectorLayoutScale', 2);
    addParameter(p, 'detectorLayoutCenterMode', 'reference');
    addParameter(p, 'detectorCenterIndex', []);

    addParameter(p, 'fitModes', {'defocus','astig_x','astig_y','coma_x','coma_y','spherical'});
    addParameter(p, 'fitXY', true);
    addParameter(p, 'fitZ', false);
    addParameter(p, 'initialXY', []);
    addParameter(p, 'estimateInitialXY', true);
    addParameter(p, 'initialXYThresholdFraction', 0.20);
    addParameter(p, 'initialCoeffs', struct());
    addParameter(p, 'maxIter', 8);
    addParameter(p, 'modelDz', 0.05);
    addParameter(p, 'modelZPadding', 0.50);
    addParameter(p, 'fdCoeff', 0.01);
    addParameter(p, 'fdZ', 0.02);
    addParameter(p, 'regCoeff', 1e-5);
    addParameter(p, 'regXY', 1e-6);
    addParameter(p, 'regZ', 1e-6);
    addParameter(p, 'maxCoeffStep', 0.05);
    addParameter(p, 'maxXYStep', 0.05);
    addParameter(p, 'maxZStep', 0.05);
    addParameter(p, 'tolStep', 1e-4);

    addParameter(p, 'hotThresholdMAD', 8);
    addParameter(p, 'hotNeighborThresholdMAD', 6);
    addParameter(p, 'hotPixelAction', 'replaceWithNeighbors');
    addParameter(p, 'boundaryWidthPx', 3);

    addParameter(p, 'diagnoseSufficiency', true);
    addParameter(p, 'sufficiencyMaxSamples', 50000);
    addParameter(p, 'sufficiencyTolerance', []);

    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;

    if ischar(opts.fitModes) || isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    end
    opts.fitModes = opts.fitModes(:).';
    opts.hotPixelAction = lower(char(opts.hotPixelAction));
    opts.frameCombine = lower(char(opts.frameCombine));
    opts.darkMode = lower(char(opts.darkMode));
    opts.squareDataMode = lower(char(opts.squareDataMode));
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
            raw4 = standardizeTwoPlaneNumeric(focusInput);
            meta.source = 'numeric two-plane input';
            meta.channelIDs = defaultModelChannelIDs(size(raw4, 3));
            return;
        end

        keepFrames = ~isempty(opts.frameIndices);
        [stack, oneMeta] = loadDetectorInput(focusInput, opts, opts.frameIndices, keepFrames);
        if ndims(stack) == 3
            error('estimateTwoPlaneISMWavefront:NeedTwoPlanes', ...
                ['A single input produced one detector plane. Provide defocusInput, ' ...
                 'or pass ''frameIndices'', [focusFrame defocusFrame] for a multi-frame PTU.']);
        end
        raw4 = standardizeTwoPlaneNumeric(stack);
        meta = oneMeta;
        meta.source = 'single two-plane input';
        if numel(meta.channelIDs) ~= size(raw4, 3)
            meta.channelIDs = defaultModelChannelIDs(size(raw4, 3));
        end
        return;
    end

    if isnumeric(focusInput) && isnumeric(defocusInput)
        focusStack = standardizeSinglePlaneNumeric(focusInput);
        defocusStack = standardizeSinglePlaneNumeric(defocusInput);
        assertSamePlaneSize(focusStack, defocusStack);
        raw4 = cat(4, focusStack, defocusStack);
        meta.source = 'numeric focus and defocus inputs';
        meta.channelIDs = defaultModelChannelIDs(size(raw4, 3));
        return;
    end

    [focusStack, focusMeta] = loadDetectorInput(focusInput, opts, opts.focusFrameIndices, false);
    [defocusStack, defocusMeta] = loadDetectorInput(defocusInput, opts, opts.defocusFrameIndices, false);
    focusStack = standardizeSinglePlaneNumeric(focusStack);
    defocusStack = standardizeSinglePlaneNumeric(defocusStack);
    assertSamePlaneSize(focusStack, defocusStack);

    raw4 = cat(4, focusStack, defocusStack);
    meta = focusMeta;
    meta.source = 'focus and defocus inputs';
    meta.focus = focusMeta;
    meta.defocus = defocusMeta;
    if ~isempty(focusMeta.channelIDs) && ~isempty(defocusMeta.channelIDs) && ...
            ~isequal(focusMeta.channelIDs(:), defocusMeta.channelIDs(:))
        warning('estimateTwoPlaneISMWavefront:ChannelIDMismatch', ...
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

function ids = defaultModelChannelIDs(nCh)
    ids = (1:nCh).';
end

function [stack, meta] = loadDetectorInput(inputValue, opts, frameIndices, keepFrames)
    if isnumeric(inputValue)
        stack = inputValue;
        meta = emptyInputMeta();
        meta.source = 'numeric input';
        meta.channelIDs = defaultModelChannelIDs(size(stack, 3));
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
            error('estimateTwoPlaneISMWavefront:BadInputFile', ...
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
            error('estimateTwoPlaneISMWavefront:NoInputFile', ...
                'No PTU or MAT file found in %s.', inputPath);
        end
        [~, newest] = max([files.datenum]);
        fileName = fullfile(files(newest).folder, files(newest).name);
        return;
    end

    if exist(inputPath, 'file') ~= 2
        error('estimateTwoPlaneISMWavefront:MissingInputFile', ...
            'Input file was not found: %s', inputPath);
    end
    fileName = inputPath;
end

function [stack, meta] = readMatDetectorStack(fileName, opts)
    S = load(fileName);
    varName = char(opts.matVariable);
    if ~isempty(varName)
        if ~isfield(S, varName)
            error('estimateTwoPlaneISMWavefront:MissingMatVariable', ...
                'Variable "%s" was not found in %s.', varName, fileName);
        end
        stack = S.(varName);
    else
        stack = chooseMatDetectorVariable(S);
    end

    meta = emptyInputMeta();
    meta.source = 'MAT file';
    meta.file = fileName;
    meta.channelIDs = defaultModelChannelIDs(size(stack, 3));
    if isfield(S, 'channelIDs')
        meta.channelIDs = double(S.channelIDs(:));
    elseif isfield(S, 'dind')
        meta.channelIDs = double(S.dind(:));
    end
    if isfield(S, 'xyPixelSizeUm')
        meta.xyPixelSizeUm = double(S.xyPixelSizeUm);
    elseif isfield(S, 'opts') && isstruct(S.opts) && isfield(S.opts, 'xyPixelSizeUm')
        meta.xyPixelSizeUm = double(S.opts.xyPixelSizeUm);
    end
end

function stack = chooseMatDetectorVariable(S)
    preferred = {'raw4','rawData','detectorVolume','detectorStack', ...
        'channelVolume','channelVolumes','spadVolume','spadStack'};
    for k = 1:numel(preferred)
        name = preferred{k};
        if isfield(S, name) && isnumeric(S.(name)) && ndims(S.(name)) >= 3
            stack = S.(name);
            return;
        end
    end

    names = fieldnames(S);
    for k = 1:numel(names)
        value = S.(names{k});
        if isnumeric(value) && ndims(value) >= 3 && any(size(value) == 23)
            stack = value;
            return;
        end
    end

    error('estimateTwoPlaneISMWavefront:NoMatDetectorData', ...
        'No numeric detector-resolved variable was found in the MAT file.');
end

function [stack, meta] = readPtuDetectorStack(fileName, opts, frameIndices, keepFrames)
    if exist('PTU_MultiFrameScanReadFast', 'file') ~= 2 || exist('PTU_Read_Head', 'file') ~= 2
        error('estimateTwoPlaneISMWavefront:MissingPtuReader', ...
            'PTU reader functions are not on the MATLAB path.');
    end

    try
        ptuOut = PTU_MultiFrameScanReadFast(fileName, opts.ptuPhotonsPerChunk, false, false);
    catch fastErr
        try
            ptuOut = PTU_MultiFrameScanRead(fileName, opts.ptuPhotonsPerChunk);
        catch slowErr
            error('estimateTwoPlaneISMWavefront:PtuReadFailed', ...
                'Could not read %s as a detector scan. Fast: %s Slow: %s', ...
                fileName, fastErr.message, slowErr.message);
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

function [stackYX, channelIDs] = ptuOutputToDetectorStack(ptuOut, opts, frameIndices, keepFrames)
    if isfield(ptuOut, 'dind')
        channelIDs = double(ptuOut.dind(:));
    else
        channelIDs = [];
    end

    if keepFrames
        if ~isfield(ptuOut, 'tag') || isempty(ptuOut.tag)
            error('estimateTwoPlaneISMWavefront:NoPtuFrames', ...
                'The PTU output does not contain per-frame tag data.');
        end
        nFrame = size(ptuOut.tag, 4);
        frameIndices = validateFrameIndices(frameIndices, nFrame);
        stack = double(ptuOut.tag(:,:,:,frameIndices));
        stackYX = permute(stack, [2 1 3 4]);
    else
        if isempty(frameIndices)
            if isfield(ptuOut, 'tags') && ~isempty(ptuOut.tags)
                stack = double(ptuOut.tags);
            elseif isfield(ptuOut, 'tag') && ~isempty(ptuOut.tag)
                stack = sum(double(ptuOut.tag), 4);
            else
                error('estimateTwoPlaneISMWavefront:NoPtuImages', ...
                    'No detector image stack was found in the PTU output.');
            end
        else
            if ~isfield(ptuOut, 'tag') || isempty(ptuOut.tag)
                error('estimateTwoPlaneISMWavefront:NoPtuFrames', ...
                    'Frame selection requires ptuOut.tag.');
            end
            frameIndices = validateFrameIndices(frameIndices, size(ptuOut.tag, 4));
            switch opts.frameCombine
                case 'sum'
                    stack = sum(double(ptuOut.tag(:,:,:,frameIndices)), 4);
                case 'mean'
                    stack = mean(double(ptuOut.tag(:,:,:,frameIndices)), 4);
                otherwise
                    error('estimateTwoPlaneISMWavefront:BadFrameCombine', ...
                        'frameCombine must be ''sum'' or ''mean''.');
            end
        end
        stackYX = permute(stack, [2 1 3]);
    end

    if isempty(channelIDs)
        channelIDs = defaultModelChannelIDs(size(stackYX, 3));
    end

    [stackYX, channelIDs] = selectAndOrderChannels(stackYX, channelIDs, opts);
end

function frameIndices = validateFrameIndices(frameIndices, nFrame)
    if isempty(frameIndices)
        frameIndices = 1:nFrame;
    end
    frameIndices = round(frameIndices(:)).';
    if any(frameIndices < 1) || any(frameIndices > nFrame)
        error('estimateTwoPlaneISMWavefront:BadFrameIndices', ...
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
        elseif size(stack, 3) ~= numel(requested)
            error('estimateTwoPlaneISMWavefront:MissingChannels', ...
                'Could not find all requested channel IDs in PTU output.');
        end
    end

    order = opts.channelOrder;
    if ~isempty(order)
        order = double(order(:));
        if numel(order) ~= size(stack, 3)
            error('estimateTwoPlaneISMWavefront:BadChannelOrder', ...
                'channelOrder must contain one entry per selected detector channel.');
        end

        if all(ismember(order, channelIDs))
            [~, loc] = ismember(order, channelIDs);
        elseif all(order >= 1 & order <= size(stack, 3))
            loc = order;
        else
            error('estimateTwoPlaneISMWavefront:BadChannelOrder', ...
                'channelOrder must be either selected channel IDs or 1-based indices.');
        end
        stack = stack(:,:,loc,:);
        channelIDs = channelIDs(loc);
    end
end

function raw4 = standardizeTwoPlaneNumeric(raw)
    raw = double(raw);
    dims = size(raw);
    if ndims(raw) ~= 4
        error('estimateTwoPlaneISMWavefront:BadTwoPlaneData', ...
            'Two-plane numeric data must be [y x 23 2] or [y x 2 23].');
    end

    if dims(3) == 23 && dims(4) == 2
        raw4 = raw;
    elseif dims(3) == 2 && dims(4) == 23
        raw4 = permute(raw, [1 2 4 3]);
    else
        error('estimateTwoPlaneISMWavefront:BadTwoPlaneData', ...
            'Could not identify detector and two-plane dimensions in numeric data.');
    end
end

function stack = standardizeSinglePlaneNumeric(stack)
    stack = double(stack);
    if ndims(stack) == 4
        if size(stack, 4) == 1
            stack = stack(:,:,:,1);
        elseif size(stack, 3) == 1
            stack = permute(stack(:,:,1,:), [1 2 4 3]);
            stack = stack(:,:,:,1);
        else
            error('estimateTwoPlaneISMWavefront:BadSinglePlaneData', ...
                'A focus or defocus input must contain one detector plane.');
        end
    end
    if ndims(stack) ~= 3 || size(stack, 3) ~= 23
        error('estimateTwoPlaneISMWavefront:BadSinglePlaneData', ...
            'A focus or defocus detector stack must have size [y x 23].');
    end
end

function assertSamePlaneSize(a, b)
    if ~isequal(size(a), size(b))
        error('estimateTwoPlaneISMWavefront:PlaneSizeMismatch', ...
            'Focus and defocus stacks must have the same [y x 23] size.');
    end
end

function assertTwoPlaneDetectorStack(raw4)
    if ndims(raw4) ~= 4 || size(raw4, 3) ~= 23 || size(raw4, 4) ~= 2
        error('estimateTwoPlaneISMWavefront:BadRaw4', ...
            'The estimator requires a [y x 23 2] detector-resolved two-plane stack.');
    end
end

function warnIfLargePlaneSeparation(planeZ)
    dz = abs(diff(planeZ));
    if dz > 10
        warning('estimateTwoPlaneISMWavefront:LargePlaneSeparation', ...
            ['planeZ is in micrometers. A %.3g um separation will create a very large ' ...
             'model grid; use [0 0.3] for a 300 nm diversity plane.'], dz);
    end
end

function sim = configureDetectorOnlySim(opts)
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

function [dark, corrected] = applyDarkAndHotPixelCorrection(raw4, channelIDs, sim, opts)
    dark = struct();
    corrected = double(raw4);

    if strcmpi(opts.darkMode, 'none')
        dark.method = 'none';
        dark.darkPerChannel = zeros(size(raw4, 3), 1);
        dark.darkPerPixel = zeros(size(raw4, 3), 1);
    elseif ~isempty(opts.darkPerChannel)
        dark = darkStructFromChannelVector(opts.darkPerChannel(:), size(raw4), 'explicit darkPerChannel option');
        corrected = subtractDarkVector(raw4, dark.darkPerPixel, opts.darkScale);
    else
        [dark, corrected] = darkCorrectionFromFileOrBoundary(raw4, channelIDs, sim, opts);
    end

    [corrected, hot] = applyHotChannelAction(corrected, dark.darkPerChannel, sim, opts);
    dark.hotChannels = hot.hotChannels;
    dark.hotScoreGlobal = hot.hotScoreGlobal;
    dark.hotScoreNeighbor = hot.hotScoreNeighbor;
    dark.neighborMedian = hot.neighborMedian;
    dark.neighborMad = hot.neighborMad;
    dark.hotPixelAction = opts.hotPixelAction;
    dark.table = detectorDarkTable(sim, dark);
end

function [dark, corrected] = darkCorrectionFromFileOrBoundary(raw4, channelIDs, sim, opts)
    darkFile = char(opts.darkFile);
    if ~isempty(darkFile) && exist(darkFile, 'file') == 2
        try
            [dark, corrected] = darkCorrectionFromPtu(raw4, channelIDs, opts);
            return;
        catch err
            if ~opts.useBoundaryDarkFallback
                rethrow(err);
            end
            warning('estimateTwoPlaneISMWavefront:DarkFileFallback', ...
                'Dark PTU correction failed (%s). Falling back to image-boundary dark.', err.message);
        end
    elseif ~isempty(darkFile) && opts.verbose
        warning('estimateTwoPlaneISMWavefront:MissingDarkFile', ...
            'Dark PTU file was not found: %s', darkFile);
    end

    if opts.useBoundaryDarkFallback
        [dark, corrected] = estimateBoundaryDark(raw4, opts);
    else
        dark = darkStructFromChannelVector(zeros(size(raw4, 3), 1), size(raw4), 'none');
        corrected = raw4;
    end
end

function [dark, corrected] = darkCorrectionFromPtu(raw4, channelIDs, opts)
    darkFile = char(opts.darkFile);
    dark = struct();
    dark.file = darkFile;
    dark.scale = opts.darkScale;

    try
        darkOpts = opts;
        darkOpts.frameCombine = 'sum';
        [darkStack, darkMeta] = readPtuDetectorStack(darkFile, darkOpts, opts.darkFrameIndices, false);
        darkStack = alignStackChannels(darkStack, darkMeta.channelIDs, channelIDs);
        dark.method = 'detector dark image PTU';
        dark.darkImage = double(darkStack);
        dark.darkPerChannel = squeeze(sum(sum(max(double(darkStack), 0), 1), 2));
        dark.darkPerChannel = dark.darkPerChannel(:);

        if size(darkStack, 1) == size(raw4, 1) && size(darkStack, 2) == size(raw4, 2)
            corrected = raw4 - opts.darkScale * reshape(double(darkStack), ...
                size(darkStack,1), size(darkStack,2), size(darkStack,3), 1);
            dark.darkPerPixel = squeeze(mean(mean(max(double(darkStack), 0), 1), 2));
            dark.darkPerPixel = dark.darkPerPixel(:);
        else
            dark.darkPerPixel = dark.darkPerChannel / max(1, size(darkStack,1) * size(darkStack,2));
            corrected = subtractDarkVector(raw4, dark.darkPerPixel, opts.darkScale);
            dark.method = 'detector dark image PTU, channel-mean subtraction';
        end
        corrected = max(corrected, 0);
        return;
    catch imageErr
        [darkPerChannel, darkIDs, darkHead] = readPtuChannelCounts(darkFile, opts);
        darkPerChannel = alignVectorChannels(darkPerChannel, darkIDs, channelIDs);
        dark.method = sprintf('PTU photon channel counts (%s)', imageErr.identifier);
        dark.head = darkHead;
        dark.darkPerChannel = darkPerChannel(:);
        dark.darkPerPixel = dark.darkPerChannel / max(1, size(raw4,1) * size(raw4,2));
        corrected = subtractDarkVector(raw4, dark.darkPerPixel, opts.darkScale);
        corrected = max(corrected, 0);
    end
end

function corrected = subtractDarkVector(raw4, darkPerPixel, darkScale)
    corrected = double(raw4) - darkScale * reshape(darkPerPixel(:), 1, 1, [], 1);
    corrected = max(corrected, 0);
end

function [darkPerChannel, channelIDs, head] = readPtuChannelCounts(fileName, opts)
    if exist('PTU_Read_Head', 'file') ~= 2 || exist('PTU_Read', 'file') ~= 2
        error('estimateTwoPlaneISMWavefront:MissingPtuReader', ...
            'PTU_Read_Head/PTU_Read are not on the MATLAB path.');
    end

    head = PTU_Read_Head(fileName);
    validPhoton = ptuPhotonValidity(head);
    maxChannel = 64;
    counts = zeros(maxChannel, 1);

    cnt = 0;
    num = 1;
    while num > 0
        [~, tcspc, chan, special, num] = PTU_Read(fileName, [cnt+1 opts.ptuPhotonsPerChunk], head);
        cnt = cnt + num;
        if num == 0
            break;
        end

        keep = special == 0 & validPhoton(chan, tcspc);
        chan = double(chan(keep));
        chan = chan(chan >= 0 & chan < maxChannel);
        if ~isempty(chan)
            counts = counts + accumarray(chan(:)+1, 1, [maxChannel 1], @sum, 0);
        end
    end

    if ~isempty(opts.channelIDs)
        channelIDs = double(opts.channelIDs(:));
        idx = channelIDs + 1;
        if any(idx < 1) || any(idx > numel(counts))
            error('estimateTwoPlaneISMWavefront:BadDarkChannels', ...
                'Requested channel IDs are outside the PTU channel range.');
        end
        darkPerChannel = counts(idx);
    else
        idx = find(counts > 0);
        channelIDs = idx - 1;
        darkPerChannel = counts(idx);
    end
end

function validPhoton = ptuPhotonValidity(head)
    anzch = 32;
    rawResolution = numericField(head, 'MeasDesc_Resolution', NaN);
    globalResolution = numericField(head, 'MeasDesc_GlobalResolution', NaN);
    if isfinite(rawResolution) && rawResolution > 0 && ...
            isfinite(globalResolution) && globalResolution > 0
        resolutionNs = max(1e9 * rawResolution, 0.128);
        chDiv = max(1, round((resolutionNs * 1e-9) / rawResolution));
        nGate = min(1024, ceil(1e9 * globalResolution / resolutionNs) + 1);
        validPhoton = @(chan, tcspc) (chan < anzch) & (tcspc < nGate * chDiv);
    else
        validPhoton = @(chan, tcspc) chan < anzch; %#ok<INUSD>
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
        error('estimateTwoPlaneISMWavefront:DarkChannelMismatch', ...
            'The dark PTU does not contain all detector channels used by the data.');
    end
    aligned = stack(:,:,loc,:);
end

function aligned = alignVectorChannels(values, sourceIDs, targetIDs)
    if isempty(sourceIDs) || isempty(targetIDs) || isequal(sourceIDs(:), targetIDs(:))
        aligned = values;
        return;
    end

    if numel(sourceIDs) == numel(targetIDs) && isequal(targetIDs(:), (1:numel(targetIDs)).')
        aligned = values;
        return;
    end

    [present, loc] = ismember(double(targetIDs(:)), double(sourceIDs(:)));
    if ~all(present)
        error('estimateTwoPlaneISMWavefront:DarkChannelMismatch', ...
            'The dark PTU does not contain all detector channels used by the data.');
    end
    aligned = values(loc);
end

function [dark, corrected] = estimateBoundaryDark(raw4, opts)
    raw4 = double(raw4);
    [ny, nx, nCh, nz] = size(raw4);
    width = max(1, min([opts.boundaryWidthPx, floor(ny/4), floor(nx/4)]));
    borderMask = false(ny, nx);
    borderMask(1:width, :) = true;
    borderMask(end-width+1:end, :) = true;
    borderMask(:, 1:width) = true;
    borderMask(:, end-width+1:end) = true;

    darkPerPixel = zeros(nCh, 1);
    noiseMad = zeros(nCh, 1);
    for k = 1:nCh
        vals = [];
        for z = 1:nz
            frame = raw4(:,:,k,z);
            vals = [vals; frame(borderMask)]; %#ok<AGROW>
        end
        vals = vals(isfinite(vals));
        if isempty(vals)
            darkPerPixel(k) = 0;
            noiseMad(k) = 0;
        else
            darkPerPixel(k) = median(vals);
            noiseMad(k) = robustMad(vals);
        end
    end

    corrected = subtractDarkVector(raw4, darkPerPixel, 1);
    dark = darkStructFromChannelVector(darkPerPixel * ny * nx, size(raw4), ...
        'image-boundary per detector channel');
    dark.darkPerPixel = darkPerPixel;
    dark.noiseMadPerChannel = noiseMad * ny * nx;
    dark.boundaryWidthPx = width;
end

function dark = darkStructFromChannelVector(darkPerChannel, rawSize, method)
    dark = struct();
    dark.method = method;
    dark.darkPerChannel = double(darkPerChannel(:));
    dark.darkPerPixel = dark.darkPerChannel / max(1, rawSize(1) * rawSize(2));
    dark.noiseMadPerChannel = nan(size(dark.darkPerChannel));
    dark.globalMedian = median(dark.darkPerChannel, 'omitnan');
    dark.globalMad = robustMad(dark.darkPerChannel);
    dark.neighborMedian = nan(size(dark.darkPerChannel));
    dark.neighborMad = nan(size(dark.darkPerChannel));
    dark.hotChannels = [];
end

function [corrected, hot] = applyHotChannelAction(corrected, darkPerChannel, sim, opts)
    neighbors = detectorNeighborListLocal(sim.detXY, sim.detPitch);
    nCh = numel(darkPerChannel);
    neighborMedian = nan(nCh, 1);
    neighborMad = nan(nCh, 1);
    for k = 1:nCh
        nb = neighbors{k};
        if isempty(nb)
            nb = setdiff(1:nCh, k);
        end
        neighborMedian(k) = median(darkPerChannel(nb), 'omitnan');
        neighborMad(k) = robustMad(darkPerChannel(nb));
    end

    globalMedian = median(darkPerChannel, 'omitnan');
    globalMad = max(robustMad(darkPerChannel), eps);
    hotScoreGlobal = (darkPerChannel(:) - globalMedian) / globalMad;
    hotScoreNeighbor = (darkPerChannel(:) - neighborMedian) ./ max(neighborMad, eps);
    hotChannels = find(hotScoreGlobal > opts.hotThresholdMAD | ...
                       hotScoreNeighbor > opts.hotNeighborThresholdMAD);

    switch opts.hotPixelAction
        case {'none','subtractonly','subtract_only'}
        case {'replace','replacewithneighbors','replace_with_neighbors'}
            for k = hotChannels(:).'
                nb = neighbors{k};
                if isempty(nb)
                    nb = setdiff(1:nCh, k);
                end
                corrected(:,:,k,:) = median(corrected(:,:,nb,:), 3);
            end
        case {'exclude','zero'}
            corrected(:,:,hotChannels,:) = 0;
        otherwise
            error('estimateTwoPlaneISMWavefront:BadHotPixelAction', ...
                'Unknown hotPixelAction "%s".', opts.hotPixelAction);
    end

    hot = struct();
    hot.hotChannels = hotChannels;
    hot.hotScoreGlobal = hotScoreGlobal;
    hot.hotScoreNeighbor = hotScoreNeighbor;
    hot.neighborMedian = neighborMedian;
    hot.neighborMad = neighborMad;
end

function T = detectorDarkTable(sim, dark)
    channel = (1:numel(dark.darkPerChannel)).';
    detectorXUm = sim.detXY(:,1);
    detectorYUm = sim.detXY(:,2);
    darkEstimateCounts = dark.darkPerChannel(:);
    darkPerPixelCounts = dark.darkPerPixel(:);
    hotScoreGlobal = dark.hotScoreGlobal(:);
    hotScoreNeighbor = dark.hotScoreNeighbor(:);
    isHot = false(numel(channel), 1);
    if isfield(dark, 'hotChannels') && ~isempty(dark.hotChannels)
        isHot(dark.hotChannels) = true;
    end
    T = table(channel, detectorXUm, detectorYUm, darkEstimateCounts, ...
        darkPerPixelCounts, hotScoreGlobal, hotScoreNeighbor, isHot);
end

function neighbors = detectorNeighborListLocal(detXY, detPitch)
    nCh = size(detXY, 1);
    D = zeros(nCh, nCh);
    for i = 1:nCh
        for j = 1:nCh
            D(i,j) = hypot(detXY(i,1) - detXY(j,1), detXY(i,2) - detXY(j,2));
        end
    end
    radius = 1.15 * detPitch;
    neighbors = cell(nCh, 1);
    for k = 1:nCh
        neighbors{k} = find(D(k,:) > 0 & D(k,:) <= radius);
    end
end

function [out, info] = makeRaw4Square(raw4, mode)
    info = struct();
    info.mode = mode;
    info.originalSize = size(raw4);
    info.applied = false;

    ny = size(raw4, 1);
    nx = size(raw4, 2);
    if ny == nx
        out = raw4;
        info.outputSize = size(out);
        info.yIndex = 1:ny;
        info.xIndex = 1:nx;
        return;
    end

    switch lower(mode)
        case {'crop','centercrop','center_crop'}
            n = min(ny, nx);
            y0 = floor((ny - n) / 2) + 1;
            x0 = floor((nx - n) / 2) + 1;
            yIdx = y0:(y0+n-1);
            xIdx = x0:(x0+n-1);
            out = raw4(yIdx, xIdx, :, :);
            info.action = 'center crop';
            info.yIndex = yIdx;
            info.xIndex = xIdx;
        case {'pad','zeropad','zero_pad'}
            n = max(ny, nx);
            y0 = floor((n - ny) / 2) + 1;
            x0 = floor((n - nx) / 2) + 1;
            out = zeros(n, n, size(raw4, 3), size(raw4, 4));
            out(y0:(y0+ny-1), x0:(x0+nx-1), :, :) = raw4;
            info.action = 'zero pad';
            info.yIndex = y0:(y0+ny-1);
            info.xIndex = x0:(x0+nx-1);
        case {'error','strict'}
            error('estimateTwoPlaneISMWavefront:NonSquareData', ...
                'The forward model expects square y-x images. Use squareDataMode crop or pad.');
        otherwise
            error('estimateTwoPlaneISMWavefront:BadSquareDataMode', ...
                'Unknown squareDataMode "%s".', mode);
    end
    info.applied = true;
    info.outputSize = size(out);
end

function sim = configureFitSim(raw4, inputMeta, opts)
    sim = configureDetectorOnlySim(opts);
    sim = applyOpticsOptions(sim, inputMeta, opts);

    if ~isempty(opts.detectorSubsamples)
        sim.detectorSubsamples = opts.detectorSubsamples;
    end

    xyPixelSizeUm = opts.xyPixelSizeUm;
    if isempty(xyPixelSizeUm) || ~isfinitePositiveScalar(xyPixelSizeUm)
        xyPixelSizeUm = inputMeta.xyPixelSizeUm;
    end
    if isfinitePositiveScalar(xyPixelSizeUm)
        sim.fovXY = xyPixelSizeUm * (size(raw4, 1) - 1);
        sim.xyPixelSizeUm = xyPixelSizeUm;
    end

    sim.nx = size(raw4, 1);
    sim.x = linspace(-sim.fovXY/2, sim.fovXY/2, sim.nx);
    sim.y = sim.x;
    if sim.nx > 1
        sim.dx = abs(sim.x(2) - sim.x(1));
    else
        sim.dx = sim.fovXY;
    end
    sim.obj = beadObject3D(sim);

    if opts.verbose
        fprintf('[estimateTwoPlaneISMWavefront] fit grid: %d x %d, dx=%.4f um.\n', ...
            sim.nx, sim.nx, sim.dx);
    end
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

function tf = isfinitePositiveScalar(v)
    tf = isnumeric(v) && isscalar(v) && isfinite(v) && v > 0;
end

function [sim, diagnostics] = overrideDetectorLayoutFromData(raw4, sim, inputMeta, opts)
    xyPixelSizeUm = sim.dx;
    if isfield(inputMeta, 'xyPixelSizeUm') && isfinitePositiveScalar(inputMeta.xyPixelSizeUm)
        xyPixelSizeUm = inputMeta.xyPixelSizeUm;
    end

    layoutOpts = struct();
    layoutOpts.positionSign = opts.detectorLayoutPositionSign;
    layoutOpts.detectorScale = opts.detectorLayoutScale;
    layoutOpts.centerMode = opts.detectorLayoutCenterMode;
    layoutOpts.centerDetectorIndex = opts.detectorCenterIndex;
    layoutOpts.planeIndex = 1;
    layoutOpts.airyUnitUm = 1.22 * sim.lamEm / sim.NA;
    [detXY, diagnostics] = estimateDetectorLayoutFromStack( ...
        raw4(:,:,:,1), xyPixelSizeUm, layoutOpts);
    oldPitch = sim.detPitch;
    oldSize = sim.detSize;
    sim.detXY = detXY;
    sim.nDet = size(detXY, 1);
    recoveredPitch = median(nearestDistanceLocal(detXY));
    if isfinitePositiveScalar(recoveredPitch)
        sim.detPitch = recoveredPitch;
        sim.detSize = (oldSize / oldPitch) * recoveredPitch;
        if strcmpi(sim.detectorPixelShape, 'hex')
            sim.detectorHexRadius = sim.detSize / sqrt(3);
        end
    end
    [~, sim.detectorIndexGrid] = detectorLayout(sim.detectorLayout, 1);
    sim.detectorGridSize = size(sim.detectorIndexGrid);
    sim.arrayN = sim.detectorGridSize(1);
end

function fitOpts = buildPhaseRetrievalOptions(sim, data4, opts)
    fitOpts = struct();
    fitOpts.sim = sim;
    fitOpts.fitModes = opts.fitModes;
    fitOpts.maxIter = opts.maxIter;
    fitOpts.fitXY = opts.fitXY;
    fitOpts.fitZ = opts.fitZ;
    fitOpts.initialCoeffs = opts.initialCoeffs;
    fitOpts.modelDz = opts.modelDz;
    fitOpts.modelZPadding = opts.modelZPadding;
    fitOpts.fdCoeff = opts.fdCoeff;
    fitOpts.fdZ = opts.fdZ;
    fitOpts.regCoeff = opts.regCoeff;
    fitOpts.regXY = opts.regXY;
    fitOpts.regZ = opts.regZ;
    fitOpts.maxCoeffStep = opts.maxCoeffStep;
    fitOpts.maxXYStep = opts.maxXYStep;
    fitOpts.maxZStep = opts.maxZStep;
    fitOpts.tolStep = opts.tolStep;
    fitOpts.verbose = opts.verbose;

    if opts.fitXY
        if ~isempty(opts.initialXY)
            fitOpts.initialXY = double(opts.initialXY(:).');
        elseif opts.estimateInitialXY
            fitOpts.initialXY = estimateInitialXYFromData(data4, sim, opts);
        end
    end
end

function xy = estimateInitialXYFromData(data4, sim, opts)
    img = squeeze(sum(sum(max(data4, 0), 3), 4));
    img = double(img);
    img = img - min(img(:));
    peak = max(img(:));
    if ~isfinite(peak) || peak <= 0
        xy = [0 0];
        return;
    end

    weights = img;
    weights(img < opts.initialXYThresholdFraction * peak) = 0;
    sw = sum(weights(:));
    if sw <= 0
        weights = img;
        sw = sum(weights(:));
    end

    [X, Y] = meshgrid(sim.x, sim.y);
    xy = [sum(X(:).*weights(:))/sw, sum(Y(:).*weights(:))/sw];
end

function diag = localIdentifiabilityDiagnostic(sim, fitModes, planeZ, p, opts, label)
    fitXY = opts.fitXY;
    fitZ = opts.fitZ;
    pLocal = restrictParamVectorForPlaneSet(p, fitModes, fitXY, fitZ);
    step = finiteDifferenceSteps(numel(fitModes), sim, opts, fitXY, fitZ);

    m0 = normalizedModelFromParams(sim, fitModes, planeZ, pLocal, fitXY, fitZ);
    sampleIdx = chooseDiagnosticSamples(m0, opts.sufficiencyMaxSamples);
    J = zeros(numel(sampleIdx), numel(pLocal));

    for q = 1:numel(pLocal)
        pp = pLocal;
        pm = pLocal;
        pp(q) = pp(q) + step(q);
        pm(q) = pm(q) - step(q);
        mp = normalizedModelFromParams(sim, fitModes, planeZ, pp, fitXY, fitZ);
        mm = normalizedModelFromParams(sim, fitModes, planeZ, pm, fitXY, fitZ);
        d = (mp(:) - mm(:)) / (2*step(q));
        J(:,q) = d(sampleIdx);
    end

    s = svd(J, 'econ');
    if isempty(opts.sufficiencyTolerance)
        tol = max(size(J)) * eps(max(s));
    else
        tol = opts.sufficiencyTolerance;
    end
    rankJ = sum(s > tol);

    diag = struct();
    diag.label = label;
    diag.nObservationsUsed = numel(sampleIdx);
    diag.nModelPixels = numel(m0);
    diag.nParameters = numel(pLocal);
    diag.parameterNames = buildParamNames(fitModes, fitXY, fitZ);
    diag.rank = rankJ;
    diag.isFullRank = rankJ == numel(pLocal);
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

function pLocal = restrictParamVectorForPlaneSet(p, fitModes, fitXY, fitZ)
    nExpected = numel(fitModes) + 2*double(fitXY) + double(fitZ);
    if numel(p) < nExpected
        error('estimateTwoPlaneISMWavefront:BadParamVector', ...
            'Fit parameter vector has fewer values than expected.');
    end
    pLocal = p(1:nExpected);
end

function names = buildParamNames(fitModes, fitXY, fitZ)
    names = fitModes;
    if fitXY
        names = [names {'x_shift','y_shift'}];
    end
    if fitZ
        names = [names {'z_offset'}];
    end
end

function step = finiteDifferenceSteps(nModes, sim, opts, fitXY, fitZ)
    step = opts.fdCoeff * ones(1, nModes);
    if fitXY
        step = [step sim.dx/4 sim.dx/4];
    end
    if fitZ
        step = [step opts.fdZ];
    end
end

function stack = normalizedModelFromParams(sim, fitModes, planeZ, p, fitXY, fitZ)
    [coeffs, xy, zOffset] = unpackParams(sim, fitModes, p, fitXY, fitZ);
    stack = normalizedStackExplicitDetectorZPlanes( ...
        sim, coeffs, planeZ, xy(1), xy(2), zOffset);
end

function [coeffs, xy, zOffset] = unpackParams(sim, fitModes, p, fitXY, fitZ)
    fullVec = zeros(1, numel(sim.modeOrder));
    for k = 1:numel(fitModes)
        idx = find(strcmp(sim.modeOrder, fitModes{k}), 1, 'first');
        if isempty(idx)
            error('estimateTwoPlaneISMWavefront:UnknownMode', ...
                'Unknown fit mode "%s".', fitModes{k});
        end
        fullVec(idx) = p(k);
    end
    coeffs = coeffStruct(sim, fullVec);
    next = numel(fitModes) + 1;
    xy = [0 0];
    if fitXY
        xy = [p(next) p(next+1)];
        next = next + 2;
    end
    zOffset = 0;
    if fitZ
        zOffset = p(next);
    end
end

function idx = chooseDiagnosticSamples(modelStack, maxSamples)
    n = numel(modelStack);
    maxSamples = max(1, min(n, round(maxSamples)));
    if n <= maxSamples
        idx = (1:n).';
        return;
    end

    uniformCount = maxSamples;
    idx = unique(round(linspace(1, n, uniformCount))).';

    nPeak = min(floor(maxSamples/4), n);
    if nPeak > 0
        [~, order] = sort(modelStack(:), 'descend');
        idx = unique([idx; order(1:nPeak)]);
        if numel(idx) > maxSamples
            idx = idx(1:maxSamples);
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
    outDir = fullfile(rootDir, 'output_matlab', 'two_plane_ism_wavefront');
end

function writeTwoPlaneOutputs(result)
    outDir = result.outputDir;
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    save(fullfile(outDir, 'two_plane_ism_wavefront_fit.mat'), 'result', '-v7.3');
    writetable(coefficientTable(result.fit), fullfile(outDir, 'estimated_wavefront_coefficients.csv'));
    writetable(fittedOffsetTable(result), fullfile(outDir, 'fitted_bead_offset.csv'));
    if ~isempty(result.detectorLayoutDiagnostics)
        writetable(detectorLayoutTable(result), fullfile(outDir, 'estimated_detector_pixel_shifts.csv'));
    end
    if isfield(result.dark, 'table')
        writetable(result.dark.table, fullfile(outDir, 'dark_counts_from_ptu.csv'));
    end
    if isfield(result.sufficiency, 'twoPlane')
        writetable(sufficiencyTable(result.sufficiency), ...
            fullfile(outDir, 'two_plane_sufficiency_diagnostic.csv'));
    end
    writeSummaryFigure(result, fullfile(outDir, 'two_plane_ism_wavefront_fit_summary.png'));
end

function T = coefficientTable(fit)
    mode = fit.sim.modeOrder(:);
    coeffWavesRMS = fit.estCoeffVector(:);
    coeffNmRMS = coeffWavesRMS * fit.sim.lamRef * 1000;
    T = table(mode, coeffWavesRMS, coeffNmRMS);
end

function T = fittedOffsetTable(result)
    fit = result.fit;
    xShiftUm = fit.estXY(1);
    yShiftUm = fit.estXY(2);
    xShiftNm = 1000 * xShiftUm;
    yShiftNm = 1000 * yShiftUm;
    xyPixelSizeUm = NaN;
    if isfield(fit.sim, 'xyPixelSizeUm') && isfinitePositiveScalar(fit.sim.xyPixelSizeUm)
        xyPixelSizeUm = fit.sim.xyPixelSizeUm;
    elseif isfield(fit.sim, 'dx') && isfinitePositiveScalar(fit.sim.dx)
        xyPixelSizeUm = fit.sim.dx;
    end
    xyPixelSizeNm = 1000 * xyPixelSizeUm;
    xShiftPixels = xShiftUm / xyPixelSizeUm;
    yShiftPixels = yShiftUm / xyPixelSizeUm;
    zOffsetUm = fit.estZOffset;
    zOffsetNm = 1000 * zOffsetUm;
    objectiveNA = fit.sim.NA;
    objectiveMagnification = NaN;
    if isfield(fit.sim, 'objectiveMagnification')
        objectiveMagnification = fit.sim.objectiveMagnification;
    end
    T = table(xShiftUm, yShiftUm, xShiftNm, yShiftNm, ...
        xShiftPixels, yShiftPixels, xyPixelSizeUm, xyPixelSizeNm, ...
        zOffsetUm, zOffsetNm, objectiveNA, objectiveMagnification);
end

function T = detectorLayoutTable(result)
    d = result.detectorLayoutDiagnostics;
    nCh = numel(d.signal);
    detectorIndex = (1:nCh).';
    channelID = detectorIndex;
    if isfield(result.inputMeta, 'channelIDs') && numel(result.inputMeta.channelIDs) == nCh
        channelID = result.inputMeta.channelIDs(:);
    end
    signal = d.signal(:);
    phaseCorrPeak = d.peakValue(:);
    shiftXPixel = d.shiftsPx(:,1);
    shiftYPixel = d.shiftsPx(:,2);
    shiftXAU = d.shiftsAU(:,1);
    shiftYAU = d.shiftsAU(:,2);
    shiftMagnitudeAU = d.shiftMagnitudeAU(:);
    detectorXAU = d.detXYAU(:,1);
    detectorYAU = d.detXYAU(:,2);
    detectorRadiusAU = d.detectorRadiusAU(:);
    isCenterDetector = false(nCh,1);
    isCenterDetector(d.centerDetectorIndex) = true;
    pixelSizeUm = repmat(d.pixelSizeUm, nCh, 1);
    pixelSizeNm = repmat(d.pixelSizeNm, nCh, 1);
    airyUnitUm = repmat(d.airyUnitUm, nCh, 1);
    T = table(detectorIndex, channelID, isCenterDetector, signal, phaseCorrPeak, ...
        shiftXPixel, shiftYPixel, shiftXAU, shiftYAU, shiftMagnitudeAU, ...
        detectorXAU, detectorYAU, detectorRadiusAU, airyUnitUm, ...
        pixelSizeUm, pixelSizeNm);
end

function T = sufficiencyTable(sufficiency)
    fields = {'focusPlaneOnly','twoPlane'};
    label = {};
    nObservationsUsed = [];
    nModelPixels = [];
    nParameters = [];
    rank = [];
    isFullRank = [];
    conditionNumber = [];
    minSingularValue = [];
    for k = 1:numel(fields)
        if ~isfield(sufficiency, fields{k})
            continue;
        end
        d = sufficiency.(fields{k});
        label{end+1,1} = d.label; %#ok<AGROW>
        nObservationsUsed(end+1,1) = d.nObservationsUsed; %#ok<AGROW>
        nModelPixels(end+1,1) = d.nModelPixels; %#ok<AGROW>
        nParameters(end+1,1) = d.nParameters; %#ok<AGROW>
        rank(end+1,1) = d.rank; %#ok<AGROW>
        isFullRank(end+1,1) = d.isFullRank; %#ok<AGROW>
        conditionNumber(end+1,1) = d.conditionNumber; %#ok<AGROW>
        minSingularValue(end+1,1) = d.minSingularValue; %#ok<AGROW>
    end
    T = table(label, nObservationsUsed, nModelPixels, nParameters, rank, ...
        isFullRank, conditionNumber, minSingularValue);
end

function writeSummaryFigure(result, outFile)
    fit = result.fit;
    dataN = fit.dataN;
    fitStack = fit.fitStack;
    residual = dataN - fitStack;

    measured = squeeze(sum(dataN, 3));
    fitted = squeeze(sum(fitStack, 3));
    resid = squeeze(sum(abs(residual), 3));

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1200 760]);
    tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    labels = {'focus', 'diversity'};

    for p = 1:2
        plotPlane(nexttile(tl), measured(:,:,p), fit.sim.x, fit.sim.y, ...
            sprintf('Measured %s', labels{p}));
        plotPlane(nexttile(tl), fitted(:,:,p), fit.sim.x, fit.sim.y, ...
            sprintf('Fitted %s', labels{p}));
        plotPlane(nexttile(tl), resid(:,:,p), fit.sim.x, fit.sim.y, ...
            sprintf('|Residual| %s', labels{p}));
    end

    ax = nexttile(tl, [1 3]);
    coeff = fit.estCoeffVector(:);
    bar(ax, coeff);
    set(ax, 'XTick', 1:numel(coeff), 'XTickLabel', fit.sim.modeOrder, 'XTickLabelRotation', 35);
    ylabel(ax, 'waves RMS');
    title(ax, 'Estimated low-order wavefront');
    grid(ax, 'on');

    title(tl, 'Two-plane ISM wavefront fit');
    saveFigurePng(fig, outFile, 180);
    close(fig);
end

function plotPlane(ax, img, x, y, titleText)
    imagesc(ax, x, y, img);
    axis(ax, 'image');
    set(ax, 'YDir', 'normal');
    colormap(ax, hot);
    colorbar(ax);
    xlabel(ax, 'x (um)');
    ylabel(ax, 'y (um)');
    title(ax, titleText);
end

function saveFigurePng(fig, outFile, dpi)
    if nargin < 3 || isempty(dpi)
        dpi = 180;
    end
    try
        exportgraphics(fig, outFile, 'Resolution', dpi);
    catch
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, outFile, '-dpng', sprintf('-r%d', dpi));
    end
end

function printTwoPlaneSummary(result)
    fit = result.fit;
    fprintf('[estimateTwoPlaneISMWavefront] fit complete.\n');
    fprintf('  residual norm: %.4e\n', fit.residualNorm);
    fprintf('  planeZ: [%g %g] um\n', result.planeZ(1), result.planeZ(2));
    xyPixelSizeUm = fit.sim.dx;
    if isfield(fit.sim, 'xyPixelSizeUm') && isfinitePositiveScalar(fit.sim.xyPixelSizeUm)
        xyPixelSizeUm = fit.sim.xyPixelSizeUm;
    end
    fprintf('  bead offset: x=%+.1f nm, y=%+.1f nm', ...
        1000*fit.estXY(1), 1000*fit.estXY(2));
    if isfinitePositiveScalar(xyPixelSizeUm)
        fprintf('  (%+.2f px, %+.2f px)', fit.estXY(1)/xyPixelSizeUm, fit.estXY(2)/xyPixelSizeUm);
    end
    fprintf('\n');
    if isfield(fit.sim, 'objectiveMagnification')
        fprintf('  objective: %.3g x, NA %.3g\n', fit.sim.objectiveMagnification, fit.sim.NA);
    else
        fprintf('  objective NA: %.3g\n', fit.sim.NA);
    end
    if ~isempty(result.detectorLayoutDiagnostics)
        dlay = result.detectorLayoutDiagnostics;
        pitchNm = median(nearestDistanceLocal(dlay.detXYUm)) * 1000;
        fprintf(['  detector shifts estimated once from focus plane: scan pixel %.2f nm, ' ...
            'median pitch %.3f AU, center detector %d\n'], ...
            dlay.pixelSizeNm, pitchNm/(1000*dlay.airyUnitUm), ...
            dlay.centerDetectorIndex);
    end
    for k = 1:numel(fit.estCoeffVector)
        fprintf('  %-10s %+8.4f waves RMS (%+7.2f nm RMS)\n', ...
            fit.sim.modeOrder{k}, fit.estCoeffVector(k), ...
            fit.estCoeffVector(k) * fit.sim.lamRef * 1000);
    end

    if isfield(result.sufficiency, 'twoPlane')
        d = result.sufficiency.twoPlane;
        fprintf('  two-plane rank: %d/%d, cond %.3g, locally sufficient: %d\n', ...
            d.rank, d.nParameters, d.conditionNumber, d.isFullRank);
    end
    if ~isempty(result.outputDir)
        fprintf('  wrote results to: %s\n', result.outputDir);
    end
end

function d = nearestDistanceLocal(xy)
    n = size(xy, 1);
    d = nan(n, 1);
    for i = 1:n
        dist = [];
        for j = i+1:n
            dist(end+1, 1) = hypot(xy(i,1)-xy(j,1), xy(i,2)-xy(j,2)); %#ok<AGROW>
        end
        for j = 1:i-1
            dist(end+1, 1) = hypot(xy(i,1)-xy(j,1), xy(i,2)-xy(j,2)); %#ok<AGROW>
        end
        dist = dist(isfinite(dist) & dist > 0);
        if ~isempty(dist)
            d(i) = min(dist);
        end
    end
    d = d(isfinite(d) & d > 0);
end

function v = numericField(s, name, defaultValue)
    v = defaultValue;
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name))
        v = double(s.(name)(1));
    end
end

function m = robustMad(vals)
    vals = vals(isfinite(vals));
    if isempty(vals)
        m = NaN;
        return;
    end
    med = median(vals);
    m = 1.4826 * median(abs(vals - med));
    if ~isfinite(m) || m == 0
        m = eps;
    end
end
