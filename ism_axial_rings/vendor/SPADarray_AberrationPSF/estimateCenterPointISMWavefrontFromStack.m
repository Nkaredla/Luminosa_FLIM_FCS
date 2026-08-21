function result = estimateCenterPointISMWavefrontFromStack(stackInput, varargin)
%ESTIMATECENTERPOINTISMWAVEFRONTFROMSTACK Select two z planes, then fit center ISM signal.
%
%   result = estimateCenterPointISMWavefrontFromStack()
%
%   Uses the current Luminosa dataset by default:
%       D:\Luminosa\Aberration_ISM\PDA23_centered_bead_point_10s_20260702-124209
%
%   The wrapper selects the focal plane from the full z-stack signal trace,
%   selects the nearest +500 nm diversity plane, and then calls
%   estimateCenterPointISMWavefront on those two detector-resolved planes.
%   The detector microimage is sampled at the bead centre fitted in the
%   detector-summed focal image. The default fit is focus-anchored: the focal
%   detector microimage sets coefficient magnitudes, while the diversity plane
%   resolves the signs of fitted even modes before a focus-dominant refinement.
%
%   For raw Luminosa scan folders, it first looks for existing batch
%   alignment outputs. If none are found and autoPreprocessRawStack is true,
%   it runs batchAnalyzeLuminosaPSFs on that scan folder and then uses the
%   generated alignment CSV.
%
%   Example:
%       res = estimateCenterPointISMWavefrontFromStack();
%
%       res = estimateCenterPointISMWavefrontFromStack([], ...
%           'scanName', '0.5uW_0.15collar_80mmlens_20260515-164337', ...
%           'defocusOffsetUm', 0.5);

    if nargin < 1 || isempty(stackInput)
        stackInput = 'D:\Luminosa\Aberration_ISM\PDA23_centered_bead_point_10s_20260702-124209';
    end

    opts = parseStackOptions(varargin{:});
    addRequiredPaths(opts);

    selection = selectCenterPointPlanes(stackInput, opts);
    coreArgs = buildCoreArgs(opts, selection);

    if ~isempty(selection.raw4)
        result = estimateCenterPointISMWavefront(selection.raw4, [], coreArgs{:});
    else
        result = estimateCenterPointISMWavefront(selection.focusFile, selection.defocusFile, ...
            coreArgs{:});
    end

    result.stackSelection = selection;

    if opts.writeOutputs
        writeStackSelectionOutputs(result, opts);
    end

    if opts.verbose
        fprintf('[estimateCenterPointISMWavefrontFromStack] focus plane %d at z=%.4f um.\n', ...
            selection.focusIndex, selection.focusZUm);
        fprintf('[estimateCenterPointISMWavefrontFromStack] defocus plane %d at z=%.4f um (relative %.4f um).\n', ...
            selection.defocusIndex, selection.defocusZUm, selection.relativePlaneZ(2));
    end
end

function opts = parseStackOptions(varargin)
    p = inputParser;
    p.FunctionName = 'estimateCenterPointISMWavefrontFromStack';
    p.KeepUnmatched = true;

    addParameter(p, 'dataRoot', 'D:\Luminosa\Aberration_ISM');
    addParameter(p, 'scanName', 'PDA23_centered_bead_point_10s_20260702-124209');
    addParameter(p, 'scanPattern', '');
    addParameter(p, 'scanIndex', []);
    addParameter(p, 'alignmentCsv', '');
    addParameter(p, 'volumeMat', '');
    addParameter(p, 'autoPreprocessRawStack', true);
    addParameter(p, 'batchResultsRoot', '');
    addParameter(p, 'batchInputSource', 'auto');
    addParameter(p, 'batchOptions', struct());
    addParameter(p, 'planeZ', []);
    addParameter(p, 'zStepUm', 0.05);
    addParameter(p, 'xyPixelSizeUm', []);
    addParameter(p, 'defocusOffsetUm', 0.5);
    addParameter(p, 'defocusSide', 'positive');
    addParameter(p, 'focusIndex', []);
    addParameter(p, 'defocusIndex', []);
    addParameter(p, 'focusMetric', 'signal');
    addParameter(p, 'fitStrategy', 'focusAnchored');
    addParameter(p, 'focusWeight', 1);
    addParameter(p, 'diversityWeight', 0.20);
    addParameter(p, 'focusSeedAmplitude', 0.05);
    addParameter(p, 'signSelectionModes', ...
        {'astig_x','astig_y','spherical'});
    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'ptuReaderFolder', defaultPtuReaderFolder());
    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;
    opts.unmatched = p.Unmatched;
    opts.defocusSide = lower(char(opts.defocusSide));
    opts.focusMetric = lower(char(opts.focusMetric));
    opts.batchInputSource = lower(char(opts.batchInputSource));
    if ~isstruct(opts.batchOptions) || ~isscalar(opts.batchOptions)
        error('estimateCenterPointISMWavefrontFromStack:BadBatchOptions', ...
            'batchOptions must be a scalar struct.');
    end
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

function selection = selectCenterPointPlanes(stackInput, opts)
    if isnumeric(stackInput) || iscell(stackInput)
        detStack = standardizeDetectorZStack(stackInput, opts.channelIDs);
        signal = zSignalTrace(detStack, opts.focusMetric);
        planeZ = resolvePlaneZ(numel(signal), [], opts);
        [focusIdx, defocusIdx] = choosePlanePair(signal, planeZ, opts);
        selection = baseSelection('numeric detector stack');
        selection.raw4 = cat(4, detStack(:,:,:,focusIdx), detStack(:,:,:,defocusIdx));
        selection.signalTrace = signal(:).';
        selection.planeZ = planeZ(:).';
        selection = fillSelectionPair(selection, focusIdx, defocusIdx, planeZ);
        return;
    end

    inputPath = char(stackInput);
    if exist(inputPath, 'dir') == 7
        scanFolder = resolveScanFolder(inputPath, opts);
        selection = selectFromScanFolder(scanFolder, opts);
        return;
    end

    if exist(inputPath, 'file') ~= 2
        error('estimateCenterPointISMWavefrontFromStack:MissingInput', ...
            'Stack input was not found: %s', inputPath);
    end

    [~,~,ext] = fileparts(inputPath);
    switch lower(ext)
        case '.mat'
            selection = selectFromMatFile(inputPath, opts);
        case '.csv'
            selection = selectFromAlignmentCsv(inputPath, opts);
        otherwise
            error('estimateCenterPointISMWavefrontFromStack:BadInputFile', ...
                'Unsupported stack input "%s". Use a folder, MAT, CSV, or numeric stack.', inputPath);
    end
end

function selection = baseSelection(label)
    selection = struct();
    selection.label = label;
    selection.scanFolder = '';
    selection.focusFile = '';
    selection.defocusFile = '';
    selection.focusFrame = [];
    selection.defocusFrame = [];
    selection.focusIndex = NaN;
    selection.defocusIndex = NaN;
    selection.focusZUm = NaN;
    selection.defocusZUm = NaN;
    selection.relativePlaneZ = [0 0.5];
    selection.xyPixelSizeUm = NaN;
    selection.signalTrace = [];
    selection.planeZ = [];
    selection.raw4 = [];
end

function selection = fillSelectionPair(selection, focusIdx, defocusIdx, planeZ)
    selection.focusIndex = focusIdx;
    selection.defocusIndex = defocusIdx;
    selection.focusZUm = planeZ(focusIdx);
    selection.defocusZUm = planeZ(defocusIdx);
    selection.relativePlaneZ = [0, planeZ(defocusIdx) - planeZ(focusIdx)];
end

function scanFolder = resolveScanFolder(inputFolder, opts)
    files = sortedSeriesFiles(inputFolder);
    if ~isempty(files)
        scanFolder = inputFolder;
        return;
    end

    dataRoot = inputFolder;
    if ~isempty(opts.scanName)
        candidate = fullfile(dataRoot, char(opts.scanName));
        if exist(candidate, 'dir') == 7
            scanFolder = candidate;
            return;
        end
    end

    folders = dir(dataRoot);
    folders = folders([folders.isdir]);
    names = {folders.name};
    folders = folders(~ismember(names, {'.','..'}));
    if ~isempty(opts.scanPattern)
        folders = folders(contains({folders.name}, char(opts.scanPattern)));
    end
    if isempty(folders)
        error('estimateCenterPointISMWavefrontFromStack:NoScanFolders', ...
            'No scan folders were found below %s.', dataRoot);
    end

    if ~isempty(opts.scanIndex)
        idx = round(opts.scanIndex);
        if idx < 1 || idx > numel(folders)
            error('estimateCenterPointISMWavefrontFromStack:BadScanIndex', ...
                'scanIndex must be between 1 and %d.', numel(folders));
        end
    else
        [~, idx] = max([folders.datenum]);
    end
    scanFolder = fullfile(folders(idx).folder, folders(idx).name);
end

function selection = selectFromScanFolder(scanFolder, opts)
    csvFile = char(opts.alignmentCsv);
    if isempty(csvFile)
        csvFile = findBatchAlignmentCsvForScanFolder(scanFolder, opts);
    end
    if ~isempty(csvFile) && exist(csvFile, 'file') == 2
        selection = selectFromAlignmentCsv(csvFile, opts);
        selection.scanFolder = scanFolder;
        return;
    end

    matFile = char(opts.volumeMat);
    if isempty(matFile)
        matFile = findBatchVolumeMatForScanFolder(scanFolder, opts);
    end
    if ~isempty(matFile) && exist(matFile, 'file') == 2
        selection = selectFromMatFile(matFile, opts);
        selection.scanFolder = scanFolder;
        return;
    end

    if opts.autoPreprocessRawStack && hasRawStackFiles(scanFolder)
        preprocessRawScanFolder(scanFolder, opts);

        csvFile = findBatchAlignmentCsvForScanFolder(scanFolder, opts);
        if ~isempty(csvFile) && exist(csvFile, 'file') == 2
            selection = selectFromAlignmentCsv(csvFile, opts);
            selection.scanFolder = scanFolder;
            return;
        end

        matFile = findBatchVolumeMatForScanFolder(scanFolder, opts);
        if ~isempty(matFile) && exist(matFile, 'file') == 2
            selection = selectFromMatFile(matFile, opts);
            selection.scanFolder = scanFolder;
            return;
        end
    end

    error('estimateCenterPointISMWavefrontFromStack:NeedStackMetadata', ...
        ['No batch alignment CSV or volume MAT was found for %s. Run ' ...
         'batchAnalyzeLuminosaPSFs first, or pass focus/defocus planes directly ' ...
         'to estimateCenterPointISMWavefront.'], scanFolder);
end

function tf = hasRawStackFiles(scanFolder)
    tf = ~isempty(sortedSeriesFiles(scanFolder)) || ...
        exist(fullfile(scanFolder, 'IntensityImage.pqdat'), 'file') == 2;
end

function preprocessRawScanFolder(scanFolder, opts)
    batchOpts = opts.batchOptions;
    batchOpts.inputSource = resolveBatchInputSource(scanFolder, opts);
    batchOpts.zStepUm = opts.zStepUm;
    batchOpts.verbose = opts.verbose;
    if ~isempty(opts.ptuReaderFolder)
        batchOpts.ptuReaderFolder = opts.ptuReaderFolder;
    end
    if isnumeric(opts.xyPixelSizeUm) && isscalar(opts.xyPixelSizeUm) && ...
            isfinite(opts.xyPixelSizeUm) && opts.xyPixelSizeUm > 0
        batchOpts.xyPixelSizeUm = opts.xyPixelSizeUm;
    end

    resultsRoot = char(opts.batchResultsRoot);
    if isempty(resultsRoot)
        batchAnalyzeLuminosaPSFs(scanFolder, [], batchOpts);
    else
        batchAnalyzeLuminosaPSFs(scanFolder, resultsRoot, batchOpts);
    end
end

function source = resolveBatchInputSource(scanFolder, opts)
    source = lower(char(opts.batchInputSource));
    if strcmp(source, 'auto')
        if ~isempty(sortedSeriesFiles(scanFolder))
            source = 'ptu';
        elseif exist(fullfile(scanFolder, 'IntensityImage.pqdat'), 'file') == 2
            source = 'pqdat';
        else
            source = 'ptu';
        end
    end
    if ~ismember(source, {'ptu','pqdat'})
        error('estimateCenterPointISMWavefrontFromStack:BadBatchInputSource', ...
            'batchInputSource must be ''auto'', ''ptu'', or ''pqdat''.');
    end
end

function selection = selectFromAlignmentCsv(csvFile, opts)
    T = readtable(csvFile);
    if ~all(ismember({'total_signal','source_file'}, T.Properties.VariableNames))
        error('estimateCenterPointISMWavefrontFromStack:BadAlignmentCsv', ...
            'Alignment CSV must contain total_signal and source_file columns.');
    end

    signal = double(T.total_signal(:)).';
    if ismember('z_um', T.Properties.VariableNames)
        planeZ = double(T.z_um(:)).';
    else
        planeZ = resolvePlaneZ(numel(signal), [], opts);
    end

    [focusIdx, defocusIdx] = choosePlanePair(signal, planeZ, opts);
    selection = baseSelection(localFileStem(csvFile));
    selection.signalTrace = signal;
    selection.planeZ = planeZ;
    selection.focusFile = tableText(T.source_file, focusIdx);
    selection.defocusFile = tableText(T.source_file, defocusIdx);
    selection.focusFrame = tableFrame(T, focusIdx);
    selection.defocusFrame = tableFrame(T, defocusIdx);
    selection = fillSelectionPair(selection, focusIdx, defocusIdx, planeZ);
end

function selection = selectFromMatFile(matFile, opts)
    S = load(matFile);
    [value, isDetectorResolved] = chooseStackVariable(S);
    if isDetectorResolved
        detStack = standardizeDetectorZStack(value, opts.channelIDs);
        signal = zSignalTrace(detStack, opts.focusMetric);
        nPlane = size(detStack, 4);
    else
        vol = double(value);
        signal = zSignalTrace(vol, opts.focusMetric);
        nPlane = size(vol, 3);
        detStack = [];
    end

    planeZ = planeZFromMatStruct(S);
    planeZ = resolvePlaneZ(nPlane, planeZ, opts);
    [focusIdx, defocusIdx] = choosePlanePair(signal, planeZ, opts);

    selection = baseSelection(localFileStem(matFile));
    selection.xyPixelSizeUm = xyPixelSizeFromMatStruct(S);
    selection.signalTrace = signal(:).';
    selection.planeZ = planeZ(:).';
    selection = fillSelectionPair(selection, focusIdx, defocusIdx, planeZ);

    if ~isempty(detStack)
        selection.raw4 = cat(4, detStack(:,:,:,focusIdx), detStack(:,:,:,defocusIdx));
        return;
    end

    [selection.focusFile, selection.focusFrame] = sourceFromFrameMeta(S, focusIdx);
    [selection.defocusFile, selection.defocusFrame] = sourceFromFrameMeta(S, defocusIdx);
    if isempty(selection.focusFile) || isempty(selection.defocusFile)
        error('estimateCenterPointISMWavefrontFromStack:NeedDetectorSource', ...
            'Scalar MAT stack has no frameMeta.sourceFile mapping to detector PTUs.');
    end
end

function [focusIdx, defocusIdx] = choosePlanePair(signal, planeZ, opts)
    nPlane = numel(signal);
    if ~isempty(opts.focusIndex)
        focusIdx = round(opts.focusIndex);
        if focusIdx < 1 || focusIdx > nPlane
            error('estimateCenterPointISMWavefrontFromStack:BadFocusIndex', ...
                'focusIndex must be between 1 and %d.', nPlane);
        end
    else
        [~, focusIdx] = max(signal);
    end

    if ~isempty(opts.defocusIndex)
        defocusIdx = round(opts.defocusIndex);
        if defocusIdx < 1 || defocusIdx > nPlane || defocusIdx == focusIdx
            error('estimateCenterPointISMWavefrontFromStack:BadDefocusIndex', ...
                'defocusIndex must be a valid plane different from focusIndex.');
        end
    else
        defocusIdx = chooseDiversityPlane(planeZ, focusIdx, opts);
    end
end

function defocusIdx = chooseDiversityPlane(planeZ, focusIdx, opts)
    offset = abs(opts.defocusOffsetUm);
    if offset <= 0
        error('estimateCenterPointISMWavefrontFromStack:BadOffset', ...
            'defocusOffsetUm must be positive.');
    end

    switch opts.defocusSide
        case {'positive','plus','+'}
            signs = [1 -1];
        case {'negative','minus','-'}
            signs = [-1 1];
        case 'auto'
            signs = [1 -1];
        otherwise
            error('estimateCenterPointISMWavefrontFromStack:BadDefocusSide', ...
                'defocusSide must be positive, negative, or auto.');
    end

    bestIdx = [];
    bestErr = inf;
    for s = signs
        target = planeZ(focusIdx) + s * offset;
        [err, idx] = min(abs(planeZ - target));
        if idx ~= focusIdx && err < bestErr
            bestIdx = idx;
            bestErr = err;
        end
        if ~strcmp(opts.defocusSide, 'auto') && idx ~= focusIdx
            break;
        end
    end

    if isempty(bestIdx)
        candidates = setdiff(1:numel(planeZ), focusIdx);
        [~, ii] = min(abs(abs(planeZ(candidates) - planeZ(focusIdx)) - offset));
        bestIdx = candidates(ii);
    end
    defocusIdx = bestIdx;
end

function args = buildCoreArgs(opts, selection)
    args = unmatchedToArgs(opts.unmatched);
    args = setNameValue(args, 'planeZ', selection.relativePlaneZ);
    args = setNameValue(args, 'channelIDs', opts.channelIDs);
    args = setNameValue(args, 'ptuReaderFolder', opts.ptuReaderFolder);
    args = setNameValue(args, 'writeOutputs', opts.writeOutputs);
    args = setNameValue(args, 'verbose', opts.verbose);
    args = setNameValue(args, 'fitStrategy', opts.fitStrategy);
    args = setNameValue(args, 'focusWeight', opts.focusWeight);
    args = setNameValue(args, 'diversityWeight', opts.diversityWeight);
    args = setNameValue(args, 'focusSeedAmplitude', opts.focusSeedAmplitude);
    args = setNameValue(args, 'signSelectionModes', opts.signSelectionModes);
    args = setNameValue(args, 'focusFrameIndices', selection.focusFrame);
    args = setNameValue(args, 'defocusFrameIndices', selection.defocusFrame);
    if isnumeric(opts.xyPixelSizeUm) && isscalar(opts.xyPixelSizeUm) && ...
            isfinite(opts.xyPixelSizeUm) && opts.xyPixelSizeUm > 0 && ...
            ~hasNameValue(args, 'xyPixelSizeUm')
        args = setNameValue(args, 'xyPixelSizeUm', opts.xyPixelSizeUm);
    elseif isfinite(selection.xyPixelSizeUm) && selection.xyPixelSizeUm > 0 && ...
            ~hasNameValue(args, 'xyPixelSizeUm')
        args = setNameValue(args, 'xyPixelSizeUm', selection.xyPixelSizeUm);
    end

    if ~isempty(opts.outputDir)
        args = setNameValue(args, 'outputDir', opts.outputDir);
    else
        args = setNameValue(args, 'outputDir', defaultOutputDir(selection, opts));
    end
end

function args = unmatchedToArgs(unmatched)
    names = fieldnames(unmatched);
    args = cell(1, 2*numel(names));
    for k = 1:numel(names)
        args{2*k-1} = names{k};
        args{2*k} = unmatched.(names{k});
    end
end

function args = setNameValue(args, name, value)
    for k = 1:2:numel(args)
        if strcmpi(args{k}, name)
            args{k+1} = value;
            return;
        end
    end
    args = [args {name, value}];
end

function tf = hasNameValue(args, name)
    tf = false;
    for k = 1:2:numel(args)
        if strcmpi(args{k}, name)
            tf = true;
            return;
        end
    end
end

function xyPixelSizeUm = xyPixelSizeFromMatStruct(S)
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

function detStack = standardizeDetectorZStack(value, channelIDs)
    if iscell(value)
        nCh = numel(value);
        first = double(value{1});
        if ndims(first) == 2
            first = reshape(first, size(first,1), size(first,2), 1);
        end
        [ny, nx, nz] = size(first);
        detStack = zeros(ny, nx, nCh, nz);
        detStack(:,:,1,:) = reshape(first, ny, nx, 1, nz);
        for k = 2:nCh
            vol = double(value{k});
            if ndims(vol) == 2
                vol = reshape(vol, size(vol,1), size(vol,2), 1);
            end
            if ~isequal(size(vol), [ny nx nz])
                error('estimateCenterPointISMWavefrontFromStack:BadCellStack', ...
                    'All detector cell volumes must have the same [y x z] size.');
            end
            detStack(:,:,k,:) = reshape(vol, ny, nx, 1, nz);
        end
        return;
    end

    value = double(value);
    nCh = numel(channelIDs);
    if ndims(value) ~= 4
        error('estimateCenterPointISMWavefrontFromStack:BadDetectorStack', ...
            'Detector-resolved full stack must be [y x 23 x z] or [y x z x 23].');
    end
    if size(value, 3) == nCh
        detStack = value;
    elseif size(value, 4) == nCh
        detStack = permute(value, [1 2 4 3]);
    elseif size(value, 3) == 23
        detStack = value;
    elseif size(value, 4) == 23
        detStack = permute(value, [1 2 4 3]);
    else
        error('estimateCenterPointISMWavefrontFromStack:BadDetectorStack', ...
            'Could not identify the detector channel dimension.');
    end
end

function [value, isDetectorResolved] = chooseStackVariable(S)
    isDetectorResolved = false;
    detectorNames = {'raw4','rawData','detectorVolume','detectorStack', ...
        'channelVolume','channelVolumes','spadVolume','spadStack'};
    for k = 1:numel(detectorNames)
        name = detectorNames{k};
        if isfield(S, name) && isDetectorStack(S.(name))
            value = S.(name);
            isDetectorResolved = true;
            return;
        end
    end

    names = fieldnames(S);
    for k = 1:numel(names)
        if isDetectorStack(S.(names{k}))
            value = S.(names{k});
            isDetectorResolved = true;
            return;
        end
    end

    scalarNames = {'volume','alignedVolume','rawVolume','rawVol','alignedVol'};
    for k = 1:numel(scalarNames)
        name = scalarNames{k};
        if isfield(S, name) && isnumeric(S.(name)) && ndims(S.(name)) == 3
            value = S.(name);
            return;
        end
    end

    for k = 1:numel(names)
        v = S.(names{k});
        if isnumeric(v) && ndims(v) == 3 && size(v, 3) ~= 23
            value = v;
            return;
        end
    end

    error('estimateCenterPointISMWavefrontFromStack:NoStackVariable', ...
        'No detector-resolved 4-D stack or scalar 3-D volume was found.');
end

function tf = isDetectorStack(value)
    tf = isnumeric(value) && ndims(value) == 4 && ...
        (size(value, 3) == 23 || size(value, 4) == 23);
    if iscell(value) && numel(value) == 23
        tf = true;
    end
end

function signal = zSignalTrace(stack, metric)
    switch metric
        case {'signal','sum','total'}
            if ndims(stack) == 4
                nz = size(stack, 4);
                signal = zeros(1, nz);
                for z = 1:nz
                    signal(z) = sum(max(reshape(stack(:,:,:,z), [], 1), 0));
                end
            else
                nz = size(stack, 3);
                signal = zeros(1, nz);
                for z = 1:nz
                    signal(z) = sum(max(reshape(stack(:,:,z), [], 1), 0));
                end
            end
        case {'peak','max'}
            if ndims(stack) == 4
                nz = size(stack, 4);
                signal = zeros(1, nz);
                for z = 1:nz
                    signal(z) = max(reshape(stack(:,:,:,z), [], 1));
                end
            else
                nz = size(stack, 3);
                signal = zeros(1, nz);
                for z = 1:nz
                    signal(z) = max(reshape(stack(:,:,z), [], 1));
                end
            end
        otherwise
            error('estimateCenterPointISMWavefrontFromStack:BadFocusMetric', ...
                'focusMetric must be ''signal'' or ''peak''.');
    end
end

function planeZ = resolvePlaneZ(nPlane, metaPlaneZ, opts)
    if ~isempty(opts.planeZ)
        planeZ = double(opts.planeZ(:)).';
    elseif ~isempty(metaPlaneZ)
        planeZ = double(metaPlaneZ(:)).';
    else
        planeZ = ((1:nPlane) - (nPlane+1)/2) * opts.zStepUm;
    end
    if numel(planeZ) ~= nPlane
        error('estimateCenterPointISMWavefrontFromStack:PlaneZMismatch', ...
            'planeZ has %d entries but the stack has %d z planes.', numel(planeZ), nPlane);
    end
end

function planeZ = planeZFromMatStruct(S)
    planeZ = [];
    if isfield(S, 'frameMeta') && isstruct(S.frameMeta) && isfield(S.frameMeta, 'zUm')
        planeZ = [S.frameMeta.zUm];
    end
end

function [fileName, frame] = sourceFromFrameMeta(S, idx)
    fileName = '';
    frame = [];
    if ~isfield(S, 'frameMeta') || ~isstruct(S.frameMeta) || numel(S.frameMeta) < idx
        return;
    end
    fm = S.frameMeta(idx);
    if isfield(fm, 'sourceFile') && ~isempty(fm.sourceFile)
        fileName = char(fm.sourceFile);
    elseif isfield(fm, 'frameFiles') && ~isempty(fm.frameFiles)
        fileName = char(fm.frameFiles);
    end
    if isfield(fm, 'sourceFrame') && ~isempty(fm.sourceFrame)
        frameValue = double(fm.sourceFrame);
        frameValue = frameValue(1);
        if isfinite(frameValue) && frameValue >= 1
            frame = frameValue;
        end
    elseif isfield(fm, 'frameFileFrames') && ~isempty(fm.frameFileFrames)
        frameValue = double(fm.frameFileFrames);
        frameValue = frameValue(1);
        if isfinite(frameValue) && frameValue >= 1
            frame = frameValue;
        end
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
    csvFile = '';
    if nargin < 2
        opts = struct();
    end
    candidates = batchOutputCandidates(scanFolder, opts, 'alignment');
    csvFile = firstExistingFile(candidates);
end

function matFile = findBatchVolumeMatForScanFolder(scanFolder, opts)
    matFile = '';
    if nargin < 2
        opts = struct();
    end
    candidates = batchOutputCandidates(scanFolder, opts, 'volume');
    matFile = newestExistingFile(candidates);
end

function candidates = batchOutputCandidates(scanFolder, opts, kind)
    scanFolder = stripTrailingFilesep(scanFolder);
    [dataRoot, scanName] = fileparts(scanFolder);
    [dataParent, datasetName] = fileparts(dataRoot);
    stem = scanOutputStem(scanName);

    parentOutputRoot = fullfile(dataParent, 'PSF_batch_outputs', ...
        sanitizeFileName(datasetName));
    directOutputRoot = fullfile(dataRoot, 'PSF_batch_outputs', ...
        sanitizeFileName(scanName));
    pqdatOutputRoot = fullfile(scanFolder, 'results_psf3d');

    if isfield(opts, 'batchResultsRoot') && ~isempty(opts.batchResultsRoot)
        explicitOutputRoot = char(opts.batchResultsRoot);
    else
        explicitOutputRoot = '';
    end

    roots = {parentOutputRoot, directOutputRoot, pqdatOutputRoot};
    if ~isempty(explicitOutputRoot)
        roots = [{explicitOutputRoot}, roots];
    end

    candidates = {};
    for r = 1:numel(roots)
        root = roots{r};
        switch kind
            case 'alignment'
                candidates{end+1} = fullfile(root, 'xz_yz_plots', ...
                    sprintf('%s_frame_alignment.csv', stem)); %#ok<AGROW>
                candidates{end+1} = fullfile(root, sanitizeFileName(scanName), ...
                    'frame_alignment.csv'); %#ok<AGROW>
                candidates = [candidates groupedBatchFiles(root, ...
                    'frame_alignment.csv', scanName)]; %#ok<AGROW>
            case 'volume'
                candidates{end+1} = fullfile(root, 'volumes_mat', ...
                    sprintf('%s_volume_raw.mat', stem)); %#ok<AGROW>
                candidates{end+1} = fullfile(root, sanitizeFileName(scanName), ...
                    'psf_volume_aligned.mat'); %#ok<AGROW>
                candidates = [candidates groupedBatchFiles(root, ...
                    'psf_volume_aligned.mat', scanName)]; %#ok<AGROW>
        end
    end
end

function files = groupedBatchFiles(root, fileName, scanName)
    files = {};
    hits = dir(fullfile(root, '*', fileName));
    if numel(hits) == 1 && ~hits(1).isdir
        files{1} = fullfile(hits(1).folder, hits(1).name);
        return;
    end
    scanKey = lower(sanitizeFileName(scanName));
    stemKey = lower(scanOutputStem(scanName));
    for k = 1:numel(hits)
        if ~hits(k).isdir
            [~, folderName] = fileparts(hits(k).folder);
            folderKey = lower(folderName);
            if contains(folderKey, scanKey) || contains(folderKey, stemKey)
                files{end+1} = fullfile(hits(k).folder, hits(k).name); %#ok<AGROW>
            end
        end
    end
end

function fileName = firstExistingFile(candidates)
    fileName = '';
    for k = 1:numel(candidates)
        if exist(candidates{k}, 'file') == 2
            fileName = candidates{k};
            return;
        end
    end
end

function fileName = newestExistingFile(candidates)
    fileName = '';
    newestTime = -inf;
    for k = 1:numel(candidates)
        if exist(candidates{k}, 'file') ~= 2
            continue;
        end
        info = dir(candidates{k});
        if ~isempty(info) && info.datenum > newestTime
            newestTime = info.datenum;
            fileName = candidates{k};
        end
    end
end

function text = tableText(value, idx)
    if iscell(value)
        text = char(value{idx});
    elseif isstring(value)
        text = char(value(idx));
    else
        text = char(value(idx));
    end
end

function frame = tableFrame(T, idx)
    frame = [];
    if ismember('source_frame', T.Properties.VariableNames)
        frameValue = double(T.source_frame(idx));
        if isfinite(frameValue) && frameValue >= 1
            frame = frameValue;
        end
    end
end

function outDir = defaultOutputDir(selection, opts)
    thisDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(thisDir);
    label = selection.label;
    if isempty(label) && ~isempty(selection.scanFolder)
        label = localFileStem(selection.scanFolder);
    end
    if isempty(label)
        label = 'stack';
    end
    dzNm = round(abs(selection.relativePlaneZ(2)) * 1000);
    label = sprintf('%s_%s_%dnm', label, char(opts.fitStrategy), dzNm);
    outDir = fullfile(rootDir, 'output_matlab', ...
        'center_point_ism_wavefront', sanitizeFileName(label));
end

function writeStackSelectionOutputs(result, opts)
    outDir = result.outputDir;
    if isempty(outDir)
        return;
    end
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    selection = result.stackSelection;
    zUm = selection.planeZ(:);
    signal = selection.signalTrace(:);
    isFocus = false(numel(zUm), 1);
    isDefocus = false(numel(zUm), 1);
    isFocus(selection.focusIndex) = true;
    isDefocus(selection.defocusIndex) = true;
    T = table((1:numel(zUm)).', zUm, signal, isFocus, isDefocus, ...
        'VariableNames', {'planeIndex','zUm','signal','isFocus','isDefocus'});
    writetable(T, fullfile(outDir, 'center_point_auto_selected_planes.csv'));
    if opts.writeOutputs
        save(fullfile(outDir, 'center_point_ism_wavefront_from_stack.mat'), 'result', '-v7.3');
    end
end

function stem = localFileStem(pathName)
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
