function result = estimateTwoPlaneISMWavefrontFromStack(stackInput, varargin)
%ESTIMATETWOPLANEISMWAVEFRONTFROMSTACK Select focus/diversity planes from a z stack.
%
%   result = estimateTwoPlaneISMWavefrontFromStack()
%
%   Uses the current Luminosa aberration dataset by default:
%       D:\Luminosa\Data\ISM_Aberation2_73
%
%   It locates a full z stack, chooses the focal plane automatically from the
%   z signal trace, chooses the plane closest to +300 nm by default, and
%   calls estimateTwoPlaneISMWavefront on the selected two detector-resolved
%   ISM planes.
%
%   Accepted stackInput values:
%       - detector-resolved numeric stack [y x 23 z] or [y x z 23]
%       - MAT file containing a detector-resolved stack
%       - saved batch scalar-volume MAT file, when frameMeta/source PTUs are
%         present so the selected detector planes can be read from PTU
%       - scan folder containing Series_*.ptu
%       - dataset root containing scan folders, e.g. ISM_Aberation2_73
%
%   Examples:
%       result = estimateTwoPlaneISMWavefrontFromStack();
%
%       result = estimateTwoPlaneISMWavefrontFromStack( ...
%           'D:\Luminosa\Data\ISM_Aberation2_73\0.4uW_0.19collar_80mmlens_20260515-155248');
%
%       result = estimateTwoPlaneISMWavefrontFromStack(raw4d, ...
%           'zStepUm', 0.05, 'defocusOffsetUm', 0.3);
%
%   The returned result is the two-plane fit result with additional
%   stackSelection and stackMeta fields.

    if nargin < 1 || isempty(stackInput)
        stackInput = 'D:\Luminosa\Data\ISM_Aberation2_73';
    end

    opts = parseStackOptions(varargin{:});
    addRequiredPaths(opts);

    [stackData, meta] = loadStackForPlaneSelection(stackInput, opts);
    [selection, twoPlaneData] = selectFocusAndDiversityPlanes(stackData, meta, opts);

    twoPlaneArgs = buildTwoPlaneArgs(opts, meta, selection);
    if isempty(twoPlaneData.raw4)
        focusSource = selection.focusSource;
        diversitySource = selection.diversitySource;
        if isempty(focusSource.file) || isempty(diversitySource.file)
            error('estimateTwoPlaneISMWavefrontFromStack:NeedDetectorSource', ...
                ['The full stack is detector-summed only and has no source PTU mapping. ' ...
                 'Provide detector-resolved [y x 23 z] data or a stack MAT with frameMeta.sourceFile.']);
        end

        twoPlaneArgs = setNameValue(twoPlaneArgs, 'focusFrameIndices', focusSource.frame);
        twoPlaneArgs = setNameValue(twoPlaneArgs, 'defocusFrameIndices', diversitySource.frame);
        result = estimateTwoPlaneISMWavefront(focusSource.file, diversitySource.file, twoPlaneArgs{:});
    else
        result = estimateTwoPlaneISMWavefront(twoPlaneData.raw4, [], twoPlaneArgs{:});
    end

    result.stackSelection = selection;
    result.stackMeta = meta;

    if opts.writeOutputs
        writeStackSelectionOutputs(result, opts);
    end

    if opts.verbose
        fprintf('[estimateTwoPlaneISMWavefrontFromStack] selected focus plane %d at z=%.4f um.\n', ...
            selection.focusIndex, selection.focusZUm);
        fprintf('[estimateTwoPlaneISMWavefrontFromStack] selected diversity plane %d at z=%.4f um (relative %.4f um).\n', ...
            selection.diversityIndex, selection.diversityZUm, selection.relativePlaneZ(2));
    end
end

function opts = parseStackOptions(varargin)
    p = inputParser;
    p.FunctionName = 'estimateTwoPlaneISMWavefrontFromStack';
    p.KeepUnmatched = true;

    addParameter(p, 'dataRoot', 'D:\Luminosa\Data\ISM_Aberation2_73');
    addParameter(p, 'scanName', '0.4uW_0.19collar_80mmlens_20260515-155248');
    addParameter(p, 'scanPattern', '');
    addParameter(p, 'scanIndex', []);
    addParameter(p, 'volumeMat', '');
    addParameter(p, 'stackVariable', '');
    addParameter(p, 'planeFiles', {});
    addParameter(p, 'planeZ', []);
    addParameter(p, 'zStepUm', 0.05);
    addParameter(p, 'defocusOffsetUm', 0.3);
    addParameter(p, 'defocusSide', 'positive');
    addParameter(p, 'focusIndex', []);
    addParameter(p, 'defocusIndex', []);
    addParameter(p, 'focusMetric', 'signal');
    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'ptuReaderFolder', defaultPtuReaderFolder());
    addParameter(p, 'ptuPhotonsPerChunk', 1e6);
    addParameter(p, 'darkFile', 'D:\Luminosa\Data\ISMdark_counts.ptu');
    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'writeSelectionFigure', true);
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;
    opts.unmatched = p.Unmatched;
    opts.defocusSide = lower(char(opts.defocusSide));
    opts.focusMetric = lower(char(opts.focusMetric));
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

function [stackData, meta] = loadStackForPlaneSelection(stackInput, opts)
    if isnumeric(stackInput) || iscell(stackInput)
        [detStack, channelIDs] = standardizeDetectorZStack(stackInput, opts.channelIDs);
        stackData = struct('detectorStack', detStack, 'signalVolume', [], ...
            'signalTrace', zSignalTrace(detStack, opts.focusMetric), ...
            'planeSources', emptyPlaneSource(0));
        meta = baseMeta('numeric detector-resolved stack');
        meta.channelIDs = channelIDs(:);
        meta.planeZ = resolvePlaneZ(size(detStack, 4), [], opts);
        return;
    end

    inputPath = char(stackInput);
    if exist(inputPath, 'dir') == 7
        [stackData, meta] = loadFolderStack(inputPath, opts);
        return;
    end

    if exist(inputPath, 'file') ~= 2
        error('estimateTwoPlaneISMWavefrontFromStack:MissingInput', ...
            'Stack input was not found: %s', inputPath);
    end

    [~,~,ext] = fileparts(inputPath);
    switch lower(ext)
        case '.mat'
            [stackData, meta] = loadMatStack(inputPath, opts);
        case '.ptu'
            [detStack, ptuMeta] = readPtuDetectorZStack(inputPath, opts);
            stackData = struct('detectorStack', detStack, 'signalVolume', [], ...
                'signalTrace', zSignalTrace(detStack, opts.focusMetric), ...
                'planeSources', makeSingleFilePlaneSources(inputPath, size(detStack, 4)));
            meta = ptuMeta;
            meta.source = 'PTU detector z stack';
            meta.planeZ = resolvePlaneZ(size(detStack, 4), [], opts);
        otherwise
            error('estimateTwoPlaneISMWavefrontFromStack:BadInputFile', ...
                'Unsupported stack file "%s". Use MAT or PTU.', inputPath);
    end
end

function meta = baseMeta(source)
    meta = struct();
    meta.source = source;
    meta.inputPath = '';
    meta.scanFolder = '';
    meta.volumeMat = '';
    meta.channelIDs = [];
    meta.xyPixelSizeUm = NaN;
    meta.planeZ = [];
end

function [stackData, meta] = loadFolderStack(folderPath, opts)
    seriesFiles = sortedSeriesFiles(folderPath);
    if isempty(seriesFiles)
        scanFolder = resolveScanFolderFromDatasetRoot(folderPath, opts);
    else
        scanFolder = folderPath;
    end

    matFile = char(opts.volumeMat);
    if isempty(matFile)
        matFile = findBatchVolumeMatForScanFolder(scanFolder);
    end

    if ~isempty(matFile) && exist(matFile, 'file') == 2
        [stackData, meta] = loadMatStack(matFile, opts);
        meta.scanFolder = scanFolder;
        if isempty(stackData.planeSources)
            stackData.planeSources = inferPlaneSourcesFromScanFolder(scanFolder, numel(stackData.signalTrace));
        end
        return;
    end

    files = sortedSeriesFiles(scanFolder);
    if isempty(files)
        error('estimateTwoPlaneISMWavefrontFromStack:NoSeriesFiles', ...
            'No Series_*.ptu files or batch volume MAT were found for %s.', scanFolder);
    end

    fileNames = fullfile({files.folder}, {files.name});
    signal = ptuFileSignalTrace(fileNames, opts);
    stackData = struct('detectorStack', [], 'signalVolume', [], ...
        'signalTrace', signal(:).', ...
        'planeSources', makeFileListPlaneSources(fileNames));

    meta = baseMeta('scan folder Series_*.ptu');
    meta.scanFolder = scanFolder;
    meta.inputPath = scanFolder;
    meta.channelIDs = opts.channelIDs(:);
    meta.planeZ = resolvePlaneZ(numel(signal), [], opts);
end

function scanFolder = resolveScanFolderFromDatasetRoot(dataRoot, opts)
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
    keep = ~ismember(names, {'.','..'});
    folders = folders(keep);

    if ~isempty(opts.scanPattern)
        match = contains({folders.name}, char(opts.scanPattern));
        folders = folders(match);
    end

    if isempty(folders)
        error('estimateTwoPlaneISMWavefrontFromStack:NoScanFolders', ...
            'No scan folders were found below %s.', dataRoot);
    end

    if ~isempty(opts.scanIndex)
        idx = round(opts.scanIndex);
        if idx < 1 || idx > numel(folders)
            error('estimateTwoPlaneISMWavefrontFromStack:BadScanIndex', ...
                'scanIndex must be between 1 and %d.', numel(folders));
        end
    else
        [~, idx] = max([folders.datenum]);
    end
    scanFolder = fullfile(folders(idx).folder, folders(idx).name);
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

function matFile = findBatchVolumeMatForScanFolder(scanFolder)
    matFile = '';
    [dataRoot, scanName] = fileparts(stripTrailingFilesep(scanFolder));
    [dataParent, datasetName] = fileparts(dataRoot);
    volumeDir = fullfile(dataParent, 'PSF_batch_outputs', sanitizeFileName(datasetName), 'volumes_mat');
    if exist(volumeDir, 'dir') ~= 7
        return;
    end

    stem = scanOutputStem(scanName);
    fallbackStem = sanitizeFileName(scanName);
    timeToken = regexp(scanName, '\d{8}[-_]\d{6}', 'match', 'once');
    if ~isempty(timeToken)
        timeToken = regexprep(timeToken, '[^0-9]+', '_');
    end
    patterns = { ...
        sprintf('%s_volume_raw.mat', stem), ...
        sprintf('*%s*volume*.mat', stem), ...
        sprintf('*%s*volume*.mat', fallbackStem)};
    if ~isempty(timeToken)
        patterns{end+1} = sprintf('*%s*volume*.mat', timeToken);
    end

    hits = [];
    for k = 1:numel(patterns)
        hits = [hits; dir(fullfile(volumeDir, patterns{k}))]; %#ok<AGROW>
    end
    if isempty(hits)
        return;
    end
    [~, newest] = max([hits.datenum]);
    matFile = fullfile(hits(newest).folder, hits(newest).name);
end

function [stackData, meta] = loadMatStack(matFile, opts)
    S = load(matFile);
    meta = baseMeta('MAT stack');
    meta.inputPath = matFile;
    meta.volumeMat = matFile;

    [value, variableName, isDetectorResolved] = chooseStackVariable(S, opts.stackVariable);
    meta.stackVariable = variableName;

    if isDetectorResolved
        [detStack, channelIDs] = standardizeDetectorZStack(value, opts.channelIDs);
        stackData.detectorStack = detStack;
        stackData.signalVolume = [];
        stackData.signalTrace = zSignalTrace(detStack, opts.focusMetric);
        meta.channelIDs = channelIDs(:);
        nPlane = size(detStack, 4);
    else
        vol = double(value);
        if ndims(vol) ~= 3
            error('estimateTwoPlaneISMWavefrontFromStack:BadScalarVolume', ...
                'Scalar stack variable "%s" must be a 3-D volume.', variableName);
        end
        stackData.detectorStack = [];
        stackData.signalVolume = vol;
        stackData.signalTrace = zSignalTrace(vol, opts.focusMetric);
        nPlane = size(vol, 3);
    end

    meta.xyPixelSizeUm = xyPixelSizeFromMatStruct(S);
    meta.planeZ = resolvePlaneZ(nPlane, planeZFromMatStruct(S), opts);
    stackData.planeSources = planeSourcesFromMatStruct(S, nPlane);
end

function [value, variableName, isDetectorResolved] = chooseStackVariable(S, variableName)
    isDetectorResolved = false;
    variableName = char(variableName);
    if ~isempty(variableName)
        if ~isfield(S, variableName)
            error('estimateTwoPlaneISMWavefrontFromStack:MissingVariable', ...
                'Variable "%s" was not found in the MAT file.', variableName);
        end
        value = S.(variableName);
        isDetectorResolved = isDetectorStack(value);
        return;
    end

    detectorNames = {'raw4','rawData','detectorVolume','detectorStack', ...
        'channelVolume','channelVolumes','spadVolume','spadStack'};
    for k = 1:numel(detectorNames)
        name = detectorNames{k};
        if isfield(S, name) && isDetectorStack(S.(name))
            value = S.(name);
            variableName = name;
            isDetectorResolved = true;
            return;
        end
    end

    names = fieldnames(S);
    for k = 1:numel(names)
        if isDetectorStack(S.(names{k}))
            value = S.(names{k});
            variableName = names{k};
            isDetectorResolved = true;
            return;
        end
    end

    scalarNames = {'volume','alignedVolume','rawVolume','rawVol','alignedVol'};
    for k = 1:numel(scalarNames)
        name = scalarNames{k};
        if isfield(S, name) && isnumeric(S.(name)) && ndims(S.(name)) == 3
            value = S.(name);
            variableName = name;
            return;
        end
    end

    for k = 1:numel(names)
        v = S.(names{k});
        if isnumeric(v) && ndims(v) == 3 && size(v, 3) ~= 23
            value = v;
            variableName = names{k};
            return;
        end
    end

    error('estimateTwoPlaneISMWavefrontFromStack:NoStackVariable', ...
        'No detector-resolved 4-D stack or scalar 3-D volume was found.');
end

function tf = isDetectorStack(value)
    tf = isnumeric(value) && ndims(value) == 4 && ...
        (size(value, 3) == 23 || size(value, 4) == 23);
    if iscell(value) && numel(value) == 23
        tf = true;
    end
end

function [detStack, channelIDs] = standardizeDetectorZStack(value, defaultChannelIDs)
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
                error('estimateTwoPlaneISMWavefrontFromStack:BadCellStack', ...
                    'All detector cell volumes must have the same [y x z] size.');
            end
            detStack(:,:,k,:) = reshape(vol, ny, nx, 1, nz);
        end
        channelIDs = defaultChannelIDs(:);
        return;
    end

    value = double(value);
    if ndims(value) ~= 4
        error('estimateTwoPlaneISMWavefrontFromStack:BadDetectorStack', ...
            'Detector-resolved full stack must be [y x 23 z] or [y x z 23].');
    end

    if size(value, 3) == 23
        detStack = value;
    elseif size(value, 4) == 23
        detStack = permute(value, [1 2 4 3]);
    else
        error('estimateTwoPlaneISMWavefrontFromStack:BadDetectorStack', ...
            'Could not identify the 23-channel detector dimension.');
    end
    channelIDs = defaultChannelIDs(:);
    if numel(channelIDs) ~= size(detStack, 3)
        channelIDs = (1:size(detStack, 3)).';
    end
end

function [detStack, meta] = readPtuDetectorZStack(fileName, opts)
    if exist('PTU_MultiFrameScanReadFast', 'file') ~= 2 || exist('PTU_Read_Head', 'file') ~= 2
        error('estimateTwoPlaneISMWavefrontFromStack:MissingPtuReader', ...
            'PTU reader functions are not on the MATLAB path.');
    end

    try
        ptuOut = PTU_MultiFrameScanReadFast(fileName, opts.ptuPhotonsPerChunk, false, false);
    catch fastErr
        try
            ptuOut = PTU_MultiFrameScanRead(fileName, opts.ptuPhotonsPerChunk);
        catch slowErr
            error('estimateTwoPlaneISMWavefrontFromStack:PtuReadFailed', ...
                'Could not read %s as a detector z stack. Fast: %s Slow: %s', ...
                fileName, fastErr.message, slowErr.message);
        end
    end
    if isfield(ptuOut, 'dind')
        channelIDs = double(ptuOut.dind(:));
    else
        channelIDs = [];
    end

    if isfield(ptuOut, 'tag') && ~isempty(ptuOut.tag)
        detStack = permute(double(ptuOut.tag), [2 1 3 4]);
    elseif isfield(ptuOut, 'tags') && ~isempty(ptuOut.tags)
        detStack = permute(double(ptuOut.tags), [2 1 3]);
        detStack = reshape(detStack, size(detStack,1), size(detStack,2), size(detStack,3), 1);
    else
        error('estimateTwoPlaneISMWavefrontFromStack:NoPtuDetectorStack', ...
            'No detector image stack was found in %s.', fileName);
    end

    [detStack, channelIDs] = selectChannels(detStack, channelIDs, opts.channelIDs);

    meta = baseMeta('PTU detector stack');
    meta.inputPath = fileName;
    meta.channelIDs = channelIDs(:);
    if isfield(ptuOut, 'head')
        meta.xyPixelSizeUm = numericField(ptuOut.head, 'ImgHdr_PixResol', NaN);
    end
end

function [stack, channelIDs] = selectChannels(stack, channelIDs, requested)
    if isempty(channelIDs)
        channelIDs = (1:size(stack, 3)).';
    end

    requested = requested(:);
    if isempty(requested)
        return;
    end

    [present, loc] = ismember(double(requested), double(channelIDs(:)));
    if all(present)
        stack = stack(:,:,loc,:);
        channelIDs = channelIDs(loc);
    elseif size(stack, 3) == numel(requested)
        channelIDs = requested;
    else
        error('estimateTwoPlaneISMWavefrontFromStack:MissingChannels', ...
            'Could not find all requested detector channel IDs in the stack.');
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
        error('estimateTwoPlaneISMWavefrontFromStack:PlaneZMismatch', ...
            'planeZ has %d entries but the stack has %d z planes.', numel(planeZ), nPlane);
    end
end

function planeZ = planeZFromMatStruct(S)
    planeZ = [];
    if isfield(S, 'frameMeta') && isstruct(S.frameMeta) && isfield(S.frameMeta, 'zUm')
        planeZ = [S.frameMeta.zUm];
    end
end

function xy = xyPixelSizeFromMatStruct(S)
    xy = NaN;
    if isfield(S, 'opts') && isstruct(S.opts) && isfield(S.opts, 'xyPixelSizeUm')
        v = S.opts.xyPixelSizeUm;
        if isnumeric(v) && isscalar(v) && isfinite(v) && v > 0
            xy = double(v);
            return;
        end
    end
    if isfield(S, 'scans') && isstruct(S.scans) && isfield(S.scans, 'xyPixelSizeUm')
        vals = [S.scans.xyPixelSizeUm];
        vals = vals(isfinite(vals) & vals > 0);
        if ~isempty(vals)
            xy = vals(1);
        end
    end
end

function sources = planeSourcesFromMatStruct(S, nPlane)
    sources = emptyPlaneSource(nPlane);
    if ~isfield(S, 'frameMeta') || ~isstruct(S.frameMeta)
        sources = emptyPlaneSource(0);
        return;
    end

    fm = S.frameMeta;
    if numel(fm) < nPlane
        sources = emptyPlaneSource(0);
        return;
    end

    haveAnyFile = false;
    for k = 1:nPlane
        fileName = '';
        frame = [];
        if isfield(fm, 'sourceFile') && ~isempty(fm(k).sourceFile)
            fileName = char(fm(k).sourceFile);
        elseif isfield(fm, 'frameFiles') && ~isempty(fm(k).frameFiles)
            fileName = char(fm(k).frameFiles);
        end
        if isfield(fm, 'sourceFrame') && ~isempty(fm(k).sourceFrame)
            frameValue = double(fm(k).sourceFrame);
            frameValue = frameValue(1);
            if isfinite(frameValue) && frameValue >= 1
                frame = frameValue;
            end
        elseif isfield(fm, 'frameFileFrames') && ~isempty(fm(k).frameFileFrames)
            frameValue = double(fm(k).frameFileFrames);
            frameValue = frameValue(1);
            if isfinite(frameValue) && frameValue >= 1
                frame = frameValue;
            end
        end
        sources(k).file = fileName;
        sources(k).frame = frame;
        haveAnyFile = haveAnyFile || ~isempty(fileName);
    end

    if ~haveAnyFile
        sources = emptyPlaneSource(0);
    end
end

function sources = inferPlaneSourcesFromScanFolder(scanFolder, nPlane)
    files = sortedSeriesFiles(scanFolder);
    if isempty(files)
        sources = emptyPlaneSource(0);
        return;
    end

    if numel(files) == nPlane
        fileNames = fullfile({files.folder}, {files.name});
        sources = makeFileListPlaneSources(fileNames);
    else
        sources = emptyPlaneSource(0);
    end
end

function sources = makeSingleFilePlaneSources(fileName, nPlane)
    sources = emptyPlaneSource(nPlane);
    for k = 1:nPlane
        sources(k).file = fileName;
        sources(k).frame = k;
    end
end

function sources = makeFileListPlaneSources(fileNames)
    sources = emptyPlaneSource(numel(fileNames));
    for k = 1:numel(fileNames)
        sources(k).file = char(fileNames{k});
        sources(k).frame = [];
    end
end

function sources = emptyPlaneSource(n)
    sources = repmat(struct('file', '', 'frame', []), n, 1);
end

function signal = zSignalTrace(stack, metric)
    switch lower(metric)
        case {'signal','sum','total'}
            if ndims(stack) == 4
                nz = size(stack, 4);
                signal = zeros(1, nz);
                for z = 1:nz
                    plane = max(stack(:,:,:,z), 0);
                    signal(z) = sum(plane(:));
                end
            elseif ndims(stack) == 3
                nz = size(stack, 3);
                signal = zeros(1, nz);
                for z = 1:nz
                    plane = max(stack(:,:,z), 0);
                    signal(z) = sum(plane(:));
                end
            else
                error('estimateTwoPlaneISMWavefrontFromStack:BadSignalStack', ...
                    'Signal stack must be 3-D or 4-D.');
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
            error('estimateTwoPlaneISMWavefrontFromStack:BadFocusMetric', ...
                'focusMetric must be ''signal'' or ''peak''.');
    end
end

function signal = ptuFileSignalTrace(fileNames, opts)
    signal = zeros(1, numel(fileNames));
    for k = 1:numel(fileNames)
        signal(k) = ptuFilePhotonSignal(fileNames{k}, opts);
        if opts.verbose
            fprintf('[estimateTwoPlaneISMWavefrontFromStack] signal %d/%d: %.0f photons\n', ...
                k, numel(fileNames), signal(k));
        end
    end
end

function signal = ptuFilePhotonSignal(fileName, opts)
    if exist('PTU_Read_Head', 'file') ~= 2 || exist('PTU_Read', 'file') ~= 2
        error('estimateTwoPlaneISMWavefrontFromStack:MissingPtuReader', ...
            'PTU_Read_Head/PTU_Read are not on the MATLAB path.');
    end

    head = PTU_Read_Head(fileName);
    validPhoton = ptuPhotonValidity(head);
    channelIDs = double(opts.channelIDs(:));

    cnt = 0;
    num = 1;
    signal = 0;
    while num > 0
        [~, tcspc, chan, special, num] = PTU_Read(fileName, [cnt+1 opts.ptuPhotonsPerChunk], head);
        cnt = cnt + num;
        if num == 0
            break;
        end
        keep = special == 0 & validPhoton(chan, tcspc);
        if ~isempty(channelIDs)
            keep = keep & ismember(double(chan), channelIDs);
        end
        signal = signal + sum(keep);
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

function [selection, twoPlaneData] = selectFocusAndDiversityPlanes(stackData, meta, opts)
    signal = stackData.signalTrace(:).';
    nPlane = numel(signal);
    planeZ = resolvePlaneZ(nPlane, meta.planeZ, opts);

    if ~isempty(opts.focusIndex)
        focusIdx = round(opts.focusIndex);
        if focusIdx < 1 || focusIdx > nPlane
            error('estimateTwoPlaneISMWavefrontFromStack:BadFocusIndex', ...
                'focusIndex must be between 1 and %d.', nPlane);
        end
    else
        [~, focusIdx] = max(signal);
    end

    if ~isempty(opts.defocusIndex)
        diversityIdx = round(opts.defocusIndex);
        if diversityIdx < 1 || diversityIdx > nPlane || diversityIdx == focusIdx
            error('estimateTwoPlaneISMWavefrontFromStack:BadDefocusIndex', ...
                'defocusIndex must be a valid plane different from focusIndex.');
        end
    else
        diversityIdx = chooseDiversityPlane(planeZ, focusIdx, opts);
    end

    relativePlaneZ = [0, planeZ(diversityIdx) - planeZ(focusIdx)];

    selection = struct();
    selection.focusIndex = focusIdx;
    selection.diversityIndex = diversityIdx;
    selection.focusZUm = planeZ(focusIdx);
    selection.diversityZUm = planeZ(diversityIdx);
    selection.relativePlaneZ = relativePlaneZ;
    selection.requestedDefocusOffsetUm = opts.defocusOffsetUm;
    selection.signalTrace = signal;
    selection.planeZ = planeZ;

    if ~isempty(stackData.planeSources)
        selection.focusSource = stackData.planeSources(focusIdx);
        selection.diversitySource = stackData.planeSources(diversityIdx);
    else
        selection.focusSource = struct('file', '', 'frame', []);
        selection.diversitySource = struct('file', '', 'frame', []);
    end

    twoPlaneData = struct();
    twoPlaneData.raw4 = [];
    if ~isempty(stackData.detectorStack)
        twoPlaneData.raw4 = cat(4, ...
            stackData.detectorStack(:,:,:,focusIdx), ...
            stackData.detectorStack(:,:,:,diversityIdx));
    end
end

function diversityIdx = chooseDiversityPlane(planeZ, focusIdx, opts)
    offset = abs(opts.defocusOffsetUm);
    if offset <= 0
        error('estimateTwoPlaneISMWavefrontFromStack:BadOffset', ...
            'defocusOffsetUm must be positive.');
    end

    side = opts.defocusSide;
    switch side
        case {'positive','plus','+'}
            signs = [1 -1];
        case {'negative','minus','-'}
            signs = [-1 1];
        case 'auto'
            signs = [1 -1];
        otherwise
            error('estimateTwoPlaneISMWavefrontFromStack:BadDefocusSide', ...
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
        if ~strcmp(side, 'auto') && idx ~= focusIdx
            break;
        end
    end

    if isempty(bestIdx)
        candidates = setdiff(1:numel(planeZ), focusIdx);
        [~, ii] = min(abs(abs(planeZ(candidates) - planeZ(focusIdx)) - offset));
        bestIdx = candidates(ii);
    end
    diversityIdx = bestIdx;
end

function args = buildTwoPlaneArgs(opts, meta, selection)
    args = unmatchedToArgs(opts.unmatched);
    args = setNameValue(args, 'planeZ', selection.relativePlaneZ);
    args = setNameValue(args, 'channelIDs', opts.channelIDs);
    args = setNameValue(args, 'darkFile', opts.darkFile);
    args = setNameValue(args, 'ptuReaderFolder', opts.ptuReaderFolder);
    args = setNameValue(args, 'ptuPhotonsPerChunk', opts.ptuPhotonsPerChunk);
    args = setNameValue(args, 'writeOutputs', opts.writeOutputs);
    args = setNameValue(args, 'verbose', opts.verbose);

    if ~isempty(opts.outputDir)
        args = setNameValue(args, 'outputDir', opts.outputDir);
    else
        args = setNameValue(args, 'outputDir', defaultOutputDir(meta));
    end

    if isfinite(meta.xyPixelSizeUm) && meta.xyPixelSizeUm > 0 && ~hasNameValue(args, 'xyPixelSizeUm')
        args = setNameValue(args, 'xyPixelSizeUm', meta.xyPixelSizeUm);
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

function tf = hasNameValue(args, name)
    tf = false;
    for k = 1:2:numel(args)
        if strcmpi(args{k}, name)
            tf = true;
            return;
        end
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

function outDir = defaultOutputDir(meta)
    thisDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(thisDir);
    label = 'stack';
    if isfield(meta, 'scanFolder') && ~isempty(meta.scanFolder)
        [~, label] = fileparts(stripTrailingFilesep(meta.scanFolder));
    elseif isfield(meta, 'inputPath') && ~isempty(meta.inputPath)
        [~, label] = fileparts(meta.inputPath);
    end
    outDir = fullfile(rootDir, 'output_matlab', 'auto_two_plane_ism_wavefront', sanitizeFileName(label));
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
    isDiversity = false(numel(zUm), 1);
    isFocus(selection.focusIndex) = true;
    isDiversity(selection.diversityIndex) = true;
    T = table((1:numel(zUm)).', zUm, signal, isFocus, isDiversity, ...
        'VariableNames', {'planeIndex','zUm','signal','isFocus','isDiversity'});
    writetable(T, fullfile(outDir, 'auto_selected_planes.csv'));
    save(fullfile(outDir, 'auto_two_plane_ism_wavefront_from_stack.mat'), 'result', '-v7.3');

    if opts.writeSelectionFigure
        writeSelectionFigure(selection, fullfile(outDir, 'auto_selected_planes.png'));
    end
end

function writeSelectionFigure(selection, outFile)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [120 120 720 420]);
    ax = axes(fig);
    plot(ax, selection.planeZ, selection.signalTrace, '-o', 'LineWidth', 1.2);
    hold(ax, 'on');
    plot(ax, selection.focusZUm, selection.signalTrace(selection.focusIndex), ...
        'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(ax, selection.diversityZUm, selection.signalTrace(selection.diversityIndex), ...
        'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    hold(ax, 'off');
    xlabel(ax, 'z (um)');
    ylabel(ax, 'signal');
    title(ax, 'Automatic focus/diversity plane selection');
    legend(ax, {'z trace','focus','diversity'}, 'Location', 'best');
    grid(ax, 'on');
    try
        exportgraphics(fig, outFile, 'Resolution', 180);
    catch
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, outFile, '-dpng', '-r180');
    end
    close(fig);
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

function v = numericField(s, name, defaultValue)
    v = defaultValue;
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name))
        v = double(s.(name)(1));
    end
end
