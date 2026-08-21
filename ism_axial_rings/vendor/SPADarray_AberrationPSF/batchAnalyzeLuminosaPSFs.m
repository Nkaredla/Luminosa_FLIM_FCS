function results = batchAnalyzeLuminosaPSFs(dataRoot, resultsRoot, opts)
%--------------------------------------------------------------------------
% batchAnalyzeLuminosaPSFs
%
% PURPOSE
%   Batch-load Luminosa PSF scans, read measured 3-D PSF volumes from each
%   folder's Series_*.ptu files by default, estimate and correct random
%   lateral frame shifts, and save orthogonal 3-D projection figures in a
%   results folder.
%
% USAGE
%   results = batchAnalyzeLuminosaPSFs();
%   results = batchAnalyzeLuminosaPSFs('D:\Luminosa\Data\ISM_Aberation2_73');
%   results = batchAnalyzeLuminosaPSFs(..., ..., struct('inputSource','pqdat'));
%
% OUTPUTS
%   resultsRoot/   (default: D:\Luminosa\Data\PSF_batch_outputs\<data folder> for PTU input)
%     scan_inventory.csv
%     batch_summary.csv
%     volumes_mat/
%       <scan>_volume_raw.mat
%     xz_yz_plots/
%       <scan>_XY.png
%       <scan>_XZ_YZ.png
%       <scan>_frame_alignment.csv
%
% NOTES
%   The alignment deliberately does not enforce mirror or radial symmetry.
%   These data have intrinsic coma, so every frame is first centered by a
%   robust high-signal support estimate, then optionally refined against an
%   empirical template made from the data itself. This preserves asymmetric
%   PSF structure while correcting random whole-frame software shifts.
%
%   XZ and YZ views use a physical z scale. The default frame spacing is
%   0.05 um, so a 60-frame stack spans 3.0 um in z.
%
%   PTU input uses PTU_Read_Head and PTU_Read from the sibling
%   Luminosa_FLIM_FCS folder. Each Series_*.ptu is read as one z source
%   file, sorted by series number, and binned directly into a [y,x,z]
%   intensity volume before alignment and projection. RawImage.ptu fallback
%   is opt-in with opts.includeRawPtuWhenNoSeries = true.
%--------------------------------------------------------------------------

    addThisFolderToPath();

    if nargin < 1 || isempty(dataRoot)
        dataRoot = 'D:\Luminosa\Data\ISM_Aberation2_73';
%         dataRoot = 'D:\Luminosa\Data\ISM_Aberation_72'; 
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    opts = fillDefaultOptions(opts);
    addPtuReaderPath(opts);

    if nargin < 2 || isempty(resultsRoot)
        resultsRoot = defaultResultsRoot(dataRoot, opts);
    end
    opts = configureOutputLayout(resultsRoot, opts);

    if ~exist(dataRoot, 'dir')
        error('batchAnalyzeLuminosaPSFs:MissingDataRoot', ...
            'Data folder not found: %s', dataRoot);
    end
    ensureDir(resultsRoot);

    scans = discoverLuminosaScans(dataRoot, opts);
    if isempty(scans)
        if strcmpi(opts.inputSource, 'ptu')
            error('batchAnalyzeLuminosaPSFs:NoScans', ...
                'No %s PTU series files found below %s.', opts.ptuSeriesPattern, dataRoot);
        else
            error('batchAnalyzeLuminosaPSFs:NoScans', ...
                'No %s files found below %s.', opts.inputFileName, dataRoot);
        end
    end

    writeScanInventory(scans, fullfile(resultsRoot, 'scan_inventory.csv'));
    groups = buildScanGroups(scans, opts);

    if opts.verbose
        fprintf('[batchAnalyzeLuminosaPSFs] Found %d scans in %d groups.\n', ...
            numel(scans), numel(groups));
        fprintf('[batchAnalyzeLuminosaPSFs] Results folder: %s\n', resultsRoot);
    end

    results = struct();
    results.dataRoot = dataRoot;
    results.resultsRoot = resultsRoot;
    results.options = opts;
    results.scans = scans;
    results.groups = repmat(emptyGroupResult(), 0, 1);

    for g = 1:numel(groups)
        group = groups(g);
        if numel(group.scanIdx) < opts.minFramesPerGroup && ~opts.writeSingleFrameGroups
            continue;
        end

        groupDir = outputGroupDir(resultsRoot, group, opts);
        ensureDir(groupDir);

        if opts.verbose
            fprintf('[batchAnalyzeLuminosaPSFs] Group %d/%d: %s (%d source file(s))\n', ...
                g, numel(groups), group.label, numel(group.scanIdx));
        end

        [rawVol, frameMeta] = loadGroupVolume(scans(group.scanIdx), opts);
        [alignedVol, alignInfo] = alignPsfVolume(rawVol, opts);

        outputFiles = writeGroupOutputs( ...
            rawVol, alignedVol, scans(group.scanIdx), group, frameMeta, ...
            alignInfo, groupDir, opts);

        out = emptyGroupResult();
        out.label = group.label;
        out.key = group.key;
        out.groupDir = groupDir;
        out.numFrames = size(rawVol, 3);
        out.imageSize = [size(rawVol, 1), size(rawVol, 2)];
        groupScale = coordinateScale(scans(group.scanIdx), size(rawVol), opts);
        out.xyPixelSizeUm = groupScale.xyPixelSizeUm;
        out.zStepUm = opts.zStepUm;
        out.zSpanUm = size(rawVol, 3) * opts.zStepUm;
        out.outputFiles = outputFiles;
        out.alignment = alignInfo;
        results.groups(end+1, 1) = out; %#ok<AGROW>
    end

    writeBatchSummary(results.groups, fullfile(resultsRoot, 'batch_summary.csv'));

    if opts.verbose
        fprintf('[batchAnalyzeLuminosaPSFs] Complete. Wrote %d group result folders.\n', ...
            numel(results.groups));
    end
end

function addThisFolderToPath()
    thisDir = fileparts(mfilename('fullpath'));
    if exist(thisDir, 'dir')
        addpath(thisDir);
    end
end

function opts = fillDefaultOptions(opts)
    opts.inputSource = lower(char(getOption(opts, 'inputSource', 'ptu'))); % ptu | pqdat
    if ~ismember(opts.inputSource, {'ptu', 'pqdat'})
        error('batchAnalyzeLuminosaPSFs:BadInputSource', ...
            'opts.inputSource must be ''ptu'' or ''pqdat''.');
    end

    if isfield(opts, 'inputFileName') && ~isempty(opts.inputFileName)
        opts.inputFileName = char(opts.inputFileName);
    elseif strcmp(opts.inputSource, 'ptu')
        opts.inputFileName = 'RawImage.ptu';
    else
        opts.inputFileName = 'IntensityImage.pqdat';
    end

    if isfield(opts, 'groupMode') && ~isempty(opts.groupMode)
        opts.groupMode = lower(char(opts.groupMode));
    elseif strcmp(opts.inputSource, 'ptu')
        opts.groupMode = 'folder';
    else
        opts.groupMode = 'metadata';
    end

    opts.ptuReaderFolder = getOption(opts, 'ptuReaderFolder', defaultPtuReaderFolder());
    opts.ptuReader = char(getOption(opts, 'ptuReader', 'localIntensity'));
    opts.ptuSeriesPattern = char(getOption(opts, 'ptuSeriesPattern', 'Series*.ptu'));
    opts.preferPtuSeries = getOption(opts, 'preferPtuSeries', true);
    opts.includeRawPtuWhenNoSeries = getOption(opts, 'includeRawPtuWhenNoSeries', false);
    opts.skipBadPtuFiles = getOption(opts, 'skipBadPtuFiles', true);
    opts.ptuPhotonsPerChunk = getOption(opts, 'ptuPhotonsPerChunk', 1e6);
    opts.ptuStoreTcspcPix = getOption(opts, 'ptuStoreTcspcPix', false);
    opts.ptuUseGPU = getOption(opts, 'ptuUseGPU', false);
    opts.minFramesPerGroup = getOption(opts, 'minFramesPerGroup', 1);
    opts.writeSingleFrameGroups = getOption(opts, 'writeSingleFrameGroups', true);
    opts.centerThresholdFractions = getOption(opts, 'centerThresholdFractions', ...
        [0.18 0.30 0.45 0.60]);
    opts.smoothSigmaPx = getOption(opts, 'smoothSigmaPx', 0.8);
    opts.targetCenter = getOption(opts, 'targetCenter', 'median'); % median | imageCenter | [x y]
    opts.smoothCenterTrajectory = getOption(opts, 'smoothCenterTrajectory', true);
    opts.centerSmoothingWindow = getOption(opts, 'centerSmoothingWindow', 9);
    opts.centerOutlierThresholdPx = getOption(opts, 'centerOutlierThresholdPx', 2.5);
    opts.minCenterScoreForTrajectory = getOption(opts, 'minCenterScoreForTrajectory', 0);
    opts.preserveSmoothedCenterTrajectory = getOption(opts, 'preserveSmoothedCenterTrajectory', true);
    opts.secondPassAlignment = getOption(opts, 'secondPassAlignment', true);
    opts.secondPassSmoothingWindow = getOption(opts, 'secondPassSmoothingWindow', 11);
    opts.maxSecondPassShiftPx = getOption(opts, 'maxSecondPassShiftPx', 1.5);
    opts.refineWithEmpiricalTemplate = getOption(opts, 'refineWithEmpiricalTemplate', false);
    opts.maxTemplateShiftPx = getOption(opts, 'maxTemplateShiftPx', 2);
    opts.minTemplatePeak = getOption(opts, 'minTemplatePeak', 0.04);
    opts.projectionMode = getOption(opts, 'projectionMode', 'max'); % max | sum
    opts.zStepUm = getOption(opts, 'zStepUm', 0.05);
    opts.xyPixelSizeUm = getOption(opts, 'xyPixelSizeUm', []);
    opts.inferXYPixelSizeFromSpqr = getOption(opts, 'inferXYPixelSizeFromSpqr', true);
    opts.makeIsosurface = getOption(opts, 'makeIsosurface', true);
    opts.isosurfaceLevel = getOption(opts, 'isosurfaceLevel', 0.25);
    opts.flatOutputLayout = getOption(opts, 'flatOutputLayout', strcmp(opts.inputSource, 'ptu'));
    opts.writeFigFiles = getOption(opts, 'writeFigFiles', true);
    opts.verbose = getOption(opts, 'verbose', true);
end

function resultsRoot = defaultResultsRoot(dataRoot, opts)
    if strcmpi(opts.inputSource, 'ptu')
        cleanRoot = regexprep(char(dataRoot), '[\\/]+$', '');
        [dataParent, dataName] = fileparts(cleanRoot);
        if isempty(dataName)
            dataName = 'psf_batch';
        end
        resultsRoot = fullfile(dataParent, 'PSF_batch_outputs', sanitizeFileName(dataName));
    else
        resultsRoot = fullfile(dataRoot, 'results_psf3d');
    end
end

function opts = configureOutputLayout(resultsRoot, opts)
    if opts.flatOutputLayout
        opts.volumeOutputDir = getOption(opts, 'volumeOutputDir', ...
            fullfile(resultsRoot, 'volumes_mat'));
        opts.projectionOutputDir = getOption(opts, 'projectionOutputDir', ...
            fullfile(resultsRoot, 'xz_yz_plots'));
        ensureDir(opts.volumeOutputDir);
        ensureDir(opts.projectionOutputDir);
    else
        opts.volumeOutputDir = '';
        opts.projectionOutputDir = '';
    end
end

function groupDir = outputGroupDir(resultsRoot, group, opts)
    if opts.flatOutputLayout
        groupDir = resultsRoot;
    else
        groupDir = fullfile(resultsRoot, sanitizeFileName(group.label));
    end
end

function folder = defaultPtuReaderFolder()
    thisDir = fileparts(mfilename('fullpath'));
    aberrationRoot = fileparts(thisDir);
    luminosaRoot = fileparts(aberrationRoot);
    folder = fullfile(luminosaRoot, 'Luminosa_FLIM_FCS');
end

function addPtuReaderPath(opts)
    if ~strcmpi(opts.inputSource, 'ptu')
        return;
    end

    readerFolder = char(opts.ptuReaderFolder);
    if exist(readerFolder, 'dir') ~= 7
        error('batchAnalyzeLuminosaPSFs:MissingPtuReaderFolder', ...
            'PTU reader folder not found: %s', readerFolder);
    end

    addpath(readerFolder);
    if ~strcmpi(opts.ptuReader, 'localIntensity') && exist(opts.ptuReader, 'file') ~= 2
        error('batchAnalyzeLuminosaPSFs:MissingPtuReader', ...
            'PTU reader function %s was not found on the MATLAB path.', opts.ptuReader);
    end
    if exist('PTU_Read_Head', 'file') ~= 2
        error('batchAnalyzeLuminosaPSFs:MissingPtuHeaderReader', ...
            'PTU_Read_Head was not found on the MATLAB path.');
    end
    if exist('PTU_Read', 'file') ~= 2
        error('batchAnalyzeLuminosaPSFs:MissingPtuRecordReader', ...
            'PTU_Read was not found on the MATLAB path.');
    end
end

function v = getOption(opts, name, defaultValue)
    if isfield(opts, name) && ~isempty(opts.(name))
        v = opts.(name);
    else
        v = defaultValue;
    end
end

function scans = discoverLuminosaScans(dataRoot, opts)
    d = dir(dataRoot);
    scans = repmat(emptyScan(), 0, 1);

    [rootFiles, rootFileMode] = discoverScanInputFiles(dataRoot, opts);
    if ~isempty(rootFiles)
        scans = scanFromDirectory(dataRoot, leafFolderName(dataRoot), ...
            rootFiles, rootFileMode, opts);
        return;
    end

    for k = 1:numel(d)
        if ~d(k).isdir || any(strcmp(d(k).name, {'.', '..'}))
            continue;
        end

        scanDir = fullfile(dataRoot, d(k).name);
        [fileNames, fileMode] = discoverScanInputFiles(scanDir, opts);
        if isempty(fileNames)
            continue;
        end
        s = scanFromDirectory(scanDir, d(k).name, fileNames, fileMode, opts);
        if ~isempty(s)
            scans(end+1, 1) = s; %#ok<AGROW>
        end
    end
end

function s = scanFromDirectory(scanDir, scanName, fileNames, fileMode, opts)
    fileName = fileNames{1};

    try
        imageMeta = readScanMetadata(fileName, scanDir, opts);
    catch err
        warning('batchAnalyzeLuminosaPSFs:BadInputFile', ...
            'Skipping %s: %s', fileName, err.message);
        s = [];
        return;
    end

    parsed = parseScanFolderName(scanName);

    s = emptyScan();
    s.name = scanName;
    s.folder = scanDir;
    s.file = fileName;
    s.files = fileNames(:).';
    s.fileMode = fileMode;
    s.numSourceFiles = numel(fileNames);
    s.inputSource = opts.inputSource;
    s.power = parsed.power;
    s.collarText = parsed.collarText;
    s.collar = parsed.collar;
    s.lens = parsed.lens;
    s.timestamp = parsed.timestamp;
    s.pixelX = imageMeta.pixelX;
    s.pixelY = imageMeta.pixelY;
    s.xyPixelSizeUm = resolveScanXYPixelSizeUm(scanDir, imageMeta, opts);
    s.numFrames = estimateScanFrameCount(imageMeta, fileNames, opts);
    s.blockCount = imageMeta.blockCount;
    s.pixelatedValueBytes = imageMeta.pixelatedValueBytes;
    s.groupKey = makeGroupKey(s);
end

function s = emptyScan()
    s = struct( ...
        'name', '', ...
        'folder', '', ...
        'file', '', ...
        'files', {{}}, ...
        'fileMode', '', ...
        'numSourceFiles', 0, ...
        'inputSource', '', ...
        'power', '', ...
        'collarText', '', ...
        'collar', NaN, ...
        'lens', '', ...
        'timestamp', '', ...
        'pixelX', NaN, ...
        'pixelY', NaN, ...
        'xyPixelSizeUm', NaN, ...
        'numFrames', NaN, ...
        'blockCount', NaN, ...
        'pixelatedValueBytes', NaN, ...
        'groupKey', '');
end

function [fileNames, fileMode] = discoverScanInputFiles(scanDir, opts)
    fileNames = {};
    fileMode = '';

    switch lower(opts.inputSource)
        case 'ptu'
            if opts.preferPtuSeries
                seriesFiles = dir(fullfile(scanDir, opts.ptuSeriesPattern));
                seriesFiles = seriesFiles(~[seriesFiles.isdir]);
                if ~isempty(seriesFiles)
                    seriesFiles = sortSeriesFiles(seriesFiles);
                    fileNames = fullfile(scanDir, {seriesFiles.name});
                    fileNames = filterReadablePtuFiles(fileNames, opts);
                    fileMode = 'ptu_series';
                    return;
                end
            end

            rawFile = fullfile(scanDir, opts.inputFileName);
            if opts.includeRawPtuWhenNoSeries && exist(rawFile, 'file') == 2
                fileNames = {rawFile};
                fileNames = filterReadablePtuFiles(fileNames, opts);
                fileMode = 'ptu_raw';
            end

        otherwise
            fileName = fullfile(scanDir, opts.inputFileName);
            if exist(fileName, 'file') == 2
                fileNames = {fileName};
                fileMode = opts.inputSource;
            end
    end
end

function fileNames = filterReadablePtuFiles(fileNames, opts)
    if isempty(fileNames) || ~opts.skipBadPtuFiles
        return;
    end

    keep = false(size(fileNames));
    for k = 1:numel(fileNames)
        try
            head = PTU_Read_Head(fileNames{k});
            keep(k) = ~isempty(head);
        catch err
            warning('batchAnalyzeLuminosaPSFs:SkippingBadPtuHeader', ...
                'Skipping unreadable PTU file %s: %s', fileNames{k}, err.message);
        end
    end

    fileNames = fileNames(keep);
end

function files = sortSeriesFiles(files)
    seriesNumber = nan(numel(files), 1);
    for k = 1:numel(files)
        tok = regexp(files(k).name, '(\d+)(?=\.ptu$)', 'tokens', 'once');
        if ~isempty(tok)
            seriesNumber(k) = str2double(tok{1});
        end
    end

    if all(isfinite(seriesNumber))
        [~, order] = sort(seriesNumber);
    else
        [~, order] = sort({files.name});
    end
    files = files(order);
end

function numFrames = estimateScanFrameCount(imageMeta, fileNames, opts)
    numFrames = imageMeta.numFrames;
    if strcmpi(opts.inputSource, 'ptu') && numel(fileNames) > 1
        framesPerFile = imageMeta.numFrames;
        if ~isfinite(framesPerFile) || framesPerFile < 1
            framesPerFile = 1;
        end
        numFrames = numel(fileNames) * framesPerFile;
    end
end

function meta = readScanMetadata(fileName, scanDir, opts)
    switch lower(opts.inputSource)
        case 'ptu'
            meta = readPtuMetadata(fileName);
        case 'pqdat'
            [~, meta] = readPqdatImage(fileName, false);
            meta.inputSource = 'pqdat';
            meta.xyPixelSizeUm = NaN;
        otherwise
            error('batchAnalyzeLuminosaPSFs:BadInputSource', ...
                'Unsupported input source: %s', opts.inputSource);
    end

    if ~isfield(meta, 'xyPixelSizeUm') || isempty(meta.xyPixelSizeUm)
        meta.xyPixelSizeUm = NaN;
    end
    if ~isfield(meta, 'blockCount') || isempty(meta.blockCount)
        meta.blockCount = NaN;
    end
    if ~isfield(meta, 'pixelatedValueBytes') || isempty(meta.pixelatedValueBytes)
        meta.pixelatedValueBytes = NaN;
    end
    if ~isfield(meta, 'numFrames') || isempty(meta.numFrames)
        meta.numFrames = NaN;
    end

    if ~isfinite(meta.xyPixelSizeUm) || meta.xyPixelSizeUm <= 0
        meta.xyPixelSizeUm = resolveXYPixelSizeUm(scanDir, opts);
    end
end

function meta = readPtuMetadata(fileName)
    head = PTU_Read_Head(fileName);
    if isempty(head)
        error('batchAnalyzeLuminosaPSFs:BadPtuHeader', ...
            'Could not read PTU header from %s.', fileName);
    end

    meta = struct();
    meta.file = fileName;
    meta.inputSource = 'ptu';
    meta.head = head;
    meta.pixelX = numericField(head, 'ImgHdr_PixX', NaN);
    meta.pixelY = numericField(head, 'ImgHdr_PixY', NaN);
    meta.numFrames = numericField(head, 'ImgHdr_MaxFrames', NaN);
    meta.xyPixelSizeUm = numericField(head, 'ImgHdr_PixResol', NaN);
    meta.blockCount = NaN;
    meta.pixelatedValueBytes = NaN;

    if ~isfinite(meta.numFrames) || meta.numFrames <= 0
        timePerPixel = numericField(head, 'ImgHdr_TimePerPixel', NaN);
        stopAfter = numericField(head, 'TTResult_StopAfter', NaN);
        if isfinite(timePerPixel) && isfinite(stopAfter) && ...
                isfinite(meta.pixelX) && isfinite(meta.pixelY) && ...
                timePerPixel > 0 && stopAfter > 0
            meta.numFrames = ceil(stopAfter / (timePerPixel * meta.pixelX * meta.pixelY));
        end
    end
end

function value = numericField(s, fieldName, defaultValue)
    value = defaultValue;
    if isfield(s, fieldName) && ~isempty(s.(fieldName)) && isnumeric(s.(fieldName))
        value = double(s.(fieldName)(1));
    end
end

function xyPixelSizeUm = resolveScanXYPixelSizeUm(scanDir, imageMeta, opts)
    if isnumeric(opts.xyPixelSizeUm) && isscalar(opts.xyPixelSizeUm) && ...
            isfinite(opts.xyPixelSizeUm) && opts.xyPixelSizeUm > 0
        xyPixelSizeUm = double(opts.xyPixelSizeUm);
        return;
    end

    xyPixelSizeUm = NaN;
    if isfield(imageMeta, 'xyPixelSizeUm') && isfinite(imageMeta.xyPixelSizeUm) && ...
            imageMeta.xyPixelSizeUm > 0
        xyPixelSizeUm = double(imageMeta.xyPixelSizeUm);
        return;
    end

    xyPixelSizeUm = resolveXYPixelSizeUm(scanDir, opts);
end

function xyPixelSizeUm = resolveXYPixelSizeUm(scanDir, opts)
    xyPixelSizeUm = NaN;

    if isnumeric(opts.xyPixelSizeUm) && isscalar(opts.xyPixelSizeUm) && ...
            isfinite(opts.xyPixelSizeUm) && opts.xyPixelSizeUm > 0
        xyPixelSizeUm = double(opts.xyPixelSizeUm);
        return;
    end

    if ~opts.inferXYPixelSizeFromSpqr
        return;
    end

    spqrFile = fullfile(scanDir, 'RawImage.spqr');
    if exist(spqrFile, 'file') ~= 2
        return;
    end

    try
        spqrMeta = readSpqrScanMetadata(spqrFile);
        if isfield(spqrMeta, 'pixelResolutionUm') && ...
                isfinite(spqrMeta.pixelResolutionUm) && spqrMeta.pixelResolutionUm > 0
            xyPixelSizeUm = spqrMeta.pixelResolutionUm;
        end
    catch err
        warning('batchAnalyzeLuminosaPSFs:BadSpqrMetadata', ...
            'Could not read pixel size from %s: %s', spqrFile, err.message);
    end
end

function meta = readSpqrScanMetadata(fileName)
    fid = fopen(fileName, 'r', 'ieee-le');
    if fid < 0
        error('batchAnalyzeLuminosaPSFs:OpenFailed', ...
            'Could not open %s.', fileName);
    end
    cleanup = onCleanup(@() fclose(fid));

    magic = fread(fid, 8, '*char')';
    fread(fid, 8, '*char');
    if numel(magic) < 6 || ~strcmp(char(magic(1:6)), 'PQSPQR')
        error('batchAnalyzeLuminosaPSFs:NotSpqr', ...
            '%s is not a PQSPQR file.', fileName);
    end

    meta = struct('pixelResolutionUm', NaN, 'pixelX', NaN, 'pixelY', NaN);
    nGuard = 0;

    while ~feof(fid)
        nGuard = nGuard + 1;
        if nGuard > 10000
            error('batchAnalyzeLuminosaPSFs:TooManyTags', ...
                'Too many tags while parsing %s.', fileName);
        end

        identRaw = fread(fid, 32, '*char')';
        if numel(identRaw) < 32
            break;
        end
        ident = trimNulls(char(identRaw));
        fread(fid, 1, 'int32=>double');
        tagType = fread(fid, 1, 'uint32=>uint32');
        valueBytes = fread(fid, 8, 'uint8=>uint8')';
        if numel(valueBytes) < 8
            break;
        end

        valueInt64 = typecast(valueBytes, 'int64');
        valueDouble = typecast(valueBytes, 'double');
        nBytes = double(valueInt64);

        switch tagType
            case uint32(hex2dec('10000008'))
                value = double(valueInt64);
            case {uint32(hex2dec('20000008')), uint32(hex2dec('21000008'))}
                value = double(valueDouble);
            case {uint32(hex2dec('4001FFFF')), uint32(hex2dec('4002FFFF')), ...
                    uint32(hex2dec('FFFFFFFF')), uint32(hex2dec('2001FFFF')), ...
                    uint32(hex2dec('1001FFFF'))}
                fseek(fid, nBytes, 'cof');
                value = [];
            otherwise
                value = [];
        end

        switch ident
            case 'ImgHdr_PixResol'
                meta.pixelResolutionUm = value;
            case 'ImgHdr_PixX'
                meta.pixelX = value;
            case 'ImgHdr_PixY'
                meta.pixelY = value;
        end

        if strcmp(ident, 'Header_End')
            break;
        end
        if isfinite(meta.pixelResolutionUm) && meta.pixelResolutionUm > 0 && ...
                isfinite(meta.pixelX) && isfinite(meta.pixelY)
            break;
        end
    end
end

function parsed = parseScanFolderName(folderName)
    parsed = struct('power', 'unknownPower', 'collarText', 'unknown', ...
        'collar', NaN, 'lens', 'unknownLens', 'timestamp', folderName);

    tok = regexp(folderName, ...
        '^([^_]+)_([0-9.]+)collar_([^_]+)_(\d{8}-\d{6})$', ...
        'tokens', 'once');

    if isempty(tok)
        return;
    end

    parsed.power = tok{1};
    parsed.collarText = tok{2};
    parsed.collar = str2double(tok{2});
    parsed.lens = tok{3};
    parsed.timestamp = tok{4};
end

function key = makeGroupKey(s)
    key = sprintf('%s|%s|%s|%dx%d', ...
        s.power, s.collarText, s.lens, s.pixelY, s.pixelX);
end

function groups = buildScanGroups(scans, opts)
    groups = repmat(struct('key', '', 'label', '', 'scanIdx', []), 0, 1);

    if strcmpi(opts.inputSource, 'ptu') || ...
            any(strcmpi(opts.groupMode, {'folder', 'scan', 'file'}))
        for k = 1:numel(scans)
            g = struct();
            g.key = scans(k).name;
            g.label = scans(k).name;
            g.scanIdx = k;
            groups(end+1, 1) = g; %#ok<AGROW>
        end
        return;
    end

    for k = 1:numel(scans)
        key = scans(k).groupKey;
        idx = find(strcmp(key, {groups.key}), 1, 'first');

        if isempty(idx)
            g = struct();
            g.key = key;
            g.label = sprintf('%s_%scollar_%s_%dx%d', ...
                scans(k).power, scans(k).collarText, scans(k).lens, ...
                scans(k).pixelY, scans(k).pixelX);
            g.scanIdx = k;
            groups(end+1, 1) = g; %#ok<AGROW>
        else
            groups(idx).scanIdx(end+1) = k;
        end
    end

    for g = 1:numel(groups)
        names = {scans(groups(g).scanIdx).timestamp};
        [~, order] = sort(names);
        groups(g).scanIdx = groups(g).scanIdx(order);
    end
end

function [vol, frameMeta] = loadGroupVolume(scans, opts)
    vol = [];
    frameMeta = repmat(struct('sourceName', '', 'sourceFile', '', ...
        'sourceFrame', NaN, 'zUm', NaN, 'totalSignal', NaN), 0, 1);

    for k = 1:numel(scans)
        [img, volumeMeta] = readScanVolume(scans(k), opts);

        for z = 1:size(img, 3)
            frame = double(img(:, :, z));
            if isempty(vol)
                vol = zeros(size(frame, 1), size(frame, 2), 0);
            elseif size(frame, 1) ~= size(vol, 1) || size(frame, 2) ~= size(vol, 2)
                error('batchAnalyzeLuminosaPSFs:GroupSizeMismatch', ...
                    'Grouped images do not have matching dimensions.');
            end

            vol(:, :, end+1) = frame; %#ok<AGROW>
            fm = struct();
            fm.sourceName = scans(k).name;
            [fm.sourceFile, fm.sourceFrame] = frameSourceInfo(volumeMeta, z, scans(k).file);
            fm.zUm = NaN;
            fm.totalSignal = sum(frame(:));
            frameMeta(end+1, 1) = fm; %#ok<AGROW>
        end

        if opts.verbose
            fprintf('  loaded %s (%d z plane(s) from %d %s file(s))\n', ...
                scans(k).name, size(img, 3), scans(k).numSourceFiles, scans(k).inputSource);
        end
    end

    zUm = zCoordinatesUm(numel(frameMeta), opts);
    for k = 1:numel(frameMeta)
        frameMeta(k).zUm = zUm(k);
    end
end

function [sourceFile, sourceFrame] = frameSourceInfo(volumeMeta, z, fallbackFile)
    sourceFile = fallbackFile;
    sourceFrame = z;
    if isfield(volumeMeta, 'frameFiles') && numel(volumeMeta.frameFiles) >= z
        sourceFile = volumeMeta.frameFiles{z};
    end
    if isfield(volumeMeta, 'frameFileFrames') && numel(volumeMeta.frameFileFrames) >= z
        sourceFrame = volumeMeta.frameFileFrames(z);
    end
end

function [img, meta] = readScanVolume(scan, opts)
    switch lower(opts.inputSource)
        case 'ptu'
            [img, meta] = readPtuScanVolume(scan, opts);
        case 'pqdat'
            [img, meta] = readPqdatImage(scan.file, true);
            meta.inputSource = 'pqdat';
            meta.frameFiles = repmat({scan.file}, 1, size(img, 3));
            meta.frameFileFrames = 1:size(img, 3);
        otherwise
            error('batchAnalyzeLuminosaPSFs:BadInputSource', ...
                'Unsupported input source: %s', opts.inputSource);
    end
end

function [vol, meta] = readPtuScanVolume(scan, opts)
    fileNames = scan.files;
    if isempty(fileNames)
        fileNames = {scan.file};
    end

    vol = [];
    frameFiles = {};
    frameFileFrames = [];
    firstMeta = struct();
    skippedFiles = {};

    for f = 1:numel(fileNames)
        try
            [thisVol, thisMeta] = readPtuVolume(fileNames{f}, opts);
        catch err
            if opts.skipBadPtuFiles
                warning('batchAnalyzeLuminosaPSFs:SkippingBadPtuRead', ...
                    'Skipping PTU file %s after read failure: %s', fileNames{f}, err.message);
                skippedFiles{end+1} = fileNames{f}; %#ok<AGROW>
                continue;
            end
            rethrow(err);
        end

        if isempty(fieldnames(firstMeta))
            firstMeta = thisMeta;
        end

        for z = 1:size(thisVol, 3)
            frame = double(thisVol(:, :, z));
            if isempty(vol)
                vol = zeros(size(frame, 1), size(frame, 2), 0);
            elseif size(frame, 1) ~= size(vol, 1) || size(frame, 2) ~= size(vol, 2)
                error('batchAnalyzeLuminosaPSFs:PtuSeriesSizeMismatch', ...
                    'PTU series image size mismatch in %s.', fileNames{f});
            end
            vol(:, :, end+1) = frame; %#ok<AGROW>
            frameFiles{end+1} = fileNames{f}; %#ok<AGROW>
            frameFileFrames(end+1) = z; %#ok<AGROW>
        end
    end

    if isempty(vol)
        error('batchAnalyzeLuminosaPSFs:NoReadablePtuPlanes', ...
            'No readable PTU image planes found in %s.', scan.folder);
    end

    meta = firstMeta;
    meta.file = scan.file;
    meta.files = fileNames;
    meta.fileMode = scan.fileMode;
    meta.numSourceFiles = numel(fileNames);
    meta.skippedFiles = skippedFiles;
    meta.numSkippedFiles = numel(skippedFiles);
    meta.numFrames = size(vol, 3);
    meta.frameFiles = frameFiles;
    meta.frameFileFrames = frameFileFrames;
end

function [vol, meta] = readPtuVolume(fileName, opts)
    reader = char(opts.ptuReader);
    switch reader
        case 'localIntensity'
            [vol, meta] = readPtuIntensityVolume(fileName, opts);
            return;
        case 'PTU_MultiFrameScanReadFast'
            try
                ptuOut = PTU_MultiFrameScanReadFast(fileName, ...
                    opts.ptuPhotonsPerChunk, opts.ptuStoreTcspcPix, opts.ptuUseGPU);
            catch err
                if contains(err.message, 'Out of range subscript') || ...
                        contains(err.message, 'sub2ind')
                    warning('batchAnalyzeLuminosaPSFs:PtuFastReaderFallback', ...
                        ['PTU_MultiFrameScanReadFast failed while binning %s. ', ...
                         'Falling back to local intensity reader.'], fileName);
                    [vol, meta] = readPtuIntensityVolume(fileName, opts);
                    return;
                end
                rethrow(err);
            end
        case 'PTU_MultiFrameScanRead'
            ptuOut = PTU_MultiFrameScanRead(fileName, opts.ptuPhotonsPerChunk);
        otherwise
            ptuOut = feval(reader, fileName, opts.ptuPhotonsPerChunk);
    end

    [vol, meta] = ptuOutputToVolume(ptuOut, fileName);
end

function [vol, meta] = readPtuIntensityVolume(fileName, opts)
    head = PTU_Read_Head(fileName);
    if isempty(head)
        error('batchAnalyzeLuminosaPSFs:BadPtuHeader', ...
            'Could not read PTU header from %s.', fileName);
    end

    nx = numericField(head, 'ImgHdr_PixX', NaN);
    ny = numericField(head, 'ImgHdr_PixY', NaN);
    if ~isfinite(nx) || ~isfinite(ny) || nx <= 0 || ny <= 0
        error('batchAnalyzeLuminosaPSFs:MissingPtuDimensions', ...
            'Missing ImgHdr_PixX/ImgHdr_PixY in %s.', fileName);
    end
    nx = round(nx);
    ny = round(ny);

    nzExpected = estimatePtuFrameCount(head, nx, ny);
    if ~isfinite(nzExpected) || nzExpected <= 0
        nzExpected = 1;
    end
    nzExpected = max(1, round(nzExpected));
    vol = zeros(ny, nx, nzExpected);

    lineStart = ptuMarkerMask(head, 'ImgHdr_LineStart', 4);
    lineStop = ptuMarkerMask(head, 'ImgHdr_LineStop', 2);
    frameMarker = ptuMarkerMask(head, 'ImgHdr_Frame', 3);
    validPhoton = ptuPhotonValidity(head);

    cnt = 0;
    tend = 0;
    num = 1;
    frameIdx = 1;
    yCarry = [];
    markerCarry = uint8([]);
    turnsStart = [];
    turnsStop = [];

    while num > 0
        [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = ...
            PTU_Read(fileName, [cnt+1 opts.ptuPhotonsPerChunk], head);
        cnt = cnt + num;
        if num == 0
            break;
        end

        if ~isempty(yCarry) || ~isempty(turnsStart) || ~isempty(turnsStop)
            tmpy = tmpy + tend;
        end

        keep = (tmpmarkers > 0) | validPhoton(tmpchan, tmptcspc);
        yCarry = [yCarry; tmpy(keep)]; %#ok<AGROW>
        markerCarry = [markerCarry; uint8(tmpmarkers(keep))]; %#ok<AGROW>

        [turnsStart, turnsStop] = appendLineMarkers( ...
            yCarry, markerCarry, turnsStart, turnsStop, lineStart, lineStop);
        frameChange = yCarry(markerCarry == frameMarker);

        isMarker = markerCarry ~= 0;
        yCarry(isMarker) = [];
        markerCarry(isMarker) = [];

        if ~isempty(yCarry)
            tend = yCarry(end) + loc;
        else
            tend = loc;
        end

        for kf = 1:numel(frameChange)
            frameEnd = frameChange(kf);
            [frameImg, yCarry, turnsStart, turnsStop] = consumePtuFrame( ...
                yCarry, turnsStart, turnsStop, frameEnd, nx, ny);
            if frameIdx > size(vol, 3)
                vol(:, :, frameIdx) = 0;
            end
            vol(:, :, frameIdx) = frameImg;
            frameIdx = frameIdx + 1;
        end
    end

    if ~isempty(yCarry) || frameIdx == 1
        frameImg = binPtuFramePhotons(yCarry, turnsStart, turnsStop, nx, ny);
        if any(frameImg(:)) || frameIdx == 1
            if frameIdx > size(vol, 3)
                vol(:, :, frameIdx) = 0;
            end
            vol(:, :, frameIdx) = frameImg;
            frameIdx = frameIdx + 1;
        end
    end

    lastFrame = frameIdx - 1;
    if lastFrame < 1
        error('batchAnalyzeLuminosaPSFs:NoPtuFrames', ...
            'No image frames could be reconstructed from %s.', fileName);
    end
    vol = vol(:, :, 1:lastFrame);

    meta = struct();
    meta.file = fileName;
    meta.inputSource = 'ptu';
    meta.head = head;
    meta.pixelY = ny;
    meta.pixelX = nx;
    meta.numFrames = size(vol, 3);
    meta.xyPixelSizeUm = numericField(head, 'ImgHdr_PixResol', NaN);
    meta.blockCount = NaN;
    meta.pixelatedValueBytes = NaN;
    meta.activeChannels = [];
    meta.frameFiles = repmat({fileName}, 1, size(vol, 3));
    meta.frameFileFrames = 1:size(vol, 3);
end

function nz = estimatePtuFrameCount(head, nx, ny)
    nz = numericField(head, 'ImgHdr_MaxFrames', NaN);
    if isfinite(nz) && nz > 0
        return;
    end

    timePerPixel = numericField(head, 'ImgHdr_TimePerPixel', NaN);
    stopAfter = numericField(head, 'TTResult_StopAfter', NaN);
    if isfinite(timePerPixel) && isfinite(stopAfter) && ...
            timePerPixel > 0 && stopAfter > 0
        nz = ceil(stopAfter / (timePerPixel * nx * ny));
    end
end

function marker = ptuMarkerMask(head, fieldName, defaultValue)
    marker = defaultValue;
    markerIndex = numericField(head, fieldName, NaN);
    if isfinite(markerIndex) && markerIndex > 0
        marker = 2^(markerIndex - 1);
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

function [turnsStart, turnsStop] = appendLineMarkers( ...
        yCarry, markerCarry, turnsStart, turnsStop, lineStart, lineStop)
    if lineStart == lineStop
        turns = yCarry(markerCarry == lineStart);
        if numel(turnsStart) > numel(turnsStop)
            turnsStart = [turnsStart; turns(2:2:end)]; %#ok<AGROW>
            turnsStop = [turnsStop; turns(1:2:end)]; %#ok<AGROW>
        else
            turnsStart = [turnsStart; turns(1:2:end)]; %#ok<AGROW>
            turnsStop = [turnsStop; turns(2:2:end)]; %#ok<AGROW>
        end
    else
        turnsStart = [turnsStart; yCarry(markerCarry == lineStart)]; %#ok<AGROW>
        turnsStop = [turnsStop; yCarry(markerCarry == lineStop)]; %#ok<AGROW>
    end
end

function [frameImg, yCarry, turnsStart, turnsStop] = consumePtuFrame( ...
        yCarry, turnsStart, turnsStop, frameEnd, nx, ny)
    inFrame = yCarry < frameEnd;
    yf = yCarry(inFrame);
    yCarry(inFrame) = [];

    startsInFrame = turnsStart < frameEnd;
    stopsInFrame = turnsStop < frameEnd;
    starts = turnsStart(startsInFrame);
    stops = turnsStop(stopsInFrame);
    turnsStart(startsInFrame) = [];
    turnsStop(stopsInFrame) = [];

    frameImg = binPtuFramePhotons(yf, starts, stops, nx, ny);
end

function frameImg = binPtuFramePhotons(yf, starts, stops, nx, ny)
    frameImg = zeros(ny, nx);
    [startsUse, stopsUse] = pairPtuStartStop(starts, stops);
    if isempty(yf) || isempty(startsUse) || isempty(stopsUse)
        return;
    end

    yf = double(yf(:));
    startsUse = double(startsUse(:));
    stopsUse = double(stopsUse(:));

    inSpan = (yf >= startsUse(1)) & (yf <= stopsUse(end));
    yf = yf(inSpan);
    if isempty(yf)
        return;
    end

    lineIdx = discretize(yf, [startsUse; inf]);
    ok = ~isnan(lineIdx) & lineIdx >= 1 & ...
        lineIdx <= numel(stopsUse) & lineIdx <= ny;
    yf = yf(ok);
    lineIdx = lineIdx(ok);
    if isempty(yf)
        return;
    end

    okStop = yf <= stopsUse(lineIdx);
    yf = yf(okStop);
    lineIdx = lineIdx(okStop);
    if isempty(yf)
        return;
    end

    dwell = max(stopsUse(lineIdx) - startsUse(lineIdx), 1);
    frac = (yf - startsUse(lineIdx)) ./ dwell;
    col = 1 + floor(nx .* frac);
    col = min(nx, max(1, col));

    pixIdx = sub2ind([ny, nx], double(lineIdx), double(col));
    counts = accumarray(pixIdx, 1, [ny * nx, 1], @sum, 0);
    frameImg = reshape(counts, [ny, nx]);
end

function [startUse, stopUse] = pairPtuStartStop(starts, stops)
    starts = starts(:);
    stops = stops(:);
    maxPairs = min(numel(starts), numel(stops));
    startUse = zeros(maxPairs, 1, class(starts));
    stopUse = zeros(maxPairs, 1, class(stops));

    i = 1;
    j = 1;
    k = 0;
    while i <= numel(starts) && j <= numel(stops)
        if stops(j) <= starts(i)
            j = j + 1;
        else
            k = k + 1;
            startUse(k) = starts(i);
            stopUse(k) = stops(j);
            i = i + 1;
            j = j + 1;
        end
    end

    startUse = startUse(1:k);
    stopUse = stopUse(1:k);
end

function [vol, meta] = ptuOutputToVolume(ptuOut, fileName)
    meta = struct();
    meta.file = fileName;
    meta.inputSource = 'ptu';
    if isfield(ptuOut, 'head')
        meta.head = ptuOut.head;
    else
        meta.head = struct();
    end
    if isfield(ptuOut, 'dind')
        meta.activeChannels = ptuOut.dind;
    else
        meta.activeChannels = [];
    end

    if isfield(ptuOut, 'tag') && ~isempty(ptuOut.tag) && size(ptuOut.tag, 4) > 0
        tag = double(ptuOut.tag);
        dims = size(tag);
        while numel(dims) < 4
            dims(end+1) = 1; %#ok<AGROW>
        end

        nx = dims(1);
        ny = dims(2);
        nz = dims(4);
        volNxNyNz = zeros(nx, ny, nz);
        for z = 1:nz
            volNxNyNz(:, :, z) = sum(tag(:, :, :, z), 3);
        end
        vol = permute(volNxNyNz, [2 1 3]);
    elseif isfield(ptuOut, 'tags') && ~isempty(ptuOut.tags)
        tags = double(ptuOut.tags);
        volNxNy = sum(tags, 3);
        vol = permute(volNxNy, [2 1]);
        vol = reshape(vol, size(vol, 1), size(vol, 2), 1);
    else
        error('batchAnalyzeLuminosaPSFs:NoPtuImage', ...
            'No tag or tags image data found after reading %s.', fileName);
    end

    meta.pixelY = size(vol, 1);
    meta.pixelX = size(vol, 2);
    meta.numFrames = size(vol, 3);
    meta.xyPixelSizeUm = numericField(meta.head, 'ImgHdr_PixResol', NaN);
    meta.blockCount = NaN;
    meta.pixelatedValueBytes = NaN;
    meta.frameFiles = repmat({fileName}, 1, size(vol, 3));
    meta.frameFileFrames = 1:size(vol, 3);
end

function [img, meta] = readPqdatImage(fileName, readPixels)
    if nargin < 2 || isempty(readPixels)
        readPixels = true;
    end

    fid = fopen(fileName, 'r', 'ieee-le');
    if fid < 0
        error('batchAnalyzeLuminosaPSFs:OpenFailed', ...
            'Could not open %s.', fileName);
    end
    cleanup = onCleanup(@() fclose(fid));

    magic = fread(fid, 8, '*char')';
    version = fread(fid, 8, '*char')';
    if numel(magic) < 6 || ~strcmp(char(magic(1:6)), 'PQDATA')
        error('batchAnalyzeLuminosaPSFs:NotPqdata', ...
            '%s is not a PQDATA file.', fileName);
    end

    meta = struct();
    meta.file = fileName;
    meta.magic = trimNulls(char(magic));
    meta.version = trimNulls(char(version));
    meta.pixelX = NaN;
    meta.pixelY = NaN;
    meta.numFrames = 1;
    meta.blockCount = 1;
    meta.pixelatedValueBytes = NaN;
    meta.pixelatedValueOffset = NaN;

    imgValues = [];
    nGuard = 0;

    while ~feof(fid)
        nGuard = nGuard + 1;
        if nGuard > 10000
            error('batchAnalyzeLuminosaPSFs:TooManyTags', ...
                'Too many tags while parsing %s.', fileName);
        end

        tagStart = ftell(fid);
        identRaw = fread(fid, 32, '*char')';
        if numel(identRaw) < 32
            break;
        end
        ident = trimNulls(char(identRaw));
        tagIdx = fread(fid, 1, 'int32=>double');
        tagType = fread(fid, 1, 'uint32=>uint32');
        valueBytes = fread(fid, 8, 'uint8=>uint8')';
        if numel(valueBytes) < 8
            break;
        end

        valueInt64 = typecast(valueBytes, 'int64');
        valueDouble = typecast(valueBytes, 'double');
        nBytes = double(valueInt64);

        switch tagType
            case uint32(hex2dec('00000008')) % bool8
                value = logical(valueInt64 ~= 0);
            case uint32(hex2dec('10000008')) % int8
                value = double(valueInt64);
            case {uint32(hex2dec('20000008')), uint32(hex2dec('21000008'))}
                value = double(valueDouble);
            case uint32(hex2dec('4001FFFF')) % ANSI string
                raw = fread(fid, nBytes, 'uint8=>uint8')';
                value = trimNulls(char(raw));
            case uint32(hex2dec('4002FFFF')) % UTF-16 string
                raw = fread(fid, nBytes, 'uint8=>uint8')';
                value = decodeUtf16LE(raw);
            case uint32(hex2dec('FFFFFFFF')) % binary blob
                value = sprintf('<%d byte binary blob>', nBytes);
                fseek(fid, nBytes, 'cof');
            case uint32(hex2dec('2001FFFF')) % float8 array
                meta.pixelatedValueBytes = nBytes;
                meta.pixelatedValueOffset = ftell(fid);
                if strcmp(ident, 'LSDPixelatedValue')
                    if readPixels
                        newValues = fread(fid, nBytes/8, 'double=>double');
                        imgValues = [imgValues; newValues]; %#ok<AGROW>
                    else
                        fseek(fid, nBytes, 'cof');
                    end
                    value = sprintf('<%d byte float8 array>', nBytes);
                else
                    fseek(fid, nBytes, 'cof');
                    value = sprintf('<%d byte float8 array>', nBytes);
                end
            otherwise
                value = double(valueInt64);
        end

        meta = captureKnownPqdatTag(meta, ident, value, tagIdx, tagStart);

        if strcmp(ident, 'Header_End')
            break;
        end
    end

    if readPixels
        if isempty(imgValues)
            error('batchAnalyzeLuminosaPSFs:NoImageArray', ...
                'No LSDPixelatedValue array found in %s.', fileName);
        end
        img = reshapePqdatImage(imgValues, meta, fileName);
    else
        img = [];
    end
end

function meta = captureKnownPqdatTag(meta, ident, value, tagIdx, tagStart)
    switch ident
        case 'LSDPixelX'
            meta.pixelX = value;
        case 'LSDPixelY'
            meta.pixelY = value;
        case 'LSDNumFrames'
            meta.numFrames = value;
        case 'LSDBlockCount'
            meta.blockCount = value;
        case 'LSDName'
            meta.name = value;
        case 'LSDType'
            meta.type = value;
        case 'LSDScaleMax'
            meta.scaleMax = value;
        case 'ImgSamplePosX'
            meta.samplePosX = value;
        case 'ImgSamplePosY'
            meta.samplePosY = value;
    end

    validIdent = matlab.lang.makeValidName(ident);
    if ~isempty(validIdent)
        fieldName = validIdent;
        if tagIdx >= 0
            fieldName = sprintf('%s_%d', validIdent, tagIdx);
        end
        if ~isfield(meta, 'tags')
            meta.tags = struct();
            meta.tagOffsets = struct();
        end
        meta.tags.(fieldName) = value;
        meta.tagOffsets.(fieldName) = tagStart;
    end
end

function img = reshapePqdatImage(imgValues, meta, fileName)
    nx = double(meta.pixelX);
    ny = double(meta.pixelY);
    if ~isfinite(nx) || ~isfinite(ny) || nx <= 0 || ny <= 0
        error('batchAnalyzeLuminosaPSFs:MissingDimensions', ...
            'Missing LSDPixelX/LSDPixelY in %s.', fileName);
    end

    nPerFrame = nx * ny;
    nFrameFromData = numel(imgValues) / nPerFrame;
    if abs(nFrameFromData - round(nFrameFromData)) > eps(max(1, nFrameFromData))
        error('batchAnalyzeLuminosaPSFs:BadImageLength', ...
            'Image array length in %s does not match %dx%d pixels.', ...
            fileName, ny, nx);
    end

    nFrame = max(1, round(nFrameFromData));
    img = reshape(imgValues, [nx, ny, nFrame]);
    img = permute(img, [2 1 3]); % file order is x-fast, MATLAB image order is y,x
end

function txt = trimNulls(txt)
    txt = txt(txt ~= char(0));
    txt = strtrim(txt);
end

function txt = decodeUtf16LE(raw)
    try
        txt = native2unicode(raw, 'UTF-16LE');
    catch
        raw = raw(:);
        raw = raw(1:2*floor(numel(raw)/2));
        txt = char(typecast(uint8(raw), 'uint16')).';
    end
    txt = trimNulls(txt);
end

function [alignedVol, alignInfo] = alignPsfVolume(rawVol, opts)
    [ny, nx, nz] = size(rawVol);
    procVol = zeros(ny, nx, nz);
    centers = zeros(nz, 2);
    centerScores = zeros(nz, 1);
    backgrounds = zeros(nz, 1);

    for z = 1:nz
        [procVol(:, :, z), backgrounds(z)] = preprocessPsfFrame(rawVol(:, :, z), opts);
        [centers(z, :), centerScores(z)] = estimatePsfCenter(procVol(:, :, z), opts);
    end

    alignmentCenters = regularizeCenterTrajectory(centers, centerScores, [ny nx], opts);
    target = resolveTargetCenter(alignmentCenters, centerScores, [ny nx], opts);
    desiredCenters = desiredCenterTrajectory(alignmentCenters, centerScores, target, opts);
    initialShiftXY = desiredCenters - centers;

    preAligned = zeros(size(rawVol));
    preAlignedProc = zeros(size(procVol));
    for z = 1:nz
        preAligned(:, :, z) = fourierShift2D(rawVol(:, :, z), ...
            initialShiftXY(z, 1), initialShiftXY(z, 2));
        preAlignedProc(:, :, z) = fourierShift2D(procVol(:, :, z), ...
            initialShiftXY(z, 1), initialShiftXY(z, 2));
    end

    templateShiftXY = zeros(nz, 2);
    templatePeak = zeros(nz, 1);

    if opts.refineWithEmpiricalTemplate && nz > 1
        template = projection(preAlignedProc, 'max');
        template = normalizeForRegistration(template);

        for z = 1:nz
            moving = normalizeForRegistration(preAlignedProc(:, :, z));
            [dxy, peakValue] = phaseCorrelationShift(moving, template);
            if norm(dxy) <= opts.maxTemplateShiftPx && peakValue >= opts.minTemplatePeak
                templateShiftXY(z, :) = dxy;
                templatePeak(z) = peakValue;
            end
        end
    end

    totalShiftXY = initialShiftXY + templateShiftXY;
    alignedVol = zeros(size(rawVol));
    for z = 1:nz
        alignedVol(:, :, z) = fourierShift2D(rawVol(:, :, z), ...
            totalShiftXY(z, 1), totalShiftXY(z, 2));
    end
    alignedVol(alignedVol < 0) = 0;

    secondPassCenters = nan(nz, 2);
    secondPassScores = zeros(nz, 1);
    secondPassTargetCenters = nan(nz, 2);
    secondPassShiftXY = zeros(nz, 2);
    if opts.secondPassAlignment && nz > 3
        [alignedVol, secondPassCenters, secondPassScores, ...
            secondPassTargetCenters, secondPassShiftXY] = ...
            refineAlignedVolumeSecondPass(alignedVol, target, opts);
        totalShiftXY = totalShiftXY + secondPassShiftXY;
    end

    alignInfo = struct();
    alignInfo.background = backgrounds;
    alignInfo.estimatedCenterXY = centers;
    alignInfo.alignmentCenterXY = alignmentCenters;
    alignInfo.desiredCenterXY = desiredCenters;
    alignInfo.centerScore = centerScores;
    alignInfo.targetCenterXY = target;
    alignInfo.initialShiftXY = initialShiftXY;
    alignInfo.templateShiftXY = templateShiftXY;
    alignInfo.templatePeak = templatePeak;
    alignInfo.secondPassCenterXY = secondPassCenters;
    alignInfo.secondPassCenterScore = secondPassScores;
    alignInfo.secondPassTargetCenterXY = secondPassTargetCenters;
    alignInfo.secondPassShiftXY = secondPassShiftXY;
    alignInfo.totalShiftXY = totalShiftXY;
end

function [refinedVol, centers, scores, targetCenters, residualShiftXY] = ...
        refineAlignedVolumeSecondPass(alignedVol, target, opts)
    [ny, nx, nz] = size(alignedVol);
    centers = zeros(nz, 2);
    scores = zeros(nz, 1);

    for z = 1:nz
        procFrame = preprocessPsfFrame(alignedVol(:, :, z), opts);
        [centers(z, :), scores(z)] = estimatePsfCenter(procFrame, opts);
    end

    secondOpts = opts;
    secondOpts.centerSmoothingWindow = opts.secondPassSmoothingWindow;
    smoothCenters = regularizeCenterTrajectory(centers, scores, [ny nx], secondOpts);

    if opts.preserveSmoothedCenterTrajectory
        targetCenters = smoothCenters;
    else
        targetCenters = repmat(target, nz, 1);
    end

    residualShiftXY = targetCenters - centers;
    residualShiftXY = limitShiftMagnitude(residualShiftXY, opts.maxSecondPassShiftPx);

    refinedVol = zeros(size(alignedVol));
    for z = 1:nz
        refinedVol(:, :, z) = fourierShift2D(alignedVol(:, :, z), ...
            residualShiftXY(z, 1), residualShiftXY(z, 2));
    end
    refinedVol(refinedVol < 0) = 0;
end

function shiftXY = limitShiftMagnitude(shiftXY, maxShiftPx)
    if ~isfinite(maxShiftPx) || maxShiftPx <= 0
        shiftXY(:) = 0;
        return;
    end

    mag = sqrt(sum(shiftXY.^2, 2));
    tooLarge = mag > maxShiftPx;
    if any(tooLarge)
        scale = maxShiftPx ./ mag(tooLarge);
        shiftXY(tooLarge, :) = shiftXY(tooLarge, :) .* scale;
    end
end

function centersOut = regularizeCenterTrajectory(centers, scores, imageSize, opts)
    centersOut = centers;
    if ~opts.smoothCenterTrajectory || size(centers, 1) < 3
        return;
    end

    ny = imageSize(1);
    nx = imageSize(2);
    valid = all(isfinite(centers), 2) & scores >= opts.minCenterScoreForTrajectory & ...
        centers(:, 1) >= 1 & centers(:, 1) <= nx & ...
        centers(:, 2) >= 1 & centers(:, 2) <= ny;

    window = oddSmoothingWindow(opts.centerSmoothingWindow, size(centers, 1));
    if window < 3
        return;
    end

    centersOut(:, 1) = smoothCenterCoordinate(centers(:, 1), valid, ...
        (nx+1)/2, window, opts.centerOutlierThresholdPx);
    centersOut(:, 2) = smoothCenterCoordinate(centers(:, 2), valid, ...
        (ny+1)/2, window, opts.centerOutlierThresholdPx);
end

function desiredCenters = desiredCenterTrajectory(alignmentCenters, scores, target, opts)
    desiredCenters = repmat(target, size(alignmentCenters, 1), 1);
    if ~opts.preserveSmoothedCenterTrajectory
        return;
    end

    valid = all(isfinite(alignmentCenters), 2) & scores >= opts.minCenterScoreForTrajectory;
    if any(valid)
        anchor = [median(alignmentCenters(valid, 1)), median(alignmentCenters(valid, 2))];
    else
        anchor = [median(alignmentCenters(:, 1)), median(alignmentCenters(:, 2))];
    end
    desiredCenters = alignmentCenters + (target - anchor);
end

function y = smoothCenterCoordinate(x, valid, fallbackValue, window, outlierThreshold)
    x = double(x(:));
    valid = valid(:) & isfinite(x);
    y = fillVectorByInterpolation(x, valid, fallbackValue);

    baseline = movingMedianVector(y, window);
    residual = y - baseline;
    sigma = robustSigma(residual(valid));
    if isfinite(outlierThreshold) && outlierThreshold > 0
        if isfinite(sigma) && sigma > 0
            outlierLimit = max(outlierThreshold, 3 * sigma);
        else
            outlierLimit = outlierThreshold;
        end
        keep = valid & abs(residual) <= outlierLimit;
        y = fillVectorByInterpolation(y, keep, fallbackValue);
    end

    y = movingMeanVector(y, window);
end

function window = oddSmoothingWindow(requestedWindow, n)
    window = max(1, round(requestedWindow));
    window = min(window, n);
    if mod(window, 2) == 0
        window = window - 1;
    end
    window = max(1, window);
end

function y = fillVectorByInterpolation(x, valid, fallbackValue)
    n = numel(x);
    idx = (1:n).';
    y = x;
    if nnz(valid) >= 2
        y(~valid) = interp1(idx(valid), x(valid), idx(~valid), 'linear', 'extrap');
    elseif nnz(valid) == 1
        y(:) = x(find(valid, 1, 'first'));
    else
        y(:) = fallbackValue;
    end
end

function sigma = robustSigma(x)
    x = x(isfinite(x));
    if isempty(x)
        sigma = NaN;
        return;
    end
    med = median(x);
    sigma = 1.4826 * median(abs(x - med));
end

function y = movingMedianVector(x, window)
    n = numel(x);
    halfWindow = floor(window / 2);
    y = zeros(size(x));
    for k = 1:n
        lo = max(1, k - halfWindow);
        hi = min(n, k + halfWindow);
        y(k) = median(x(lo:hi));
    end
end

function y = movingMeanVector(x, window)
    n = numel(x);
    halfWindow = floor(window / 2);
    y = zeros(size(x));
    for k = 1:n
        lo = max(1, k - halfWindow);
        hi = min(n, k + halfWindow);
        y(k) = mean(x(lo:hi));
    end
end

function [img, background] = preprocessPsfFrame(frame, opts)
    frame = double(frame);
    border = [frame(1, :).'; frame(end, :).'; frame(:, 1); frame(:, end)];
    finiteBorder = border(isfinite(border));
    if isempty(finiteBorder)
        background = 0;
    else
        background = max(median(finiteBorder), 0);
    end
    img = frame - background;
    img(~isfinite(img)) = 0;
    img(img < 0) = 0;

    if opts.smoothSigmaPx > 0
        img = gaussianSmooth2D(img, opts.smoothSigmaPx);
    end

    s = sum(img(:));
    if s > 0
        img = img / s;
    end
end

function img = gaussianSmooth2D(img, sigma)
    radius = max(1, ceil(3*sigma));
    x = -radius:radius;
    kernel = exp(-(x.^2) / (2*sigma^2));
    kernel = kernel / sum(kernel);
    img = conv2(conv2(img, kernel, 'same'), kernel.', 'same');
end

function [centerXY, score] = estimatePsfCenter(img, opts)
    [ny, nx] = size(img);
    fallback = [(nx+1)/2, (ny+1)/2];
    peak = max(img(:));

    if ~isfinite(peak) || peak <= 0
        centerXY = fallback;
        score = 0;
        return;
    end

    centers = zeros(0, 2);
    weights = opts.centerThresholdFractions(:).';

    [xGrid, yGrid] = meshgrid(1:nx, 1:ny);
    for k = 1:numel(weights)
        w = max(img - weights(k)*peak, 0);
        mass = sum(w(:));
        if mass <= 0
            continue;
        end
        centers(end+1, :) = [sum(xGrid(:).*w(:))/mass, ...
                             sum(yGrid(:).*w(:))/mass]; %#ok<AGROW>
    end

    if isempty(centers)
        mass = sum(img(:));
        if mass > 0
            centers = [sum(xGrid(:).*img(:))/mass, sum(yGrid(:).*img(:))/mass];
        else
            centers = fallback;
        end
    end

    centerXY = [median(centers(:, 1)), median(centers(:, 2))];
    score = peak / max(eps, sum(img(:)));
end

function target = resolveTargetCenter(centers, scores, imageSize, opts)
    ny = imageSize(1);
    nx = imageSize(2);
    imageCenter = [(nx+1)/2, (ny+1)/2];

    if isnumeric(opts.targetCenter) && numel(opts.targetCenter) == 2
        target = double(opts.targetCenter(:)).';
        return;
    end

    if ischar(opts.targetCenter) || isstring(opts.targetCenter)
        mode = char(opts.targetCenter);
    else
        mode = 'median';
    end

    switch lower(mode)
        case {'imagecenter', 'image_center', 'center'}
            target = imageCenter;
        otherwise
            valid = all(isfinite(centers), 2) & scores > 0;
            if any(valid)
                target = [median(centers(valid, 1)), median(centers(valid, 2))];
            else
                target = imageCenter;
            end
    end
end

function img = normalizeForRegistration(img)
    img = double(img);
    img(~isfinite(img)) = 0;
    img = img - median(img(:));
    s = norm(img(:));
    if s > 0
        img = img / s;
    end
    [ny, nx] = size(img);
    img = img .* hannWindow2D(ny, nx);
end

function w = hannWindow2D(ny, nx)
    if ny > 1
        wy = 0.5 - 0.5*cos(2*pi*(0:ny-1)'/(ny-1));
    else
        wy = 1;
    end
    if nx > 1
        wx = 0.5 - 0.5*cos(2*pi*(0:nx-1)/(nx-1));
    else
        wx = 1;
    end
    w = wy * wx;
end

function [shiftXY, peakValue] = phaseCorrelationShift(moving, fixedRef)
    [ny, nx] = size(moving);
    A = fft2(moving);
    B = fft2(fixedRef);
    cps = B .* conj(A);
    denom = abs(cps);
    cps(denom > 0) = cps(denom > 0) ./ denom(denom > 0);
    cps(denom == 0) = 0;

    corr = real(ifft2(cps));
    [peakValue, idx] = max(corr(:));
    [py, px] = ind2sub(size(corr), idx);

    dy = py - 1;
    dx = px - 1;
    if dy > floor(ny/2)
        dy = dy - ny;
    end
    if dx > floor(nx/2)
        dx = dx - nx;
    end

    shiftXY = [dx, dy];
end

function outputFiles = writeGroupOutputs(rawVol, alignedVol, scans, group, ...
        frameMeta, alignInfo, groupDir, opts)
    outputFiles = struct();

    if opts.flatOutputLayout
        outputFiles = writeFlatScanOutputs(rawVol, alignedVol, scans, group, ...
            frameMeta, alignInfo, opts);
        return;
    end

    outputFiles.projectionFigure = fullfile(groupDir, 'aligned_psf_projections.png');
    outputFiles.xyMIP = fullfile(groupDir, 'aligned_xy_mip.png');
    outputFiles.xzMIP = fullfile(groupDir, 'aligned_xz_mip.png');
    outputFiles.yzMIP = fullfile(groupDir, 'aligned_yz_mip.png');
    outputFiles.alignmentCsv = fullfile(groupDir, 'frame_alignment.csv');
    outputFiles.volumeMat = fullfile(groupDir, 'psf_volume_aligned.mat');

    writeProjectionFigure(rawVol, alignedVol, scans, group, frameMeta, ...
        alignInfo, outputFiles.projectionFigure, opts);
    writeProjectionPngs(alignedVol, outputFiles, scans, group, opts);
    writeAlignmentCsv(frameMeta, alignInfo, outputFiles.alignmentCsv);

    if opts.makeIsosurface && size(alignedVol, 3) >= 3
        outputFiles.isosurface = fullfile(groupDir, 'aligned_psf_isosurface.png');
        writeIsosurfaceFigure(alignedVol, scans, group, alignInfo, ...
            outputFiles.isosurface, opts);
    else
        outputFiles.isosurface = '';
    end

    rawVolume = rawVol; %#ok<NASGU>
    alignedVolume = alignedVol; %#ok<NASGU>
    save(outputFiles.volumeMat, 'rawVolume', 'alignedVolume', ...
        'scans', 'group', 'frameMeta', 'alignInfo', 'opts');
end

function outputFiles = writeFlatScanOutputs(rawVol, alignedVol, scans, group, ...
        frameMeta, alignInfo, opts)
    stem = scanOutputStem(group.label);
    outputFiles = struct();
    outputFiles.volumeMat = fullfile(opts.volumeOutputDir, ...
        sprintf('%s_volume_raw.mat', stem));
    outputFiles.xyPlot = fullfile(opts.projectionOutputDir, ...
        sprintf('%s_XY.png', stem));
    outputFiles.xzYzPlot = fullfile(opts.projectionOutputDir, ...
        sprintf('%s_XZ_YZ.png', stem));
    outputFiles.alignmentCsv = fullfile(opts.projectionOutputDir, ...
        sprintf('%s_frame_alignment.csv', stem));

    writeFlatProjectionFigures(alignedVol, scans, group, outputFiles, opts);
    writeAlignmentCsv(frameMeta, alignInfo, outputFiles.alignmentCsv);

    volume = alignedVol; %#ok<NASGU>
    rawVolume = rawVol; %#ok<NASGU>
    alignedVolume = alignedVol; %#ok<NASGU>
    save(outputFiles.volumeMat, 'volume', 'rawVolume', 'alignedVolume', ...
        'scans', 'group', 'frameMeta', 'alignInfo', 'opts');
end

function writeFlatProjectionFigures(vol, scans, group, outputFiles, opts)
    xy = projection(vol, opts.projectionMode);
    xz = xzProjection(vol, opts.projectionMode);
    yz = yzProjection(vol, opts.projectionMode);
    scale = coordinateScale(scans, size(vol), opts);
    titleLabel = strrep(group.label, '_', '\_');

    writeSingleProjectionFigure(xy, scale.x, scale.y, scale.xLabel, scale.yLabel, ...
        scale.xyPixelSizeUm, scale.xyPixelSizeUm, ...
        sprintf('Aligned XY projection: %s', titleLabel), outputFiles.xyPlot, ...
        opts.writeFigFiles);

    writeXzYzProjectionFigure(xz, yz, scale, ...
        sprintf('Aligned XZ / YZ projections: %s', titleLabel), ...
        outputFiles.xzYzPlot, opts);
end

function writeXzYzProjectionFigure(xz, yz, scale, titleText, outFile, opts)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [120 120 980 460]);
    subplot(1, 2, 1);
    showImageScaled(xz, scale.x, scale.z, scale.xLabel, scale.zLabel, ...
        scale.xyPixelSizeUm, scale.zStepUm);
    title('Aligned XZ projection');

    subplot(1, 2, 2);
    showImageScaled(yz, scale.y, scale.z, scale.yLabel, scale.zLabel, ...
        scale.xyPixelSizeUm, scale.zStepUm);
    title('Aligned YZ projection');

    addSuperTitle(fig, titleText);
    colormap(fig, hot);
    saveFigure(fig, outFile, 180);
    if opts.writeFigFiles
        savefig(fig, replaceFileExtension(outFile, '.fig'));
    end
    close(fig);
end

function stem = scanOutputStem(label)
    stem = regexprep(label, '[^A-Za-z0-9]+', '_');
    stem = regexprep(stem, '_+', '_');
    stem = regexprep(stem, '^_+|_+$', '');
    if isempty(stem)
        stem = 'scan';
    end
    if ~isletter(stem(1))
        stem = ['x', stem];
    end
end

function outFile = replaceFileExtension(fileName, newExt)
    [folder, name] = fileparts(fileName);
    outFile = fullfile(folder, [name, newExt]);
end

function writeProjectionFigure(rawVol, alignedVol, scans, group, frameMeta, ...
        alignInfo, outFile, opts)
    rawXY = projection(rawVol, opts.projectionMode);
    alignedXY = projection(alignedVol, opts.projectionMode);
    alignedXZ = xzProjection(alignedVol, opts.projectionMode);
    alignedYZ = yzProjection(alignedVol, opts.projectionMode);
    scale = coordinateScale(scans, size(alignedVol), opts);

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 850]);
    titleText = sprintf('%s  (%d frames, %dx%d, dz = %.3g um)', ...
        strrep(group.label, '_', '\_'), size(alignedVol, 3), ...
        size(alignedVol, 1), size(alignedVol, 2), opts.zStepUm);
    addSuperTitle(fig, titleText);

    subplot(2, 3, 1);
    showImageScaled(rawXY, scale.x, scale.y, scale.xLabel, scale.yLabel, ...
        scale.xyPixelSizeUm, scale.xyPixelSizeUm);
    title('Raw XY projection');

    subplot(2, 3, 2);
    showImageScaled(alignedXY, scale.x, scale.y, scale.xLabel, scale.yLabel, ...
        scale.xyPixelSizeUm, scale.xyPixelSizeUm);
    title('Aligned XY projection');

    subplot(2, 3, 3);
    showImageScaled(alignedXZ, scale.x, scale.z, scale.xLabel, scale.zLabel, ...
        scale.xyPixelSizeUm, scale.zStepUm);
    title('Aligned XZ projection');

    subplot(2, 3, 4);
    showImageScaled(alignedYZ, scale.y, scale.z, scale.yLabel, scale.zLabel, ...
        scale.xyPixelSizeUm, scale.zStepUm);
    title('Aligned YZ projection');

    subplot(2, 3, 5);
    z = scale.z;
    plot(z, alignInfo.estimatedCenterXY(:, 1), 'o-', 'LineWidth', 1.2);
    hold on;
    plot(z, alignInfo.estimatedCenterXY(:, 2), 's-', 'LineWidth', 1.2);
    yline(alignInfo.targetCenterXY(1), '--', 'target x');
    yline(alignInfo.targetCenterXY(2), ':', 'target y');
    hold off;
    grid on;
    xlabel(scale.zLabel);
    ylabel('center position (px)');
    legend({'x', 'y'}, 'Location', 'best', 'Box', 'off');
    title('Estimated frame centers');

    subplot(2, 3, 6);
    plot(z, alignInfo.totalShiftXY(:, 1), 'o-', 'LineWidth', 1.2);
    hold on;
    plot(z, alignInfo.totalShiftXY(:, 2), 's-', 'LineWidth', 1.2);
    signalTrace = [frameMeta.totalSignal].';
    signalScale = max(signalTrace);
    if signalScale > 0
        signalTrace = signalTrace ./ signalScale;
    end
    plot(z, signalTrace, 'k.-', 'LineWidth', 1.0);
    hold off;
    grid on;
    xlabel(scale.zLabel);
    ylabel('shift (px) / normalized signal');
    legend({'x shift', 'y shift', 'signal'}, 'Location', 'best', 'Box', 'off');
    title(sprintf('Frame shifts, first: %s', scans(1).timestamp));

    colormap(fig, hot);
    saveFigure(fig, outFile, 180);
    close(fig);
end

function addSuperTitle(fig, titleText)
    if exist('sgtitle', 'file') == 2
        sgtitle(titleText, 'FontSize', 13, 'FontWeight', 'bold');
    else
        annotation(fig, 'textbox', [0.02 0.955 0.96 0.035], ...
            'String', titleText, 'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, ...
            'FontWeight', 'bold', 'Interpreter', 'tex');
    end
end

function showImageScaled(img, xCoord, yCoord, xLabel, yLabel, xStep, yStep)
    if nargin < 6, xStep = []; end
    if nargin < 7, yStep = []; end

    imagesc(axisLimitsFromCenters(xCoord, xStep), ...
        axisLimitsFromCenters(yCoord, yStep), img);
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel(xLabel);
    ylabel(yLabel);
    colorbar;
end

function writeProjectionPngs(vol, outputFiles, scans, group, opts)
    xy = projection(vol, opts.projectionMode);
    xz = xzProjection(vol, opts.projectionMode);
    yz = yzProjection(vol, opts.projectionMode);
    scale = coordinateScale(scans, size(vol), opts);

    imwrite(normalizeImageForWrite(xy), outputFiles.xyMIP);
    writeSingleProjectionFigure(xz, scale.x, scale.z, scale.xLabel, scale.zLabel, ...
        scale.xyPixelSizeUm, scale.zStepUm, ...
        sprintf('Aligned XZ projection: %s', strrep(group.label, '_', '\_')), ...
        outputFiles.xzMIP);
    writeSingleProjectionFigure(yz, scale.y, scale.z, scale.yLabel, scale.zLabel, ...
        scale.xyPixelSizeUm, scale.zStepUm, ...
        sprintf('Aligned YZ projection: %s', strrep(group.label, '_', '\_')), ...
        outputFiles.yzMIP);
end

function writeSingleProjectionFigure(img, xCoord, yCoord, xLabel, yLabel, ...
        xStep, yStep, titleText, outFile, writeFigFile)
    if nargin < 9 || isempty(writeFigFile)
        writeFigFile = false;
    end
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [120 120 720 560]);
    showImageScaled(img, xCoord, yCoord, xLabel, yLabel, xStep, yStep);
    title(titleText, 'Interpreter', 'tex');
    colormap(fig, hot);
    saveFigure(fig, outFile, 180);
    if writeFigFile
        savefig(fig, replaceFileExtension(outFile, '.fig'));
    end
    close(fig);
end

function writeIsosurfaceFigure(vol, scans, group, alignInfo, outFile, opts)
    v = double(vol);
    v = v - min(v(:));
    if max(v(:)) > 0
        v = v / max(v(:));
    end

    scale = coordinateScale(scans, size(v), opts);
    [X, Y, Z] = meshgrid(scale.x, scale.y, scale.z);

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [120 120 850 720]);
    p = patch(isosurface(X, Y, Z, v, opts.isosurfaceLevel));
    set(p, 'FaceColor', [0.95 0.45 0.10], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.55);
    try
        isonormals(X, Y, Z, v, p);
    catch
    end
    hold on;
    targetX = coordinateAtPixel(alignInfo.targetCenterXY(1), scale.x);
    targetY = coordinateAtPixel(alignInfo.targetCenterXY(2), scale.y);
    plot3(targetX, targetY, 0, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
    hold off;
    axis image;
    axis vis3d;
    grid on;
    view(3);
    xlabel(scale.xLabel);
    ylabel(scale.yLabel);
    zlabel(scale.zLabel);
    title(sprintf('Aligned PSF isosurface: %s', strrep(group.label, '_', '\_')), ...
        'Interpreter', 'tex');
    light('Style', 'infinite', 'Position', [0 1 1]);
    light('Style', 'infinite', 'Position', [1 -1 0.5]);
    lighting gouraud;
    saveFigure(fig, outFile, 180);
    close(fig);
end

function scale = coordinateScale(scans, volSize, opts)
    if numel(volSize) < 3
        volSize(3) = 1;
    end

    ny = volSize(1);
    nx = volSize(2);
    nz = volSize(3);

    xyPixelSizeUm = firstFinite([scans.xyPixelSizeUm]);
    if isnumeric(opts.xyPixelSizeUm) && isscalar(opts.xyPixelSizeUm) && ...
            isfinite(opts.xyPixelSizeUm) && opts.xyPixelSizeUm > 0
        xyPixelSizeUm = double(opts.xyPixelSizeUm);
    end

    if isfinite(xyPixelSizeUm) && xyPixelSizeUm > 0
        scale.x = ((1:nx) - (nx+1)/2) * xyPixelSizeUm;
        scale.y = ((1:ny) - (ny+1)/2) * xyPixelSizeUm;
        scale.xLabel = 'x (um)';
        scale.yLabel = 'y (um)';
        scale.xyPixelSizeUm = xyPixelSizeUm;
    else
        scale.x = 1:nx;
        scale.y = 1:ny;
        scale.xLabel = 'x pixel';
        scale.yLabel = 'y pixel';
        scale.xyPixelSizeUm = NaN;
    end

    scale.z = zCoordinatesUm(nz, opts);
    if isfinite(opts.zStepUm) && opts.zStepUm > 0
        scale.zLabel = 'z (um)';
    else
        scale.zLabel = 'frame index';
    end
    scale.zStepUm = opts.zStepUm;
end

function zUm = zCoordinatesUm(nz, opts)
    if isfield(opts, 'zStepUm') && isfinite(opts.zStepUm) && opts.zStepUm > 0
        zUm = ((1:nz) - (nz+1)/2) * opts.zStepUm;
    else
        zUm = 1:nz;
    end
end

function v = firstFinite(values)
    values = values(isfinite(values));
    if isempty(values)
        v = NaN;
    else
        v = values(1);
    end
end

function lim = axisLimitsFromCenters(coord, defaultStep)
    if nargin < 2 || isempty(defaultStep) || ~isfinite(defaultStep) || defaultStep <= 0
        defaultStep = 1;
    end

    coord = double(coord(:)).';
    if numel(coord) > 1
        step = median(diff(coord));
        if ~isfinite(step) || step == 0
            step = defaultStep;
        end
    else
        step = defaultStep;
    end
    lim = [coord(1) - step/2, coord(end) + step/2];
end

function value = coordinateAtPixel(pixelPosition, coord)
    coord = double(coord(:)).';
    if numel(coord) <= 1
        value = coord(1);
    else
        value = interp1(1:numel(coord), coord, pixelPosition, 'linear', 'extrap');
    end
end

function out = projection(vol, mode)
    out = projectionOverDim(vol, 3, mode);
end

function img = xzProjection(vol, mode)
    p = projectionOverDim(vol, 1, mode);       % 1 x nx x nz
    img = reshape(p, size(vol, 2), size(vol, 3)).';
end

function img = yzProjection(vol, mode)
    p = projectionOverDim(vol, 2, mode);       % ny x 1 x nz
    img = reshape(p, size(vol, 1), size(vol, 3)).';
end

function out = projectionOverDim(vol, dim, mode)
    switch lower(mode)
        case 'sum'
            out = sum(vol, dim);
        otherwise
            out = max(vol, [], dim);
    end
end

function img = normalizeImageForWrite(img)
    img = double(img);
    img(~isfinite(img)) = 0;
    img = img - min(img(:));
    m = max(img(:));
    if m > 0
        img = img / m;
    end
    img = uint16(round(65535 * img));
end

function saveFigure(fig, outFile, resolution)
    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, outFile, 'Resolution', resolution);
    else
        print(fig, outFile, '-dpng', sprintf('-r%d', resolution));
    end
end

function writeScanInventory(scans, outFile)
    rows = cell(numel(scans), 17);
    for k = 1:numel(scans)
        rows(k, :) = { ...
            scans(k).name, scans(k).inputSource, scans(k).fileMode, ...
            scans(k).numSourceFiles, scans(k).power, scans(k).collarText, ...
            scans(k).collar, scans(k).lens, scans(k).timestamp, ...
            scans(k).pixelY, scans(k).pixelX, scans(k).xyPixelSizeUm, ...
            scans(k).numFrames, scans(k).blockCount, ...
            scans(k).pixelatedValueBytes, scans(k).groupKey, scans(k).file};
    end
    writeCsv(outFile, {'name','input_source','file_mode','num_source_files', ...
        'power','collar_text','collar','lens','timestamp','pixel_y','pixel_x', ...
        'xy_pixel_size_um','num_frames','block_count','pixelated_value_bytes', ...
        'group_key','file'}, rows);
end

function writeAlignmentCsv(frameMeta, alignInfo, outFile)
    n = numel(frameMeta);
    rows = cell(n, 21);
    for k = 1:n
        if isfield(alignInfo, 'alignmentCenterXY')
            alignmentCenter = alignInfo.alignmentCenterXY(k, :);
        else
            alignmentCenter = alignInfo.estimatedCenterXY(k, :);
        end
        if isfield(alignInfo, 'secondPassCenterXY') && all(isfinite(alignInfo.secondPassCenterXY(k, :)))
            secondPassCenter = alignInfo.secondPassCenterXY(k, :);
        else
            secondPassCenter = [NaN NaN];
        end
        if isfield(alignInfo, 'secondPassShiftXY')
            secondPassShift = alignInfo.secondPassShiftXY(k, :);
        else
            secondPassShift = [0 0];
        end
        rows(k, :) = { ...
            k, frameMeta(k).sourceName, frameMeta(k).sourceFile, frameMeta(k).sourceFrame, ...
            frameMeta(k).zUm, frameMeta(k).totalSignal, alignInfo.background(k), ...
            alignInfo.estimatedCenterXY(k, 1), alignInfo.estimatedCenterXY(k, 2), ...
            alignmentCenter(1), alignmentCenter(2), ...
            secondPassCenter(1), secondPassCenter(2), ...
            alignInfo.targetCenterXY(1), alignInfo.targetCenterXY(2), ...
            alignInfo.initialShiftXY(k, 1), alignInfo.initialShiftXY(k, 2), ...
            secondPassShift(1), secondPassShift(2), ...
            alignInfo.totalShiftXY(k, 1), alignInfo.totalShiftXY(k, 2)};
    end
    writeCsv(outFile, {'frame_index','source_name','source_file','source_frame', ...
        'z_um','total_signal','background','estimated_center_x','estimated_center_y', ...
        'alignment_center_x','alignment_center_y','second_pass_center_x','second_pass_center_y', ...
        'target_center_x','target_center_y','initial_shift_x','initial_shift_y', ...
        'second_pass_shift_x','second_pass_shift_y','total_shift_x','total_shift_y'}, rows);
end

function writeBatchSummary(groupResults, outFile)
    rows = cell(numel(groupResults), 12);
    for k = 1:numel(groupResults)
        projectionFigure = outputField(groupResults(k).outputFiles, 'projectionFigure', '');
        xyPlot = outputField(groupResults(k).outputFiles, 'xyPlot', '');
        xzYzPlot = outputField(groupResults(k).outputFiles, 'xzYzPlot', '');
        volumeMat = outputField(groupResults(k).outputFiles, 'volumeMat', '');
        rows(k, :) = { ...
            groupResults(k).label, groupResults(k).numFrames, ...
            groupResults(k).imageSize(1), groupResults(k).imageSize(2), ...
            groupResults(k).xyPixelSizeUm, groupResults(k).zStepUm, ...
            groupResults(k).zSpanUm, groupResults(k).groupDir, ...
            projectionFigure, xyPlot, xzYzPlot, volumeMat};
    end
    writeCsv(outFile, {'group_label','num_frames','pixel_y','pixel_x', ...
        'xy_pixel_size_um','z_step_um','z_span_um','group_dir', ...
        'projection_figure','xy_plot','xz_yz_plot','volume_mat'}, rows);
end

function value = outputField(s, name, defaultValue)
    value = defaultValue;
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    end
end

function writeCsv(outFile, headers, rows)
    fid = fopen(outFile, 'w');
    if fid < 0
        error('batchAnalyzeLuminosaPSFs:CsvOpenFailed', ...
            'Could not open CSV for writing: %s', outFile);
    end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '%s\n', joinCsvRow(headers));
    for r = 1:size(rows, 1)
        fprintf(fid, '%s\n', joinCsvRow(rows(r, :)));
    end
end

function line = joinCsvRow(values)
    parts = cell(1, numel(values));
    for k = 1:numel(values)
        parts{k} = csvValue(values{k});
    end
    line = strjoin(parts, ',');
end

function txt = csvValue(v)
    if isnumeric(v)
        if isempty(v) || (isscalar(v) && isnan(v))
            txt = '';
        elseif isscalar(v)
            txt = sprintf('%.15g', v);
        else
            txt = sprintf('%.15g ', v);
            txt = strtrim(txt);
        end
    elseif islogical(v)
        txt = sprintf('%d', v);
    else
        txt = char(v);
    end

    txt = strrep(txt, '"', '""');
    if any(txt == ',') || any(txt == '"') || any(txt == sprintf('\n')) || any(txt == sprintf('\r'))
        txt = ['"', txt, '"'];
    end
end

function ensureDir(pathName)
    if exist(pathName, 'dir') ~= 7
        mkdir(pathName);
    end
end

function clean = sanitizeFileName(txt)
    clean = regexprep(txt, '[<>:"/\\|?*]', '_');
    clean = regexprep(clean, '\s+', '_');
    clean = regexprep(clean, '_+', '_');
    clean = regexprep(clean, '^_+|_+$', '');
end

function name = leafFolderName(pathName)
    pathName = char(pathName);
    while ~isempty(pathName) && (pathName(end) == filesep || pathName(end) == '/' || pathName(end) == '\')
        pathName(end) = [];
    end
    [~, name] = fileparts(pathName);
    if isempty(name)
        name = 'scan';
    end
end

function out = emptyGroupResult()
    out = struct('label', '', 'key', '', 'groupDir', '', ...
        'numFrames', 0, 'imageSize', [0 0], 'xyPixelSizeUm', NaN, ...
        'zStepUm', NaN, 'zSpanUm', NaN, 'outputFiles', struct(), ...
        'alignment', struct());
end
