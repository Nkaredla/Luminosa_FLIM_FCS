function data = loadFullStackISMData(inputValue, varargin)
%LOADFULLSTACKISMDATA Load raw detector-resolved ISM planes for Poisson fits.
%
%   data = loadFullStackISMData(alignmentCsv)
%
%   The loader deliberately does not interpolate, align, background-subtract,
%   or flat-field the photon counts. Saved scan-alignment shifts, measured
%   backgrounds, and independent detector gains are returned separately so
%   they can be applied to the forward model without destroying Poisson count
%   statistics.

    opts = parseOptions(varargin{:});
    addRequiredPaths(opts);

    if isnumeric(inputValue)
        data = numericInputData(inputValue, opts);
        return;
    end

    inputPath = char(inputValue);
    if exist(inputPath, 'dir') == 7
        inputPath = findAlignmentCsv(inputPath);
    end
    if exist(inputPath, 'file') ~= 2
        error('loadFullStackISMData:MissingInput', ...
            'Input was not found: %s', inputPath);
    end

    [~,~,ext] = fileparts(inputPath);
    switch lower(ext)
        case '.csv'
            data = loadFromAlignmentCsv(inputPath, opts);
        case '.mat'
            data = loadFromMat(inputPath, opts);
        otherwise
            error('loadFullStackISMData:BadInput', ...
                'Use an alignment CSV, detector-stack MAT file, folder, or numeric stack.');
    end
end

function opts = parseOptions(varargin)
    p = inputParser;
    p.FunctionName = 'loadFullStackISMData';
    addParameter(p, 'stageZUm', []);
    addParameter(p, 'xyPixelSizeUm', []);
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
    addParameter(p, 'verbose', true);
    parse(p, varargin{:});
    opts = p.Results;
    opts.channelIDs = double(opts.channelIDs(:)).';
    opts.channelOrder = double(opts.channelOrder(:)).';
    opts.ptuReaderFolder = char(opts.ptuReaderFolder);
    opts.cacheFile = char(opts.cacheFile);
    opts.flatFieldVariable = char(opts.flatFieldVariable);
    opts.darkFile = char(opts.darkFile);
    opts.darkMeasurementMode = lower(char(opts.darkMeasurementMode));
    opts.backgroundMode = lower(char(opts.backgroundMode));
    if ~ismember(opts.darkMeasurementMode,{'tttr','image','auto'})
        error('loadFullStackISMData:BadDarkMeasurementMode', ...
            'darkMeasurementMode must be tttr, image, or auto.');
    end
    if ~ismember(opts.backgroundMode,{'auto','dark','boundary'})
        error('loadFullStackISMData:BadBackgroundMode', ...
            'backgroundMode must be auto, dark, or boundary.');
    end
    if ~isscalar(opts.darkScale) || ~isfinite(opts.darkScale) || opts.darkScale < 0
        error('loadFullStackISMData:BadDarkScale', ...
            'darkScale must be a finite nonnegative scalar.');
    end
end

function folder = defaultPtuReaderFolder()
    thisDir = fileparts(mfilename('fullpath'));
    folder = fullfile(fileparts(fileparts(thisDir)), 'Luminosa_FLIM_FCS');
end

function addRequiredPaths(opts)
    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    addpath(fileparts(thisDir));
    if ~isempty(opts.ptuReaderFolder) && exist(opts.ptuReaderFolder, 'dir') == 7
        addpath(opts.ptuReaderFolder);
    end
end

function csvFile = findAlignmentCsv(folder)
    files = dir(fullfile(folder, '*frame_alignment.csv'));
    if isempty(files)
        files = dir(fullfile(folder, '**', '*frame_alignment.csv'));
    end
    if isempty(files)
        error('loadFullStackISMData:MissingAlignmentCsv', ...
            'No frame-alignment CSV was found below %s.', folder);
    end
    [~, idx] = max([files.datenum]);
    csvFile = fullfile(files(idx).folder, files(idx).name);
end

function data = loadFromAlignmentCsv(csvFile, opts)
    T = readtable(csvFile);
    required = {'source_file','z_um'};
    if ~all(ismember(required, T.Properties.VariableNames))
        error('loadFullStackISMData:BadAlignmentCsv', ...
            'Alignment CSV must contain source_file and z_um.');
    end

    cacheFile = opts.cacheFile;
    if isempty(cacheFile)
        [folder, stem] = fileparts(csvFile);
        cacheFile = fullfile(folder, [stem '_raw_detector_cache.mat']);
    end
    if opts.reuseCache && exist(cacheFile, 'file') == 2
        S = load(cacheFile);
        if isfield(S, 'data') && isstruct(S.data) && ...
                isfield(S.data, 'rawCounts')
            data = S.data;
            data = attachBackgroundModel(data, opts);
            data = attachFlatField(data, opts);
            if opts.verbose
                fprintf('[loadFullStackISMData] Loaded cache: %s\n', cacheFile);
            end
            return;
        end
    end

    nPlane = height(T);
    if nPlane < 1
        error('loadFullStackISMData:EmptyAlignmentCsv', ...
            'Alignment CSV contains no planes: %s', csvFile);
    end

    fileName = tableText(T.source_file, 1);
    frameIndex = tableFrame(T, 1);
    [firstPlane, channelIDs, firstHead] = readPtuDetectorPlane( ...
        fileName, frameIndex, opts);
    rawCounts = zeros(size(firstPlane,1), size(firstPlane,2), ...
        size(firstPlane,3), nPlane, 'like', firstPlane);
    rawCounts(:,:,:,1) = firstPlane;
    if opts.verbose
        fprintf('[loadFullStackISMData] PTU plane 1/%d\n', nPlane);
    end

    for iz = 2:nPlane
        fileName = tableText(T.source_file, iz);
        frameIndex = tableFrame(T, iz);
        [plane, planeIDs, ~] = readPtuDetectorPlane( ...
            fileName, frameIndex, opts);
        if ~isequal(size(plane), size(firstPlane))
            error('loadFullStackISMData:PlaneSizeMismatch', ...
                'Plane %d has size [%s], expected [%s].', iz, ...
                num2str(size(plane)), num2str(size(firstPlane)));
        end
        if ~isequal(planeIDs(:), channelIDs(:))
            error('loadFullStackISMData:ChannelMismatch', ...
                'Plane %d has a different detector-channel order.', iz);
        end
        rawCounts(:,:,:,iz) = plane;
        if opts.verbose
            fprintf('[loadFullStackISMData] PTU plane %d/%d\n', iz, nPlane);
        end
    end

    data = assembleData(rawCounts, double(T.z_um(:)).', channelIDs, ...
        firstHead, opts);
    data.alignmentCsv = csvFile;
    data.sourceFiles = tableTextColumn(T.source_file);
    data.sourceFrames = tableFrameColumn(T);
    data.alignmentShiftXPx = tableNumericOrZeros(T, 'total_shift_x');
    data.alignmentShiftYPx = tableNumericOrZeros(T, 'total_shift_y');
    data.targetCenterXY = resolveTargetCenter(T, size(rawCounts));
    data.modelShiftToRawXPx = data.targetCenterXY(1) - ...
        (size(rawCounts,2)+1)/2 - data.alignmentShiftXPx;
    data.modelShiftToRawYPx = data.targetCenterXY(2) - ...
        (size(rawCounts,1)+1)/2 - data.alignmentShiftYPx;
    data.cacheFile = cacheFile;
    cacheData = data;
    cacheData.flatFieldGain = ones(size(data.flatFieldGain));
    cacheData.flatFieldCalibrated = false;
    cacheData.flatFieldSource = 'not stored in raw cache';
    cacheData = useBoundaryBackground(cacheData);
    data = cacheData;
    save(cacheFile, 'data', '-v7.3');
    data = attachBackgroundModel(data, opts);
    data = attachFlatField(data, opts);
end

function data = loadFromMat(matFile, opts)
    S = load(matFile);
    raw = findDetectorStack(S);
    stageZ = opts.stageZUm;
    if isempty(stageZ)
        stageZ = stageZFromStruct(S);
    end
    data = numericInputData(raw, setfieldLocal(opts, 'stageZUm', stageZ));
    data.sourceMat = matFile;
    if isfield(S, 'alignmentShiftXPx')
        data.alignmentShiftXPx = double(S.alignmentShiftXPx(:)).';
    end
    if isfield(S, 'alignmentShiftYPx')
        data.alignmentShiftYPx = double(S.alignmentShiftYPx(:)).';
    end
    data.modelShiftToRawXPx = -data.alignmentShiftXPx;
    data.modelShiftToRawYPx = -data.alignmentShiftYPx;
end

function data = numericInputData(raw, opts)
    raw = double(raw);
    nCh = numel(opts.channelIDs);
    if ndims(raw) ~= 4
        error('loadFullStackISMData:BadNumericStack', ...
            'Numeric input must be [y x detector z] or [y x z detector].');
    end
    if size(raw,3) == nCh
        raw4 = raw;
    elseif size(raw,4) == nCh
        raw4 = permute(raw, [1 2 4 3]);
    elseif size(raw,3) == 23
        raw4 = raw;
    elseif size(raw,4) == 23
        raw4 = permute(raw, [1 2 4 3]);
    else
        error('loadFullStackISMData:BadNumericStack', ...
            'Could not identify the detector dimension.');
    end

    stageZ = double(opts.stageZUm(:)).';
    if numel(stageZ) ~= size(raw4,4)
        error('loadFullStackISMData:MissingStageZ', ...
            'stageZUm must contain one recorded position per plane.');
    end
    head = struct();
    if ~isempty(opts.xyPixelSizeUm)
        head.ImgHdr_PixResol = opts.xyPixelSizeUm;
    end
    data = assembleData(raw4, stageZ, opts.channelIDs, head, opts);
    data.alignmentShiftXPx = zeros(1,size(raw4,4));
    data.alignmentShiftYPx = zeros(1,size(raw4,4));
    data.modelShiftToRawXPx = zeros(1,size(raw4,4));
    data.modelShiftToRawYPx = zeros(1,size(raw4,4));
    data.targetCenterXY = [(size(raw4,2)+1)/2, (size(raw4,1)+1)/2];
    data = attachBackgroundModel(data, opts);
    data = attachFlatField(data, opts);
end

function data = assembleData(rawCounts, stageZ, channelIDs, head, opts)
    if any(rawCounts(:) < 0) || any(~isfinite(rawCounts(:)))
        error('loadFullStackISMData:BadCounts', ...
            'Raw counts must be finite and nonnegative.');
    end
    data = struct();
    data.rawCounts = rawCounts;
    data.stageZUm = double(stageZ(:)).';
    data.channelIDs = double(channelIDs(:)).';
    data.head = head;
    data.xyPixelSizeUm = resolvePixelSize(head, opts);
    data.xUm = centeredAxis(size(rawCounts,2), data.xyPixelSizeUm);
    data.yUm = centeredAxis(size(rawCounts,1), data.xyPixelSizeUm);
    data.boundaryBackgroundPerPixel = measuredBoundaryBackground( ...
        rawCounts, opts.boundaryWidthPx);
    data = useBoundaryBackground(data);
    data.axialRawCounts = squeeze(sum(sum(sum(rawCounts,1),2),3)).';
    data.imageSize = size(rawCounts);
    data.flatFieldGain = ones(1,size(rawCounts,3));
    data.flatFieldCalibrated = false;
    data.flatFieldSource = 'unity';
end

function data = attachBackgroundModel(data, opts)
    if ~isfield(data,'boundaryBackgroundPerPixel') || ...
            isempty(data.boundaryBackgroundPerPixel)
        data.boundaryBackgroundPerPixel = measuredBoundaryBackground( ...
            data.rawCounts,opts.boundaryWidthPx);
    end

    if strcmp(opts.backgroundMode,'boundary')
        data = useBoundaryBackground(data);
        return;
    end

    try
        if ~isempty(opts.darkPerPixel)
            darkPerPixel = double(opts.darkPerPixel(:));
            source = 'explicit darkPerPixel option';
            exposureCalibrated = true;
            exposureScale = 1;
            darkInfo = struct('method','explicit darkPerPixel option', ...
                'expectedCountsPerPixel',darkPerPixel);
        elseif ~isempty(opts.darkFile) && exist(opts.darkFile,'file') == 2
            [darkPerPixel,source,exposureCalibrated,exposureScale,darkInfo] = ...
                readIndependentDarkPerPixel( ...
                opts.darkFile,data,opts);
        else
            error('loadFullStackISMData:MissingDarkFile', ...
                'Independent dark PTU was not found: %s',opts.darkFile);
        end

        nCh = size(data.rawCounts,3);
        if numel(darkPerPixel) ~= nCh || any(~isfinite(darkPerPixel)) || ...
                any(darkPerPixel < 0)
            error('loadFullStackISMData:BadDarkVector', ...
                'Independent dark background must contain %d nonnegative values.',nCh);
        end
        nPlane = size(data.rawCounts,4);
        data.backgroundPerPixel = repmat(reshape( ...
            opts.darkScale*darkPerPixel,1,1,nCh,1),1,1,1,nPlane);
        data.backgroundIndependent = true;
        data.backgroundExposureCalibrated = exposureCalibrated;
        data.backgroundExposureScale = exposureScale;
        data.backgroundSource = source;
        data.darkScale = opts.darkScale;
        darkInfo.appliedCountsPerPixel = opts.darkScale*darkPerPixel(:);
        darkInfo.channelIDs = data.channelIDs(:);
        data.darkDiagnostics = darkInfo;
        data = updateBackgroundTotals(data);
        if ~exposureCalibrated
            message = ['Dark counts were loaded, but acquisition-time ' ...
                'normalization could not be established from the PTU headers.'];
            if strcmp(opts.backgroundMode,'dark')
                error('loadFullStackISMData:UncalibratedDarkExposure',message);
            end
            warning('loadFullStackISMData:UncalibratedDarkExposure', ...
                '%s Acceptance will fail.',message);
        end
    catch darkError
        if strcmp(opts.backgroundMode,'dark')
            rethrow(darkError);
        end
        warning('loadFullStackISMData:DarkFallback', ...
            ['Independent dark background could not be loaded (%s). ' ...
            'Using boundary estimates; acceptance will fail.'],darkError.message);
        data = useBoundaryBackground(data);
        data.backgroundFallbackError = darkError.message;
    end
end

function data = useBoundaryBackground(data)
    if ~isfield(data,'boundaryBackgroundPerPixel')
        if isfield(data,'backgroundPerPixel')
            data.boundaryBackgroundPerPixel = data.backgroundPerPixel;
        else
            error('loadFullStackISMData:MissingBoundaryBackground', ...
                'No boundary background is stored in the stack data.');
        end
    end
    data.backgroundPerPixel = data.boundaryBackgroundPerPixel;
    data.backgroundIndependent = false;
    data.backgroundExposureCalibrated = false;
    data.backgroundExposureScale = NaN;
    data.backgroundSource = 'per-plane image-boundary estimate';
    data.darkScale = NaN;
    data.darkDiagnostics = struct('method','boundary fallback');
    if isfield(data,'axialRawCounts')
        data = updateBackgroundTotals(data);
    end
end

function data = updateBackgroundTotals(data)
    bgTotals = squeeze(sum(data.backgroundPerPixel,3)) * ...
        size(data.rawCounts,1) * size(data.rawCounts,2);
    data.axialBackgroundCounts = double(bgTotals(:)).';
    data.axialSignalCounts = max(data.axialRawCounts - ...
        data.axialBackgroundCounts,0);
end

function [darkPerPixel,source,exposureCalibrated,exposureScale,darkInfo] = ...
        readIndependentDarkPerPixel(fileName,data,opts)
    mode=opts.darkMeasurementMode;
    if strcmp(mode,'auto')
        head=PTU_Read_Head(fileName);
        isScan=isstruct(head) && isfield(head,'ImgHdr_BiDirect') && ...
            isfield(head,'ImgHdr_PixX') && isfield(head,'ImgHdr_PixY');
        if isScan
            mode='image';
        else
            mode='tttr';
        end
    end

    if strcmp(mode,'tttr')
        [darkPerPixel,exposureCalibrated,exposureScale,darkInfo] = ...
            readDarkTttrPerPixel(fileName,data,opts);
        source = sprintf(['independent TTTR point dark counts: %s ' ...
            '(%.6g counts/pixel exposure scale)'],fileName,exposureScale);
        if opts.verbose
            fprintf(['[loadFullStackISMData] Read point dark PTU directly ' ...
                'with PTU_Read: %.3f s, scan dwell %.6g ms/pixel\n'], ...
                darkInfo.darkDurationSeconds,darkInfo.scanPixelDwellMs);
        end
        return;
    end

    try
        darkOpts = opts;
        darkOpts.channelIDs = data.channelIDs;
        darkOpts.channelOrder = [];
        [darkPlane,darkIDs,darkHead,darkReadInfo] = ...
            readPtuDetectorPlane(fileName,[],darkOpts);
        if ~isequal(darkIDs(:),data.channelIDs(:))
            [present,loc] = ismember(data.channelIDs(:),darkIDs(:));
            if ~all(present)
                error('loadFullStackISMData:DarkChannelMismatch', ...
                    'Dark PTU is missing requested detector channels.');
            end
            darkPlane = darkPlane(:,:,loc);
        end
        [exposureScale,exposureCalibrated] = imageDarkExposureScale( ...
            data.head,darkHead,darkReadInfo.nFramesCombined);
        darkPerPixel = exposureScale * ...
            squeeze(mean(mean(max(double(darkPlane),0),1),2));
        source = sprintf(['independent detector dark image: %s ' ...
            '(automatic exposure scale %.6g)'],fileName,exposureScale);
        darkInfo = struct('method','detector dark image', ...
            'expectedCountsPerPixel',darkPerPixel(:), ...
            'nFramesCombined',darkReadInfo.nFramesCombined);
    catch scanError
        if strcmp(opts.darkMeasurementMode,'image')
            rethrow(scanError);
        end
        [darkPerPixel,exposureCalibrated,exposureScale,darkInfo] = ...
            readDarkTttrPerPixel(fileName,data,opts);
        source = sprintf(['independent TTTR dark counts: %s ' ...
            '(automatic exposure scale %.6g)'],fileName,exposureScale);
        if opts.verbose
            fprintf(['[loadFullStackISMData] Dark PTU is not an image scan; ' ...
                'using TTTR channel counts. Scan-reader error: %s\n'], ...
                scanError.message);
        end
    end
end

function [darkPerPixel,exposureCalibrated,exposureScale,darkInfo] = ...
        readDarkTttrPerPixel(fileName,data,opts)
    if exist('PTU_Read_Head','file') ~= 2 || exist('PTU_Read','file') ~= 2
        error('loadFullStackISMData:MissingTttrReader', ...
            'PTU_Read_Head/PTU_Read is required for a non-image dark PTU.');
    end
    head = PTU_Read_Head(fileName);
    if isempty(head) || ~isfield(head,'TTResult_NumberOfRecords')
        error('loadFullStackISMData:BadDarkPtu', ...
            'Dark PTU has no TTResult_NumberOfRecords field.');
    end

    channelIDs = double(data.channelIDs(:));
    nBins = max([64;channelIDs+1]);
    counts = zeros(nBins,1);
    validPhoton = ptuPhotonValidity(head);
    nRecords = double(head.TTResult_NumberOfRecords);
    cursor = 0;
    while cursor < nRecords
        nRead = min(double(opts.ptuPhotonsPerChunk),nRecords-cursor);
        [~,tcspc,chan,special,num] = PTU_Read( ...
            fileName,[cursor+1,nRead],head);
        cursor = cursor+double(num);
        if num <= 0
            break;
        end
        keep = special == 0 & validPhoton(chan,tcspc);
        chan = double(chan(keep));
        chan = chan(isfinite(chan) & chan >= 0 & chan < nBins);
        if ~isempty(chan)
            counts = counts+accumarray(chan(:)+1,1,[nBins,1],@sum,0);
        end
    end
    if any(channelIDs < 0 | channelIDs+1 > nBins)
        error('loadFullStackISMData:DarkChannelMismatch', ...
            'Dark counts do not cover all requested detector channels.');
    end
    dataDwellMs = numericField(data.head,'ImgHdr_TimePerPixel',NaN);
    darkDurationMs = numericField(head,'TTResult_StopAfter',NaN);
    exposureCalibrated = isfinite(dataDwellMs) && dataDwellMs > 0 && ...
        isfinite(darkDurationMs) && darkDurationMs > 0;
    if exposureCalibrated
        exposureScale = dataDwellMs/darkDurationMs;
    else
        exposureScale = 1;
    end
    selectedCounts=counts(channelIDs+1);
    darkPerPixel = exposureScale*selectedCounts;
    darkInfo=struct();
    darkInfo.method='TTTR point measurement';
    darkInfo.totalPhotonCounts=selectedCounts(:);
    darkInfo.darkDurationMs=darkDurationMs;
    darkInfo.darkDurationSeconds=darkDurationMs/1000;
    darkInfo.scanPixelDwellMs=dataDwellMs;
    darkInfo.countRateHz=selectedCounts(:)/(darkDurationMs/1000);
    darkInfo.expectedCountsPerPixel=darkPerPixel(:);
    darkInfo.totalRecords=nRecords;
    darkInfo.recordsRead=cursor;
end

function [scale,calibrated] = imageDarkExposureScale( ...
        dataHead,darkHead,nCombinedFrames)
    nCombinedFrames = max(1,double(nCombinedFrames));
    dataDwellMs = numericField(dataHead,'ImgHdr_TimePerPixel',NaN);
    darkDwellMs = numericField(darkHead,'ImgHdr_TimePerPixel',NaN);
    calibrated = isfinite(dataDwellMs) && dataDwellMs > 0 && ...
        isfinite(darkDwellMs) && darkDwellMs > 0;
    if calibrated
        scale = dataDwellMs/darkDwellMs/nCombinedFrames;
    else
        scale = 1/nCombinedFrames;
    end
end

function validPhoton = ptuPhotonValidity(head)
    rawResolution = numericField(head,'MeasDesc_Resolution',NaN);
    globalResolution = numericField(head,'MeasDesc_GlobalResolution',NaN);
    if isfinite(rawResolution) && rawResolution > 0 && ...
            isfinite(globalResolution) && globalResolution > 0
        resolutionNs = max(1e9*rawResolution,0.128);
        chDiv = max(1,round((resolutionNs*1e-9)/rawResolution));
        nGate = min(1024,ceil(1e9*globalResolution/resolutionNs)+1);
        validPhoton = @(chan,tcspc) (chan < 32) & (tcspc < nGate*chDiv);
    else
        validPhoton = @(chan,tcspc) chan < 32;
    end
end

function bg = measuredBoundaryBackground(rawCounts, width)
    width = max(1, min(round(width), ...
        floor(min(size(rawCounts,1),size(rawCounts,2))/2)));
    nCh = size(rawCounts,3);
    nPlane = size(rawCounts,4);
    bg = zeros(1,1,nCh,nPlane);
    for iz = 1:nPlane
        for ch = 1:nCh
            img = rawCounts(:,:,ch,iz);
            border = [reshape(img(1:width,:),[],1); ...
                reshape(img(end-width+1:end,:),[],1); ...
                reshape(img(:,1:width),[],1); ...
                reshape(img(:,end-width+1:end),[],1)];
            bg(1,1,ch,iz) = max(median(border, 'omitnan'), 0);
        end
    end
end

function data = attachFlatField(data, opts)
    [gain, calibrated, source] = resolveFlatField( ...
        opts.flatField, data.channelIDs, opts);
    data.flatFieldGain = gain(:).';
    data.flatFieldCalibrated = calibrated;
    data.flatFieldSource = source;
end

function [gain, calibrated, source] = resolveFlatField(inputValue, channelIDs, opts)
    nCh = numel(channelIDs);
    gain = ones(nCh,1);
    calibrated = false;
    source = 'unity; no independent flat-field supplied';
    if isempty(inputValue)
        return;
    end
    if isnumeric(inputValue)
        values = double(inputValue(:));
        source = 'numeric flatField option';
    else
        fileName = char(inputValue);
        if exist(fileName, 'file') ~= 2
            error('loadFullStackISMData:MissingFlatField', ...
                'Flat-field input was not found: %s', fileName);
        end
        [~,~,ext] = fileparts(fileName);
        switch lower(ext)
            case '.csv'
                T = readtable(fileName);
                values = flatFieldFromTable(T, channelIDs);
            case '.mat'
                S = load(fileName);
                values = flatFieldFromStruct(S, opts.flatFieldVariable);
            case '.ptu'
                [plane, ids] = readPtuDetectorPlane(fileName, [], opts);
                [present, loc] = ismember(channelIDs(:), ids(:));
                if ~all(present)
                    error('loadFullStackISMData:FlatFieldChannels', ...
                        'Flat-field PTU is missing requested detector channels.');
                end
                plane = plane(:,:,loc);
                bg = measuredBoundaryBackground(reshape(plane, ...
                    size(plane,1),size(plane,2),size(plane,3),1), ...
                    opts.boundaryWidthPx);
                bg = reshape(bg,1,1,[]);
                signal = squeeze(mean(mean(max(plane-bg,0),1),2));
                values = signal(:);
            otherwise
                error('loadFullStackISMData:BadFlatField', ...
                    'flatField must be numeric, CSV, MAT, or PTU.');
        end
        source = fileName;
    end
    if numel(values) ~= nCh || any(~isfinite(values)) || any(values <= 0)
        error('loadFullStackISMData:BadFlatField', ...
            'Flat-field gain must contain %d finite positive values.', nCh);
    end
    gain = values(:) / mean(values);
    calibrated = true;
end

function values = flatFieldFromTable(T, channelIDs)
    names = T.Properties.VariableNames;
    gainName = firstPresent(names, {'gain','qe','efficiency','relativeGain'});
    if isempty(gainName)
        error('loadFullStackISMData:BadFlatFieldCsv', ...
            'Flat-field CSV requires gain, qe, efficiency, or relativeGain.');
    end
    values = double(T.(gainName));
    channelName = firstPresent(names, {'channelID','channel','ChannelID'});
    if ~isempty(channelName)
        [present, loc] = ismember(channelIDs(:), double(T.(channelName)));
        if ~all(present)
            error('loadFullStackISMData:FlatFieldChannels', ...
                'Flat-field CSV is missing requested channels.');
        end
        values = values(loc);
    end
end

function values = flatFieldFromStruct(S, requested)
    if ~isempty(requested)
        if ~isfield(S, requested)
            error('loadFullStackISMData:BadFlatFieldMat', ...
                'Variable %s was not found.', requested);
        end
        values = double(S.(requested));
        return;
    end
    names = {'gain','gains','detectorGain','flatFieldGain','qe','efficiency'};
    for k = 1:numel(names)
        if isfield(S,names{k}) && isnumeric(S.(names{k}))
            values = double(S.(names{k}));
            return;
        end
    end
    error('loadFullStackISMData:BadFlatFieldMat', ...
        'No gain/QE vector was found in the MAT file.');
end

function [plane, channelIDs, head, readInfo] = ...
        readPtuDetectorPlane(fileName, frameIndex, opts)
    if exist(fileName, 'file') ~= 2
        error('loadFullStackISMData:MissingPtu', ...
            'PTU file was not found: %s', fileName);
    end
    waitbarCleanup = suppressPtuWaitbars(); %#ok<NASGU>
    ptuOut = [];
    if exist('PTU_MultiFrameScanReadFast','file') == 2
        try
            ptuOut = PTU_MultiFrameScanReadFast( ...
                fileName, opts.ptuPhotonsPerChunk, false, false);
        catch
            ptuOut = [];
        end
    end
    if isempty(ptuOut)
        if exist('PTU_MultiFrameScanRead','file') ~= 2
            error('loadFullStackISMData:MissingReader', ...
                'No PTU scan reader is available.');
        end
        ptuOut = PTU_MultiFrameScanRead(fileName, opts.ptuPhotonsPerChunk);
    end

    readInfo = struct('nFramesCombined',1);
    if ~isempty(frameIndex) && isfield(ptuOut,'tag') && ~isempty(ptuOut.tag)
        frameIndex = min(max(round(frameIndex),1),size(ptuOut.tag,4));
        raw = double(ptuOut.tag(:,:,:,frameIndex));
    elseif isfield(ptuOut,'tags') && ~isempty(ptuOut.tags)
        raw = double(ptuOut.tags);
    elseif isfield(ptuOut,'tag') && ~isempty(ptuOut.tag)
        readInfo.nFramesCombined = size(ptuOut.tag,4);
        raw = sum(double(ptuOut.tag),4);
    else
        error('loadFullStackISMData:NoImage', ...
            'No detector image was found in %s.', fileName);
    end
    if isfield(ptuOut,'dind')
        channelIDs = double(ptuOut.dind(:));
    else
        channelIDs = (1:size(raw,3)).';
    end
    plane = permute(raw,[2 1 3]);
    [plane, channelIDs] = selectChannels(plane, channelIDs, opts);
    head = struct();
    if isfield(ptuOut,'head'), head = ptuOut.head; end
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

function [plane, channelIDs] = selectChannels(plane, channelIDs, opts)
    requested = opts.channelIDs(:);
    if all(ismember(requested,channelIDs))
        [~,loc] = ismember(requested,channelIDs);
        plane = plane(:,:,loc);
        channelIDs = channelIDs(loc);
    elseif size(plane,3) == numel(requested)
        channelIDs = requested;
    else
        error('loadFullStackISMData:MissingChannels', ...
            'The PTU does not contain all requested detector channels.');
    end
    if ~isempty(opts.channelOrder)
        order = opts.channelOrder(:);
        if all(ismember(order,channelIDs))
            [~,loc] = ismember(order,channelIDs);
        elseif all(order >= 1 & order <= numel(channelIDs))
            loc = order;
        else
            error('loadFullStackISMData:BadChannelOrder', ...
                'channelOrder must contain channel IDs or detector indices.');
        end
        plane = plane(:,:,loc);
        channelIDs = channelIDs(loc);
    end
end

function pixelSize = resolvePixelSize(head, opts)
    pixelSize = opts.xyPixelSizeUm;
    if isempty(pixelSize) || ~isfinite(pixelSize) || pixelSize <= 0
        pixelSize = numericField(head,'ImgHdr_PixResol',NaN);
    end
    if ~isfinite(pixelSize) || pixelSize <= 0
        error('loadFullStackISMData:MissingPixelSize', ...
            'xyPixelSizeUm is required when the PTU header lacks ImgHdr_PixResol.');
    end
end

function center = resolveTargetCenter(T, rawSize)
    center = [(rawSize(2)+1)/2, (rawSize(1)+1)/2];
    if all(ismember({'target_center_x','target_center_y'}, ...
            T.Properties.VariableNames))
        x = double(T.target_center_x);
        y = double(T.target_center_y);
        x = x(isfinite(x));
        y = y(isfinite(y));
        if ~isempty(x) && ~isempty(y)
            center = [median(x), median(y)];
        end
    end
end

function raw = findDetectorStack(S)
    names = {'rawCounts','raw4','rawData','detectorVolume','detectorStack', ...
        'channelVolume','channelVolumes','spadVolume','spadStack'};
    for k = 1:numel(names)
        if isfield(S,names{k}) && isnumeric(S.(names{k})) && ...
                ndims(S.(names{k})) == 4
            raw = S.(names{k});
            return;
        end
    end
    fields = fieldnames(S);
    for k = 1:numel(fields)
        value = S.(fields{k});
        if isnumeric(value) && ndims(value) == 4 && ...
                (size(value,3)==23 || size(value,4)==23)
            raw = value;
            return;
        end
    end
    error('loadFullStackISMData:NoDetectorStack', ...
        'No detector-resolved 4-D stack was found.');
end

function stageZ = stageZFromStruct(S)
    stageZ = [];
    candidates = {'stageZUm','planeZ','zUm','z_um'};
    for k = 1:numel(candidates)
        if isfield(S,candidates{k}) && isnumeric(S.(candidates{k}))
            stageZ = double(S.(candidates{k})(:)).';
            return;
        end
    end
    if isfield(S,'scans') && isstruct(S.scans) && isfield(S.scans,'zUm')
        stageZ = double([S.scans.zUm]);
    end
end

function text = tableText(value, row)
    if iscell(value), text = char(value{row});
    elseif isstring(value), text = char(value(row));
    elseif iscategorical(value), text = char(value(row));
    else, text = char(value(row,:));
    end
end

function values = tableTextColumn(value)
    if iscell(value), values = value;
    elseif isstring(value), values = cellstr(value);
    elseif iscategorical(value), values = cellstr(value);
    else, values = cellstr(value);
    end
end

function frame = tableFrame(T,row)
    frame = [];
    if ismember('source_frame',T.Properties.VariableNames)
        value = double(T.source_frame(row));
        if isfinite(value) && value >= 1, frame = round(value); end
    end
end

function frames = tableFrameColumn(T)
    frames = nan(1,height(T));
    if ismember('source_frame',T.Properties.VariableNames)
        frames = double(T.source_frame(:)).';
    end
end

function values = tableNumericOrZeros(T,name)
    values = zeros(1,height(T));
    if ismember(name,T.Properties.VariableNames)
        values = double(T.(name)(:)).';
        values(~isfinite(values)) = 0;
    end
end

function name = firstPresent(names,candidates)
    name = '';
    for k = 1:numel(candidates)
        idx = find(strcmpi(names,candidates{k}),1);
        if ~isempty(idx)
            name = names{idx};
            return;
        end
    end
end

function axisValues = centeredAxis(n,pixelSize)
    axisValues = ((1:n)-(n+1)/2)*pixelSize;
end

function value = numericField(s,name,defaultValue)
    value = defaultValue;
    if isstruct(s) && isfield(s,name) && isnumeric(s.(name)) && ~isempty(s.(name))
        value = double(s.(name)(1));
    end
end

function out = setfieldLocal(s,name,value)
    out = s;
    out.(name) = value;
end
