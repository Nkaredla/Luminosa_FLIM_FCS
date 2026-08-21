function result = estimateExperimentalWavefrontFromPSF(dataPath, varargin)
%ESTIMATEEXPERIMENTALWAVEFRONTFROMPSF Fit effective aberrations from a measured 3-D PSF.
%
%   result = estimateExperimentalWavefrontFromPSF(dataPath)
%
%   dataPath can be either a scan folder or a saved MAT file. The preferred
%   input is detector-resolved data with size [y x 23 z]. If only the
%   existing batch-output scalar volume is available, the function fits the
%   detector-summed 3-D PSF and reports this as a summed-volume fit.
%
%   For detector-resolved input, per-channel dark count is estimated from
%   image-boundary pixels over all z planes. A hot detector channel is
%   detected by comparing that boundary level with both the global detector
%   population and neighbouring honeycomb detector pixels.

    if nargin < 1 || isempty(dataPath)
        dataPath = 'D:\Luminosa\Data\ISM_Aberation2_73\0.4uW_0.19collar_80mmlens_20260515-155248';
    end
    dataPath = char(dataPath);

    addLocalPaths();
    opts = parseExperimentalFitOptions(varargin{:});

    [data, meta] = loadExperimentalPsfInput(dataPath, opts);
    [data, squareInfo] = makePsfDataSquare(data, opts);
    meta.squareInfo = squareInfo;
    if opts.verbose && squareInfo.applied
        fprintf('[estimateExperimentalWavefrontFromPSF] %s data from %dx%d to %dx%d.\n', ...
            squareInfo.action, squareInfo.originalSize(1), squareInfo.originalSize(2), ...
            squareInfo.outputSize(1), squareInfo.outputSize(2));
    end
    sim = configureExperimentalSim(data, meta, opts);

    planeZ = resolvePlaneZ(data, meta, opts);
    [fitData, fitPlaneZ, fitPlaneIndex] = selectFitPlanes(data, planeZ, opts);

    if meta.detectorResolved
        [dark, correctedData] = estimateDetectorDarkAndHotPixels(fitData, sim, opts);
        fitOpts = buildPhaseRetrievalOptions(sim, opts);
        [initialXY, fitOpts] = setInitialXYFromData(correctedData, sim, fitOpts, opts);
        fit = phaseRetrieval3DBead(true, correctedData, fitPlaneZ, fitOpts);
        fitType = 'detector-resolved';
        compare = detectorFitComparisonVolumes(correctedData, fit.fitStack);
    else
        if ~opts.projectionFallback
            error('estimateExperimentalWavefrontFromPSF:NeedDetectorData', ...
                ['The loaded data are a detector-summed 3-D volume. Provide a ' ...
                 '[y x 23 z] detector stack, or set projectionFallback=true.']);
        end
        [dark, correctedData] = estimateScalarVolumeBackground(fitData, opts);
        tcspcDark = estimateTcspcDetectorDarkFromFolder(dataPath, sim, opts);
        if ~isempty(tcspcDark)
            dark.tcspcDetector = tcspcDark;
            dark.hotChannels = tcspcDark.hotChannels;
        end
        fitRunOpts = opts;
        [initialXY, fitRunOpts] = setInitialXYFromData(correctedData, sim, fitRunOpts, opts);
        fit = fitSummedVolumeAberration(correctedData, fitPlaneZ, sim, fitRunOpts);
        fitType = 'summed-volume';
        compare = scalarFitComparisonVolumes(correctedData, fit.fitVolume);
    end

    result = struct();
    result.dataPath = dataPath;
    result.inputMeta = meta;
    result.fitType = fitType;
    result.fitPlaneIndex = fitPlaneIndex;
    result.planeZ = fitPlaneZ;
    result.dark = dark;
    result.initialXY = initialXY;
    result.rawFitData = fitData;
    result.correctedFitData = correctedData;
    result.fit = fit;
    result.compare = compare;
    result.options = opts;

    result.outputDir = resolveOutputDir(dataPath, meta, opts);
    if ~isempty(result.outputDir)
        writeExperimentalFitOutputs(result);
    end

    if opts.verbose
        fprintf('[estimateExperimentalWavefrontFromPSF] %s fit complete.\n', fitType);
        fprintf('  residual norm: %.4e\n', fit.residualNorm);
        printCoefficientSummary(fit);
        if isfield(dark, 'hotChannels') && ~isempty(dark.hotChannels)
            fprintf('  hot detector channel(s): %s\n', mat2str(dark.hotChannels(:).'));
        end
        if ~isempty(result.outputDir)
            fprintf('  wrote results to: %s\n', result.outputDir);
        end
    end
end

function addLocalPaths()
    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    parentDir = fileparts(thisDir);
    if exist(fullfile(parentDir, 'SPADarray_AberrationPSF'), 'dir') == 7
        addpath(parentDir);
    end
end

function opts = parseExperimentalFitOptions(varargin)
    p = inputParser;
    p.FunctionName = 'estimateExperimentalWavefrontFromPSF';
    addParameter(p, 'matFile', '');
    addParameter(p, 'dataVariable', '');
    addParameter(p, 'rawData', []);
    addParameter(p, 'planeZ', []);
    addParameter(p, 'outputRoot', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'zStepUm', 0.05);
    addParameter(p, 'xyPixelSizeUm', []);
    addParameter(p, 'objectiveNA', []);
    addParameter(p, 'emissionWavelengthUm', []);
    addParameter(p, 'excitationWavelengthUm', []);
    addParameter(p, 'mediumRefractiveIndex', []);
    addParameter(p, 'inferOpticsFromPtuHeader', true);
    addParameter(p, 'fitModes', {'defocus','astig_x','astig_y','coma_x','coma_y','spherical'});
    addParameter(p, 'maxIter', 8);
    addParameter(p, 'maxFitPlanes', 31);
    addParameter(p, 'fitXY', true);
    addParameter(p, 'fitZ', false);
    addParameter(p, 'initialXY', []);
    addParameter(p, 'estimateInitialXY', true);
    addParameter(p, 'initialXYThresholdFraction', 0.20);
    addParameter(p, 'projectionFallback', true);
    addParameter(p, 'useAlignedVolume', true);
    addParameter(p, 'squareDataMode', 'crop');
    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detectorPixelShape', 'hex');
    addParameter(p, 'boundaryWidthPx', 3);
    addParameter(p, 'estimateTcspcDark', true);
    addParameter(p, 'tcspcDarkPlaneFraction', 0.20);
    addParameter(p, 'hotThresholdMAD', 8);
    addParameter(p, 'hotNeighborThresholdMAD', 6);
    addParameter(p, 'hotPixelAction', 'replaceWithNeighbors');
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
    addParameter(p, 'verbose', true);
    parse(p, varargin{:});
    opts = p.Results;

    if ischar(opts.fitModes) || isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    end
    opts.fitModes = opts.fitModes(:).';
    opts.hotPixelAction = lower(char(opts.hotPixelAction));
end

function [data, meta] = loadExperimentalPsfInput(dataPath, opts)
    meta = struct();
    meta.sourcePath = dataPath;
    meta.sourceMatFile = '';
    meta.sourceVariable = '';
    meta.detectorResolved = false;
    meta.scanName = scanNameFromPath(dataPath);
    meta.xyPixelSizeUm = NaN;
    meta.planeZ = [];
    meta.headerOptics = readOpticsFromPtuHeader(dataPath, opts);

    if ~isempty(opts.rawData)
        sim0 = defaultParams();
        data = opts.rawData;
        [data, meta.detectorResolved] = standardizeLoadedPsfData(data, sim0.nDet, opts.planeZ);
        meta.sourceVariable = 'rawData option';
        return;
    end

    matFile = resolveInputMatFile(dataPath, opts);
    if isempty(matFile) || exist(matFile, 'file') ~= 2
        error('estimateExperimentalWavefrontFromPSF:NoInputMat', ...
            ['Could not find a saved PSF MAT file for "%s". Run batchAnalyzeLuminosaPSFs ' ...
             'first, or pass ''matFile''/''rawData'' with detector-resolved [y x 23 z] data.'], ...
            dataPath);
    end

    S = load(matFile);
    meta.sourceMatFile = matFile;
    [data, variableName] = choosePsfDataVariable(S, opts);
    meta.sourceVariable = variableName;

    sim0 = defaultParams();
    [data, meta.detectorResolved] = standardizeLoadedPsfData(data, sim0.nDet, opts.planeZ);

    if isfield(S, 'frameMeta') && isfield(S.frameMeta, 'zUm')
        meta.planeZ = [S.frameMeta.zUm];
    end

    meta.xyPixelSizeUm = xyPixelSizeFromLoadedStruct(S);
    if isfield(S, 'group') && isfield(S.group, 'label') && ~isempty(S.group.label)
        meta.scanName = char(S.group.label);
    end
end

function matFile = resolveInputMatFile(dataPath, opts)
    matFile = char(opts.matFile);
    if ~isempty(matFile)
        return;
    end

    if exist(dataPath, 'file') == 2
        [~,~,ext] = fileparts(dataPath);
        if strcmpi(ext, '.mat')
            matFile = dataPath;
            return;
        end
    end

    if exist(dataPath, 'dir') ~= 7
        matFile = '';
        return;
    end

    [scanParent, scanName] = fileparts(stripTrailingFilesep(dataPath));
    [dataParent, datasetName] = fileparts(scanParent);
    stem = scanOutputStem(scanName);
    datasetStem = sanitizeFileName(datasetName);
    timeToken = regexp(stem, '\d{8}_\d{6}', 'match', 'once');

    candidates = {};
    candidates{end+1} = fullfile(dataParent, 'PSF_batch_outputs', datasetStem, ...
        'volumes_mat', sprintf('%s_volume_raw.mat', stem));
    candidates{end+1} = fullfile(dataPath, 'psf_volume_aligned.mat');
    candidates{end+1} = fullfile(dataPath, sprintf('%s_volume_raw.mat', stem));

    localMat = dir(fullfile(dataPath, '*.mat'));
    for k = 1:numel(localMat)
        candidates{end+1} = fullfile(localMat(k).folder, localMat(k).name); %#ok<AGROW>
    end

    matFile = '';
    for k = 1:numel(candidates)
        if exist(candidates{k}, 'file') == 2
            matFile = candidates{k};
            return;
        end
    end

    searchDirs = { ...
        fullfile(dataParent, 'PSF_batch_outputs', datasetStem, 'volumes_mat'), ...
        fullfile(dataParent, 'PSF_batch_outputs', datasetStem), ...
        fullfile(scanParent, 'results_psf3d')};
    searchPatterns = {sprintf('*%s*volume*.mat', stem)};
    if ~isempty(timeToken)
        searchPatterns{end+1} = sprintf('*%s*volume*.mat', timeToken);
    end

    for d = 1:numel(searchDirs)
        if exist(searchDirs{d}, 'dir') ~= 7
            continue;
        end
        for q = 1:numel(searchPatterns)
            hits = dir(fullfile(searchDirs{d}, searchPatterns{q}));
            if ~isempty(hits)
                [~, newest] = max([hits.datenum]);
                matFile = fullfile(hits(newest).folder, hits(newest).name);
                return;
            end
        end
    end
end

function [data, variableName] = choosePsfDataVariable(S, opts)
    variableName = char(opts.dataVariable);
    if ~isempty(variableName)
        if ~isfield(S, variableName)
            error('estimateExperimentalWavefrontFromPSF:MissingVariable', ...
                'Variable "%s" was not found in the MAT file.', variableName);
        end
        data = S.(variableName);
        return;
    end

    detectorNames = {'raw4','rawData','detectorVolume','detectorStack', ...
        'channelVolume','channelVolumes','spadVolume','spadStack'};
    for k = 1:numel(detectorNames)
        name = detectorNames{k};
        if isfield(S, name) && isnumeric(S.(name)) && ndims(S.(name)) >= 3
            data = S.(name);
            variableName = name;
            return;
        end
    end

    names = fieldnames(S);
    for k = 1:numel(names)
        value = S.(names{k});
        if isnumeric(value) && ndims(value) == 4 && any(size(value) == 23)
            data = value;
            variableName = names{k};
            return;
        end
    end

    scalarPriority = {};
    if opts.useAlignedVolume
        scalarPriority = {'volume','alignedVolume','rawVolume'};
    else
        scalarPriority = {'rawVolume','volume','alignedVolume'};
    end

    for k = 1:numel(scalarPriority)
        name = scalarPriority{k};
        if isfield(S, name) && isnumeric(S.(name)) && ndims(S.(name)) == 3
            data = S.(name);
            variableName = name;
            return;
        end
    end

    for k = 1:numel(names)
        value = S.(names{k});
        if isnumeric(value) && ndims(value) == 3
            data = value;
            variableName = names{k};
            return;
        end
    end

    error('estimateExperimentalWavefrontFromPSF:NoImageVariable', ...
        'No numeric 3-D volume or 4-D detector stack was found in the MAT file.');
end

function optics = readOpticsFromPtuHeader(dataPath, opts)
    optics = struct('source', '', 'headerFile', '', 'NA', NaN, ...
        'lamExc', NaN, 'lamEm', NaN, 'matchedFields', struct());

    if ~opts.inferOpticsFromPtuHeader || exist('PTU_Read_Head', 'file') ~= 2 || ...
            exist(dataPath, 'dir') ~= 7
        return;
    end

    files = dir(fullfile(dataPath, 'Series*.ptu'));
    if isempty(files)
        files = dir(fullfile(dataPath, '*.ptu'));
    end
    if isempty(files)
        return;
    end

    idx = zeros(numel(files), 1);
    for k = 1:numel(files)
        idx(k) = ptuSeriesFileIndex(files(k).name);
    end
    [~, order] = sort(idx);
    fileName = fullfile(files(order(1)).folder, files(order(1)).name);

    try
        head = PTU_Read_Head(fileName);
    catch
        return;
    end

    [optics.NA, optics.matchedFields.NA] = firstHeaderNumericByRegex(head, ...
        {'(^|_)NA($|_)', 'Objective.*NA', 'Obj.*NA', 'Numerical.*Aperture'}, ...
        'na');
    [optics.lamExc, optics.matchedFields.lamExc] = firstHeaderNumericByRegex(head, ...
        {'Exc.*Wave', 'Excitation.*Wave', 'Laser.*Wave', 'Lambda.*Exc', 'Wavelength.*Exc'}, ...
        'wavelength');
    [optics.lamEm, optics.matchedFields.lamEm] = firstHeaderNumericByRegex(head, ...
        {'Em.*Wave', 'Emission.*Wave', 'Detection.*Wave', 'Lambda.*Em', 'Wavelength.*Em'}, ...
        'wavelength');

    optics.source = 'PTU header';
    optics.headerFile = fileName;
end

function idx = ptuSeriesFileIndex(name)
    tok = regexp(name, '^Series_(\d+)\.ptu$', 'tokens', 'once');
    if isempty(tok)
        idx = inf;
    else
        idx = str2double(tok{1});
    end
end

function [value, fieldName] = firstHeaderNumericByRegex(head, patterns, kind)
    value = NaN;
    fieldName = '';
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
                fieldName = names{k};
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
    tf = isfinite(v) && v > 0;
    if ~tf
        return;
    end
    switch lower(kind)
        case 'na'
            tf = v >= 0.1 && v <= 2.0;
        case 'wavelength'
            vu = normalizeHeaderValue(v, kind);
            tf = vu >= 0.3 && vu <= 1.1;
    end
end

function v = normalizeHeaderValue(v, kind)
    switch lower(kind)
        case 'wavelength'
            if v > 10
                v = v / 1000; % nm -> um
            end
    end
end

function [data, detectorResolved] = standardizeLoadedPsfData(data, nDet, planeZ)
    detectorResolved = false;

    if iscell(data)
        if numel(data) ~= nDet
            error('estimateExperimentalWavefrontFromPSF:BadCellData', ...
                'Cell-array detector data must contain %d channel volumes.', nDet);
        end
        first = double(data{1});
        if ndims(first) == 2
            first = reshape(first, size(first,1), size(first,2), 1);
        end
        [ny, nx, nz] = size(first);
        raw4 = zeros(ny, nx, nDet, nz);
        raw4(:,:,1,:) = reshape(first, ny, nx, 1, nz);
        for k = 2:nDet
            vol = double(data{k});
            if ndims(vol) == 2
                vol = reshape(vol, size(vol,1), size(vol,2), 1);
            end
            if size(vol,1) ~= ny || size(vol,2) ~= nx || size(vol,3) ~= nz
                error('estimateExperimentalWavefrontFromPSF:BadCellVolume', ...
                    'All detector channel volumes must have the same [y x z] size.');
            end
            raw4(:,:,k,:) = reshape(vol, ny, nx, 1, nz);
        end
        data = raw4;
        detectorResolved = true;
        return;
    end

    data = double(data);

    dims = size(data);
    if ndims(data) == 4
        if dims(3) == nDet
            detectorResolved = true;
        elseif dims(4) == nDet
            data = permute(data, [1 2 4 3]);
            detectorResolved = true;
        else
            error('estimateExperimentalWavefrontFromPSF:BadDetectorStack', ...
                '4-D data must contain a %d-channel detector dimension.', nDet);
        end
        return;
    end

    if ndims(data) == 3
        if dims(3) == nDet && (~isempty(planeZ) && numel(planeZ) == 1)
            data = reshape(data, dims(1), dims(2), nDet, 1);
            detectorResolved = true;
        else
            detectorResolved = false;
        end
        return;
    end

    error('estimateExperimentalWavefrontFromPSF:BadDataShape', ...
        'Expected a 3-D scalar volume or a 4-D [y x 23 z] detector stack.');
end

function [data, info] = makePsfDataSquare(data, opts)
    info = struct();
    info.mode = char(opts.squareDataMode);
    info.originalSize = size(data);
    info.applied = false;

    ny = size(data, 1);
    nx = size(data, 2);
    if ny == nx
        info.outputSize = size(data);
        info.yIndex = 1:ny;
        info.xIndex = 1:nx;
        return;
    end

    mode = lower(char(opts.squareDataMode));
    switch mode
        case {'crop','centercrop','center_crop'}
            n = min(ny, nx);
            y0 = floor((ny - n) / 2) + 1;
            x0 = floor((nx - n) / 2) + 1;
            yIdx = y0:(y0+n-1);
            xIdx = x0:(x0+n-1);
            if ndims(data) == 4
                data = data(yIdx, xIdx, :, :);
            else
                data = data(yIdx, xIdx, :);
            end
            info.applied = true;
            info.action = 'center crop';
            info.yIndex = yIdx;
            info.xIndex = xIdx;

        case {'pad','zeropad','zero_pad'}
            n = max(ny, nx);
            y0 = floor((n - ny) / 2) + 1;
            x0 = floor((n - nx) / 2) + 1;
            if ndims(data) == 4
                padded = zeros(n, n, size(data,3), size(data,4));
                padded(y0:(y0+ny-1), x0:(x0+nx-1), :, :) = data;
            else
                padded = zeros(n, n, size(data,3));
                padded(y0:(y0+ny-1), x0:(x0+nx-1), :) = data;
            end
            data = padded;
            info.applied = true;
            info.action = 'zero pad';
            info.yIndex = y0:(y0+ny-1);
            info.xIndex = x0:(x0+nx-1);

        case {'error','strict'}
            error('estimateExperimentalWavefrontFromPSF:NonSquareData', ...
                ['Loaded PSF data are %d x %d. The current forward model expects ' ...
                 'square y-x images. Set squareDataMode to ''crop'' or ''pad''.'], ...
                ny, nx);

        otherwise
            error('estimateExperimentalWavefrontFromPSF:BadSquareDataMode', ...
                'Unknown squareDataMode "%s". Use ''crop'', ''pad'', or ''error''.', ...
                opts.squareDataMode);
    end

    info.outputSize = size(data);
end

function xyPixelSizeUm = xyPixelSizeFromLoadedStruct(S)
    xyPixelSizeUm = NaN;
    if isfield(S, 'opts') && isstruct(S.opts) && isfield(S.opts, 'xyPixelSizeUm')
        v = S.opts.xyPixelSizeUm;
        if isnumeric(v) && isscalar(v) && isfinite(v) && v > 0
            xyPixelSizeUm = double(v);
            return;
        end
    end
    if isfield(S, 'scans') && isstruct(S.scans) && isfield(S.scans, 'xyPixelSizeUm')
        vals = [S.scans.xyPixelSizeUm];
        vals = vals(isfinite(vals) & vals > 0);
        if ~isempty(vals)
            xyPixelSizeUm = double(vals(1));
        end
    end
end

function sim = configureExperimentalSim(data, meta, opts)
    sim = defaultParams();
    sim = applyOpticsOptions(sim, meta, opts);
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

    xyPixelSizeUm = opts.xyPixelSizeUm;
    if isempty(xyPixelSizeUm) || ~isnumeric(xyPixelSizeUm) || ...
            ~isscalar(xyPixelSizeUm) || ~isfinite(xyPixelSizeUm) || xyPixelSizeUm <= 0
        xyPixelSizeUm = meta.xyPixelSizeUm;
    end

    nxy = size(data, 1);
    if size(data, 2) ~= nxy
        error('estimateExperimentalWavefrontFromPSF:NonSquareData', ...
            'The current forward model expects square y-x PSF images.');
    end

    if ~isfinite(xyPixelSizeUm) || xyPixelSizeUm <= 0
        error('estimateExperimentalWavefrontFromPSF:MissingLateralScale', ...
            ['Measured PSF fitting requires xyPixelSizeUm from the PTU/MAT ' ...
            'metadata or an explicit option.']);
    end
    sim.xyPixelSizeUm = xyPixelSizeUm;
    sim.fovXY = xyPixelSizeUm * (nxy - 1);
    sim = setSimLateralGrid(sim, nxy);
end

function sim = setSimLateralGrid(sim, nxy)
    sim.nx = nxy;
    sim.x = linspace(-sim.fovXY/2, sim.fovXY/2, nxy);
    sim.y = sim.x;
    if nxy > 1
        sim.dx = abs(sim.x(2) - sim.x(1));
    else
        sim.dx = sim.fovXY;
    end
    sim.obj = beadObject3D(sim);
end

function sim = applyOpticsOptions(sim, meta, opts)
    optics = struct();
    if isfield(meta, 'headerOptics') && ~isempty(meta.headerOptics)
        optics = meta.headerOptics;
    end

    sim.opticsSource = 'defaultParams';
    if isfield(optics, 'NA') && isfinite(optics.NA) && optics.NA > 0
        sim.NA = optics.NA;
        sim.opticsSource = 'PTU header';
    end
    if isfield(optics, 'lamEm') && isfinite(optics.lamEm) && optics.lamEm > 0
        sim.lamEm = optics.lamEm;
        sim.lamRef = optics.lamEm;
        sim.opticsSource = 'PTU header';
    end
    if isfield(optics, 'lamExc') && isfinite(optics.lamExc) && optics.lamExc > 0
        sim.lamExc = optics.lamExc;
        sim.opticsSource = 'PTU header';
    end

    if isPositiveScalar(opts.objectiveNA)
        sim.NA = double(opts.objectiveNA);
        sim.opticsSource = 'explicit options';
    end
    if isPositiveScalar(opts.emissionWavelengthUm)
        sim.lamEm = double(opts.emissionWavelengthUm);
        sim.lamRef = double(opts.emissionWavelengthUm);
        sim.opticsSource = 'explicit options';
    end
    if isPositiveScalar(opts.excitationWavelengthUm)
        sim.lamExc = double(opts.excitationWavelengthUm);
        sim.opticsSource = 'explicit options';
    end
    if isPositiveScalar(opts.mediumRefractiveIndex)
        sim.nMedium = double(opts.mediumRefractiveIndex);
        sim.opticsSource = 'explicit options';
    end

    if opts.verbose
        fprintf(['[estimateExperimentalWavefrontFromPSF] optics: NA=%.3f, ' ...
            'lambda_exc=%.3f um, lambda_em/ref=%.3f um, n=%.3f (%s).\n'], ...
            sim.NA, sim.lamExc, sim.lamEm, sim.nMedium, sim.opticsSource);
    end
end

function tf = isPositiveScalar(v)
    tf = isnumeric(v) && isscalar(v) && isfinite(v) && v > 0;
end

function [initialXY, runOpts] = setInitialXYFromData(data, sim, runOpts, opts)
    if ~opts.fitXY
        initialXY = [0 0];
        return;
    end

    if isPositiveVector2(opts.initialXY)
        initialXY = double(opts.initialXY(:).');
    elseif opts.estimateInitialXY
        initialXY = estimateInitialXYFromData(data, sim, opts);
    else
        initialXY = [0 0];
    end

    runOpts.initialXY = initialXY;
    if opts.verbose
        fprintf('[estimateExperimentalWavefrontFromPSF] initial XY shift: x=%+.4f um, y=%+.4f um.\n', ...
            initialXY(1), initialXY(2));
    end
end

function tf = isPositiveVector2(v)
    tf = isnumeric(v) && numel(v) == 2 && all(isfinite(v(:)));
end

function xy = estimateInitialXYFromData(data, sim, opts)
    if ndims(data) == 4
        img = squeeze(sum(sum(max(data, 0), 3), 4));
    else
        img = sum(max(data, 0), 3);
    end

    img = double(img);
    img = img - min(img(:));
    peak = max(img(:));
    if ~isfinite(peak) || peak <= 0
        xy = [0 0];
        return;
    end

    thresh = opts.initialXYThresholdFraction * peak;
    weights = img;
    weights(img < thresh) = 0;
    sw = sum(weights(:));
    if sw <= 0
        weights = img;
        sw = sum(weights(:));
    end

    [X, Y] = meshgrid(sim.x, sim.y);
    xy = [sum(X(:).*weights(:))/sw, sum(Y(:).*weights(:))/sw];
end

function planeZ = resolvePlaneZ(data, meta, opts)
    if ~isempty(opts.planeZ)
        planeZ = opts.planeZ(:).';
    elseif isfield(meta, 'planeZ') && ~isempty(meta.planeZ)
        planeZ = meta.planeZ(:).';
    else
        nPlane = numberOfPlanes(data);
        planeZ = ((1:nPlane) - (nPlane+1)/2) * opts.zStepUm;
    end

    if numel(planeZ) ~= numberOfPlanes(data)
        error('estimateExperimentalWavefrontFromPSF:PlaneZMismatch', ...
            'planeZ has %d entries but data have %d z planes.', ...
            numel(planeZ), numberOfPlanes(data));
    end
end

function nPlane = numberOfPlanes(data)
    if ndims(data) == 4
        nPlane = size(data, 4);
    else
        nPlane = size(data, 3);
    end
end

function [fitData, fitPlaneZ, fitPlaneIndex] = selectFitPlanes(data, planeZ, opts)
    nPlane = numberOfPlanes(data);
    maxFitPlanes = opts.maxFitPlanes;
    if isempty(maxFitPlanes) || ~isfinite(maxFitPlanes) || maxFitPlanes <= 0 || maxFitPlanes >= nPlane
        fitPlaneIndex = 1:nPlane;
    else
        maxFitPlanes = max(3, round(maxFitPlanes));
        sig = zSignalTrace(data);
        [~, z0] = max(sig);
        first = max(1, min(z0 - floor((maxFitPlanes-1)/2), nPlane - maxFitPlanes + 1));
        fitPlaneIndex = first:(first + maxFitPlanes - 1);
    end

    fitPlaneZ = planeZ(fitPlaneIndex);
    if ndims(data) == 4
        fitData = data(:,:,:,fitPlaneIndex);
    else
        fitData = data(:,:,fitPlaneIndex);
    end
end

function sig = zSignalTrace(data)
    if ndims(data) == 4
        sig = squeeze(sum(sum(sum(max(data, 0), 1), 2), 3));
    else
        sig = squeeze(sum(sum(max(data, 0), 1), 2));
    end
    sig = sig(:).';
end

function fitOpts = buildPhaseRetrievalOptions(sim, opts)
    fitOpts = struct();
    fitOpts.sim = sim;
    fitOpts.fitModes = opts.fitModes;
    fitOpts.maxIter = opts.maxIter;
    fitOpts.fitXY = opts.fitXY;
    fitOpts.fitZ = opts.fitZ;
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
end

function [dark, corrected] = estimateDetectorDarkAndHotPixels(raw4, sim, opts)
    raw4 = double(raw4);
    [ny, nx, nCh, nz] = size(raw4);
    width = max(1, min([opts.boundaryWidthPx, floor(ny/4), floor(nx/4)]));
    borderMask = false(ny, nx);
    borderMask(1:width, :) = true;
    borderMask(end-width+1:end, :) = true;
    borderMask(:, 1:width) = true;
    borderMask(:, end-width+1:end) = true;

    boundaryMedian = zeros(nCh, 1);
    boundaryMad = zeros(nCh, 1);
    for k = 1:nCh
        vals = [];
        for z = 1:nz
            frame = raw4(:,:,k,z);
            vals = [vals; frame(borderMask)]; %#ok<AGROW>
        end
        vals = vals(isfinite(vals));
        if isempty(vals)
            boundaryMedian(k) = 0;
            boundaryMad(k) = 0;
        else
            boundaryMedian(k) = median(vals);
            boundaryMad(k) = robustMad(vals);
        end
    end

    neighbors = detectorNeighborList(sim.detXY, sim.detPitch);
    neighborMedian = nan(nCh, 1);
    neighborMad = nan(nCh, 1);
    for k = 1:nCh
        nb = neighbors{k};
        if isempty(nb)
            nb = setdiff(1:nCh, k);
        end
        neighborMedian(k) = median(boundaryMedian(nb));
        neighborMad(k) = robustMad(boundaryMedian(nb));
    end

    globalMedian = median(boundaryMedian);
    globalMad = max(robustMad(boundaryMedian), eps);
    hotScoreGlobal = (boundaryMedian - globalMedian) / globalMad;
    hotScoreNeighbor = (boundaryMedian - neighborMedian) ./ max(neighborMad, eps);
    hotMask = hotScoreGlobal > opts.hotThresholdMAD | ...
              hotScoreNeighbor > opts.hotNeighborThresholdMAD;
    hotChannels = find(hotMask);

    corrected = raw4 - reshape(boundaryMedian, 1, 1, nCh, 1);
    corrected = max(corrected, 0);

    switch opts.hotPixelAction
        case {'none','subtractonly','subtract_only'}
            % Keep the dark-corrected hot channel in the fit.
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
            error('estimateExperimentalWavefrontFromPSF:BadHotPixelAction', ...
                'Unknown hotPixelAction "%s".', opts.hotPixelAction);
    end

    dark = struct();
    dark.method = 'boundary per detector channel';
    dark.boundaryWidthPx = width;
    dark.darkPerChannel = boundaryMedian;
    dark.noiseMadPerChannel = boundaryMad;
    dark.globalMedian = globalMedian;
    dark.globalMad = globalMad;
    dark.neighborMedian = neighborMedian;
    dark.neighborMad = neighborMad;
    dark.hotScoreGlobal = hotScoreGlobal;
    dark.hotScoreNeighbor = hotScoreNeighbor;
    dark.hotChannels = hotChannels;
    dark.hotPixelAction = opts.hotPixelAction;
    dark.table = detectorDarkTable(sim, dark);
end

function [dark, corrected] = estimateScalarVolumeBackground(vol, opts)
    vol = double(vol);
    [ny, nx, nz] = size(vol);
    width = max(1, min([opts.boundaryWidthPx, floor(ny/4), floor(nx/4)]));
    borderMask = false(ny, nx);
    borderMask(1:width, :) = true;
    borderMask(end-width+1:end, :) = true;
    borderMask(:, 1:width) = true;
    borderMask(:, end-width+1:end) = true;

    vals = [];
    for z = 1:nz
        frame = vol(:,:,z);
        vals = [vals; frame(borderMask)]; %#ok<AGROW>
    end
    vals = vals(isfinite(vals));
    bg = max(median(vals), 0);
    corrected = max(vol - bg, 0);

    dark = struct();
    dark.method = 'scalar volume boundary';
    dark.boundaryWidthPx = width;
    dark.background = bg;
    dark.backgroundMad = robustMad(vals);
    dark.hotChannels = [];
end

function dark = estimateTcspcDetectorDarkFromFolder(scanFolder, sim, opts)
    dark = [];
    if ~opts.estimateTcspcDark || exist(scanFolder, 'dir') ~= 7
        return;
    end

    files = dir(fullfile(scanFolder, 'TCSPCCurve*.pqdat'));
    nCh = size(sim.detXY, 1);
    if numel(files) < nCh
        return;
    end

    fileIndex = zeros(numel(files), 1);
    for k = 1:numel(files)
        fileIndex(k) = tcspcCurveFileIndex(files(k).name);
    end
    [~, order] = sort(fileIndex);
    files = files(order);

    nPlane = floor(numel(files) / nCh);
    if nPlane < 1
        return;
    end
    files = files(1:(nPlane*nCh));

    counts = nan(nCh, nPlane);
    for p = 1:nPlane
        for k = 1:nCh
            idx = (p-1)*nCh + k;
            fileName = fullfile(files(idx).folder, files(idx).name);
            try
                counts(k,p) = readPqdatCurveYSum(fileName);
            catch err
                warning('estimateExperimentalWavefrontFromPSF:BadTcspcCurve', ...
                    'Could not read %s: %s', fileName, err.message);
            end
        end
    end

    planeSignal = sum(counts, 1, 'omitnan');
    nDarkPlane = max(3, ceil(opts.tcspcDarkPlaneFraction * nPlane));
    nDarkPlane = min(nDarkPlane, nPlane);
    [~, planeOrder] = sort(planeSignal);
    darkPlanes = sort(planeOrder(1:nDarkPlane));

    darkPerChannel = nan(nCh, 1);
    noiseMadPerChannel = nan(nCh, 1);
    for k = 1:nCh
        vals = counts(k, darkPlanes);
        vals = vals(isfinite(vals));
        if isempty(vals)
            darkPerChannel(k) = NaN;
            noiseMadPerChannel(k) = NaN;
        else
            darkPerChannel(k) = median(vals);
            noiseMadPerChannel(k) = robustMad(vals);
        end
    end

    neighbors = detectorNeighborList(sim.detXY, sim.detPitch);
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

    finiteDark = darkPerChannel(isfinite(darkPerChannel));
    globalMedian = median(finiteDark);
    globalMad = max(robustMad(finiteDark), eps);
    hotScoreGlobal = (darkPerChannel - globalMedian) / globalMad;
    hotScoreNeighbor = (darkPerChannel - neighborMedian) ./ max(neighborMad, eps);
    hotMask = hotScoreGlobal > opts.hotThresholdMAD | ...
              hotScoreNeighbor > opts.hotNeighborThresholdMAD;

    dark = struct();
    dark.method = 'TCSPC curve integral on low-signal z planes';
    dark.darkPerChannel = darkPerChannel;
    dark.noiseMadPerChannel = noiseMadPerChannel;
    dark.globalMedian = globalMedian;
    dark.globalMad = globalMad;
    dark.neighborMedian = neighborMedian;
    dark.neighborMad = neighborMad;
    dark.hotScoreGlobal = hotScoreGlobal;
    dark.hotScoreNeighbor = hotScoreNeighbor;
    dark.hotChannels = find(hotMask);
    dark.darkPlaneIndex = darkPlanes(:);
    dark.planeSignal = planeSignal(:);
    dark.curveCounts = counts;
    dark.table = detectorDarkTable(sim, dark);
end

function idx = tcspcCurveFileIndex(name)
    if strcmpi(name, 'TCSPCCurve.pqdat')
        idx = 0;
        return;
    end
    tok = regexp(name, '^TCSPCCurve_(\d+)\.pqdat$', 'tokens', 'once');
    if isempty(tok)
        idx = inf;
    else
        idx = str2double(tok{1});
    end
end

function ySum = readPqdatCurveYSum(fileName)
    fid = fopen(fileName, 'r', 'ieee-le');
    if fid < 0
        error('Could not open %s.', fileName);
    end
    cleanup = onCleanup(@() fclose(fid));

    magic = fread(fid, 8, '*char')';
    fread(fid, 8, '*char');
    if numel(magic) < 6 || ~strcmp(char(magic(1:6)), 'PQDATA')
        error('%s is not a PQDATA file.', fileName);
    end

    ySum = NaN;
    while ~feof(fid)
        identRaw = fread(fid, 32, '*char')';
        if numel(identRaw) < 32
            break;
        end
        ident = trimNulls(char(identRaw));
        fread(fid, 1, 'int32');
        tagType = fread(fid, 1, 'uint32=>uint32');
        valueBytes = fread(fid, 8, 'uint8=>uint8')';
        if numel(valueBytes) < 8
            break;
        end
        nBytes = double(typecast(uint8(valueBytes), 'int64'));

        switch tagType
            case uint32(hex2dec('2001FFFF')) % float8 array
                if strcmp(ident, 'LSDCurveY')
                    y = fread(fid, nBytes/8, 'double=>double');
                    ySum = sum(y(isfinite(y)));
                else
                    fseek(fid, nBytes, 'cof');
                end
            case {uint32(hex2dec('4001FFFF')), uint32(hex2dec('4002FFFF')), ...
                    uint32(hex2dec('FFFFFFFF'))}
                fseek(fid, nBytes, 'cof');
            otherwise
                if strcmp(ident, 'Header_End')
                    break;
                end
        end

        if strcmp(ident, 'Header_End')
            break;
        end
    end

    if ~isfinite(ySum)
        error('No LSDCurveY array found in %s.', fileName);
    end
end

function T = detectorDarkTable(sim, dark)
    channel = (1:numel(dark.darkPerChannel)).';
    detectorXUm = sim.detXY(:,1);
    detectorYUm = sim.detXY(:,2);
    darkEstimateCounts = dark.darkPerChannel(:);
    darkNoiseMadCounts = dark.noiseMadPerChannel(:);
    neighborMedianCounts = dark.neighborMedian(:);
    hotScoreGlobal = dark.hotScoreGlobal(:);
    hotScoreNeighbor = dark.hotScoreNeighbor(:);
    isHot = false(numel(channel), 1);
    isHot(dark.hotChannels) = true;
    T = table(channel, detectorXUm, detectorYUm, darkEstimateCounts, ...
        darkNoiseMadCounts, neighborMedianCounts, hotScoreGlobal, ...
        hotScoreNeighbor, isHot);
end

function neighbors = detectorNeighborList(detXY, detPitch)
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

function fit = fitSummedVolumeAberration(data, planeZ, sim, opts)
    sim = prepareSummedVolumeSim(sim, planeZ, opts);
    dataN = normalizeScalarVolume(data);
    paramNames = buildExperimentalParamNames(opts.fitModes, opts.fitXY, opts.fitZ);
    p = initialExperimentalParams(sim, opts);
    step = experimentalFiniteDifferenceSteps(numel(opts.fitModes), sim, opts);
    reg = experimentalRegularization(numel(opts.fitModes), opts);
    maxStep = experimentalMaxStep(numel(opts.fitModes), opts);

    history = zeros(opts.maxIter, 2);
    for it = 1:opts.maxIter
        m0 = summedModelFromParams(sim, opts.fitModes, planeZ, p, opts.fitXY, opts.fitZ);
        r = dataN(:) - m0(:);
        J = zeros(numel(m0), numel(p));

        for q = 1:numel(p)
            pp = p;
            pm = p;
            pp(q) = pp(q) + step(q);
            pm(q) = pm(q) - step(q);
            mp = summedModelFromParams(sim, opts.fitModes, planeZ, pp, opts.fitXY, opts.fitZ);
            mm = summedModelFromParams(sim, opts.fitModes, planeZ, pm, opts.fitXY, opts.fitZ);
            J(:,q) = (mp(:) - mm(:)) / (2*step(q));
        end

        H = J.'*J + diag(reg);
        g = J.'*r;
        if rcond(H) < 1e-12
            delta = pinv(H) * g;
        else
            delta = H \ g;
        end
        delta = max(min(delta(:).', maxStep), -maxStep);
        p = p + delta;

        history(it,1) = norm(r);
        history(it,2) = norm(delta);
        if opts.verbose
            fprintf('[summedVolumeFit] iter %02d  residual %.4e  step %.4e\n', ...
                it, history(it,1), history(it,2));
        end
        if norm(delta) < opts.tolStep
            history = history(1:it,:);
            break;
        end
    end

    fitVolume = summedModelFromParams(sim, opts.fitModes, planeZ, p, opts.fitXY, opts.fitZ);
    [estCoeffs, estXY, estZOffset] = unpackExperimentalParams(sim, opts.fitModes, p, opts.fitXY, opts.fitZ);

    fit = struct();
    fit.dataN = dataN;
    fit.background = 0;
    fit.planeZ = planeZ;
    fit.fitModes = opts.fitModes;
    fit.paramNames = paramNames;
    fit.paramVector = p;
    fit.estCoeffs = estCoeffs;
    fit.estCoeffVector = coeffsToFullVector(sim, estCoeffs);
    fit.estXY = estXY;
    fit.estZOffset = estZOffset;
    fit.fitVolume = fitVolume;
    fit.fitStack = [];
    fit.residual = dataN - fitVolume;
    fit.residualNorm = norm(fit.residual(:));
    fit.history = history;
    fit.phase = zernikePhaseMap(sim, estCoeffs, sim.lamEm);
    fit.sim = sim;
    fit.options = opts;
end

function sim = prepareSummedVolumeSim(sim, planeZ, opts)
    dzTarget = opts.modelDz;
    zPadding = opts.modelZPadding;
    zMin = min([planeZ(:); 0]) - zPadding;
    zMax = max([planeZ(:); 0]) + zPadding;
    nZ = ceil((zMax - zMin) / dzTarget) + 1;
    if mod(nZ, 2) == 0
        nZ = nZ + 1;
    end
    sim.nzRange = (nZ - 1) * dzTarget;
    sim.nz = nZ;
    sim.z = linspace(-sim.nzRange/2, sim.nzRange/2, sim.nz);
    sim.obj = beadObject3D(sim);
end

function volN = normalizeScalarVolume(vol)
    vol = max(double(vol), 0);
    s = sum(vol(:));
    if s > 0
        volN = vol / s;
    else
        volN = vol;
    end
end

function vol = summedModelFromParams(sim, fitModes, planeZ, p, fitXY, fitZ)
    [coeffs, xy, zOffset] = unpackExperimentalParams(sim, fitModes, p, fitXY, fitZ);
    stack = normalizedStackExplicitDetectorZPlanes(sim, coeffs, planeZ, xy(1), xy(2), zOffset);
    summed = sum(stack, 3);
    vol = reshape(summed, size(summed,1), size(summed,2), size(summed,4));
    vol = normalizeScalarVolume(vol);
end

function names = buildExperimentalParamNames(fitModes, fitXY, fitZ)
    names = fitModes;
    if fitXY
        names = [names {'x_shift','y_shift'}];
    end
    if fitZ
        names = [names {'z_offset'}];
    end
end

function p = initialExperimentalParams(sim, opts)
    coeffs0 = struct();
    if isfield(opts, 'initialCoeffs') && ~isempty(opts.initialCoeffs)
        coeffs0 = opts.initialCoeffs;
    end
    coeffs0 = coeffStruct(sim, coeffs0);
    p = zeros(1, numel(opts.fitModes));
    for k = 1:numel(opts.fitModes)
        if isfield(coeffs0, opts.fitModes{k})
            p(k) = coeffs0.(opts.fitModes{k});
        end
    end
    if opts.fitXY
        xy0 = [0 0];
        if isPositiveVector2(opts.initialXY)
            xy0 = double(opts.initialXY(:).');
        end
        p = [p xy0(1) xy0(2)];
    end
    if opts.fitZ
        p = [p 0];
    end
end

function step = experimentalFiniteDifferenceSteps(nModes, sim, opts)
    step = opts.fdCoeff * ones(1, nModes);
    if opts.fitXY
        step = [step sim.dx/4 sim.dx/4];
    end
    if opts.fitZ
        step = [step opts.fdZ];
    end
end

function reg = experimentalRegularization(nModes, opts)
    reg = opts.regCoeff * ones(1, nModes);
    if opts.fitXY
        reg = [reg opts.regXY opts.regXY];
    end
    if opts.fitZ
        reg = [reg opts.regZ];
    end
end

function maxStep = experimentalMaxStep(nModes, opts)
    maxStep = opts.maxCoeffStep * ones(1, nModes);
    if opts.fitXY
        maxStep = [maxStep opts.maxXYStep opts.maxXYStep];
    end
    if opts.fitZ
        maxStep = [maxStep opts.maxZStep];
    end
end

function [coeffs, xy, zOffset] = unpackExperimentalParams(sim, fitModes, p, fitXY, fitZ)
    fullVec = zeros(1, numel(sim.modeOrder));
    for k = 1:numel(fitModes)
        idx = find(strcmp(sim.modeOrder, fitModes{k}), 1, 'first');
        if isempty(idx)
            error('estimateExperimentalWavefrontFromPSF:UnknownMode', ...
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

function v = coeffsToFullVector(sim, coeffs)
    coeffs = coeffStruct(sim, coeffs);
    v = zeros(1, numel(sim.modeOrder));
    for k = 1:numel(sim.modeOrder)
        name = sim.modeOrder{k};
        if isfield(coeffs, name)
            v(k) = coeffs.(name);
        end
    end
end

function compare = detectorFitComparisonVolumes(data4, fitStack)
    measuredSum = sum(data4, 3);
    fittedSum = sum(fitStack, 3);
    measured = reshape(measuredSum, size(measuredSum,1), size(measuredSum,2), size(measuredSum,4));
    fitted = reshape(fittedSum, size(fittedSum,1), size(fittedSum,2), size(fittedSum,4));
    compare = comparisonStruct(measured, fitted);
end

function compare = scalarFitComparisonVolumes(measured, fitted)
    compare = comparisonStruct(measured, fitted);
end

function compare = comparisonStruct(measured, fitted)
    measured = normalizeScalarVolume(measured);
    fitted = normalizeScalarVolume(fitted);
    compare = struct();
    compare.measuredVolume = measured;
    compare.fittedVolume = fitted;
    compare.residualVolume = measured - fitted;
end

function outDir = resolveOutputDir(dataPath, meta, opts)
    if ~opts.writeOutputs
        outDir = '';
        return;
    end

    if ~isempty(opts.outputRoot)
        outRoot = opts.outputRoot;
    elseif ~isempty(meta.sourceMatFile)
        outRoot = fullfile(fileparts(fileparts(meta.sourceMatFile)), 'aberration_fits');
    elseif exist(dataPath, 'dir') == 7
        [scanParent, ~] = fileparts(stripTrailingFilesep(dataPath));
        [dataParent, datasetName] = fileparts(scanParent);
        outRoot = fullfile(dataParent, 'PSF_batch_outputs', ...
            sanitizeFileName(datasetName), 'aberration_fits');
    else
        outRoot = fullfile(pwd, 'aberration_fits');
    end

    outDir = fullfile(outRoot, scanOutputStem(meta.scanName));
end

function writeExperimentalFitOutputs(result)
    outDir = result.outputDir;
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    save(fullfile(outDir, 'experimental_wavefront_fit.mat'), 'result', '-v7.3');
    coeffTable = coefficientTable(result.fit);
    writetable(coeffTable, fullfile(outDir, 'estimated_wavefront_coefficients.csv'));

    if isfield(result.dark, 'table')
        writetable(result.dark.table, fullfile(outDir, 'dark_hot_pixel_estimate.csv'));
        writeDarkDetectorMap(result, fullfile(outDir, 'dark_detector_map.png'));
    else
        scalarBg = table(result.dark.background, result.dark.backgroundMad, ...
            'VariableNames', {'boundaryBackgroundCounts','boundaryMadCounts'});
        writetable(scalarBg, fullfile(outDir, 'scalar_background_estimate.csv'));
        if isfield(result.dark, 'tcspcDetector') && ~isempty(result.dark.tcspcDetector)
            writetable(result.dark.tcspcDetector.table, ...
                fullfile(outDir, 'tcspc_dark_hot_pixel_estimate.csv'));
            writeDetectorMetricMap(result.fit.sim, result.dark.tcspcDetector.darkPerChannel, ...
                'TCSPC low-signal detector counts', ...
                fullfile(outDir, 'tcspc_dark_detector_map.png'));
        end
    end

    writeFitSummaryFigure(result, fullfile(outDir, 'experimental_wavefront_fit_summary.png'));
end

function T = coefficientTable(fit)
    mode = fit.sim.modeOrder(:);
    coeffWavesRMS = fit.estCoeffVector(:);
    coeffNmRMS = coeffWavesRMS * fit.sim.lamRef * 1000;
    T = table(mode, coeffWavesRMS, coeffNmRMS);

    pName = fit.paramNames(:);
    pValue = fit.paramVector(:);
    P = table(pName, pValue, 'VariableNames', {'fitParameter','fitValue'});
    T = [T; table(repmat({''}, height(P), 1), nan(height(P),1), nan(height(P),1), ...
        'VariableNames', T.Properties.VariableNames)];
    T.fitParameter = [repmat({''}, numel(mode), 1); P.fitParameter];
    T.fitValue = [nan(numel(mode), 1); P.fitValue];
end

function writeFitSummaryFigure(result, outFile)
    fit = result.fit;
    sim = fit.sim;
    measured = result.compare.measuredVolume;
    fitted = result.compare.fittedVolume;
    residual = result.compare.residualVolume;

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1350 900]);
    tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    plotVolumeProjection(nexttile(tl), measured, sim.x, sim.y, result.planeZ, 'xy', 'Measured XY');
    plotVolumeProjection(nexttile(tl), fitted, sim.x, sim.y, result.planeZ, 'xy', 'Fitted XY');
    plotVolumeProjection(nexttile(tl), abs(residual), sim.x, sim.y, result.planeZ, 'xy', '|Residual| XY');

    plotVolumeProjection(nexttile(tl), measured, sim.x, sim.y, result.planeZ, 'xz', 'Measured XZ');
    plotVolumeProjection(nexttile(tl), fitted, sim.x, sim.y, result.planeZ, 'xz', 'Fitted XZ');
    plotVolumeProjection(nexttile(tl), abs(residual), sim.x, sim.y, result.planeZ, 'xz', '|Residual| XZ');

    plotVolumeProjection(nexttile(tl), measured, sim.x, sim.y, result.planeZ, 'yz', 'Measured YZ');
    plotVolumeProjection(nexttile(tl), fitted, sim.x, sim.y, result.planeZ, 'yz', 'Fitted YZ');

    ax = nexttile(tl);
    coeff = fit.estCoeffVector(:);
    bar(ax, coeff);
    set(ax, 'XTick', 1:numel(coeff), 'XTickLabel', fit.sim.modeOrder, 'XTickLabelRotation', 35);
    ylabel(ax, 'waves RMS');
    title(ax, sprintf('Estimated modes (%s)', result.fitType), 'Interpreter', 'none');
    grid(ax, 'on');

    title(tl, sprintf('Experimental wavefront fit: %s', result.inputMeta.scanName), ...
        'Interpreter', 'none');
    saveFigurePng(fig, outFile, 180);
    close(fig);
end

function plotVolumeProjection(ax, vol, x, y, z, plane, titleText)
    switch lower(plane)
        case 'xy'
            img = max(vol, [], 3);
            imagesc(ax, x, y, img);
            xlabel(ax, 'x (um)');
            ylabel(ax, 'y (um)');
            axis(ax, 'image');
        case 'xz'
            img = squeeze(max(vol, [], 1)).';
            imagesc(ax, x, z, img);
            xlabel(ax, 'x (um)');
            ylabel(ax, 'z (um)');
            axis(ax, 'image');
        case 'yz'
            img = squeeze(max(vol, [], 2)).';
            imagesc(ax, y, z, img);
            xlabel(ax, 'y (um)');
            ylabel(ax, 'z (um)');
            axis(ax, 'image');
    end
    set(ax, 'YDir', 'normal');
    title(ax, titleText);
    colormap(ax, hot);
    colorbar(ax);
end

function writeDarkDetectorMap(result, outFile)
    writeDetectorMetricMap(result.fit.sim, result.dark.darkPerChannel, ...
        'Estimated detector dark count from boundary pixels', outFile);
end

function writeDetectorMetricMap(sim, values, titleText, outFile)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [120 120 700 520]);
    ax = axes(fig);
    plotDetectorHexMap(sim.detXY, values, 'Parent', ax, 'CellScale', 1.01);
    hold(ax, 'on');
    for k = 1:size(sim.detXY, 1)
        text(ax, sim.detXY(k,1), sim.detXY(k,2), sprintf('%d', k), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold');
    end
    hold(ax, 'off');
    title(ax, titleText);
    colormap(ax, hot);
    colorbar(ax);
    saveFigurePng(fig, outFile, 180);
    close(fig);
end

function printCoefficientSummary(fit)
    coeff = fit.estCoeffVector(:);
    for k = 1:numel(coeff)
        fprintf('  %-10s %+8.4f waves RMS (%+7.2f nm RMS)\n', ...
            fit.sim.modeOrder{k}, coeff(k), coeff(k) * fit.sim.lamRef * 1000);
    end
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

function scanName = scanNameFromPath(dataPath)
    p = stripTrailingFilesep(char(dataPath));
    [~, scanName, ext] = fileparts(p);
    if strcmpi(ext, '.mat')
        scanName = regexprep(scanName, '_volume_raw$', '');
    end
end

function s = stripTrailingFilesep(s)
    while ~isempty(s) && (s(end) == filesep || s(end) == '/' || s(end) == '\')
        s(end) = [];
    end
end

function stem = scanOutputStem(name)
    stem = sanitizeFileName(name);
    if isempty(stem)
        stem = 'scan';
    end
    if ~isletter(stem(1))
        stem = ['x' stem];
    end
end

function clean = sanitizeFileName(name)
    clean = regexprep(char(name), '[^A-Za-z0-9]+', '_');
    clean = regexprep(clean, '^_+|_+$', '');
    if isempty(clean)
        clean = 'unnamed';
    end
end

function m = robustMad(vals)
    vals = vals(:);
    vals = vals(isfinite(vals));
    if isempty(vals)
        m = 0;
        return;
    end
    med = median(vals);
    m = 1.4826 * median(abs(vals - med));
end

function txt = trimNulls(txt)
    txt = txt(txt ~= char(0));
    txt = strtrim(txt);
end
