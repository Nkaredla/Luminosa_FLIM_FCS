function results = compareMeasuredAndFitted3DPSFs(varargin)
%COMPAREMEASUREDANDFITTED3DPSFS Compare measured and fitted 23-channel PSFs.
%
%   results = compareMeasuredAndFitted3DPSFs()
%
%   This workflow:
%     1. Loads the coefficients from the center-point wavefront fit.
%     2. Reconstructs the full detector-resolved experimental z stack from
%        the PTU files listed in the batch alignment CSV.
%     3. Applies the saved subpixel shifts to every detector channel.
%     4. Computes the fitted optical effective PSF h_exc .* eta_k.
%     5. Convolves each calculated PSF with the finite bead for a fair
%        comparison against the measured bead response.
%     6. Exports XY, XZ, and YZ projections for every detector pixel.
%
%   The measured volume is cached after the first run because reading all
%   PTU planes is substantially slower than regenerating the figures.

    opts = parseOptions(varargin{:});
    addRequiredPaths(opts);

    fitResult = loadFitResult(opts.fitResultFile);
    opts = inheritFitChannelOptions(opts, fitResult);
    alignmentCsv = resolveAlignmentCsv(opts, fitResult);
    alignmentTable = readtable(alignmentCsv);
    validateAlignmentTable(alignmentTable);

    outputDir = resolveOutputDir(opts);
    ensureFolder(outputDir);

    [measuredVolume, dataInfo] = loadMeasuredDetectorVolume( ...
        alignmentTable, fitResult, opts, outputDir);

    [modelPsf, modelBead, modelInfo] = calculateFittedVolumes( ...
        fitResult, measuredVolume, dataInfo, opts);

    [measuredN, modelBeadN, modelPsfN] = normalizeComparisonVolumes( ...
        measuredVolume, modelBead, modelPsf);
    metrics = calculateMetrics(measuredN, modelBeadN, modelPsfN, ...
        measuredVolume, modelBead, dataInfo, fitResult, opts.supportThreshold);

    results = struct();
    results.fitResultFile = opts.fitResultFile;
    results.alignmentCsv = alignmentCsv;
    results.outputDir = outputDir;
    results.options = opts;
    results.fitResult = fitResult;
    results.dataInfo = dataInfo;
    results.modelInfo = modelInfo;
    results.measuredVolume = measuredVolume;
    results.modelBeadVolume = modelBead;
    results.modelEffectivePsf = modelPsf;
    results.metrics = metrics;

    writeComparisonOutputs(results, measuredN, modelBeadN, modelPsfN);

    if opts.verbose
        printSummary(results);
    end
end

function opts = inheritFitChannelOptions(opts, fitResult)
    if ~isfield(fitResult, 'options') || ~isstruct(fitResult.options)
        return;
    end
    if isempty(opts.channelOrder) && isfield(fitResult.options, 'channelOrder') && ...
            ~isempty(fitResult.options.channelOrder)
        opts.channelOrder = double(fitResult.options.channelOrder(:)).';
    end
    if isequal(opts.channelIDs, 9:31) && isfield(fitResult.options, 'channelIDs') && ...
            ~isempty(fitResult.options.channelIDs)
        opts.channelIDs = double(fitResult.options.channelIDs(:)).';
    end
end

function opts = parseOptions(varargin)
    p = inputParser;
    p.FunctionName = 'compareMeasuredAndFitted3DPSFs';

    defaultFitDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'output_matlab', 'center_point_ism_wavefront', ...
        'x0_4uW_0_19collar_80mmlens_20260515_155744_volume_raw_focusAnchored_500nm');
    addParameter(p, 'fitResultFile', ...
        fullfile(defaultFitDir, 'center_point_ism_wavefront_from_stack.mat'));
    addParameter(p, 'alignmentCsv', ...
        'D:\Luminosa\Data\PSF_batch_outputs\ISM_Aberation2_73\xz_yz_plots\x0_4uW_0_19collar_80mmlens_20260515_155744_frame_alignment.csv');
    addParameter(p, 'darkFile', 'D:\Luminosa\Data\ISMdark_counts.ptu');
    addParameter(p, 'darkScale', 1);
    addParameter(p, 'subtractBoundary', true);
    addParameter(p, 'boundaryWidthPx', 3);
    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'channelOrder', []);

    addParameter(p, 'ptuReaderFolder', defaultPtuReaderFolder());
    addParameter(p, 'ptuPhotonsPerChunk', 1e6);
    addParameter(p, 'applySavedAlignment', false);
    addParameter(p, 'alignmentShiftMinPx', 4);
    addParameter(p, 'recenterMeasuredVolume', false);
    addParameter(p, 'registerAxialShift', true);
    addParameter(p, 'maxAxialRegistrationUm', 0.30);
    addParameter(p, 'axialRegistrationStepUm', 0.01);

    addParameter(p, 'modelZStepUm', []);
    addParameter(p, 'modelZPaddingUm', 0.20);
    addParameter(p, 'detectorSubsamples', []);
    addParameter(p, 'projectionMode', 'max');
    addParameter(p, 'supportThreshold', 0.01);

    addParameter(p, 'outputDir', '');
    addParameter(p, 'reuseMeasuredCache', false);
    addParameter(p, 'writeVolumeMat', true);
    addParameter(p, 'writePerDetectorFigures', true);
    addParameter(p, 'writeArrayFigures', true);
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;
    opts.fitResultFile = char(opts.fitResultFile);
    opts.alignmentCsv = char(opts.alignmentCsv);
    opts.darkFile = char(opts.darkFile);
    opts.ptuReaderFolder = char(opts.ptuReaderFolder);
    opts.outputDir = char(opts.outputDir);
    opts.projectionMode = lower(char(opts.projectionMode));
    opts.channelIDs = double(opts.channelIDs(:)).';
    opts.channelOrder = double(opts.channelOrder(:)).';
    opts.alignmentShiftMinPx = double(opts.alignmentShiftMinPx);
    if ~isscalar(opts.alignmentShiftMinPx) || ...
            ~isfinite(opts.alignmentShiftMinPx) || opts.alignmentShiftMinPx < 0
        error('compareMeasuredAndFitted3DPSFs:BadAlignmentThreshold', ...
            'alignmentShiftMinPx must be a finite nonnegative scalar.');
    end

    if ~ismember(opts.projectionMode, {'max','sum'})
        error('compareMeasuredAndFitted3DPSFs:BadProjectionMode', ...
            'projectionMode must be ''max'' or ''sum''.');
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

function fitResult = loadFitResult(fileName)
    if exist(fileName, 'file') ~= 2
        error('compareMeasuredAndFitted3DPSFs:MissingFitResult', ...
            'Fit result file was not found: %s', fileName);
    end
    S = load(fileName);
    if isfield(S, 'result')
        fitResult = S.result;
    elseif isfield(S, 'results')
        fitResult = S.results;
    else
        names = fieldnames(S);
        fitResult = [];
        for k = 1:numel(names)
            value = S.(names{k});
            if isstruct(value) && isfield(value, 'fit')
                fitResult = value;
                break;
            end
        end
    end
    if isempty(fitResult) || ~isfield(fitResult, 'fit') || ...
            ~isfield(fitResult.fit, 'estCoeffs')
        error('compareMeasuredAndFitted3DPSFs:BadFitResult', ...
            'The MAT file does not contain a center-point wavefront result.');
    end
end

function csvFile = resolveAlignmentCsv(opts, fitResult)
    csvFile = opts.alignmentCsv;
    if ~isempty(csvFile) && exist(csvFile, 'file') == 2
        return;
    end

    scanName = '';
    if isfield(fitResult, 'stackSelection') && ...
            isfield(fitResult.stackSelection, 'focusFile') && ...
            ~isempty(fitResult.stackSelection.focusFile)
        scanFolder = fileparts(char(fitResult.stackSelection.focusFile));
        [~, scanName] = fileparts(scanFolder);
    end
    if isempty(scanName)
        error('compareMeasuredAndFitted3DPSFs:MissingAlignmentCsv', ...
            'Pass alignmentCsv explicitly; it could not be inferred from the fit result.');
    end

    searchRoot = 'D:\Luminosa\Data\PSF_batch_outputs';
    files = dir(fullfile(searchRoot, '**', '*frame_alignment.csv'));
    scanToken = regexprep(scanName, '[^A-Za-z0-9]+', '_');
    keep = contains({files.name}, scanToken, 'IgnoreCase', true);
    files = files(keep);
    if isempty(files)
        error('compareMeasuredAndFitted3DPSFs:MissingAlignmentCsv', ...
            'No frame-alignment CSV was found for scan %s.', scanName);
    end
    csvFile = fullfile(files(1).folder, files(1).name);
end

function validateAlignmentTable(T)
    required = {'source_file','z_um','total_shift_x','total_shift_y'};
    if ~all(ismember(required, T.Properties.VariableNames))
        error('compareMeasuredAndFitted3DPSFs:BadAlignmentCsv', ...
            'Alignment CSV must contain source_file, z_um, total_shift_x, and total_shift_y.');
    end
end

function outDir = resolveOutputDir(opts)
    if ~isempty(opts.outputDir)
        outDir = opts.outputDir;
        return;
    end
    outDir = fullfile(fileparts(opts.fitResultFile), 'psf3d_detector_comparison');
end

function [volume, info] = loadMeasuredDetectorVolume(T, fitResult, opts, outDir)
    cacheFile = fullfile(outDir, 'measured_detector_volume_cache.mat');
    if opts.reuseMeasuredCache && exist(cacheFile, 'file') == 2
        S = load(cacheFile);
        if isfield(S, 'measuredVolume') && isfield(S, 'dataInfo')
            volume = S.measuredVolume;
            info = S.dataInfo;
            if opts.verbose
                fprintf('[compareMeasuredAndFitted3DPSFs] Loaded measured cache: %s\n', cacheFile);
            end
            return;
        end
    end

    nPlane = height(T);
    volume = [];
    sourceChannelIDs = [];
    darkImage = [];
    darkVector = [];

    for iz = 1:nPlane
        fileName = tableText(T.source_file, iz);
        frameIndex = tableFrame(T, iz);
        [plane, channelIDs, head] = readPtuDetectorPlane(fileName, frameIndex, opts);

        if isempty(volume)
            volume = zeros(size(plane,1), size(plane,2), nPlane, size(plane,3));
            sourceChannelIDs = channelIDs(:).';
            [darkImage, darkVector] = loadDarkCorrection( ...
                opts, size(plane), sourceChannelIDs);
            firstHead = head;
        elseif ~isequal(size(plane), [size(volume,1), size(volume,2), size(volume,4)])
            error('compareMeasuredAndFitted3DPSFs:PlaneSizeMismatch', ...
                'Detector plane %d has a different size.', iz);
        end

        plane = applyCorrections(plane, darkImage, darkVector, opts);
        if opts.applySavedAlignment
            shiftX = double(T.total_shift_x(iz));
            shiftY = double(T.total_shift_y(iz));
            if hypot(shiftX, shiftY) >= opts.alignmentShiftMinPx
                for ch = 1:size(plane,3)
                    plane(:,:,ch) = max(fourierShift2D( ...
                        plane(:,:,ch), shiftX, shiftY), 0);
                end
            end
        end
        volume(:,:,iz,:) = permute(plane, [1 2 4 3]);

        if opts.verbose
            fprintf('[compareMeasuredAndFitted3DPSFs] PTU plane %d/%d\n', iz, nPlane);
        end
    end

    focusIndex = resolveFocusIndex(T, fitResult);
    focusZUm = double(T.z_um(focusIndex));
    zUm = double(T.z_um(:)).';
    relativeZUm = zUm - focusZUm;

    recenterShiftXY = [0 0];
    if opts.recenterMeasuredVolume
        focusImage = sum(squeeze(volume(:,:,focusIndex,:)), 3);
        centerXY = intensityCentroid(focusImage);
        targetXY = [(size(volume,2)+1)/2, (size(volume,1)+1)/2];
        recenterShiftXY = targetXY - centerXY;
        for iz = 1:nPlane
            for ch = 1:size(volume,4)
                volume(:,:,iz,ch) = max(fourierShift2D( ...
                    volume(:,:,iz,ch), recenterShiftXY(1), recenterShiftXY(2)), 0);
            end
        end
    end

    info = struct();
    info.channelIDs = sourceChannelIDs;
    info.zUm = zUm;
    info.relativeZUm = relativeZUm;
    info.focusIndex = focusIndex;
    info.focusZUm = focusZUm;
    info.recenterShiftXY = recenterShiftXY;
    info.xyPixelSizeUm = resolvePixelSize(fitResult, firstHead);
    info.xUm = centeredAxis(size(volume,2), info.xyPixelSizeUm);
    info.yUm = centeredAxis(size(volume,1), info.xyPixelSizeUm);
    info.alignmentApplied = opts.applySavedAlignment;
    info.alignmentShiftMinPx = opts.alignmentShiftMinPx;
    if ismember('total_shift_x', T.Properties.VariableNames) && ...
            ismember('total_shift_y', T.Properties.VariableNames)
        savedShiftMagnitudePx = hypot(double(T.total_shift_x(:)), ...
            double(T.total_shift_y(:)));
        info.savedShiftMagnitudePx = savedShiftMagnitudePx(:).';
        info.nSavedShiftsAboveThreshold = nnz(savedShiftMagnitudePx >= opts.alignmentShiftMinPx);
        info.nSavedShiftsBelowThreshold = nnz(savedShiftMagnitudePx < opts.alignmentShiftMinPx);
    end
    info.darkFile = opts.darkFile;

    measuredVolume = volume; %#ok<NASGU>
    dataInfo = info; %#ok<NASGU>
    save(cacheFile, 'measuredVolume', 'dataInfo', '-v7.3');
end

function [plane, channelIDs, head] = readPtuDetectorPlane(fileName, frameIndex, opts)
    if exist(fileName, 'file') ~= 2
        error('compareMeasuredAndFitted3DPSFs:MissingPtu', ...
            'PTU file was not found: %s', fileName);
    end
    hasFast = exist('PTU_MultiFrameScanReadFast', 'file') == 2;
    hasSlow = exist('PTU_MultiFrameScanRead', 'file') == 2;
    if ~hasFast && ~hasSlow
        error('compareMeasuredAndFitted3DPSFs:MissingReader', ...
            'PTU_MultiFrameScanReadFast/PTU_MultiFrameScanRead is not on the path.');
    end

    waitbarCleanup = suppressPtuWaitbars(); %#ok<NASGU>
    ptuOut = [];
    if hasFast
        try
            ptuOut = PTU_MultiFrameScanReadFast( ...
                fileName, opts.ptuPhotonsPerChunk, false, false);
        catch
            ptuOut = [];
        end
    end
    if isempty(ptuOut)
        ptuOut = PTU_MultiFrameScanRead(fileName, opts.ptuPhotonsPerChunk);
    end

    if ~isempty(frameIndex) && isfield(ptuOut, 'tag') && ~isempty(ptuOut.tag)
        frameIndex = min(max(round(frameIndex), 1), size(ptuOut.tag,4));
        raw = double(ptuOut.tag(:,:,:,frameIndex));
    elseif isfield(ptuOut, 'tags') && ~isempty(ptuOut.tags)
        raw = double(ptuOut.tags);
    elseif isfield(ptuOut, 'tag') && ~isempty(ptuOut.tag)
        raw = sum(double(ptuOut.tag), 4);
    else
        error('compareMeasuredAndFitted3DPSFs:NoDetectorImage', ...
            'No detector image was found in %s.', fileName);
    end
    if isfield(ptuOut, 'dind')
        channelIDs = double(ptuOut.dind(:));
    else
        channelIDs = (1:size(raw,3)).';
    end
    plane = permute(raw, [2 1 3]);
    [plane, channelIDs] = selectChannels(plane, channelIDs, opts);

    head = struct();
    if isfield(ptuOut, 'head')
        head = ptuOut.head;
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

function [plane, channelIDs] = selectChannels(plane, channelIDs, opts)
    requested = opts.channelIDs(:);
    if ~isempty(channelIDs) && all(ismember(requested, channelIDs))
        [~, loc] = ismember(requested, channelIDs);
        plane = plane(:,:,loc);
        channelIDs = channelIDs(loc);
    elseif size(plane,3) == numel(requested)
        channelIDs = requested;
    else
        error('compareMeasuredAndFitted3DPSFs:MissingChannels', ...
            'The PTU does not contain all requested detector channels.');
    end

    if ~isempty(opts.channelOrder)
        order = opts.channelOrder(:);
        if all(ismember(order, channelIDs))
            [~, loc] = ismember(order, channelIDs);
        elseif all(order >= 1 & order <= numel(channelIDs))
            loc = order;
        else
            error('compareMeasuredAndFitted3DPSFs:BadChannelOrder', ...
                'channelOrder must contain channel IDs or detector indices.');
        end
        plane = plane(:,:,loc);
        channelIDs = channelIDs(loc);
    end
end

function [darkImage, darkVector] = loadDarkCorrection(opts, planeSize, channelIDs)
    darkImage = [];
    darkVector = [];
    if isempty(opts.darkFile) || exist(opts.darkFile, 'file') ~= 2
        return;
    end
    try
        darkOpts = opts;
        darkOpts.channelIDs = channelIDs;
        darkOpts.channelOrder = [];
        [darkPlane, darkIDs] = readPtuDetectorPlane(opts.darkFile, [], darkOpts);
        if ~isequal(darkIDs(:), channelIDs(:))
            [present, loc] = ismember(channelIDs(:), darkIDs(:));
            if all(present)
                darkPlane = darkPlane(:,:,loc);
            end
        end
        if isequal(size(darkPlane), planeSize)
            darkImage = darkPlane;
        else
            darkVector = squeeze(median(median(max(darkPlane,0),1),2));
            darkVector = double(darkVector(:));
        end
    catch scanError
        [darkVector, nPhotons] = readDarkVectorFromTttr( ...
            opts.darkFile, opts.ptuPhotonsPerChunk, planeSize, channelIDs);
        if opts.verbose
            fprintf(['[compareMeasuredAndFitted3DPSFs] Dark PTU is not an image scan; ' ...
                'used TTTR channel counts (%d photons). Scan-reader error: %s\n'], ...
                nPhotons, scanError.message);
        end
    end
end

function [darkVector, nPhotons] = readDarkVectorFromTttr( ...
        darkFile, chunkSize, planeSize, channelIDs)
    if exist('PTU_Read_Head', 'file') ~= 2 || exist('PTU_Read', 'file') ~= 2
        error('compareMeasuredAndFitted3DPSFs:MissingTttrReader', ...
            'PTU_Read_Head/PTU_Read is required for the non-image dark PTU.');
    end

    head = PTU_Read_Head(darkFile);
    if isempty(head) || ~isfield(head, 'TTResult_NumberOfRecords')
        error('compareMeasuredAndFitted3DPSFs:BadDarkPtu', ...
            'Dark-count PTU has no TTResult_NumberOfRecords field: %s', darkFile);
    end

    targetIDs = double(channelIDs(:));
    nBins = max([64; targetIDs + 1]);
    channelCounts = zeros(nBins,1);
    nRecords = double(head.TTResult_NumberOfRecords);
    chunkSize = max(1,double(chunkSize));
    cnt = 0;
    nPhotons = 0;

    while cnt < nRecords
        nRead = min(chunkSize,nRecords-cnt);
        [~,~,chan,special,num] = PTU_Read( ...
            darkFile,[cnt+1,nRead],head);
        cnt = cnt + double(num);
        if num <= 0
            break;
        end

        isPhoton = special == 0;
        chan = double(chan(isPhoton));
        chan = chan(isfinite(chan) & chan >= 0 & chan < nBins);
        nPhotons = nPhotons + numel(chan);
        if ~isempty(chan)
            channelCounts = channelCounts + ...
                accumarray(chan(:)+1,1,[nBins,1]);
        end
    end

    if any(targetIDs < 0 | targetIDs+1 > nBins)
        error('compareMeasuredAndFitted3DPSFs:DarkChannelMismatch', ...
            'Dark-count vector cannot cover detector channels: %s', ...
            num2str(targetIDs(:).'));
    end
    nScanPixels = max(1,prod(double(planeSize(1:2))));
    darkVector = channelCounts(targetIDs+1) ./ nScanPixels;
end

function plane = applyCorrections(plane, darkImage, darkVector, opts)
    plane = double(plane);
    if ~isempty(darkImage)
        plane = max(plane - opts.darkScale * darkImage, 0);
    elseif ~isempty(darkVector)
        plane = max(plane - opts.darkScale * reshape(darkVector,1,1,[]), 0);
    end
    if opts.subtractBoundary
        width = max(1, min(round(opts.boundaryWidthPx), ...
            floor(min(size(plane,1),size(plane,2))/2)));
        for ch = 1:size(plane,3)
            img = plane(:,:,ch);
            border = [reshape(img(1:width,:),[],1); ...
                reshape(img(end-width+1:end,:),[],1); ...
                reshape(img(:,1:width),[],1); ...
                reshape(img(:,end-width+1:end),[],1)];
            bg = median(border(:), 'omitnan');
            plane(:,:,ch) = max(img - bg, 0);
        end
    end
end

function focusIndex = resolveFocusIndex(T, fitResult)
    focusIndex = [];
    if isfield(fitResult, 'stackSelection') && ...
            isfield(fitResult.stackSelection, 'focusZUm')
        focusZ = double(fitResult.stackSelection.focusZUm);
        [~, focusIndex] = min(abs(double(T.z_um) - focusZ));
    end
    if isempty(focusIndex)
        if ismember('total_signal', T.Properties.VariableNames)
            [~, focusIndex] = max(double(T.total_signal));
        else
            focusIndex = round((height(T)+1)/2);
        end
    end
end

function [modelPsf, modelBead, info] = calculateFittedVolumes( ...
        fitResult, measuredVolume, dataInfo, opts)
    sim = fitResult.fit.sim;
    coeffs = fitResult.fit.estCoeffs;

    pixelSize = dataInfo.xyPixelSizeUm;
    sim.nx = size(measuredVolume,2);
    sim.fovXY = pixelSize * (sim.nx - 1);
    sim.x = centeredAxis(sim.nx, pixelSize);
    sim.y = centeredAxis(size(measuredVolume,1), pixelSize);
    sim.dx = pixelSize;

    if ~isempty(opts.detectorSubsamples)
        sim.detectorSubsamples = opts.detectorSubsamples;
    end

    relativeZ = dataInfo.relativeZUm;
    dz = opts.modelZStepUm;
    if isempty(dz) || ~isfinite(dz) || dz <= 0
        dz = median(abs(diff(relativeZ)));
    end
    maxZ = max(abs(relativeZ)) + opts.modelZPaddingUm;
    nHalf = ceil(maxZ / dz);
    sim.z = (-nHalf:nHalf) * dz;
    sim.nz = numel(sim.z);
    sim.nzRange = sim.z(end) - sim.z(1);
    sim.obj = beadObject3D(sim);

    if opts.verbose
        fprintf('[compareMeasuredAndFitted3DPSFs] Calculating fitted 3D PSFs for %d detectors.\n', ...
            size(sim.detXY,1));
    end
    hEff = spadEffectivePSFArray(sim, coeffs, 'none');

    beadKernel = cropNonzeroKernel(sim.obj);
    beadResponse = zeros(size(hEff));
    for ch = 1:size(hEff,4)
        beadResponse(:,:,:,ch) = convn(hEff(:,:,:,ch), beadKernel, 'same');
    end

    axialShiftUm = 0;
    if opts.registerAxialShift
        axialShiftUm = estimateAxialShift(measuredVolume, beadResponse, ...
            sim.z, relativeZ, opts);
    end
    queryZ = relativeZ - axialShiftUm;
    modelPsf = interpolateVolumeZ(hEff, sim.z, queryZ);
    modelBead = interpolateVolumeZ(beadResponse, sim.z, queryZ);

    info = struct();
    info.sim = sim;
    info.coeffs = coeffs;
    info.modelZUm = sim.z;
    info.dataRelativeZUm = relativeZ;
    info.axialRegistrationShiftUm = axialShiftUm;
    info.queryZUm = queryZ;
    info.beadKernelSize = size(beadKernel);
end

function kernel = cropNonzeroKernel(kernel)
    mask = kernel > max(kernel(:)) * 1e-10;
    [iy, ix, iz] = ind2sub(size(kernel), find(mask));
    if isempty(ix)
        return;
    end
    kernel = kernel(min(iy):max(iy), min(ix):max(ix), min(iz):max(iz));
    s = sum(kernel(:));
    if s > 0
        kernel = kernel / s;
    end
end

function shiftUm = estimateAxialShift(measured, modelFull, modelZ, dataZ, opts)
    measuredProfile = squeeze(sum(measured, [1 2 4]));
    measuredProfile = normalizeVector(measuredProfile);
    shifts = -opts.maxAxialRegistrationUm:opts.axialRegistrationStepUm: ...
        opts.maxAxialRegistrationUm;
    score = nan(size(shifts));
    modelProfileFull = squeeze(sum(modelFull, [1 2 4]));
    for k = 1:numel(shifts)
        profile = interp1(modelZ, modelProfileFull, ...
            dataZ - shifts(k), 'linear', 0);
        profile = normalizeVector(profile);
        score(k) = dot(measuredProfile(:), profile(:));
    end
    [~, idx] = max(score);
    shiftUm = shifts(idx);
end

function out = interpolateVolumeZ(vol, sourceZ, queryZ)
    nz = size(vol,3);
    flat = reshape(permute(vol, [3 1 2 4]), nz, []);
    flatQ = interp1(sourceZ(:), flat, queryZ(:), 'linear', 0);
    out = permute(reshape(flatQ, numel(queryZ), size(vol,1), ...
        size(vol,2), size(vol,4)), [2 3 1 4]);
end

function [measuredN, modelBeadN, modelPsfN] = normalizeComparisonVolumes( ...
        measured, modelBead, modelPsf)
    measuredN = normalizeEachChannel(measured, 'sum');
    modelBeadN = normalizeEachChannel(modelBead, 'sum');
    modelPsfN = normalizeEachChannel(modelPsf, 'sum');
end

function out = normalizeEachChannel(vol, mode)
    out = max(double(vol), 0);
    for ch = 1:size(out,4)
        one = out(:,:,:,ch);
        if strcmp(mode, 'peak')
            scale = max(one(:));
        else
            scale = sum(one(:));
        end
        if scale > 0
            out(:,:,:,ch) = one / scale;
        end
    end
end

function T = calculateMetrics(measured, modelBead, modelPsf, ...
        measuredRaw, modelBeadRaw, dataInfo, fitResult, supportThreshold)
    nCh = size(measured,4);
    detectorIndex = (1:nCh).';
    channelID = dataInfo.channelIDs(:);
    detectorXUm = fitResult.fit.sim.detXY(:,1);
    detectorYUm = fitResult.fit.sim.detXY(:,2);
    airyUnitUm = resolveAiryUnitUm(fitResult.fit.sim);
    detectorXAU = detectorXUm / airyUnitUm;
    detectorYAU = detectorYUm / airyUnitUm;
    detectorRadiusAU = hypot(detectorXAU, detectorYAU);
    corr3D = nan(nCh,1);
    cosine3D = nan(nCh,1);
    nrmse3D = nan(nCh,1);
    corrXY = nan(nCh,1);
    corrXZ = nan(nCh,1);
    corrYZ = nan(nCh,1);
    opticalToBeadCosine = nan(nCh,1);
    measuredDetectorFraction = zeros(nCh,1);
    modeledDetectorFraction = zeros(nCh,1);
    measuredTotal = sum(measuredRaw(:));
    modeledTotal = sum(modelBeadRaw(:));

    for ch = 1:nCh
        A = measured(:,:,:,ch);
        B = modelBead(:,:,:,ch);
        P = modelPsf(:,:,:,ch);
        corr3D(ch) = maskedCorrelation(A, B, supportThreshold);
        cosine3D(ch) = cosineSimilarity(A, B);
        nrmse3D(ch) = norm(A(:)-B(:)) / max(norm(A(:)), eps);
        corrXY(ch) = maskedCorrelation(projectOne(A,'xy','max'), ...
            projectOne(B,'xy','max'), supportThreshold);
        corrXZ(ch) = maskedCorrelation(projectOne(A,'xz','max'), ...
            projectOne(B,'xz','max'), supportThreshold);
        corrYZ(ch) = maskedCorrelation(projectOne(A,'yz','max'), ...
            projectOne(B,'yz','max'), supportThreshold);
        opticalToBeadCosine(ch) = cosineSimilarity(P, B);
        measuredDetectorFraction(ch) = sum(measuredRaw(:,:,:,ch),'all') / ...
            max(measuredTotal,eps);
        modeledDetectorFraction(ch) = sum(modelBeadRaw(:,:,:,ch),'all') / ...
            max(modeledTotal,eps);
    end

    T = table(detectorIndex, channelID, detectorXAU, detectorYAU, ...
        detectorRadiusAU, ...
        corr3D, cosine3D, nrmse3D, corrXY, corrXZ, corrYZ, ...
        opticalToBeadCosine, measuredDetectorFraction, ...
        modeledDetectorFraction);
end

function writeComparisonOutputs(results, measuredN, modelBeadN, modelPsfN)
    outDir = results.outputDir;
    ensureFolder(outDir);
    writetable(results.metrics, fullfile(outDir, 'detector_3d_psf_metrics.csv'));

    if results.options.writePerDetectorFigures
        detectorDir = fullfile(outDir, 'per_detector');
        ensureFolder(detectorDir);
        for ch = 1:size(measuredN,4)
            outFile = fullfile(detectorDir, sprintf( ...
                'detector_%02d_channel_%02d_xy_xz_yz.png', ...
                ch, results.dataInfo.channelIDs(ch)));
            writeDetectorFigure(results, measuredN, modelBeadN, modelPsfN, ...
                ch, outFile);
        end
    end

    if results.options.writeArrayFigures
        arrayDir = fullfile(outDir, 'array_overview');
        ensureFolder(arrayDir);
        planes = {'xy','xz','yz'};
        for k = 1:numel(planes)
            outFile = fullfile(arrayDir, sprintf( ...
                'all_detectors_%s_comparison.png', upper(planes{k})));
            writeArrayProjectionFigure(results, measuredN, modelBeadN, ...
                planes{k}, outFile);
        end
        writeMetricFigure(results, fullfile(arrayDir, ...
            'detector_3d_comparison_metrics.png'));
    end

    if results.options.writeVolumeMat
        measuredVolume = results.measuredVolume; %#ok<NASGU>
        modelBeadVolume = results.modelBeadVolume; %#ok<NASGU>
        modelEffectivePsf = results.modelEffectivePsf; %#ok<NASGU>
        dataInfo = results.dataInfo; %#ok<NASGU>
        modelInfo = results.modelInfo; %#ok<NASGU>
        metrics = results.metrics; %#ok<NASGU>
        save(fullfile(outDir, 'measured_and_fitted_3d_psfs.mat'), ...
            'measuredVolume', 'modelBeadVolume', 'modelEffectivePsf', ...
            'dataInfo', 'modelInfo', 'metrics', '-v7.3');
    end
end

function writeDetectorFigure(results, measured, modelBead, modelPsf, ch, outFile)
    A = measured(:,:,:,ch);
    B = modelBead(:,:,:,ch);
    P = modelPsf(:,:,:,ch);
    R = A - B;
    projections = {'xy','xz','yz'};
    rowVolumes = {A, B, P, R};
    rowLabels = {'Measured bead PSF estimate', ...
        'Fitted finite-bead response', 'Fitted optical effective PSF', ...
        'Measured - fitted bead'};

    fig = figure('Visible','off','Color','w','Position',[40 40 1250 1180]);
    tl = tiledlayout(fig, 4, 3, 'Padding','compact','TileSpacing','compact');
    for row = 1:4
        for col = 1:3
            ax = nexttile(tl);
            image = projectOne(rowVolumes{row}, projections{col}, ...
                results.options.projectionMode);
            [axisX, axisY, xLabel, yLabel] = projectionAxes( ...
                results.dataInfo, projections{col});
            imagesc(ax, axisX, axisY, image);
            axis(ax, 'image');
            set(ax, 'YDir','normal');
            xlabel(ax, xLabel);
            ylabel(ax, yLabel);
            title(ax, sprintf('%s: %s', rowLabels{row}, upper(projections{col})));
            if row == 4
                lim = max(abs(image(:)));
                if lim > 0
                    caxis(ax, [-lim lim]);
                end
                colormap(ax, redBlueMap(256));
            else
                colormap(ax, hot);
            end
            colorbar(ax);
        end
    end
    metric = results.metrics(ch,:);
    airyUnitUm = resolveAiryUnitUm(results.fitResult.fit.sim);
    title(tl, sprintf(['Detector %d, TCSPC channel %d, detXY=(%.3f, %.3f) AU; ' ...
        '3D corr %.3f, cosine %.3f, NRMSE %.3f'], ...
        ch, results.dataInfo.channelIDs(ch), ...
        results.fitResult.fit.sim.detXY(ch,1) / airyUnitUm, ...
        results.fitResult.fit.sim.detXY(ch,2) / airyUnitUm, ...
        metric.corr3D, metric.cosine3D, metric.nrmse3D), ...
        'FontWeight','bold');
    exportFigure(fig, outFile, 170);
end

function writeArrayProjectionFigure(results, measured, modelBead, planeName, outFile)
    measuredCanvas = detectorProjectionCanvas( ...
        results.fitResult.fit.sim, measured, planeName, ...
        results.options.projectionMode, false);
    modelCanvas = detectorProjectionCanvas( ...
        results.fitResult.fit.sim, modelBead, planeName, ...
        results.options.projectionMode, false);
    residualCanvas = measuredCanvas - modelCanvas;

    fig = figure('Visible','off','Color','w','Position',[50 80 1500 780]);
    tl = tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');
    ax = nexttile(tl);
    imagesc(ax, measuredCanvas); axis(ax,'image'); axis(ax,'off'); colormap(ax,'hot');
    title(ax, sprintf('Measured %s projections', upper(planeName)));
    ax = nexttile(tl);
    imagesc(ax, modelCanvas); axis(ax,'image'); axis(ax,'off'); colormap(ax,'hot');
    title(ax, sprintf('Fitted bead %s projections', upper(planeName)));
    ax = nexttile(tl);
    imagesc(ax, residualCanvas); axis(ax,'image'); axis(ax,'off');
    lim = max(abs(residualCanvas(:)), [], 'omitnan');
    if lim > 0
        caxis(ax,[-lim lim]);
    end
    colormap(ax,redBlueMap(256));
    title(ax, 'Measured - fitted');
    title(tl, sprintf('All 23 detector pixels: %s projection comparison', ...
        upper(planeName)), 'FontWeight','bold');
    exportFigure(fig, outFile, 180);
end

function canvas = detectorProjectionCanvas(sim, vol, planeName, mode, signedMap)
    idxGrid = regularDetectorIndexGrid(sim, size(vol,4));
    sample = projectOne(vol(:,:,:,1), planeName, mode);
    tileH = size(sample,1);
    tileW = size(sample,2);
    gap = 2;
    canvas = nan(size(idxGrid,1)*tileH + (size(idxGrid,1)-1)*gap, ...
        size(idxGrid,2)*tileW + (size(idxGrid,2)-1)*gap);
    for row = 1:size(idxGrid,1)
        for col = 1:size(idxGrid,2)
            ch = idxGrid(row,col);
            if ~isfinite(ch)
                continue;
            end
            image = projectOne(vol(:,:,:,ch), planeName, mode);
            if ~signedMap
                pk = max(image(:));
                if pk > 0
                    image = image / pk;
                end
            end
            rr = (row-1)*(tileH+gap) + (1:tileH);
            cc = (col-1)*(tileW+gap) + (1:tileW);
            canvas(rr,cc) = image;
        end
    end
end

function idxGrid = regularDetectorIndexGrid(sim, nCh)
    layoutName = 'honeycomb23';
    if isfield(sim, 'detectorLayout') && ~isempty(sim.detectorLayout)
        layoutName = sim.detectorLayout;
    end
    try
        [~, idxGrid] = detectorLayout(layoutName, 1);
    catch
        idxGrid = [];
    end
    if nnz(isfinite(idxGrid)) ~= nCh
        nCol = ceil(sqrt(nCh));
        nRow = ceil(nCh/nCol);
        idxGrid = nan(nRow,nCol);
        idxGrid(1:nCh) = 1:nCh;
    end
end

function airyUnitUm = resolveAiryUnitUm(sim)
    airyUnitUm = NaN;
    if isfield(sim, 'airyUnitUm') && isfinite(sim.airyUnitUm) && sim.airyUnitUm > 0
        airyUnitUm = double(sim.airyUnitUm);
    elseif isfield(sim, 'lamEm') && isfield(sim, 'NA') && ...
            isfinite(sim.lamEm) && sim.lamEm > 0 && isfinite(sim.NA) && sim.NA > 0
        airyUnitUm = 1.22 * double(sim.lamEm) / double(sim.NA);
    end
    if ~isfinite(airyUnitUm) || airyUnitUm <= 0
        error('compareMeasuredAndFitted3DPSFs:MissingAiryScale', ...
            'Cannot convert detector distances to AU because lambda_em or NA is unavailable.');
    end
end

function writeMetricFigure(results, outFile)
    T = results.metrics;
    fig = figure('Visible','off','Color','w','Position',[80 80 1150 680]);
    tl = tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
    ax = nexttile(tl);
    bar(ax,T.corr3D); ylim(ax,[-1 1]); ylabel(ax,'Pearson correlation');
    xlabel(ax,'detector index'); title(ax,'3D correlation'); grid(ax,'on');
    ax = nexttile(tl);
    bar(ax,T.cosine3D); ylim(ax,[0 1]); ylabel(ax,'cosine similarity');
    xlabel(ax,'detector index'); title(ax,'3D cosine similarity'); grid(ax,'on');
    ax = nexttile(tl);
    bar(ax,T.nrmse3D); ylabel(ax,'NRMSE');
    xlabel(ax,'detector index'); title(ax,'3D normalized error'); grid(ax,'on');
    ax = nexttile(tl);
    plot(ax,T.detectorIndex,T.corrXY,'-o', ...
        T.detectorIndex,T.corrXZ,'-s', ...
        T.detectorIndex,T.corrYZ,'-^','LineWidth',1.1);
    ylim(ax,[-1 1]); xlabel(ax,'detector index');
    ylabel(ax,'projection correlation'); legend(ax,{'XY','XZ','YZ'});
    title(ax,'Projection correlations'); grid(ax,'on');
    title(tl,'Measured versus fitted finite-bead 3D response');
    exportFigure(fig,outFile,180);
end

function image = projectOne(vol, planeName, mode)
    switch lower(planeName)
        case 'xy'
            if strcmp(mode,'sum')
                image = squeeze(sum(vol,3));
            else
                image = squeeze(max(vol,[],3));
            end
        case 'xz'
            if strcmp(mode,'sum')
                image = squeeze(sum(vol,1)).';
            else
                image = squeeze(max(vol,[],1)).';
            end
        case 'yz'
            if strcmp(mode,'sum')
                image = squeeze(sum(vol,2)).';
            else
                image = squeeze(max(vol,[],2)).';
            end
        otherwise
            error('compareMeasuredAndFitted3DPSFs:BadProjection', ...
                'Projection must be XY, XZ, or YZ.');
    end
end

function [xAxis,yAxis,xLabel,yLabel] = projectionAxes(info, planeName)
    switch lower(planeName)
        case 'xy'
            xAxis = info.xUm; yAxis = info.yUm;
            xLabel = 'x (um)'; yLabel = 'y (um)';
        case 'xz'
            xAxis = info.xUm; yAxis = info.relativeZUm;
            xLabel = 'x (um)'; yLabel = 'z from focus (um)';
        case 'yz'
            xAxis = info.yUm; yAxis = info.relativeZUm;
            xLabel = 'y (um)'; yLabel = 'z from focus (um)';
    end
end

function value = maskedCorrelation(a,b,threshold)
    a = double(a); b = double(b);
    mask = isfinite(a) & isfinite(b);
    limit = threshold * max([a(:);b(:)]);
    mask = mask & (a >= limit | b >= limit);
    if nnz(mask) < 3
        value = NaN;
        return;
    end
    aa = a(mask); bb = b(mask);
    aa = aa - mean(aa); bb = bb - mean(bb);
    value = dot(aa,bb) / max(norm(aa)*norm(bb),eps);
end

function value = cosineSimilarity(a,b)
    a = double(a(:)); b = double(b(:));
    value = dot(a,b) / max(norm(a)*norm(b),eps);
end

function v = normalizeVector(v)
    v = max(double(v(:)),0);
    n = norm(v);
    if n > 0
        v = v/n;
    end
end

function centerXY = intensityCentroid(img)
    img = max(double(img),0);
    threshold = 0.20 * max(img(:));
    weights = img;
    weights(weights < threshold) = 0;
    if sum(weights(:)) <= 0
        [~,idx] = max(img(:));
        [iy,ix] = ind2sub(size(img),idx);
        centerXY = [ix iy];
        return;
    end
    [yy,xx] = ndgrid(1:size(img,1),1:size(img,2));
    centerXY = [sum(xx(:).*weights(:)),sum(yy(:).*weights(:))] / sum(weights(:));
end

function pixelSize = resolvePixelSize(fitResult, head)
    pixelSize = NaN;
    if isfield(fitResult.fit,'sim') && isfield(fitResult.fit.sim,'dx')
        pixelSize = double(fitResult.fit.sim.dx);
    end
    if (~isfinite(pixelSize) || pixelSize <= 0) && isstruct(head) && ...
            isfield(head,'ImgHdr_PixResol')
        pixelSize = headerNumeric(head.ImgHdr_PixResol);
    end
    if ~isfinite(pixelSize) || pixelSize <= 0
        error('compareMeasuredAndFitted3DPSFs:MissingPixelSize', ...
            'Could not determine the scan pixel size.');
    end
end

function value = headerNumeric(raw)
    value = NaN;
    if isnumeric(raw) && ~isempty(raw)
        value = double(raw(1));
    elseif isstruct(raw)
        names = {'value','Value','data','Data','TagValue','tagValue'};
        for k = 1:numel(names)
            if isfield(raw,names{k}) && isnumeric(raw.(names{k})) && ...
                    ~isempty(raw.(names{k}))
                value = double(raw.(names{k})(1));
                return;
            end
        end
    end
end

function axisValues = centeredAxis(n,pixelSize)
    axisValues = ((1:n) - (n+1)/2) * pixelSize;
end

function text = tableText(column,index)
    if iscell(column)
        text = char(column{index});
    elseif isstring(column)
        text = char(column(index));
    elseif iscategorical(column)
        text = char(string(column(index)));
    else
        text = char(column(index,:));
    end
end

function frame = tableFrame(T,index)
    frame = [];
    if ismember('source_frame',T.Properties.VariableNames)
        value = double(T.source_frame(index));
        if isfinite(value) && value >= 1
            frame = value;
        end
    end
end

function ensureFolder(folder)
    if exist(folder,'dir') ~= 7
        mkdir(folder);
    end
end

function exportFigure(fig,fileName,resolution)
    try
        exportgraphics(fig,fileName,'Resolution',resolution);
    catch
        set(fig,'PaperPositionMode','auto');
        print(fig,fileName,'-dpng',sprintf('-r%d',resolution));
    end
    close(fig);
end

function cmap = redBlueMap(n)
    if nargin < 1
        n = 256;
    end
    half = ceil(n/2);
    blue = [linspace(0.1,1,half).', linspace(0.2,1,half).', ones(half,1)];
    red = [ones(n-half,1), linspace(1,0.2,n-half).', linspace(1,0.1,n-half).'];
    cmap = [blue; flipud(red)];
end

function stem = sanitizeFileName(name)
    stem = regexprep(char(name),'[^A-Za-z0-9._-]+','_');
    stem = regexprep(stem,'_+','_');
end

function printSummary(results)
    T = results.metrics;
    fprintf('\nMeasured versus fitted 3D detector PSFs\n');
    fprintf('  detectors: %d\n',height(T));
    fprintf('  measured planes: %d, z %.3f to %.3f um from focus\n', ...
        numel(results.dataInfo.relativeZUm), ...
        min(results.dataInfo.relativeZUm),max(results.dataInfo.relativeZUm));
    fprintf('  fitted model axial registration: %+.3f um\n', ...
        results.modelInfo.axialRegistrationShiftUm);
    fprintf('  median 3D correlation: %.3f\n',median(T.corr3D,'omitnan'));
    fprintf('  median cosine similarity: %.3f\n',median(T.cosine3D,'omitnan'));
    fprintf('  median NRMSE: %.3f\n',median(T.nrmse3D,'omitnan'));
    fprintf('  outputs: %s\n\n',results.outputDir);
end
