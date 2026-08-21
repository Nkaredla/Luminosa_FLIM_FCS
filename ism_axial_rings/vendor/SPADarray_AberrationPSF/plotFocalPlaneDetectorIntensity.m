function out = plotFocalPlaneDetectorIntensity(stackInput, varargin)
%PLOTFOCALPLANEDETECTORINTENSITY Visualize bead focus and 23-detector signal.
%
%   out = plotFocalPlaneDetectorIntensity()
%
%   Uses the current Luminosa aberration dataset by default:
%       D:\Luminosa\Data\ISM_Aberation2_73
%
%   The function selects the focal z plane from the stack signal trace,
%   reads the detector-resolved focal plane, finds the bead centre in the
%   summed ISM image, and plots:
%       1) the bead intensity image at focus
%       2) the 23-pixel detector intensity distribution at the bead centre
%       3) a channel-by-channel intensity bar plot
%       4) a scan-position grid of detector honeycomb maps around the bead
%
%   Example:
%       out = plotFocalPlaneDetectorIntensity();
%
%       out = plotFocalPlaneDetectorIntensity([], ...
%           'scanName', '0.5uW_0.15collar_80mmlens_20260515-164337');
%
%       out = plotFocalPlaneDetectorIntensity(raw4d, 'planeZ', zUm);
%
%       out = plotFocalPlaneDetectorIntensity([], 'scanGridSize', 10, ...
%           'scanStepPx', 1);

    if nargin < 1 || isempty(stackInput)
        stackInput = 'D:\Luminosa\Data\ISM_Aberation2_73';
    end

    opts = parseOptions(varargin{:});
    addRequiredPaths(opts);

    source = resolveFocusPlaneSource(stackInput, opts);
    rawFocus = loadFocusDetectorPlane(source, opts);
    correctedFocus = correctFocusPlaneForDisplay(rawFocus, opts);

    focusImage = sum(correctedFocus, 3);
    [centerXY, centerYX] = estimateBeadCenter(focusImage, opts);
    detectorValues = squeeze(correctedFocus(centerYX(1), centerYX(2), :));
    detectorValues = double(detectorValues(:));
    normalizedValues = detectorValues / max(sum(detectorValues), eps);
    scanGrid = sampleScanGridDetectorValues(correctedFocus, centerXY, opts);

    out = struct();
    out.source = source;
    out.rawFocus = rawFocus;
    out.correctedFocus = correctedFocus;
    out.focusImage = focusImage;
    out.centerXY = centerXY;
    out.centerPixelYX = centerYX;
    out.channelIDs = source.channelIDs(:);
    out.detectorValues = detectorValues;
    out.normalizedDetectorValues = normalizedValues;
    out.scanGrid = scanGrid;
    out.outputDir = resolveOutputDir(source, opts);

    if opts.writeOutputs
        writeOutputs(out, opts);
    end

    if opts.verbose
        fprintf('[plotFocalPlaneDetectorIntensity] focus plane %d at z=%.4f um.\n', ...
            source.focusIndex, source.focusZUm);
        fprintf('[plotFocalPlaneDetectorIntensity] bead centre x=%.2f, y=%.2f; sampled pixel row=%d, col=%d.\n', ...
            centerXY(1), centerXY(2), centerYX(1), centerYX(2));
        if ~isempty(out.scanGrid)
            fprintf('[plotFocalPlaneDetectorIntensity] sampled %dx%d scan positions around the bead centre, step %.3g px.\n', ...
                out.scanGrid.gridSize, out.scanGrid.gridSize, out.scanGrid.stepPx);
        end
        if ~isempty(out.outputDir)
            fprintf('[plotFocalPlaneDetectorIntensity] wrote outputs to %s\n', out.outputDir);
        end
    end
end

function opts = parseOptions(varargin)
    p = inputParser;
    p.FunctionName = 'plotFocalPlaneDetectorIntensity';

    addParameter(p, 'dataRoot', 'D:\Luminosa\Data\ISM_Aberation2_73');
    addParameter(p, 'scanName', '0.4uW_0.19collar_80mmlens_20260515-155248');
    addParameter(p, 'scanPattern', '');
    addParameter(p, 'scanIndex', []);
    addParameter(p, 'volumeMat', '');
    addParameter(p, 'alignmentCsv', '');
    addParameter(p, 'planeZ', []);
    addParameter(p, 'zStepUm', 0.05);
    addParameter(p, 'focusIndex', []);
    addParameter(p, 'focusMetric', 'signal');
    addParameter(p, 'centerXY', []);
    addParameter(p, 'centerMode', 'centroid');
    addParameter(p, 'centerThresholdFraction', 0.20);
    addParameter(p, 'scanGridSize', 5);
    addParameter(p, 'scanStepPx', 1);
    addParameter(p, 'subtractBoundary', true);
    addParameter(p, 'boundaryWidthPx', 3);
    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'ptuReaderFolder', defaultPtuReaderFolder());
    addParameter(p, 'ptuPhotonsPerChunk', 1e6);
    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;
    opts.focusMetric = lower(char(opts.focusMetric));
    opts.centerMode = lower(char(opts.centerMode));
    opts.scanGridSize = normalizeScanGridSize(opts.scanGridSize);
    opts.scanStepPx = normalizePositiveScalar(opts.scanStepPx, 'scanStepPx');
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

function source = resolveFocusPlaneSource(stackInput, opts)
    if isnumeric(stackInput) || iscell(stackInput)
        detStack = standardizeDetectorZStack(stackInput, opts.channelIDs);
        signal = zSignalTrace(detStack, opts.focusMetric);
        planeZ = resolvePlaneZ(numel(signal), [], opts);
        focusIdx = chooseFocusIndex(signal, opts);
        source = emptySource();
        source.kind = 'numeric';
        source.label = 'numeric_detector_stack';
        source.detectorStack = detStack;
        source.focusIndex = focusIdx;
        source.focusZUm = planeZ(focusIdx);
        source.signalTrace = signal(:).';
        source.planeZ = planeZ(:).';
        source.channelIDs = resolveChannelIDs(opts.channelIDs, size(detStack, 3));
        return;
    end

    inputPath = char(stackInput);
    if exist(inputPath, 'dir') == 7
        scanFolder = resolveScanFolder(inputPath, opts);
        source = resolveFocusFromScanFolder(scanFolder, opts);
        return;
    end

    if exist(inputPath, 'file') ~= 2
        error('plotFocalPlaneDetectorIntensity:MissingInput', ...
            'Input was not found: %s', inputPath);
    end

    [~,~,ext] = fileparts(inputPath);
    switch lower(ext)
        case '.mat'
            source = resolveFocusFromMat(inputPath, opts);
        case '.ptu'
            detStack = readPtuDetectorStack(inputPath, opts, []);
            signal = zSignalTrace(detStack, opts.focusMetric);
            planeZ = resolvePlaneZ(numel(signal), [], opts);
            focusIdx = chooseFocusIndex(signal, opts);
            source = emptySource();
            source.kind = 'ptu_stack';
            source.label = localFileStem(inputPath);
            source.file = inputPath;
            source.frame = focusIdx;
            source.detectorStack = detStack;
            source.focusIndex = focusIdx;
            source.focusZUm = planeZ(focusIdx);
            source.signalTrace = signal(:).';
            source.planeZ = planeZ(:).';
            source.channelIDs = resolveChannelIDs(opts.channelIDs, size(detStack, 3));
        otherwise
            error('plotFocalPlaneDetectorIntensity:BadInputFile', ...
                'Unsupported input file "%s". Use MAT or PTU.', inputPath);
    end
end

function source = emptySource()
    source = struct();
    source.kind = '';
    source.label = '';
    source.scanFolder = '';
    source.file = '';
    source.frame = [];
    source.detectorStack = [];
    source.focusIndex = NaN;
    source.focusZUm = NaN;
    source.signalTrace = [];
    source.planeZ = [];
    source.channelIDs = [];
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
        error('plotFocalPlaneDetectorIntensity:NoScanFolders', ...
            'No scan folders were found below %s.', dataRoot);
    end

    if ~isempty(opts.scanIndex)
        idx = round(opts.scanIndex);
        if idx < 1 || idx > numel(folders)
            error('plotFocalPlaneDetectorIntensity:BadScanIndex', ...
                'scanIndex must be between 1 and %d.', numel(folders));
        end
    else
        [~, idx] = max([folders.datenum]);
    end
    scanFolder = fullfile(folders(idx).folder, folders(idx).name);
end

function source = resolveFocusFromScanFolder(scanFolder, opts)
    csvFile = char(opts.alignmentCsv);
    if isempty(csvFile)
        csvFile = findBatchAlignmentCsvForScanFolder(scanFolder);
    end
    if ~isempty(csvFile) && exist(csvFile, 'file') == 2
        source = resolveFocusFromAlignmentCsv(csvFile, opts);
        source.scanFolder = scanFolder;
        return;
    end

    matFile = char(opts.volumeMat);
    if isempty(matFile)
        matFile = findBatchVolumeMatForScanFolder(scanFolder);
    end
    if ~isempty(matFile) && exist(matFile, 'file') == 2
        source = resolveFocusFromMat(matFile, opts);
        source.scanFolder = scanFolder;
        return;
    end

    files = sortedSeriesFiles(scanFolder);
    if isempty(files)
        error('plotFocalPlaneDetectorIntensity:NoSeriesFiles', ...
            'No Series_*.ptu files or batch stack outputs were found for %s.', scanFolder);
    end

    fileNames = fullfile({files.folder}, {files.name});
    signal = zeros(1, numel(fileNames));
    for k = 1:numel(fileNames)
        detPlane = readPtuDetectorStack(fileNames{k}, opts, []);
        signal(k) = sum(max(detPlane(:), 0));
    end
    planeZ = resolvePlaneZ(numel(signal), [], opts);
    focusIdx = chooseFocusIndex(signal, opts);

    source = emptySource();
    source.kind = 'scan_folder';
    source.scanFolder = scanFolder;
    source.label = localFileStem(scanFolder);
    source.file = fileNames{focusIdx};
    source.frame = [];
    source.focusIndex = focusIdx;
    source.focusZUm = planeZ(focusIdx);
    source.signalTrace = signal(:).';
    source.planeZ = planeZ(:).';
    source.channelIDs = opts.channelIDs(:);
end

function source = resolveFocusFromAlignmentCsv(csvFile, opts)
    T = readtable(csvFile);
    if ~all(ismember({'total_signal','source_file'}, T.Properties.VariableNames))
        error('plotFocalPlaneDetectorIntensity:BadAlignmentCsv', ...
            'Alignment CSV must contain total_signal and source_file columns.');
    end

    signal = double(T.total_signal(:)).';
    focusIdx = chooseFocusIndex(signal, opts);
    zUm = nan(1, numel(signal));
    if ismember('z_um', T.Properties.VariableNames)
        zUm = double(T.z_um(:)).';
    end
    if any(~isfinite(zUm))
        zUm = resolvePlaneZ(numel(signal), [], opts);
    end

    frame = [];
    if ismember('source_frame', T.Properties.VariableNames)
        frameValue = double(T.source_frame(focusIdx));
        if isfinite(frameValue) && frameValue >= 1
            frame = frameValue;
        end
    end

    source = emptySource();
    source.kind = 'alignment_csv';
    source.label = localFileStem(csvFile);
    source.file = tableText(T.source_file, focusIdx);
    source.frame = frame;
    source.focusIndex = focusIdx;
    source.focusZUm = zUm(focusIdx);
    source.signalTrace = signal;
    source.planeZ = zUm;
    source.channelIDs = opts.channelIDs(:);
end

function source = resolveFocusFromMat(matFile, opts)
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
    focusIdx = chooseFocusIndex(signal, opts);
    [fileName, frame] = sourceFromFrameMeta(S, focusIdx);

    source = emptySource();
    source.kind = 'mat_stack';
    source.label = localFileStem(matFile);
    source.file = fileName;
    source.frame = frame;
    source.detectorStack = detStack;
    source.focusIndex = focusIdx;
    source.focusZUm = planeZ(focusIdx);
    source.signalTrace = signal(:).';
    source.planeZ = planeZ(:).';
    source.channelIDs = opts.channelIDs(:);
end

function rawFocus = loadFocusDetectorPlane(source, opts)
    if ~isempty(source.detectorStack)
        rawFocus = source.detectorStack(:,:,:,source.focusIndex);
        return;
    end
    if isempty(source.file)
        error('plotFocalPlaneDetectorIntensity:NeedDetectorSource', ...
            'The selected focal plane has no detector-resolved source file.');
    end
    rawFocus = readPtuDetectorStack(source.file, opts, source.frame);
end

function corrected = correctFocusPlaneForDisplay(rawFocus, opts)
    corrected = double(rawFocus);
    if ~opts.subtractBoundary
        return;
    end

    nCh = size(corrected, 3);
    bg = zeros(1, 1, nCh);
    w = max(1, min(round(opts.boundaryWidthPx), floor(min(size(corrected,1), size(corrected,2))/2)));
    for ch = 1:nCh
        img = corrected(:,:,ch);
        border = [reshape(img(1:w,:), [], 1); ...
                  reshape(img(end-w+1:end,:), [], 1); ...
                  reshape(img(:,1:w), [], 1); ...
                  reshape(img(:,end-w+1:end), [], 1)];
        bg(1,1,ch) = median(border(isfinite(border)));
    end
    corrected = max(corrected - bg, 0);
end

function [centerXY, centerYX] = estimateBeadCenter(img, opts)
    img = double(img);
    if ~isempty(opts.centerXY)
        centerXY = double(opts.centerXY(:)).';
        if numel(centerXY) ~= 2
            error('plotFocalPlaneDetectorIntensity:BadCenterXY', ...
                'centerXY must be [x y].');
        end
    else
        switch opts.centerMode
            case 'peak'
                [~, idx] = max(img(:));
                [cy, cx] = ind2sub(size(img), idx);
                centerXY = [cx, cy];
            case {'centroid','weighted'}
                positive = max(img - min(img(:)), 0);
                threshold = opts.centerThresholdFraction * max(positive(:));
                weights = positive;
                weights(weights < threshold) = 0;
                if sum(weights(:)) <= 0
                    [~, idx] = max(img(:));
                    [cy, cx] = ind2sub(size(img), idx);
                    centerXY = [cx, cy];
                else
                    [yy, xx] = ndgrid(1:size(img,1), 1:size(img,2));
                    centerXY = [sum(xx(:).*weights(:)), sum(yy(:).*weights(:))] / sum(weights(:));
                end
            otherwise
                error('plotFocalPlaneDetectorIntensity:BadCenterMode', ...
                    'centerMode must be ''centroid'' or ''peak''.');
        end
    end

    ix = min(max(round(centerXY(1)), 1), size(img, 2));
    iy = min(max(round(centerXY(2)), 1), size(img, 1));
    centerYX = [iy, ix];
end

function scanGrid = sampleScanGridDetectorValues(focusStack, centerXY, opts)
    nGrid = opts.scanGridSize;
    if isempty(nGrid) || nGrid == 0
        scanGrid = [];
        return;
    end

    offsets = ((1:nGrid) - (nGrid+1)/2) * opts.scanStepPx;
    [dxGrid, dyGrid] = meshgrid(offsets, offsets);
    positionsXY = [centerXY(1) + reshape(dxGrid.', [], 1), ...
        centerXY(2) + reshape(dyGrid.', [], 1)];
    scanRowCol = zeros(numel(dxGrid), 2);
    scanRowCol(:,1) = reshape(repmat(1:nGrid, nGrid, 1), [], 1);
    scanRowCol(:,2) = reshape(repmat((1:nGrid).', 1, nGrid), [], 1);

    values = interpolateDetectorStack(focusStack, positionsXY);
    scanTotal = sum(values, 2);
    denom = scanTotal;
    bad = ~isfinite(denom) | denom <= 0;
    denom(bad) = 1;
    normalizedValues = bsxfun(@rdivide, values, denom);
    normalizedValues(bad,:) = NaN;

    scanGrid = struct();
    scanGrid.gridSize = nGrid;
    scanGrid.stepPx = opts.scanStepPx;
    scanGrid.offsetAxisPx = offsets(:);
    scanGrid.offsetXY = [reshape(dxGrid.', [], 1), reshape(dyGrid.', [], 1)];
    scanGrid.positionXY = positionsXY;
    scanGrid.scanRowCol = scanRowCol;
    scanGrid.detectorValues = values;
    scanGrid.normalizedDetectorValues = normalizedValues;
    scanGrid.totalIntensity = scanTotal;
end

function values = interpolateDetectorStack(focusStack, positionsXY)
    focusStack = double(focusStack);
    [ny, nx, nCh] = size(focusStack);
    [xx, yy] = meshgrid(1:nx, 1:ny);
    nPos = size(positionsXY, 1);
    values = nan(nPos, nCh);

    for ch = 1:nCh
        img = focusStack(:,:,ch);
        values(:,ch) = interp2(xx, yy, img, ...
            positionsXY(:,1), positionsXY(:,2), 'linear', NaN);
    end

    values(values < 0) = 0;
end

function nGrid = normalizeScanGridSize(value)
    if isempty(value)
        nGrid = [];
        return;
    end

    nGrid = round(double(value(1)));
    if nGrid == 0
        return;
    end
    if ~ismember(nGrid, [5 10])
        error('plotFocalPlaneDetectorIntensity:BadScanGridSize', ...
            'scanGridSize must be 5 or 10. Use 0 or [] to disable the scan grid.');
    end
end

function value = normalizePositiveScalar(value, name)
    value = double(value(1));
    if ~isfinite(value) || value <= 0
        error('plotFocalPlaneDetectorIntensity:BadPositiveScalar', ...
            '%s must be a finite positive scalar.', name);
    end
end

function writeOutputs(out, opts)
    outDir = out.outputDir;
    if isempty(outDir)
        return;
    end
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    writeSummaryFigure(out, opts, fullfile(outDir, 'focal_plane_detector_intensity_summary.png'));
    writeFocalPlaneImage(out, fullfile(outDir, 'focal_plane_bead_image.png'));
    writeDetectorDistributionFigure(out, opts, fullfile(outDir, 'center_detector_distribution.png'));
    if ~isempty(out.scanGrid)
        writeScanGridOverlayFigure(out, fullfile(outDir, 'scan_position_grid_on_bead.png'));
        writeScanGridDetectorFigure(out, opts, fullfile(outDir, 'scan_position_detector_hexmaps.png'));
        writeScanGridDetectorTable(out, fullfile(outDir, 'scan_position_detector_distribution.csv'));
        writeScanGridTotalTable(out, fullfile(outDir, 'scan_position_total_intensity.csv'));
    end

    T = table((1:numel(out.detectorValues)).', out.channelIDs(:), ...
        out.detectorValues(:), out.normalizedDetectorValues(:), ...
        'VariableNames', {'detectorIndex','channelID','intensity','normalizedIntensity'});
    writetable(T, fullfile(outDir, 'center_detector_distribution.csv'));
    save(fullfile(outDir, 'focal_plane_detector_intensity.mat'), 'out', '-v7.3');
end

function writeSummaryFigure(out, opts, outFile)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1350 430]);
    tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile(tl, 1);
    imagesc(ax1, out.focusImage);
    axis(ax1, 'image');
    colormap(ax1, 'hot');
    colorbar(ax1);
    hold(ax1, 'on');
    plot(ax1, out.centerXY(1), out.centerXY(2), 'co', 'MarkerSize', 8, 'LineWidth', 1.2);
    hold(ax1, 'off');
    title(ax1, sprintf('Focal bead image, z = %.3f um', out.source.focusZUm));
    xlabel(ax1, 'x pixel');
    ylabel(ax1, 'y pixel');

    ax2 = nexttile(tl, 2);
    plotDetectorMap(ax2, out.normalizedDetectorValues, opts);
    title(ax2, 'Detector distribution at bead centre');
    cb = colorbar(ax2);
    setColorbarLabel(cb, 'fraction of centre signal');

    ax3 = nexttile(tl, 3);
    bar(ax3, out.channelIDs, out.normalizedDetectorValues, 0.85);
    xlabel(ax3, 'PTU channel ID');
    ylabel(ax3, 'fraction of centre signal');
    title(ax3, 'Centre-pixel detector intensities');
    grid(ax3, 'on');

    exportFigure(fig, outFile);
end

function writeFocalPlaneImage(out, outFile)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 560 500]);
    ax = axes(fig);
    imagesc(ax, out.focusImage);
    axis(ax, 'image');
    colormap(ax, 'hot');
    colorbar(ax);
    hold(ax, 'on');
    plot(ax, out.centerXY(1), out.centerXY(2), 'co', 'MarkerSize', 8, 'LineWidth', 1.2);
    hold(ax, 'off');
    title(ax, sprintf('Focal bead image, z = %.3f um', out.source.focusZUm));
    xlabel(ax, 'x pixel');
    ylabel(ax, 'y pixel');
    exportFigure(fig, outFile);
end

function writeDetectorDistributionFigure(out, opts, outFile)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 560 500]);
    ax = axes(fig);
    plotDetectorMap(ax, out.normalizedDetectorValues, opts);
    title(ax, '23-pixel detector distribution at bead centre');
    cb = colorbar(ax);
    setColorbarLabel(cb, 'fraction of centre signal');
    exportFigure(fig, outFile);
end

function writeScanGridOverlayFigure(out, outFile)
    sg = out.scanGrid;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 620 560]);
    ax = axes(fig);
    imagesc(ax, out.focusImage);
    axis(ax, 'image');
    colormap(ax, 'hot');
    colorbar(ax);
    hold(ax, 'on');
    scatter(ax, sg.positionXY(:,1), sg.positionXY(:,2), 34, 'c', 'filled', ...
        'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    plot(ax, out.centerXY(1), out.centerXY(2), 'wo', ...
        'MarkerSize', 9, 'LineWidth', 1.2);
    hold(ax, 'off');
    title(ax, sprintf('%dx%d scan positions at focal plane, z = %.3f um', ...
        sg.gridSize, sg.gridSize, out.source.focusZUm));
    xlabel(ax, 'x pixel');
    ylabel(ax, 'y pixel');
    exportFigure(fig, outFile);
end

function writeScanGridDetectorFigure(out, opts, outFile)
    sg = out.scanGrid;
    nGrid = sg.gridSize;
    values = sg.normalizedDetectorValues;
    vmax = max(values(isfinite(values)));
    if isempty(vmax) || ~isfinite(vmax) || vmax <= 0
        vmax = 1;
    end
    clim = [0 vmax];

    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', scanGridFigurePosition(nGrid));
    tl = tiledlayout(fig, nGrid, nGrid, ...
        'Padding', 'compact', 'TileSpacing', 'compact');

    lastAx = [];
    for pos = 1:size(values, 1)
        ax = nexttile(tl);
        plotDetectorMap(ax, values(pos,:), opts, clim);
        labelScanGridTile(ax, sg, pos);
        lastAx = ax;
    end

    if ~isempty(lastAx)
        cb = colorbar(lastAx);
        cb.Layout.Tile = 'east';
        setColorbarLabel(cb, 'fraction of local scan-position signal');
    end

    title(tl, sprintf('%dx%d scan-position detector distributions, step %.3g px', ...
        nGrid, nGrid, sg.stepPx), 'FontWeight', 'bold');
    exportFigure(fig, outFile);
end

function writeScanGridDetectorTable(out, outFile)
    sg = out.scanGrid;
    nPos = size(sg.detectorValues, 1);
    nDet = numel(out.channelIDs);
    perDetector = ones(nDet, 1);

    scanIndex = kron((1:nPos).', perDetector);
    scanRow = kron(sg.scanRowCol(:,1), perDetector);
    scanCol = kron(sg.scanRowCol(:,2), perDetector);
    xPixel = kron(sg.positionXY(:,1), perDetector);
    yPixel = kron(sg.positionXY(:,2), perDetector);
    dxPixel = kron(sg.offsetXY(:,1), perDetector);
    dyPixel = kron(sg.offsetXY(:,2), perDetector);
    detectorIndex = repmat((1:nDet).', nPos, 1);
    channelID = repmat(out.channelIDs(:), nPos, 1);
    intensity = reshape(sg.detectorValues.', [], 1);
    normalizedIntensity = reshape(sg.normalizedDetectorValues.', [], 1);
    scanTotalIntensity = kron(sg.totalIntensity(:), perDetector);

    T = table(scanIndex, scanRow, scanCol, xPixel, yPixel, ...
        dxPixel, dyPixel, detectorIndex, channelID, intensity, ...
        normalizedIntensity, scanTotalIntensity);
    writetable(T, outFile);
end

function writeScanGridTotalTable(out, outFile)
    sg = out.scanGrid;
    scanIndex = (1:size(sg.positionXY, 1)).';
    scanRow = sg.scanRowCol(:,1);
    scanCol = sg.scanRowCol(:,2);
    xPixel = sg.positionXY(:,1);
    yPixel = sg.positionXY(:,2);
    dxPixel = sg.offsetXY(:,1);
    dyPixel = sg.offsetXY(:,2);
    totalIntensity = sg.totalIntensity(:);

    T = table(scanIndex, scanRow, scanCol, xPixel, yPixel, ...
        dxPixel, dyPixel, totalIntensity);
    writetable(T, outFile);
end

function labelScanGridTile(ax, sg, pos)
    set(ax, 'Visible', 'on', 'Box', 'off', 'Color', 'none', ...
        'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', []);

    row = sg.scanRowCol(pos, 1);
    col = sg.scanRowCol(pos, 2);
    offset = sg.offsetXY(pos, :);
    if row == 1
        title(ax, sprintf('x %+0.1f', offset(1)), ...
            'FontSize', 7, 'FontWeight', 'normal');
    end
    if col == 1
        ylabel(ax, sprintf('y %+0.1f', offset(2)), ...
            'Rotation', 0, 'Color', 'k', 'FontSize', 7, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
end

function pos = scanGridFigurePosition(nGrid)
    panelW = 120;
    panelH = 105;
    if nGrid <= 5
        panelW = 145;
        panelH = 125;
    end
    pos = [80 80 panelW*nGrid + 190 panelH*nGrid + 125];
end

function setColorbarLabel(cb, labelText)
    try
        drawnow;
        cb.Label.String = labelText;
    catch
        try
            ylabel(cb, labelText);
        catch
        end
    end
end

function plotDetectorMap(ax, values, opts, clim)
    if nargin < 4
        clim = [];
    end
    [detXY, ~] = detectorLayout(opts.detectorLayout, 1);
    if exist('plotDetectorHexMap', 'file') == 2
        plotDetectorHexMap(detXY, values, 'Parent', ax, 'CLim', clim);
    else
        scatter(ax, detXY(:,1), detXY(:,2), 320, values, 'filled');
        if ~isempty(clim)
            caxis(ax, clim);
        end
        axis(ax, 'equal');
        axis(ax, 'off');
    end
    colormap(ax, 'parula');
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
                error('plotFocalPlaneDetectorIntensity:BadCellStack', ...
                    'All detector cell volumes must have the same [y x z] size.');
            end
            detStack(:,:,k,:) = reshape(vol, ny, nx, 1, nz);
        end
        return;
    end

    value = double(value);
    if ndims(value) ~= 4
        error('plotFocalPlaneDetectorIntensity:BadDetectorStack', ...
            'Detector-resolved stack must be [y x 23 z] or [y x z 23].');
    end
    if size(value, 3) == numel(channelIDs)
        detStack = value;
    elseif size(value, 4) == numel(channelIDs)
        detStack = permute(value, [1 2 4 3]);
    elseif size(value, 3) == 23
        detStack = value;
    elseif size(value, 4) == 23
        detStack = permute(value, [1 2 4 3]);
    else
        error('plotFocalPlaneDetectorIntensity:BadDetectorStack', ...
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

    error('plotFocalPlaneDetectorIntensity:NoStackVariable', ...
        'No detector-resolved 4-D stack or scalar 3-D volume was found.');
end

function tf = isDetectorStack(value)
    tf = isnumeric(value) && ndims(value) == 4 && ...
        (size(value, 3) == 23 || size(value, 4) == 23);
    if iscell(value) && numel(value) == 23
        tf = true;
    end
end

function stack = readPtuDetectorStack(fileName, opts, frame)
    hasFastReader = exist('PTU_MultiFrameScanReadFast', 'file') == 2;
    hasSlowReader = exist('PTU_MultiFrameScanRead', 'file') == 2;
    if ~hasFastReader && ~hasSlowReader
        error('plotFocalPlaneDetectorIntensity:MissingPtuReader', ...
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
            error('plotFocalPlaneDetectorIntensity:PtuReadFailed', ...
                'Could not read %s as a detector scan. Fast: %s Slow: slow reader not available', ...
                fileName, fastMessage);
        end
        try
            ptuOut = PTU_MultiFrameScanRead(fileName, opts.ptuPhotonsPerChunk);
        catch slowErr
            error('plotFocalPlaneDetectorIntensity:PtuReadFailed', ...
                'Could not read %s as a detector scan. Fast: %s Slow: %s', ...
                fileName, fastMessage, slowErr.message);
        end
    end

    if isfield(ptuOut, 'dind')
        channelIDs = double(ptuOut.dind(:));
    else
        channelIDs = [];
    end

    if isfield(ptuOut, 'tag') && ~isempty(ptuOut.tag)
        if isempty(frame)
            raw = sum(double(ptuOut.tag), 4);
        else
            frame = round(frame(1));
            if frame < 1 || frame > size(ptuOut.tag, 4)
                error('plotFocalPlaneDetectorIntensity:BadFrame', ...
                    'Frame index must be within 1:%d for %s.', size(ptuOut.tag, 4), fileName);
            end
            raw = double(ptuOut.tag(:,:,:,frame));
        end
        stack = permute(raw, [2 1 3]);
    elseif isfield(ptuOut, 'tags') && ~isempty(ptuOut.tags)
        stack = permute(double(ptuOut.tags), [2 1 3]);
    else
        error('plotFocalPlaneDetectorIntensity:NoPtuDetectorData', ...
            'No detector image stack was found in %s.', fileName);
    end

    [stack, ~] = selectChannels(stack, channelIDs, opts.channelIDs);
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
        error('plotFocalPlaneDetectorIntensity:MissingChannels', ...
            'Could not find all requested detector channel IDs in the stack.');
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
            error('plotFocalPlaneDetectorIntensity:BadFocusMetric', ...
                'focusMetric must be ''signal'' or ''peak''.');
    end
end

function focusIdx = chooseFocusIndex(signal, opts)
    if ~isempty(opts.focusIndex)
        focusIdx = round(opts.focusIndex);
        if focusIdx < 1 || focusIdx > numel(signal)
            error('plotFocalPlaneDetectorIntensity:BadFocusIndex', ...
                'focusIndex must be between 1 and %d.', numel(signal));
        end
    else
        [~, focusIdx] = max(signal);
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
        error('plotFocalPlaneDetectorIntensity:PlaneZMismatch', ...
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

function csvFile = findBatchAlignmentCsvForScanFolder(scanFolder)
    csvFile = '';
    [dataRoot, scanName] = fileparts(stripTrailingFilesep(scanFolder));
    [dataParent, datasetName] = fileparts(dataRoot);
    plotDir = fullfile(dataParent, 'PSF_batch_outputs', sanitizeFileName(datasetName), 'xz_yz_plots');
    if exist(plotDir, 'dir') ~= 7
        return;
    end
    stem = scanOutputStem(scanName);
    candidate = fullfile(plotDir, sprintf('%s_frame_alignment.csv', stem));
    if exist(candidate, 'file') == 2
        csvFile = candidate;
    end
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
    hits = dir(fullfile(volumeDir, sprintf('%s_volume_raw.mat', stem)));
    if isempty(hits)
        hits = dir(fullfile(volumeDir, sprintf('*%s*volume*.mat', stem)));
    end
    if isempty(hits)
        return;
    end
    [~, newest] = max([hits.datenum]);
    matFile = fullfile(hits(newest).folder, hits(newest).name);
end

function ids = resolveChannelIDs(requested, nCh)
    ids = requested(:);
    if numel(ids) ~= nCh
        ids = (1:nCh).';
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

function outDir = resolveOutputDir(source, opts)
    if ~isempty(opts.outputDir)
        outDir = char(opts.outputDir);
        return;
    end
    thisDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(thisDir);
    label = source.label;
    if isempty(label) && ~isempty(source.scanFolder)
        label = localFileStem(source.scanFolder);
    end
    if isempty(label)
        label = 'focal_plane';
    end
    outDir = fullfile(rootDir, 'output_matlab', 'focal_plane_detector_intensity', sanitizeFileName(label));
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
