function out = plotIdealBeadScanDetectorIntensity(varargin)
%PLOTIDEALBEADSCANDETECTORINTENSITY Ideal bead scan maps for aberration modes.
%
%   out = plotIdealBeadScanDetectorIntensity()
%
%   Generates ideal detector-resolved scans of a fluorescent bead for no
%   aberration plus all supported single Zernike aberration modes. Each
%   scenario samples a 10x10 lateral scan-position grid around the bead
%   centre with a 20 nm step and writes a separate output folder.
%
%   Example:
%       out = plotIdealBeadScanDetectorIntensity();
%
%       out = plotIdealBeadScanDetectorIntensity('scanGridSize', 10, ...
%           'scanStepUm', 0.020, 'amplitudeWaves', 0.15);
%
%       out = plotIdealBeadScanDetectorIntensity('modes', {'coma_x'});

    addpath(fileparts(mfilename('fullpath')));
    opts = parseInputs(varargin{:});

    sim = buildIdealSimulation(opts);
    scenarios = buildScenarios(opts);
    outputRoot = resolveOutputRoot(opts);

    out = struct();
    out.sim = sim;
    out.options = opts;
    out.detXY = sim.detXY;
    out.baselineExpectedDetectorScanShiftUm = 0.5 * sim.detXY;
    out.channelIDs = (1:sim.nDet).';
    out.outputRoot = outputRoot;
    out.scenarios = repmat(emptyScenarioOutput(), 0, 1);

    for s = 1:numel(scenarios)
        scenario = scenarios(s);
        scanStack = normalizedStackExplicitDetector(sim, scenario.coeffs, 0, 0);
        scanGrid = sampleScanGrid(scanStack, sim, opts);
        scenarioOut = buildScenarioOutput(sim, opts, scenario, scanStack, scanGrid, outputRoot);
        out.scenarios(end+1, 1) = scenarioOut; %#ok<AGROW>

        if opts.writeOutputs
            writeOutputs(scenarioOut, opts);
        end
    end

    if opts.writeOutputs
        if exist(outputRoot, 'dir') ~= 7
            mkdir(outputRoot);
        end
        writeScenarioSummary(out, fullfile(outputRoot, 'scenario_summary.csv'));
        save(fullfile(outputRoot, 'ideal_bead_scan_all_scenarios.mat'), 'out', '-v7.3');
    end

    if opts.verbose
        fprintf('[plotIdealBeadScanDetectorIntensity] generated %d scenario folders, %dx%d positions, %.1f nm step.\n', ...
            numel(out.scenarios), opts.scanGridSize, opts.scanGridSize, 1000*opts.scanStepUm);
        fprintf('[plotIdealBeadScanDetectorIntensity] simulation grid: %dx%d, dx %.2f nm, detector layout %s.\n', ...
            sim.nx, sim.ny, 1000*sim.dx, sim.detectorLayout);
        if opts.writeOutputs
            fprintf('[plotIdealBeadScanDetectorIntensity] wrote outputs below %s\n', out.outputRoot);
        end
    end
end

function opts = parseInputs(varargin)
    p = inputParser;
    p.FunctionName = 'plotIdealBeadScanDetectorIntensity';

    addParameter(p, 'scanGridSize', 10);
    addParameter(p, 'scanStepUm', 0.020);
    addParameter(p, 'modes', 'all');
    addParameter(p, 'amplitudeWaves', 0.15);
    addParameter(p, 'includeNoAberration', true);
    addParameter(p, 'maxZernikeOrder', 6);
    addParameter(p, 'fovXY', 1.8);
    addParameter(p, 'nx', []);
    addParameter(p, 'nzRange', []);
    addParameter(p, 'nz', []);
    addParameter(p, 'beadRadiusUm', []);
    addParameter(p, 'beadSubsamples', []);
    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detPitchUm', []);
    addParameter(p, 'detFillRatio', 1.0);
    addParameter(p, 'detectorSubsamples', []);
    addParameter(p, 'colormap', 'parula');
    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;
    opts.scanGridSize = normalizePositiveInteger(opts.scanGridSize, 'scanGridSize');
    opts.scanStepUm = normalizePositiveScalar(opts.scanStepUm, 'scanStepUm');
    opts.amplitudeWaves = normalizeFiniteScalar(opts.amplitudeWaves, 'amplitudeWaves');
    opts.maxZernikeOrder = normalizePositiveInteger(opts.maxZernikeOrder, 'maxZernikeOrder');
    opts.modes = normalizeModes(opts.modes, opts.maxZernikeOrder);
    opts.includeNoAberration = logical(opts.includeNoAberration);
    opts.fovXY = normalizePositiveScalar(opts.fovXY, 'fovXY');
    opts.detFillRatio = normalizePositiveScalar(opts.detFillRatio, 'detFillRatio');
end

function sim = buildIdealSimulation(opts)
    sim = defaultParams();
    sim.sampleGeometry = 'homogeneous';

    sim.fovXY = opts.fovXY;
    sim.nx = resolveNx(opts);
    sim.ny = sim.nx;
    sim.x = linspace(-sim.fovXY/2, sim.fovXY/2, sim.nx);
    sim.y = sim.x;
    sim.dx = abs(sim.x(2)-sim.x(1));

    if ~isempty(opts.nzRange)
        sim.nzRange = normalizePositiveScalar(opts.nzRange, 'nzRange');
    end
    if ~isempty(opts.nz)
        sim.nz = normalizePositiveInteger(opts.nz, 'nz');
    end
    sim.z = linspace(-sim.nzRange/2, sim.nzRange/2, sim.nz);

    if ~isempty(opts.beadRadiusUm)
        sim.beadRadius = normalizePositiveScalar(opts.beadRadiusUm, 'beadRadiusUm');
    end
    if ~isempty(opts.beadSubsamples)
        sim.beadSubsamples = opts.beadSubsamples;
    end

    sim.detectorLayout = char(opts.detectorLayout);
    if ~isempty(opts.detPitchUm)
        sim.detPitch = normalizePositiveScalar(opts.detPitchUm, 'detPitchUm');
    end
    sim.detSize = opts.detFillRatio * sim.detPitch;
    sim.detectorPixelShape = 'hex';
    sim.detectorHexRadius = sim.detSize / sqrt(3);
    [sim.detXY, sim.detectorIndexGrid, sim.detectorLayoutInfo] = ...
        detectorLayout(sim.detectorLayout, sim.detPitch);
    sim.nDet = size(sim.detXY, 1);
    sim.detectorGridSize = size(sim.detectorIndexGrid);
    sim.arrayN = sim.detectorGridSize(1);
    if ~isempty(opts.detectorSubsamples)
        sim.detectorSubsamples = normalizePositiveInteger(opts.detectorSubsamples, 'detectorSubsamples');
    end
    sim.detectorImageInverted = true;

    sim.obj = beadObject3D(sim);
end

function nx = resolveNx(opts)
    if ~isempty(opts.nx)
        nx = normalizePositiveInteger(opts.nx, 'nx');
        if mod(nx, 2) == 0
            nx = nx + 1;
        end
        return;
    end

    nIntervals = ceil(opts.fovXY / opts.scanStepUm);
    if mod(nIntervals, 2) == 1
        nIntervals = nIntervals + 1;
    end
    nx = nIntervals + 1;
end

function scenarios = buildScenarios(opts)
    scenarios = repmat(emptyScenario(), 0, 1);

    if opts.includeNoAberration
        sc = emptyScenario();
        sc.name = 'No aberration';
        sc.modeName = 'none';
        sc.folderName = '00_no_aberration';
        sc.coeffs = struct();
        sc.amplitudeWaves = 0;
        scenarios(end+1, 1) = sc; %#ok<AGROW>
    end

    modeFolderIndex = 0;
    for k = 1:numel(opts.modes)
        modeFolderIndex = modeFolderIndex + 1;
        modeName = opts.modes{k};
        sc = emptyScenario();
        sc.name = prettyModeName(modeName);
        sc.modeName = modeName;
        sc.folderName = sanitizeFileName(sprintf('%02d_%s_%gwaves', ...
            modeFolderIndex, modeName, opts.amplitudeWaves));
        sc.coeffs = struct(modeName, opts.amplitudeWaves);
        sc.amplitudeWaves = opts.amplitudeWaves;
        scenarios(end+1, 1) = sc; %#ok<AGROW>
    end

    if isempty(scenarios)
        error('plotIdealBeadScanDetectorIntensity:NoScenarios', ...
            'No scenarios were requested. Enable includeNoAberration or provide at least one mode.');
    end
end

function sc = emptyScenario()
    sc = struct();
    sc.name = '';
    sc.modeName = '';
    sc.folderName = '';
    sc.coeffs = struct();
    sc.amplitudeWaves = 0;
end

function out = buildScenarioOutput(sim, opts, scenario, scanStack, scanGrid, outputRoot)
    out = emptyScenarioOutput();
    out.name = scenario.name;
    out.modeName = scenario.modeName;
    out.folderName = scenario.folderName;
    out.amplitudeWaves = scenario.amplitudeWaves;
    out.sim = sim;
    out.options = opts;
    out.coeffs = scenario.coeffs;
    out.scanStack = scanStack;
    out.scanGrid = scanGrid;
    out.detXY = sim.detXY;
    out.expectedDetectorScanShiftUm = 0.5 * sim.detXY;
    out.channelIDs = (1:sim.nDet).';
    out.outputRoot = outputRoot;
    out.outputDir = fullfile(outputRoot, scenario.folderName);
end

function out = emptyScenarioOutput()
    out = struct();
    out.name = '';
    out.modeName = '';
    out.folderName = '';
    out.amplitudeWaves = [];
    out.sim = [];
    out.options = [];
    out.coeffs = struct();
    out.scanStack = [];
    out.scanGrid = [];
    out.detXY = [];
    out.expectedDetectorScanShiftUm = [];
    out.channelIDs = [];
    out.outputRoot = '';
    out.outputDir = '';
end

function scanGrid = sampleScanGrid(scanStack, sim, opts)
    nGrid = opts.scanGridSize;
    offsets = ((1:nGrid) - (nGrid+1)/2) * opts.scanStepUm;
    [dxGrid, dyGrid] = meshgrid(offsets, offsets);
    offsetXY = [reshape(dxGrid.', [], 1), reshape(dyGrid.', [], 1)];

    values = interpolateDetectorStack(scanStack, sim.x, sim.y, offsetXY);
    total = sum(values, 2);
    denom = total;
    bad = ~isfinite(denom) | denom <= 0;
    denom(bad) = 1;
    normalizedValues = bsxfun(@rdivide, values, denom);
    normalizedValues(bad,:) = NaN;

    totalMap = reshape(total, nGrid, nGrid).';
    totalMax = max(totalMap(:));
    normalizedTotalMap = totalMap / max(totalMax, eps);

    scanGrid = struct();
    scanGrid.gridSize = nGrid;
    scanGrid.stepUm = opts.scanStepUm;
    scanGrid.stepNm = 1000 * opts.scanStepUm;
    scanGrid.offsetAxisUm = offsets(:);
    scanGrid.offsetXYUm = offsetXY;
    scanGrid.scanRowCol = scanRowCol(nGrid);
    scanGrid.detectorValues = values;
    scanGrid.normalizedDetectorValues = normalizedValues;
    scanGrid.totalIntensity = total;
    scanGrid.totalIntensityMap = totalMap;
    scanGrid.normalizedTotalIntensityMap = normalizedTotalMap;
end

function rc = scanRowCol(nGrid)
    rc = zeros(nGrid*nGrid, 2);
    rc(:,1) = reshape(repmat(1:nGrid, nGrid, 1), [], 1);
    rc(:,2) = reshape(repmat((1:nGrid).', 1, nGrid), [], 1);
end

function values = interpolateDetectorStack(scanStack, xAxis, yAxis, positionsXY)
    [xx, yy] = meshgrid(xAxis, yAxis);
    nPos = size(positionsXY, 1);
    nDet = size(scanStack, 3);
    values = nan(nPos, nDet);

    for k = 1:nDet
        values(:,k) = interp2(xx, yy, scanStack(:,:,k), ...
            positionsXY(:,1), positionsXY(:,2), 'linear', NaN);
    end
    values(values < 0) = 0;
end

function writeOutputs(out, opts)
    outDir = out.outputDir;
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    writeTotalScanMap(out, opts, fullfile(outDir, 'ideal_scan_total_intensity_map.png'));
    writeScanGridDetectorFigure(out, opts, fullfile(outDir, 'ideal_scan_position_detector_hexmaps.png'));
    writeDetectorResolvedScanImages(out, opts, fullfile(outDir, 'ideal_detector_resolved_scan_images.png'));
    writeDetectorResolvedScanImagesHoneycomb(out, opts, fullfile(outDir, 'ideal_detector_resolved_scan_images_honeycomb.png'));
    writeScanGridDetectorTable(out, fullfile(outDir, 'ideal_scan_position_detector_distribution.csv'));
    writeScanGridTotalTable(out, fullfile(outDir, 'ideal_scan_position_total_intensity.csv'));
    writeExpectedShiftTable(out, fullfile(outDir, 'ideal_expected_detector_scan_shift.csv'));
    save(fullfile(outDir, 'ideal_bead_scan_detector_intensity.mat'), 'out', '-v7.3');
end

function writeTotalScanMap(out, opts, outFile)
    sg = out.scanGrid;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 580 520]);
    ax = axes(fig);
    imagesc(ax, sg.offsetAxisUm*1000, sg.offsetAxisUm*1000, sg.normalizedTotalIntensityMap);
    axis(ax, 'image');
    set(ax, 'YDir', 'normal');
    colormap(ax, opts.colormap);
    cb = colorbar(ax);
    setColorbarLabel(cb, 'normalized total intensity');
    title(ax, sprintf('Ideal bead scan: %s, %.1f nm step', out.name, sg.stepNm));
    xlabel(ax, 'scan x offset [nm]');
    ylabel(ax, 'scan y offset [nm]');
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

    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', scanGridFigurePosition(nGrid));
    tl = tiledlayout(fig, nGrid, nGrid, ...
        'Padding', 'compact', 'TileSpacing', 'compact');

    lastAx = [];
    for pos = 1:size(values, 1)
        ax = nexttile(tl);
        plotDetectorMap(ax, out.detXY, values(pos,:), opts, [0 vmax]);
        labelScanGridTile(ax, sg, pos);
        lastAx = ax;
    end

    if ~isempty(lastAx)
        cb = colorbar(lastAx);
        cb.Layout.Tile = 'east';
        setColorbarLabel(cb, 'fraction of local scan-position signal');
    end

    title(tl, sprintf('Ideal detector distributions: %s, %dx%d scan, %.1f nm step', ...
        out.name, nGrid, nGrid, sg.stepNm), 'FontWeight', 'bold');
    exportFigure(fig, outFile);
end

function writeDetectorResolvedScanImages(out, opts, outFile)
    sg = out.scanGrid;
    nGrid = sg.gridSize;
    nDet = size(sg.detectorValues, 2);
    values = sg.detectorValues;
    vmax = max(values(:));
    if ~isfinite(vmax) || vmax <= 0
        vmax = 1;
    end

    nCols = 5;
    nRows = ceil(nDet / nCols);
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [80 80 980 760]);
    tl = tiledlayout(fig, nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');

    lastAx = [];
    for k = 1:nDet
        ax = nexttile(tl);
        map = reshape(values(:,k), nGrid, nGrid).';
        imagesc(ax, sg.offsetAxisUm*1000, sg.offsetAxisUm*1000, map, [0 vmax]);
        axis(ax, 'image');
        set(ax, 'YDir', 'normal');
        colormap(ax, opts.colormap);
        title(ax, sprintf('detector %02d', k), 'FontSize', 8, 'FontWeight', 'normal');
        set(ax, 'FontSize', 7);
        lastAx = ax;
    end

    if ~isempty(lastAx)
        cb = colorbar(lastAx);
        cb.Layout.Tile = 'east';
        setColorbarLabel(cb, 'ideal intensity');
    end
    title(tl, sprintf('Ideal detector-resolved scan images: %s', out.name), ...
        'FontWeight', 'bold');
    xlabel(tl, 'scan x offset [nm]');
    ylabel(tl, 'scan y offset [nm]');
    exportFigure(fig, outFile);
end

function writeDetectorResolvedScanImagesHoneycomb(out, opts, outFile)
    sg = out.scanGrid;
    idxGrid = out.sim.detectorIndexGrid;
    nGrid = sg.gridSize;
    values = sg.detectorValues;
    vmax = max(values(:));
    if ~isfinite(vmax) || vmax <= 0
        vmax = 1;
    end

    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [80 80 1260 760]);
    tl = tiledlayout(fig, size(idxGrid,1), size(idxGrid,2), ...
        'Padding', 'compact', 'TileSpacing', 'compact');

    scanAxisNm = 1000 * sg.offsetAxisUm;
    scanLimits = [min(scanAxisNm), max(scanAxisNm)];
    lastAx = [];

    for r = 1:size(idxGrid,1)
        for c = 1:size(idxGrid,2)
            ax = nexttile(tl);
            k = idxGrid(r,c);
            if ~isfinite(k)
                axis(ax, 'off');
                continue;
            end

            k = round(k);
            map = reshape(values(:,k), nGrid, nGrid).';
            imagesc(ax, scanAxisNm, scanAxisNm, map, [0 vmax]);
            axis(ax, 'image');
            set(ax, 'YDir', 'normal', 'XTick', [], 'YTick', []);
            xlim(ax, scanLimits);
            ylim(ax, scanLimits);
            colormap(ax, opts.colormap);
            title(ax, sprintf('%02d', k), 'FontSize', 8, 'FontWeight', 'normal');
            markExpectedShiftIfVisible(ax, out.expectedDetectorScanShiftUm(k,:), scanLimits);
            lastAx = ax;
        end
    end

    if ~isempty(lastAx)
        cb = colorbar(lastAx);
        cb.Layout.Tile = 'east';
        setColorbarLabel(cb, 'ideal intensity');
    end
    title(tl, sprintf(['Ideal detector-resolved scan images in honeycomb layout: %s ' ...
        '(white circle = no-aberration half-pitch peak if inside %.0f nm crop)'], ...
        out.name, max(abs(scanLimits))), ...
        'FontWeight', 'bold');
    exportFigure(fig, outFile);
end

function markExpectedShiftIfVisible(ax, shiftUm, scanLimitsNm)
    shiftNm = 1000 * shiftUm;
    if shiftNm(1) < scanLimitsNm(1) || shiftNm(1) > scanLimitsNm(2) || ...
            shiftNm(2) < scanLimitsNm(1) || shiftNm(2) > scanLimitsNm(2)
        return;
    end

    hold(ax, 'on');
    plot(ax, shiftNm(1), shiftNm(2), 'wo', ...
        'MarkerSize', 4, 'LineWidth', 0.8);
    hold(ax, 'off');
end

function writeScanGridDetectorTable(out, outFile)
    sg = out.scanGrid;
    nPos = size(sg.detectorValues, 1);
    nDet = size(sg.detectorValues, 2);
    perDetector = ones(nDet, 1);

    scanIndex = kron((1:nPos).', perDetector);
    scanRow = kron(sg.scanRowCol(:,1), perDetector);
    scanCol = kron(sg.scanRowCol(:,2), perDetector);
    scanXOffsetNm = kron(1000*sg.offsetXYUm(:,1), perDetector);
    scanYOffsetNm = kron(1000*sg.offsetXYUm(:,2), perDetector);
    detectorIndex = repmat((1:nDet).', nPos, 1);
    intensity = reshape(sg.detectorValues.', [], 1);
    normalizedIntensity = reshape(sg.normalizedDetectorValues.', [], 1);
    scanTotalIntensity = kron(sg.totalIntensity(:), perDetector);

    T = table(scanIndex, scanRow, scanCol, scanXOffsetNm, scanYOffsetNm, ...
        detectorIndex, intensity, normalizedIntensity, scanTotalIntensity);
    writetable(T, outFile);
end

function writeScanGridTotalTable(out, outFile)
    sg = out.scanGrid;
    scanIndex = (1:size(sg.offsetXYUm, 1)).';
    scanRow = sg.scanRowCol(:,1);
    scanCol = sg.scanRowCol(:,2);
    scanXOffsetNm = 1000 * sg.offsetXYUm(:,1);
    scanYOffsetNm = 1000 * sg.offsetXYUm(:,2);
    totalIntensity = sg.totalIntensity(:);
    totalMax = max(totalIntensity(:));
    normalizedTotalIntensity = totalIntensity / max(totalMax, eps);

    T = table(scanIndex, scanRow, scanCol, scanXOffsetNm, scanYOffsetNm, ...
        totalIntensity, normalizedTotalIntensity);
    writetable(T, outFile);
end

function writeExpectedShiftTable(out, outFile)
    detectorIndex = (1:size(out.detXY, 1)).';
    detectorXNm = 1000 * out.detXY(:,1);
    detectorYNm = 1000 * out.detXY(:,2);
    expectedScanShiftXNm = 1000 * out.expectedDetectorScanShiftUm(:,1);
    expectedScanShiftYNm = 1000 * out.expectedDetectorScanShiftUm(:,2);

    T = table(detectorIndex, detectorXNm, detectorYNm, ...
        expectedScanShiftXNm, expectedScanShiftYNm);
    writetable(T, outFile);
end

function writeScenarioSummary(out, outFile)
    n = numel(out.scenarios);
    scenarioIndex = (1:n).';
    name = cell(n, 1);
    modeName = cell(n, 1);
    amplitudeWaves = zeros(n, 1);
    outputDir = cell(n, 1);

    for k = 1:n
        name{k} = out.scenarios(k).name;
        modeName{k} = out.scenarios(k).modeName;
        amplitudeWaves(k) = out.scenarios(k).amplitudeWaves;
        outputDir{k} = out.scenarios(k).outputDir;
    end

    T = table(scenarioIndex, name, modeName, amplitudeWaves, outputDir);
    writetable(T, outFile);
end

function labelScanGridTile(ax, sg, pos)
    set(ax, 'Visible', 'on', 'Box', 'off', 'Color', 'none', ...
        'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', []);

    row = sg.scanRowCol(pos, 1);
    col = sg.scanRowCol(pos, 2);
    offsetNm = 1000 * sg.offsetXYUm(pos, :);
    if row == 1
        title(ax, sprintf('x %+0.0f', offsetNm(1)), ...
            'FontSize', 7, 'FontWeight', 'normal');
    end
    if col == 1
        ylabel(ax, sprintf('y %+0.0f', offsetNm(2)), ...
            'Rotation', 0, 'Color', 'k', 'FontSize', 7, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
end

function plotDetectorMap(ax, detXY, values, opts, clim)
    if exist('plotDetectorHexMap', 'file') == 2
        plotDetectorHexMap(detXY, values, 'Parent', ax, 'CLim', clim);
    else
        scatter(ax, detXY(:,1), detXY(:,2), 320, values, 'filled');
        caxis(ax, clim);
        axis(ax, 'equal');
        axis(ax, 'off');
    end
    colormap(ax, opts.colormap);
end

function pos = scanGridFigurePosition(nGrid)
    pos = [80 80 120*nGrid + 190 105*nGrid + 125];
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

function exportFigure(fig, outFile)
    try
        exportgraphics(fig, outFile, 'Resolution', 180);
    catch
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, outFile, '-dpng', '-r180');
    end
    close(fig);
end

function outDir = resolveOutputRoot(opts)
    if ~isempty(opts.outputDir)
        outDir = char(opts.outputDir);
        return;
    end

    rootDir = fileparts(fileparts(mfilename('fullpath')));
    label = sprintf('all_aberrations_%dx%d_%gnm_step_%gwaves', ...
        opts.scanGridSize, opts.scanGridSize, 1000*opts.scanStepUm, opts.amplitudeWaves);
    outDir = fullfile(rootDir, 'output_matlab', ...
        'ideal_bead_scan_detector_intensity', sanitizeFileName(label));
end

function modes = normalizeModes(value, maxOrder)
    validModes = zernikeModeNames(maxOrder);

    if isstring(value) && numel(value) > 1
        modes = cellstr(value(:).');
    elseif ischar(value) || isstring(value)
        text = strtrim(char(value));
        key = lower(strrep(strrep(text, '-', '_'), ' ', '_'));
        switch key
            case {'all', '*'}
                modes = validModes;
                return;
            case {'none', 'no_aberration', 'noaberration', 'baseline'}
                modes = {};
                return;
            otherwise
                modes = cellstr(value);
        end
    elseif iscell(value)
        modes = value(:).';
    else
        error('plotIdealBeadScanDetectorIntensity:BadModes', ...
            'modes must be ''all'', ''none'', a string mode name, or a cell array of mode names.');
    end

    modes = cellfun(@char, modes, 'UniformOutput', false);
    modes = cellfun(@strtrim, modes, 'UniformOutput', false);
    modes = modes(~cellfun(@isempty, modes));

    [tf, loc] = ismember(lower(modes), lower(validModes));
    if ~all(tf)
        bad = strjoin(modes(~tf), ', ');
        error('plotIdealBeadScanDetectorIntensity:UnknownMode', ...
            'Unknown Zernike mode(s): %s', bad);
    end
    modes = validModes(loc);
end

function name = prettyModeName(modeName)
    name = strrep(char(modeName), '_', ' ');
    parts = regexp(name, '\s+', 'split');
    for k = 1:numel(parts)
        if ~isempty(parts{k})
            parts{k} = [upper(parts{k}(1)) parts{k}(2:end)];
        end
    end
    name = strjoin(parts, ' ');
end

function value = normalizePositiveScalar(value, name)
    value = double(value(1));
    if ~isfinite(value) || value <= 0
        error('plotIdealBeadScanDetectorIntensity:BadPositiveScalar', ...
            '%s must be a finite positive scalar.', name);
    end
end

function value = normalizeFiniteScalar(value, name)
    value = double(value(1));
    if ~isfinite(value)
        error('plotIdealBeadScanDetectorIntensity:BadFiniteScalar', ...
            '%s must be a finite scalar.', name);
    end
end

function value = normalizePositiveInteger(value, name)
    value = round(double(value(1)));
    if ~isfinite(value) || value <= 0
        error('plotIdealBeadScanDetectorIntensity:BadPositiveInteger', ...
            '%s must be a finite positive integer.', name);
    end
end

function stem = sanitizeFileName(name)
    stem = regexprep(char(name), '[^A-Za-z0-9._-]+', '_');
    stem = regexprep(stem, '_+', '_');
    stem = regexprep(stem, '^_+|_+$', '');
    if isempty(stem)
        stem = 'ideal_scan';
    end
end
