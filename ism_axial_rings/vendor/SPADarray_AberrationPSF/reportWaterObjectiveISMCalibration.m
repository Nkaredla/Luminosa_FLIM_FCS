function report = reportWaterObjectiveISMCalibration(calibrationFile, varargin)
%REPORTWATEROBJECTIVEISMCALIBRATION Print wavefront and export detector PSFs.
%
%   report = reportWaterObjectiveISMCalibration(calibrationFile)
%
%   Loads a calibration MAT produced by calibrateWaterObjectiveISMFromBead3D,
%   prints the fitted wavefront coefficients, and writes XY/XZ/YZ projection
%   figures for the measured and fitted 3-D PSFs of all 23 detectors.
%
%   The function does not rerun the fit; it uses finalResult saved in the
%   calibration MAT.

    if nargin < 1 || isempty(calibrationFile)
        calibrationFile = defaultCalibrationFile();
    end

    p = inputParser;
    p.FunctionName = 'reportWaterObjectiveISMCalibration';
    addParameter(p, 'outputDir', '');
    addParameter(p, 'projectionMode', 'max');
    addParameter(p, 'figureVisible', 'off');
    addParameter(p, 'normalizeEachDetector', true);
    addParameter(p, 'writePerDetectorFigures', true);
    parse(p, varargin{:});
    opts = p.Results;
    opts.outputDir = char(opts.outputDir);
    opts.projectionMode = lower(char(opts.projectionMode));
    opts.figureVisible = char(opts.figureVisible);
    if ~ismember(opts.projectionMode, {'max','sum'})
        error('reportWaterObjectiveISMCalibration:BadProjectionMode', ...
            'projectionMode must be max or sum.');
    end

    [calibration, fitResult] = loadCalibrationAndFit(calibrationFile);
    if isempty(opts.outputDir)
        opts.outputDir = fullfile(calibration.outputDir, 'calibration_psf_report');
    end
    if exist(opts.outputDir, 'dir') ~= 7
        mkdir(opts.outputDir);
    end

    coeffs = wavefrontTable(calibration, fitResult);
    disp(coeffs);
    writetable(coeffs, fullfile(opts.outputDir, ...
        'calibration_wavefront_coefficients.csv'));

    [measured, fitted, residual, info] = calibrationVolumes(fitResult);
    if opts.normalizeEachDetector
        measuredPlot = normalizeEachDetectorVolume(measured);
        fittedPlot = normalizeEachDetectorVolume(fitted);
        residualPlot = measuredPlot - fittedPlot;
        normalizationNote = 'per-detector normalized';
    else
        measuredPlot = measured;
        fittedPlot = fitted;
        residualPlot = residual;
        normalizationNote = 'absolute fitted-count scale';
    end

    projections = {'xy','xz','yz'};
    for k = 1:numel(projections)
        planeName = projections{k};
        measuredCanvas = detectorProjectionCanvas(measuredPlot, planeName, ...
            fitResult.fit.sim, opts.projectionMode);
        fittedCanvas = detectorProjectionCanvas(fittedPlot, planeName, ...
            fitResult.fit.sim, opts.projectionMode);
        residualCanvas = detectorProjectionCanvas(residualPlot, planeName, ...
            fitResult.fit.sim, opts.projectionMode);
        outFile = fullfile(opts.outputDir, sprintf( ...
            'calibration_all_detectors_%s_psf_projection.png', upper(planeName)));
        writeAllDetectorProjectionFigure(measuredCanvas, fittedCanvas, ...
            residualCanvas, planeName, normalizationNote, opts, outFile);
    end

    if opts.writePerDetectorFigures
        perDir = fullfile(opts.outputDir, 'per_detector');
        if exist(perDir, 'dir') ~= 7
            mkdir(perDir);
        end
        for ch = 1:size(measuredPlot,4)
            outFile = fullfile(perDir, sprintf( ...
                'detector_%02d_channel_%02d_calibration_psf_xy_xz_yz.png', ...
                ch, fitResult.data.channelIDs(ch)));
            writePerDetectorFigure(measuredPlot(:,:,:,ch), ...
                fittedPlot(:,:,:,ch), residualPlot(:,:,:,ch), ...
                ch, fitResult.data.channelIDs(ch), info, normalizationNote, ...
                opts, outFile);
        end
    end

    qeTable = detectorQETable(calibration);
    writetable(qeTable, fullfile(opts.outputDir, ...
        'calibration_detector_relative_qe.csv'));

    report = struct();
    report.calibration = calibration;
    report.fitResult = fitResult;
    report.wavefrontTable = coeffs;
    report.detectorQETable = qeTable;
    report.outputDir = opts.outputDir;

    fprintf('\nCalibration PSF report written to:\n  %s\n', opts.outputDir);
    fprintf('  all-detector figures: calibration_all_detectors_XY/XZ/YZ_psf_projection.png\n');
    if opts.writePerDetectorFigures
        fprintf('  per-detector figures: %s\n', fullfile(opts.outputDir, 'per_detector'));
    end
end

function [calibration, fitResult] = loadCalibrationAndFit(fileName)
    fileName = char(fileName);
    if exist(fileName, 'file') ~= 2
        error('reportWaterObjectiveISMCalibration:MissingCalibration', ...
            'Calibration file was not found: %s', fileName);
    end
    S = load(fileName);
    if ~isfield(S, 'calibration') || ~isstruct(S.calibration)
        error('reportWaterObjectiveISMCalibration:MissingCalibrationStruct', ...
            'No calibration struct was found in %s.', fileName);
    end
    calibration = S.calibration;
    if isfield(S, 'finalResult') && isstruct(S.finalResult)
        fitResult = S.finalResult;
    elseif isfield(S, 'noQEResult') && isstruct(S.noQEResult)
        fitResult = S.noQEResult;
    else
        error('reportWaterObjectiveISMCalibration:MissingFitResult', ...
            'Calibration MAT does not contain finalResult or noQEResult.');
    end
end

function [measured, fitted, residual, info] = calibrationVolumes(result)
    raw = double(result.data.rawCounts);
    bg = double(result.data.backgroundPerPixel);
    bg = repmat(bg, size(raw,1), size(raw,2), 1, 1);
    measuredYXChZ = max(raw - bg, 0);
    fittedYXChZ = result.fit.globalPhotonScale * max(double(result.fit.model), 0);
    measured = permute(measuredYXChZ, [1 2 4 3]); % [y x z detector]
    fitted = permute(fittedYXChZ, [1 2 4 3]);
    residual = measured - fitted;
    info = struct();
    info.xUm = result.data.xUm;
    info.yUm = result.data.yUm;
    info.stageZUm = result.data.stageZUm;
    info.fitZ0Um = result.fit.estZ0Um;
    info.relativeZUm = result.data.stageZUm - result.fit.estZ0Um;
end

function Vn = normalizeEachDetectorVolume(V)
    Vn = zeros(size(V));
    for ch = 1:size(V,4)
        s = sum(V(:,:,:,ch), 'all');
        if isfinite(s) && s > 0
            Vn(:,:,:,ch) = V(:,:,:,ch) / s;
        end
    end
end

function canvas = detectorProjectionCanvas(volume, planeName, sim, mode)
    nCh = size(volume, 4);
    if isfield(sim, 'detectorIndexGrid') && nnz(isfinite(sim.detectorIndexGrid)) == nCh
        idxGrid = sim.detectorIndexGrid;
    else
        try
            [~, idxGrid] = detectorLayout('honeycomb23', 1);
        catch
            idxGrid = detectorIndexGridFallback(nCh);
        end
    end
    if nnz(isfinite(idxGrid)) ~= nCh
        idxGrid = detectorIndexGridFallback(nCh);
    end

    firstImage = projectOne(volume(:,:,:,1), planeName, mode);
    tileSize = size(firstImage);
    gap = max(2, round(min(tileSize) * 0.08));
    [nRow, nCol] = size(idxGrid);
    canvas = nan(nRow * tileSize(1) + (nRow-1) * gap, ...
        nCol * tileSize(2) + (nCol-1) * gap);

    for r = 1:nRow
        for c = 1:nCol
            ch = idxGrid(r,c);
            if ~isfinite(ch) || ch < 1 || ch > nCh
                continue;
            end
            img = projectOne(volume(:,:,:,ch), planeName, mode);
            y0 = (r-1) * (tileSize(1) + gap) + 1;
            x0 = (c-1) * (tileSize(2) + gap) + 1;
            canvas(y0:y0+tileSize(1)-1, x0:x0+tileSize(2)-1) = img;
        end
    end
end

function image = projectOne(vol, planeName, mode)
    switch lower(planeName)
        case 'xy'
            if strcmp(mode, 'sum')
                image = squeeze(sum(vol, 3));
            else
                image = squeeze(max(vol, [], 3));
            end
        case 'xz'
            if strcmp(mode, 'sum')
                image = squeeze(sum(vol, 1)).';
            else
                image = squeeze(max(vol, [], 1)).';
            end
        case 'yz'
            if strcmp(mode, 'sum')
                image = squeeze(sum(vol, 2)).';
            else
                image = squeeze(max(vol, [], 2)).';
            end
        otherwise
            error('reportWaterObjectiveISMCalibration:BadProjection', ...
                'Projection must be XY, XZ, or YZ.');
    end
end

function writeAllDetectorProjectionFigure(measuredCanvas, fittedCanvas, ...
        residualCanvas, planeName, normalizationNote, opts, outFile)
    fig = figure('Visible', opts.figureVisible, 'Color', 'w', ...
        'Position', [60 60 1280 920]);
    tl = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    canvases = {measuredCanvas, fittedCanvas, residualCanvas};
    titles = {'measured calibration bead', 'fitted calibration model', ...
        'measured - fitted'};
    for row = 1:3
        ax = nexttile(tl);
        imagesc(ax, canvases{row});
        axis(ax, 'image');
        axis(ax, 'off');
        title(ax, titles{row}, 'Interpreter', 'none');
        if row == 3
            lim = max(abs(canvases{row}(:)), [], 'omitnan');
            if ~isfinite(lim) || lim <= 0, lim = 1; end
            caxis(ax, [-lim lim]);
            colormap(ax, redBlueMap(256));
        else
            setPositiveImageLimits(ax, canvases{row});
            colormap(ax, hot);
        end
        colorbar(ax);
    end
    title(tl, sprintf('Calibration all detectors: %s PSF projection (%s)', ...
        upper(planeName), normalizationNote), 'Interpreter', 'none');
    exportgraphics(fig, outFile, 'Resolution', 180);
    close(fig);
end

function writePerDetectorFigure(A, B, R, detIdx, channelID, info, ...
        normalizationNote, opts, outFile)
    fig = figure('Visible', opts.figureVisible, 'Color', 'w', ...
        'Position', [80 80 1350 820]);
    tl = tiledlayout(fig, 3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    volumes = {A, B, R};
    rows = {'measured', 'fitted', 'measured - fitted'};
    projections = {'xy','xz','yz'};
    for row = 1:3
        for col = 1:3
            ax = nexttile(tl);
            img = projectOne(volumes{row}, projections{col}, opts.projectionMode);
            imagesc(ax, img);
            axis(ax, 'image');
            axis(ax, 'off');
            title(ax, sprintf('%s %s', rows{row}, upper(projections{col})), ...
                'Interpreter', 'none');
            if row == 3
                lim = max(abs(img(:)), [], 'omitnan');
                if ~isfinite(lim) || lim <= 0, lim = 1; end
                caxis(ax, [-lim lim]);
                colormap(ax, redBlueMap(256));
            else
                setPositiveImageLimits(ax, img);
                colormap(ax, hot);
            end
            colorbar(ax);
        end
    end
    title(tl, sprintf(['Calibration detector %02d / channel %02d; ' ...
        'z0 %.4f um; %s'], detIdx, channelID, info.fitZ0Um, normalizationNote), ...
        'Interpreter', 'none');
    exportgraphics(fig, outFile, 'Resolution', 180);
    close(fig);
end

function setPositiveImageLimits(ax, img)
    finiteValues = img(isfinite(img));
    if isempty(finiteValues)
        return;
    end
    lo = min(finiteValues);
    hi = max(finiteValues);
    if hi > lo
        caxis(ax, [lo hi]);
    end
end

function T = wavefrontTable(calibration, result)
    modes = calibration.wavefront.fitModes(:);
    waves = zeros(numel(modes),1);
    for k = 1:numel(modes)
        if isfield(calibration.wavefront.coeffs, modes{k})
            waves(k) = calibration.wavefront.coeffs.(modes{k});
        elseif isfield(result.fit.estCoeffs, modes{k})
            waves(k) = result.fit.estCoeffs.(modes{k});
        end
    end
    lambdaRefUm = calibration.wavefront.lambdaRefUm;
    opdNmAtLambdaRef = waves * lambdaRefUm * 1000;
    T = table(modes, waves, opdNmAtLambdaRef);
end

function T = detectorQETable(calibration)
    detectorIndex = (1:numel(calibration.detector.channelIDs)).';
    channelID = calibration.detector.channelIDs(:);
    relativeQE = calibration.detector.relativeQE(:);
    T = table(detectorIndex, channelID, relativeQE);
end

function idxGrid = detectorIndexGridFallback(nCh)
    nCol = ceil(sqrt(nCh));
    nRow = ceil(nCh / nCol);
    idxGrid = nan(nRow, nCol);
    for k = 1:nCh
        r = floor((k-1) / nCol) + 1;
        c = mod(k-1, nCol) + 1;
        idxGrid(r,c) = k;
    end
end

function cmap = redBlueMap(n)
    if nargin < 1, n = 256; end
    x = linspace(-1,1,n).';
    cmap = [max(0,-x), 1-abs(x), max(0,x)];
    cmap = max(min(cmap,1),0);
end

function fileName = defaultCalibrationFile()
    root = fileparts(fileparts(mfilename('fullpath')));
    fileName = fullfile(root, 'output_matlab', ...
        'water_objective_calibration', 'test_20260515_144001', ...
        'water_objective_ism_calibration.mat');
end
