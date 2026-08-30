function heightMaps = immune_cell_MIET_height_maps(biasCorrectionMat, calibrationFile, cfg)
%IMMUNE_CELL_MIET_HEIGHT_MAPS Convert bias-corrected lifetimes into heights.
%
% heightMaps = immune_cell_MIET_height_maps(biasCorrectionMat, calibrationFile)
% heightMaps = immune_cell_MIET_height_maps(biasCorrectionMat, calibrationFile, cfg)
%
% Takes the output of immune_cell_MIET_apply_long_lifetime_bias_correction and
% maps the membrane component's lifetime through a MIET calibration curve to
% get height above the metal-plus-spacer surface. Every available pixel
% binning - native, 2x2 sliding and 4x4 sliding - is converted. Neither the
% source analysis nor the bias correction is modified; results are written to a
% new MAT and one PNG per binning beside the bias correction.
%
% WHICH LIFETIME IS INVERTED
%
% The LONG component - the cell membrane - is the one that carries height
% information. Its bias-corrected value is used, not the raw posterior mean,
% because the raw estimate is biased by the SLB-amplitude prior and the
% correction is exactly what removes that bias. Raw heights are computed too,
% but only as a diagnostic panel, so the size of the correction is visible in
% nanometres rather than only in nanoseconds.
%
% WHY PIXELS COME BACK BLANK
%
% A MIET curve is invertible only between its most-quenched value at the
% surface and the first lifetime that recurs at a greater height. A measured
% lifetime outside that window has NO height, and inventing one by clamping to
% the nearest end would manufacture a flat plateau exactly where the data has
% run out of information. Such pixels are therefore left NaN and counted
% separately as belowFloor (shorter than the curve reaches; more quenched than
% the model allows) and aboveCeiling (longer than the last unique value; the
% membrane is beyond MIET's reach, or the stack/quantum-yield assumptions make
% the curve too shallow).
%
% HEIGHT UNCERTAINTY
%
% The bias correction predicts a per-pixel lifetime standard deviation. It is
% propagated as dz = sigma_tau * |dz/dtau| using the local slope of the
% calibration curve. This matters because the curve is steep near the surface
% and flat far from it: the same lifetime precision is worth a few nanometres
% low down and tens of nanometres high up, and a height map without that
% context invites over-reading its own far field.
%
% INPUTS
%   biasCorrectionMat  path to immune_cell_MIET_long_lifetime_bias_correction.mat
%   calibrationFile    MAT holding a calibration curve. A file written by
%                      make_rt_miet_calibration is used with its full
%                      provenance; any file readable by
%                      load_lifetime_height_calibration also works.
%   cfg .outputDir     defaults to the bias correction's folder
%       .outputMatFile defaults to immune_cell_MIET_height_maps<suffix>.mat
%       .nameSuffix    appended to the MAT and figure names, so a second
%                      calibration variant can be written beside the first
%                      instead of silently replacing it
%       .writeFigures  true
%       .showFigures   false
%       .colormap      'viridis'
%       .heightClipQuantiles  [0.02 0.98] for the display colour limits
%
% See also MAKE_RT_MIET_CALIBRATION, MIET_CALIBRATION_CURVE,
% IMMUNE_CELL_MIET_APPLY_LONG_LIFETIME_BIAS_CORRECTION.

    if nargin < 3 || isempty(cfg); cfg = struct(); end
    biasCorrectionMat = char(biasCorrectionMat);
    if ~isfile(biasCorrectionMat)
        error('immune_cell_MIET_height_maps:BiasFile', ...
            'Bias-correction MAT does not exist: %s', biasCorrectionMat);
    end
    cfg = fillDefaults(cfg, biasCorrectionMat);
    if ~isfolder(cfg.outputDir); mkdir(cfg.outputDir); end

    loaded = load(biasCorrectionMat, 'biasCorrection');
    if ~isfield(loaded, 'biasCorrection') || ~isstruct(loaded.biasCorrection)
        error('immune_cell_MIET_height_maps:BiasSchema', ...
            'MAT must contain a biasCorrection struct: %s', biasCorrectionMat);
    end
    biasCorrection = loaded.biasCorrection;
    calibration = readCalibration(calibrationFile);

    stageNames = fieldnames(biasCorrection.stages);
    stages = struct();
    for index = 1:numel(stageNames)
        name = stageNames{index};
        source = biasCorrection.stages.(name);
        if ~isfield(source, 'available') || ~source.available
            stages.(name) = struct('available', false, ...
                'label', getFieldOr(source, 'label', name), ...
                'reason', 'the binning was not available in the source analysis');
            continue;
        end
        stages.(name) = convertStage(source, calibration);
    end

    heightMaps = struct();
    heightMaps.method = ['MIET height from the bias-corrected long-lifetime ' ...
        'component; no clamping outside the invertible window'];
    heightMaps.biasCorrectionMat = biasCorrectionMat;
    heightMaps.sourceResultMat = getFieldOr(biasCorrection, 'sourceResultMat', '');
    heightMaps.fixedSlbLifetimeNs = getFieldOr(biasCorrection, 'fixedSlbLifetimeNs', NaN);
    heightMaps.calibration = calibration;
    heightMaps.config = cfg;
    heightMaps.stages = stages;
    % The consistency check must use whichever lifetime the calibration was
    % actually anchored on. Under the H2 reading the pipeline's fixed short
    % component is substrate emission, not bilayer fluorescence, so checking it
    % against the curve would report a spurious failure.
    slbCheckNs = heightMaps.fixedSlbLifetimeNs;
    if ~isempty(cfg.slbLifetimeNsOverride)
        slbCheckNs = double(cfg.slbLifetimeNsOverride);
    end
    heightMaps.slbReferenceLifetimeNs = slbCheckNs;
    heightMaps.slbConsistency = slbConsistency(slbCheckNs, calibration);
    heightMaps.outputFiles = struct('mat', cfg.outputMatFile, 'figures', {{}});

    % Numbers first: a headless rendering failure must not lose the maps.
    save(cfg.outputMatFile, 'heightMaps', '-v7.3');

    if cfg.writeFigures
        figureFiles = {};
        try
            for index = 1:numel(stageNames)
                name = stageNames{index};
                if ~stages.(name).available; continue; end
                figureFile = fullfile(cfg.outputDir, sprintf( ...
                    'immune_cell_MIET_%s_height%s.png', name, cfg.nameSuffix));
                writeHeightFigure(stages.(name), calibration, ...
                    slbCheckNs, figureFile, cfg);
                figureFiles{end + 1} = figureFile; %#ok<AGROW>
            end
            heightMaps.outputFiles.figures = figureFiles;
            save(cfg.outputMatFile, 'heightMaps', '-v7.3');
        catch figureError
            warning('immune_cell_MIET_height_maps:FigureExport', ...
                ['Height maps were saved, but figure export failed: %s'], ...
                figureError.message);
        end
    end

    fprintf('immune_cell_MIET_height_maps: saved %s\n', cfg.outputMatFile);
    printStageSummary(stages);
end

% ---------------------------------------------------------------- conversion

function stage = convertStage(source, calibration)
    inverse = source.inverse;
    reliableMask = logical(inverse.reliableDisplayMask);
    displayMask = logical(inverse.inputDisplayMask);
    corrected = double(inverse.displayCorrectedLifetimeNs);
    raw = double(inverse.rawLifetimeNs);

    [heightNm, statusCode] = invertLifetime(corrected, reliableMask, calibration);
    rawHeightNm = invertLifetime(raw, displayMask, calibration);

    slopeNmPerNs = calibrationSlope(corrected, calibration);
    slopeNmPerNs(~isfinite(heightNm)) = NaN;
    heightStdNm = nan(size(heightNm));
    if isfield(inverse, 'predictedEmpiricalStdNs')
        sigmaTau = double(inverse.predictedEmpiricalStdNs);
        heightStdNm = abs(sigmaTau) .* slopeNmPerNs;
    end

    heightCorrectionNm = rawHeightNm - heightNm;

    counted = nnz(reliableMask);
    stage = struct();
    stage.available = true;
    stage.label = getFieldOr(source, 'label', '');
    stage.windowSize = getFieldOr(source, 'windowSize', [NaN NaN]);
    stage.heightNm = single(heightNm);
    stage.heightStdNm = single(heightStdNm);
    stage.rawHeightNm = single(rawHeightNm);
    stage.heightCorrectionNm = single(heightCorrectionNm);
    stage.slopeNmPerNs = single(slopeNmPerNs);
    stage.correctedLifetimeNs = single(corrected);
    stage.rawLifetimeNs = single(raw);
    stage.displayMask = displayMask;
    stage.reliableLifetimeMask = reliableMask;
    stage.heightMask = isfinite(heightNm);
    stage.statusCode = statusCode;
    stage.statusDefinitions = statusDefinitions();

    stage.summary = struct( ...
        'reliableLifetimePixelCount', counted, ...
        'heightPixelCount', nnz(stage.heightMask), ...
        'heightFractionOfReliable', nnz(stage.heightMask) / max(counted, 1), ...
        'belowCalibrationFloorCount', nnz(statusCode == 2), ...
        'aboveCalibrationCeilingCount', nnz(statusCode == 3), ...
        'medianHeightNm', finiteMedian(heightNm(stage.heightMask)), ...
        'p05HeightNm', finiteQuantile(heightNm(stage.heightMask), 0.05), ...
        'p95HeightNm', finiteQuantile(heightNm(stage.heightMask), 0.95), ...
        'medianHeightStdNm', finiteMedian(heightStdNm(stage.heightMask)), ...
        'medianHeightCorrectionNm', ...
            finiteMedian(heightCorrectionNm(stage.heightMask)));
end

function [heightNm, statusCode] = invertLifetime(lifetimeNs, mask, calibration)
% statusCode: 0 not evaluated, 1 inverted, 2 below the curve's floor,
%             3 above the last unique lifetime.
    heightNm = nan(size(lifetimeNs));
    statusCode = zeros(size(lifetimeNs), 'uint8');
    evaluate = mask & isfinite(lifetimeNs);
    if ~any(evaluate(:)); return; end

    lo = calibration.lifetimeNs(1);
    hi = calibration.lifetimeNs(end);
    values = lifetimeNs(evaluate);

    below = values < lo;
    above = values > hi;
    inside = ~below & ~above;

    heightsHere = nan(size(values));
    if any(inside)
        heightsHere(inside) = interp1(calibration.lifetimeNs, ...
            calibration.heightNm, values(inside), 'pchip', NaN);
    end
    codes = zeros(size(values), 'uint8');
    codes(inside) = 1;
    codes(below) = 2;
    codes(above) = 3;

    heightNm(evaluate) = heightsHere;
    statusCode(evaluate) = codes;
end

function slopeNmPerNs = calibrationSlope(lifetimeNs, calibration)
% dz/dtau at the measured lifetime, from the calibration curve itself. Central
% differences on the tabulated curve, then sampled at the pixel's lifetime.
    tau = calibration.lifetimeNs;
    z = calibration.heightNm;
    dzdtau = gradient(z, tau);
    slopeNmPerNs = nan(size(lifetimeNs));
    valid = isfinite(lifetimeNs) & lifetimeNs >= tau(1) & lifetimeNs <= tau(end);
    if any(valid(:))
        slopeNmPerNs(valid) = interp1(tau, dzdtau, lifetimeNs(valid), 'linear', NaN);
    end
end

function report = slbConsistency(tauSlbNs, calibration)
% The bare SLB sits a few nanometres above the spacer, so the model's floor is
% an upper bound on how short the SLB lifetime can be. If the measured SLB
% lifetime is below the floor, the curve cannot describe the SLB at all and the
% stack or the quantum yield is wrong. Reporting it is the cheapest available
% check on the calibration, so it is stored with every result rather than left
% to be noticed by eye.
    report = struct();
    report.measuredSlbLifetimeNs = tauSlbNs;
    report.calibrationFloorLifetimeNs = calibration.lifetimeNs(1);
    report.calibrationFloorHeightNm = calibration.heightNm(1);
    report.slbLifetimeIsBelowFloor = isfinite(tauSlbNs) && ...
        tauSlbNs < calibration.lifetimeNs(1);
    if report.slbLifetimeIsBelowFloor
        report.note = ['The measured SLB lifetime is SHORTER than anything the ' ...
            'calibration can produce, so the stack or the quantum yield ' ...
            'underestimates the quenching. Heights are still reported for the ' ...
            'membrane component, but their absolute zero is unverified.'];
    else
        report.impliedSlbHeightNm = interp1(calibration.lifetimeNs, ...
            calibration.heightNm, tauSlbNs, 'pchip', NaN);
        report.note = ['The measured SLB lifetime falls inside the calibration, ' ...
            'so the curve''s surface end is at least self-consistent.'];
    end
end

% ------------------------------------------------------------------ graphics

function writeHeightFigure(stage, calibration, tauSlbNs, outputFile, cfg)
    visibility = 'off';
    if cfg.showFigures; visibility = 'on'; end

    height = double(stage.heightNm);
    heightStd = double(stage.heightStdNm);
    correction = double(stage.heightCorrectionNm);
    lifetime = double(stage.correctedLifetimeNs);
    lifetime(~stage.reliableLifetimeMask) = NaN;

    heightLimits = robustLimits(height(stage.heightMask), cfg.heightClipQuantiles);
    stdLimits = robustLimits(heightStd(stage.heightMask), [0 0.98]);
    correctionLimit = max(finiteQuantile(abs(correction(stage.heightMask)), 0.98), 1);
    lifetimeLimits = robustLimits(lifetime(stage.reliableLifetimeMask), [0.02 0.98]);

    h = figure('Color', 'w', 'Visible', visibility, 'Position', [50 50 1560 900]);
    layout = tiledlayout(h, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    showPanel(nexttile(layout), height, perceptualColormap(cfg.colormap), ...
        heightLimits, 'MIET height', 'height [nm]');
    showPanel(nexttile(layout), heightStd, perceptualColormap('magma'), ...
        stdLimits, 'propagated height uncertainty', '\sigma_z [nm]');
    showPanel(nexttile(layout), correction, perceptualColormap('coolwarm'), ...
        [-correctionLimit correctionLimit], ...
        'raw minus corrected height', '\Deltaz [nm]');
    showPanel(nexttile(layout), lifetime, perceptualColormap('cividis'), ...
        lifetimeLimits, 'bias-corrected membrane lifetime', '\tau [ns]');

    status = double(stage.statusCode);
    status(status == 0) = NaN;
    ax = nexttile(layout);
    showPanel(ax, status, [0.16 0.44 0.71; 0.85 0.37 0.01; 0.55 0.13 0.42], ...
        [0.5 3.5], 'inversion status', '');
    colorbar(ax, 'Ticks', 1:3, 'TickLabels', ...
        {'inverted', 'below floor', 'above ceiling'});

    ax = nexttile(layout);
    histogram(ax, height(stage.heightMask), 60, 'FaceColor', [0.20 0.40 0.65], ...
        'EdgeColor', 'none');
    xlabel(ax, 'height [nm]'); ylabel(ax, 'pixels');
    title(ax, 'height distribution');
    grid(ax, 'on'); box(ax, 'off');

    total = stage.summary.reliableLifetimePixelCount;
    title(layout, sprintf('MIET height | %s | %d/%d reliable pixels inverted (%.1f%%)', ...
        stage.label, stage.summary.heightPixelCount, total, ...
        100 * stage.summary.heightFractionOfReliable), 'FontWeight', 'bold');
    subtitle(layout, sprintf([ ...
        'calibration %.3f-%.3f ns over %.1f-%.1f nm | above ceiling %d px, ' ...
        'below floor %d px | measured tau_SLB %.3f ns | median z %.1f nm ' ...
        '(5-95%%: %.1f-%.1f nm)'], ...
        calibration.lifetimeNs(1), calibration.lifetimeNs(end), ...
        calibration.heightNm(1), calibration.heightNm(end), ...
        stage.summary.aboveCalibrationCeilingCount, ...
        stage.summary.belowCalibrationFloorCount, tauSlbNs, ...
        stage.summary.medianHeightNm, stage.summary.p05HeightNm, ...
        stage.summary.p95HeightNm), 'Interpreter', 'none');

    exportgraphics(h, outputFile, 'Resolution', 250);
    if ~cfg.showFigures; close(h); end
end

function showPanel(ax, data, map, limits, titleText, barLabel)
    handle = imagesc(ax, data);
    set(handle, 'AlphaData', isfinite(data));
    colormap(ax, map);
    if all(isfinite(limits)) && limits(2) > limits(1)
        caxis(ax, limits); %#ok<CAXIS>
    end
    [nRow, nCol] = size(data);
    axis(ax, 'image');
    set(ax, 'PlotBoxAspectRatio', [nCol nRow 1]);
    set(ax, 'XTick', [], 'YTick', [], 'Color', [0.15 0.15 0.15]);
    title(ax, titleText);
    if ~isempty(barLabel)
        bar = colorbar(ax);
        bar.Label.String = barLabel;
    end
end

% ------------------------------------------------------------------- helpers

function calibration = readCalibration(calibrationFile)
    calibrationFile = char(calibrationFile);
    if ~isfile(calibrationFile)
        error('immune_cell_MIET_height_maps:CalibrationFile', ...
            'Calibration MAT does not exist: %s', calibrationFile);
    end
    stored = load(calibrationFile);
    if isfield(stored, 'calibration') && isstruct(stored.calibration) && ...
            isfield(stored.calibration, 'heightNm')
        calibration = stored.calibration;
        calibration.file = calibrationFile;
    else
        calibration = load_lifetime_height_calibration(calibrationFile);
    end

    lifetimeNs = double(calibration.lifetimeNs(:));
    heightNm = double(calibration.heightNm(:));
    keep = isfinite(lifetimeNs) & isfinite(heightNm);
    lifetimeNs = lifetimeNs(keep);
    heightNm = heightNm(keep);
    [lifetimeNs, order] = sort(lifetimeNs, 'ascend');
    heightNm = heightNm(order);
    [lifetimeNs, unique_idx] = unique(lifetimeNs, 'stable');
    heightNm = heightNm(unique_idx);
    if numel(lifetimeNs) < 3
        error('immune_cell_MIET_height_maps:CalibrationTooShort', ...
            'The calibration needs at least three unique lifetimes.');
    end
    if any(diff(heightNm) <= 0)
        error('immune_cell_MIET_height_maps:CalibrationNotMonotonic', ...
            ['The calibration is not monotonic after sorting by lifetime, so ' ...
             'the inverse is not a function. Regenerate it with curveType ' ...
             '''minimum''.']);
    end
    calibration.lifetimeNs = lifetimeNs;
    calibration.heightNm = heightNm;
end

function cfg = fillDefaults(cfg, biasCorrectionMat)
    if ~isfield(cfg, 'outputDir') || isempty(cfg.outputDir)
        cfg.outputDir = fileparts(biasCorrectionMat);
    end
    if ~isfield(cfg, 'nameSuffix') || isempty(cfg.nameSuffix)
        cfg.nameSuffix = '';
    end
    if ~isfield(cfg, 'slbLifetimeNsOverride')
        cfg.slbLifetimeNsOverride = [];
    end
    if ~isfield(cfg, 'outputMatFile') || isempty(cfg.outputMatFile)
        cfg.outputMatFile = fullfile(cfg.outputDir, sprintf( ...
            'immune_cell_MIET_height_maps%s.mat', cfg.nameSuffix));
    end
    if ~isfield(cfg, 'writeFigures') || isempty(cfg.writeFigures)
        cfg.writeFigures = true;
    end
    if ~isfield(cfg, 'showFigures') || isempty(cfg.showFigures)
        cfg.showFigures = false;
    end
    if ~isfield(cfg, 'colormap') || isempty(cfg.colormap)
        cfg.colormap = 'viridis';
    end
    if ~isfield(cfg, 'heightClipQuantiles') || isempty(cfg.heightClipQuantiles)
        cfg.heightClipQuantiles = [0.02 0.98];
    end
end

function definitions = statusDefinitions()
    definitions = struct( ...
        'code0', 'not evaluated (outside the display mask or not reliable)', ...
        'code1', 'inverted through the calibration curve', ...
        'code2', ['lifetime shorter than the calibration floor: more quenched ' ...
                  'than the modelled stack allows'], ...
        'code3', ['lifetime longer than the last uniquely invertible value: ' ...
                  'beyond MIET sensitivity for this curve']);
end

function limits = robustLimits(values, quantiles)
    values = double(values(:));
    values = values(isfinite(values));
    if isempty(values); limits = [0 1]; return; end
    lo = finiteQuantile(values, quantiles(1));
    hi = finiteQuantile(values, quantiles(2));
    if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
        lo = min(values); hi = max(values);
    end
    if hi <= lo; hi = lo + eps(lo) + 1e-6; end
    limits = [lo hi];
end

function value = finiteQuantile(values, probability)
    values = double(values(:));
    values = values(isfinite(values));
    if isempty(values); value = NaN; return; end
    value = quantile(values, probability);
end

function value = finiteMedian(values)
    values = double(values(:));
    values = values(isfinite(values));
    if isempty(values); value = NaN; return; end
    value = median(values);
end

function value = getFieldOr(structure, name, fallback)
    if isstruct(structure) && isfield(structure, name)
        value = structure.(name);
    else
        value = fallback;
    end
end

function printStageSummary(stages)
    names = fieldnames(stages);
    for index = 1:numel(names)
        stage = stages.(names{index});
        if ~stage.available
            fprintf('  %-11s unavailable (%s)\n', names{index}, ...
                getFieldOr(stage, 'reason', 'no reason recorded'));
            continue;
        end
        fprintf(['  %-11s %6d/%6d px inverted (%.1f%%)  median z %.1f nm  ' ...
                 'sigma_z %.1f nm  above ceiling %d  below floor %d\n'], ...
            names{index}, stage.summary.heightPixelCount, ...
            stage.summary.reliableLifetimePixelCount, ...
            100 * stage.summary.heightFractionOfReliable, ...
            stage.summary.medianHeightNm, stage.summary.medianHeightStdNm, ...
            stage.summary.aboveCalibrationCeilingCount, ...
            stage.summary.belowCalibrationFloorCount);
    end
end
