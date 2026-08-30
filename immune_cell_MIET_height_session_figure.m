function out = immune_cell_MIET_height_session_figure(summarySource, cfg)
%IMMUNE_CELL_MIET_HEIGHT_SESSION_FIGURE Session overview of MIET height maps.
%
% out = immune_cell_MIET_height_session_figure(summaryMatOrDir)
% out = immune_cell_MIET_height_session_figure(summarySource, cfg)
%
% Tiles every acquisition's height map for one binning on a SHARED colour
% scale, with the pooled height distribution and the per-acquisition inversion
% status beside it. One figure per binning.
%
% WHY A SHARED SCALE
%
% Each per-acquisition figure autoscales to its own data, which is right for
% inspecting one field and wrong for comparing six: a 20 nm spread and a 120 nm
% spread then look identical. Fixing one scale across the session is what makes
% a real difference between fields visible, and it is also what exposes the
% failure mode here - if every field saturates at the same value, the maps are
% reporting the calibration's ceiling rather than the sample.
%
% THE STATUS BAR IS PART OF THE RESULT
%
% The fraction of pixels that could NOT be inverted is not a diagnostic
% footnote; it is the headline when it is large. It is therefore plotted
% alongside the maps rather than buried in a CSV, so a reader cannot take a
% median height from the figure without also seeing how little of the field it
% came from.
%
% INPUTS
%   summarySource  the batch summary MAT written by
%                  run_batch_immune_cell_MIET_height_maps, or the directory
%                  holding it
%   cfg .binnings      {'native','sliding2x2','sliding4x4'}
%       .outputDir     defaults beside the summary
%       .colormap      'viridis'
%       .heightLimits  [] auto from the pooled 2-98% range, or [lo hi]
%       .showFigures   false
%
% See also RUN_BATCH_IMMUNE_CELL_MIET_HEIGHT_MAPS, IMMUNE_CELL_MIET_HEIGHT_MAPS.

    if nargin < 2 || isempty(cfg); cfg = struct(); end
    summaryMat = resolveSummary(summarySource);
    cfg = fillDefaults(cfg, summaryMat);

    loaded = load(summaryMat, 'summary', 'info');
    summary = loaded.summary;
    info = loaded.info;
    if isempty(summary)
        error('immune_cell_MIET_height_session_figure:EmptySummary', ...
            'The batch summary holds no rows: %s', summaryMat);
    end

    out = struct();
    out.summaryMat = summaryMat;
    out.figures = {};
    out.binnings = {};

    for index = 1:numel(cfg.binnings)
        binning = cfg.binnings{index};
        rows = summary(summary.binning == string(binning) & summary.available, :);
        if isempty(rows)
            fprintf('  %s: nothing available, skipped\n', binning);
            continue;
        end
        panels = loadPanels(rows, binning);
        if isempty(panels)
            fprintf('  %s: no height maps could be loaded, skipped\n', binning);
            continue;
        end
        figureFile = fullfile(cfg.outputDir, sprintf( ...
            'immune_cell_MIET_height_session_%s.png', binning));
        stats = writeSessionFigure(panels, rows, binning, info, figureFile, cfg);
        out.figures{end + 1} = figureFile; %#ok<AGROW>
        out.binnings{end + 1} = binning; %#ok<AGROW>
        out.stats.(binning) = stats;
        fprintf('  %s: %d fields, pooled median %.1f nm, inverted %.1f%% -> %s\n', ...
            binning, numel(panels), stats.pooledMedianNm, ...
            100 * stats.invertedFraction, figureFile);
    end
end

% ------------------------------------------------------------------- loading

function panels = loadPanels(rows, binning)
    panels = struct('label', {}, 'heightNm', {}, 'inverted', {}, ...
        'aboveCeiling', {}, 'belowFloor', {}, 'total', {});
    for index = 1:height(rows)
        matFile = char(rows.heightMapsMat(index));
        if ~isfile(matFile)
            warning('immune_cell_MIET_height_session_figure:MissingMat', ...
                'Height MAT is missing, field skipped: %s', matFile);
            continue;
        end
        loaded = load(matFile, 'heightMaps');
        stage = loaded.heightMaps.stages.(binning);
        if ~stage.available; continue; end
        entry = struct();
        entry.label = shortLabel(char(rows.acquisition(index)));
        entry.heightNm = double(stage.heightNm);
        entry.heightNm(~stage.heightMask) = NaN;
        entry.inverted = stage.summary.heightPixelCount;
        entry.aboveCeiling = stage.summary.aboveCalibrationCeilingCount;
        entry.belowFloor = stage.summary.belowCalibrationFloorCount;
        entry.total = stage.summary.reliableLifetimePixelCount;
        panels(end + 1) = entry; %#ok<AGROW>
    end
end

function label = shortLabel(acquisition)
    parts = split(string(acquisition), '|');
    label = char(regexprep(parts(1), '^_?\d{8}-', ''));
    if isempty(label); label = char(parts(1)); end
end

% ------------------------------------------------------------------ graphics

function stats = writeSessionFigure(panels, rows, binning, info, figureFile, cfg)
    pooled = [];
    for index = 1:numel(panels)
        values = panels(index).heightNm(:);
        pooled = [pooled; values(isfinite(values))]; %#ok<AGROW>
    end
    limits = cfg.heightLimits;
    if isempty(limits)
        limits = [quantile(pooled, 0.02) quantile(pooled, 0.98)];
        if ~(limits(2) > limits(1)); limits = [min(pooled) max(pooled) + 1e-6]; end
    end

    inverted = sum([panels.inverted]);
    total = sum([panels.total]);
    stats = struct( ...
        'fieldCount', numel(panels), ...
        'pooledMedianNm', median(pooled), ...
        'pooledP05Nm', quantile(pooled, 0.05), ...
        'pooledP95Nm', quantile(pooled, 0.95), ...
        'invertedFraction', inverted / max(total, 1), ...
        'aboveCeilingFraction', sum([panels.aboveCeiling]) / max(total, 1), ...
        'belowFloorFraction', sum([panels.belowFloor]) / max(total, 1), ...
        'displayLimitsNm', limits);

    visibility = 'off';
    if cfg.showFigures; visibility = 'on'; end
    fieldCount = numel(panels);
    % Beyond about six fields a single strip of maps is too thin to read
    % anything from, so the strip wraps. The sessions differ in field count
    % (6 at room temperature, 11 at 4 C) and both have to stay legible.
    mapCols = min(fieldCount, 6);
    mapRows = ceil(fieldCount / mapCols);
    h = figure('Color', 'w', 'Visible', visibility, ...
        'Position', [40 40 max(1400, 260 * mapCols), 480 + 230 * mapRows]);
    % A flat grid gives the maps whole tall tiles, so they shrink and leave a
    % band of empty white. Nesting lets the map block keep only what it needs.
    outerRows = 2 * mapRows + 3;
    outer = tiledlayout(h, outerRows, 1, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    mapRow = tiledlayout(outer, mapRows, mapCols, 'TileSpacing', 'compact', ...
        'Padding', 'none');
    mapRow.Layout.Tile = 1;
    mapRow.Layout.TileSpan = [2 * mapRows, 1];
    bottomRow = tiledlayout(outer, 1, 2, 'TileSpacing', 'compact', ...
        'Padding', 'none');
    bottomRow.Layout.Tile = 2 * mapRows + 1;
    bottomRow.Layout.TileSpan = [3 1];

    map = perceptualColormap(cfg.colormap);
    for index = 1:fieldCount
        ax = nexttile(mapRow);
        handle = imagesc(ax, panels(index).heightNm);
        set(handle, 'AlphaData', isfinite(panels(index).heightNm));
        colormap(ax, map);
        caxis(ax, limits); %#ok<CAXIS>
        [nRow, nCol] = size(panels(index).heightNm);
        axis(ax, 'image');
        set(ax, 'PlotBoxAspectRatio', [nCol nRow 1], 'XTick', [], 'YTick', [], ...
            'Color', [0.15 0.15 0.15]);
        title(ax, sprintf('%s  (%.0f%%)', panels(index).label, ...
            100 * panels(index).inverted / max(panels(index).total, 1)), ...
            'FontSize', 9);
        if index == fieldCount
            % Named colourBar, not bar: a local variable called bar shadows the
            % bar() function for the whole function body, which broke the
            % stacked status plot below.
            colourBar = colorbar(ax);
            colourBar.Label.String = 'height [nm]';
        end
    end

    ax = nexttile(bottomRow);
    histogram(ax, pooled, 80, 'FaceColor', [0.20 0.40 0.65], 'EdgeColor', 'none');
    hold(ax, 'on');
    xline(ax, stats.pooledMedianNm, 'r-', sprintf('median %.0f nm', ...
        stats.pooledMedianNm), 'LineWidth', 1.5);
    ceilingNm = info.calibration.heightNm(end);
    xline(ax, ceilingNm, 'k--', sprintf('ceiling %.0f nm', ceilingNm), ...
        'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
    xlabel(ax, 'height [nm]'); ylabel(ax, 'pixels');
    title(ax, 'pooled height distribution, all fields');
    grid(ax, 'on'); box(ax, 'off');

    ax = nexttile(bottomRow);
    fractions = [[panels.inverted]; [panels.aboveCeiling]; [panels.belowFloor]].';
    fractions = 100 * fractions ./ max([panels.total].', 1);
    barHandle = bar(ax, fractions, 'stacked', 'EdgeColor', 'none');
    barHandle(1).FaceColor = [0.16 0.44 0.71];
    barHandle(2).FaceColor = [0.55 0.13 0.42];
    barHandle(3).FaceColor = [0.85 0.37 0.01];
    set(ax, 'XTick', 1:numel(panels), 'XTickLabel', {panels.label}, ...
        'XTickLabelRotation', 45);
    ylabel(ax, '% of reliable pixels'); ylim(ax, [0 100]);
    legend(ax, {'inverted', 'above ceiling', 'below floor'}, ...
        'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);
    title(ax, 'inversion status per field');
    box(ax, 'off');

    calibrationLabel = '';
    if isfield(info.calibration, 'label'); calibrationLabel = info.calibration.label; end
    title(outer, sprintf('MIET height  |  %s  |  %s', ...
        binningLabel(binning, rows), sessionName(info)), ...
        'FontWeight', 'bold', 'Interpreter', 'none', 'FontSize', 12);
    subtitle(outer, { ...
        sprintf(['pooled median %.1f nm (5-95%%: %.1f-%.1f nm)  |  %.1f%% of ' ...
            'reliable pixels inverted, %.1f%% above the ceiling'], ...
            stats.pooledMedianNm, stats.pooledP05Nm, stats.pooledP95Nm, ...
            100 * stats.invertedFraction, 100 * stats.aboveCeilingFraction), ...
        calibrationLabel}, 'Interpreter', 'none', 'FontSize', 8);

    exportgraphics(h, figureFile, 'Resolution', 220);
    if ~cfg.showFigures; close(h); end
end

function text = binningLabel(binning, rows)
    switch binning
        case 'native';     text = 'native pixels';
        case 'sliding2x2'; text = '2x2 sliding TCSPC';
        case 'sliding4x4'; text = '4x4 sliding TCSPC';
        otherwise;         text = binning;
    end
    text = sprintf('%s (%d fields)', text, height(rows));
end

function name = sessionName(info)
    [~, name] = fileparts(char(info.dataRoot));
    if isempty(name); name = char(info.dataRoot); end
end

% ------------------------------------------------------------------- helpers

function summaryMat = resolveSummary(summarySource)
    summarySource = char(summarySource);
    if isfolder(summarySource)
        summaryMat = fullfile(summarySource, 'immune_cell_MIET_height_summary.mat');
    else
        summaryMat = summarySource;
    end
    if ~isfile(summaryMat)
        error('immune_cell_MIET_height_session_figure:SummaryMissing', ...
            'Batch summary MAT not found: %s', summaryMat);
    end
end

function cfg = fillDefaults(cfg, summaryMat)
    if ~isfield(cfg, 'binnings') || isempty(cfg.binnings)
        cfg.binnings = {'native', 'sliding2x2', 'sliding4x4'};
    end
    if ~isfield(cfg, 'outputDir') || isempty(cfg.outputDir)
        cfg.outputDir = fileparts(summaryMat);
    end
    if ~isfield(cfg, 'colormap') || isempty(cfg.colormap)
        cfg.colormap = 'viridis';
    end
    if ~isfield(cfg, 'heightLimits'); cfg.heightLimits = []; end
    if ~isfield(cfg, 'showFigures') || isempty(cfg.showFigures)
        cfg.showFigures = false;
    end
end
