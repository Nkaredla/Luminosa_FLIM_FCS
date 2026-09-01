function out = immune_cell_MIET_profile_figure(heightMapsMat, cfg)
%IMMUNE_CELL_MIET_PROFILE_FIGURE Cross-section profiles of height and CD58.
%
% out = immune_cell_MIET_profile_figure(heightMapsMat)
% out = immune_cell_MIET_profile_figure(heightMapsMat, cfg)
%
% Cuts one line along x and one along y through a cell and plots MIET height
% and CD58-Atto488 intensity against distance on a twin-axis plot, with both
% directions in the same axes. Two small maps above show where the cuts run.
%
% WHY BOTH QUANTITIES ON ONE PLOT
%
% The question these profiles answer is whether membrane height tracks ligand
% density - whether the cell sits lower where there is more CD58. Plotting the
% two on separate figures makes that comparison an act of memory. On one pair
% of axes the correlation, or its absence, is immediate. Height is on the left
% axis in nanometres, CD58 on the right in photons; they share only the
% distance axis, so neither scale distorts the other.
%
% ENCODING: COLOUR IS DIRECTION, STYLE IS QUANTITY
%
% Four curves in one axes needs a scheme that does not rely on the reader
% tracking a legend entry per line. Colour separates the x cut from the y cut,
% and line style separates height (solid) from CD58 (dashed), so any curve can
% be identified from two independent cues.
%
% yyaxis is used rather than plotyy: plotyy has been deprecated since R2016a
% and returns handles that do not behave like normal axes.
%
% BAND AVERAGING
%
% Height noise here is not ordinary shot noise. The calibration curve flattens
% near tau_0, where dz/dtau reaches several hundred nm per ns, so a pixel whose
% lifetime is a little high is thrown far up in height. The resulting scatter is
% heavy-tailed, and a MEAN over the band is dragged around by those few
% outliers. Both the band statistic and the along-cut smoothing are therefore
% MEDIANS, which ignore the tail instead of averaging it in.
%
% Two windows are applied: bandPixels perpendicular to the cut, and
% smoothAlongPixels along it. Both are printed on the figure, because together
% they set the lateral resolution of the trace.
%
% The average ignores NaNs rather than filling them, so the trace is built only
% from pixels that actually inverted. That is the opposite choice from the
% topography rendering, which fills holes to get a continuous surface: a
% profile is read quantitatively, so an interpolated value in it would be
% indistinguishable from a measured one.
%
% cfg fields
%   binning       'sliding4x4'
%   lineXUm       [] cut along x at this y; default the footprint centroid
%   lineYUm       [] cut along y at this x; default the footprint centroid
%   bandPixels    11 - perpendicular averaging band, in pixels. The raw maps
%                 are noisy pixel to pixel; a narrow band gives an unreadable
%                 trace. 11 px is about 0.27 um here
%   smoothAlongPixels 7 - moving-median window applied ALONG the cut
%   displayClip   [0.02 0.98] quantiles for the map panels' colour limits
%   colormap      'turbo' for the height map panel
%   cd58Colormap  'bone'
%   outputDir     '' writes beside the height MAT
%   namePrefix    '' derived from the acquisition folder
%   fontSize      13
%   showFigures   false
%
% See also IMMUNE_CELL_MIET_TOPOGRAPHY_FIGURE, IMMUNE_CELL_MIET_CD58_IMAGE.

    if nargin < 2 || isempty(cfg); cfg = struct(); end
    cfg = withDefaults(cfg, struct( ...
        'binning', 'sliding4x4', 'lineXUm', [], 'lineYUm', [], ...
        'bandPixels', 15, 'smoothAlongPixels', 7, 'displayClip', [0.02 0.98], ...
        'colormap', 'turbo', 'cd58Colormap', 'bone', ...
        'outputDir', '', 'namePrefix', '', 'fontSize', 13, ...
        'showFigures', false));

    heightMapsMat = char(heightMapsMat);
    loaded = load(heightMapsMat, 'heightMaps');
    heightMaps = loaded.heightMaps;
    if ~isfield(heightMaps.stages, cfg.binning)
        error('immune_cell_MIET_profile_figure:NoBinning', ...
            'Binning %s is not in %s', cfg.binning, heightMapsMat);
    end
    stage = heightMaps.stages.(cfg.binning);
    if ~stage.available
        error('immune_cell_MIET_profile_figure:StageUnavailable', ...
            'Binning %s is unavailable in %s', cfg.binning, heightMapsMat);
    end

    acquisitionDir = acquisitionFolder(heightMapsMat);
    [~, acquisition] = fileparts(acquisitionDir);
    if isempty(cfg.outputDir); cfg.outputDir = fileparts(heightMapsMat); end
    if ~isfolder(cfg.outputDir); mkdir(cfg.outputDir); end
    if isempty(cfg.namePrefix)
        cfg.namePrefix = sprintf('immune_cell_MIET_profiles_%s_%s', ...
            regexprep(acquisition, '^_', ''), cfg.binning);
    end

    heightNm = double(stage.heightNm);
    heightNm(~stage.heightMask) = NaN;
    footprint = logical(stage.displayMask);
    if ~any(footprint(:)); footprint = isfinite(heightNm); end
    cd58 = loadCd58(acquisitionDir, size(heightNm));
    pxUm = pixelSizeFromHeader(acquisitionDir);

    [rows, cols] = croppedRange(footprint, 6);
    heightCrop = heightNm(rows, cols);
    maskCrop = footprint(rows, cols);
    cd58Crop = cd58(rows, cols);
    if nnz(isfinite(heightCrop) & maskCrop) < 50
        error('immune_cell_MIET_profile_figure:TooFewPixels', ...
            'Only %d usable height pixels in %s', ...
            nnz(isfinite(heightCrop) & maskCrop), heightMapsMat);
    end

    [nRow, nCol] = size(heightCrop);
    if isempty(pxUm) || ~isfinite(pxUm) || pxUm <= 0
        pxUm = 1; lateralLabel = 'pixels';
    else
        lateralLabel = '\mum';
    end
    xUm = (0:nCol - 1) * pxUm;
    yUm = (0:nRow - 1) * pxUm;

    % Cuts default to the footprint centroid, which is the one position that
    % needs no justification and lands inside the cell for any shape that is
    % roughly convex.
    [colGrid, rowGrid] = meshgrid(1:nCol, 1:nRow);
    rowCentre = round(mean(rowGrid(maskCrop)));
    colCentre = round(mean(colGrid(maskCrop)));
    if ~isempty(cfg.lineXUm); rowCentre = clampIndex(cfg.lineXUm / pxUm + 1, nRow); end
    if ~isempty(cfg.lineYUm); colCentre = clampIndex(cfg.lineYUm / pxUm + 1, nCol); end

    band = max(1, round(cfg.bandPixels));
    rowBand = clampIndex(rowCentre - floor(band/2), nRow):clampIndex(rowCentre + floor(band/2), nRow);
    colBand = clampIndex(colCentre - floor(band/2), nCol):clampIndex(colCentre + floor(band/2), nCol);

    profile = struct();
    profile.xUm = xUm;
    profile.yUm = yUm;
    along = max(1, round(cfg.smoothAlongPixels));
    profile.heightAlongX = smoothAlong(bandMedian(heightCrop(rowBand, :), 1), along);
    profile.cd58AlongX = smoothAlong(bandMedian(cd58Crop(rowBand, :), 1), along);
    profile.heightAlongY = smoothAlong(bandMedian(heightCrop(:, colBand), 2), along);
    profile.cd58AlongY = smoothAlong(bandMedian(cd58Crop(:, colBand), 2), along);
    profile.smoothAlongPixels = along;
    profile.rowCentre = rowCentre;
    profile.colCentre = colCentre;
    profile.cutYUm = yUm(rowCentre);
    profile.cutXUm = xUm(colCentre);

    out = struct();
    out.heightMapsMat = heightMapsMat;
    out.acquisition = acquisition;
    out.binning = cfg.binning;
    out.profile = profile;
    out.figure = fullfile(cfg.outputDir, sprintf('%s.png', cfg.namePrefix));

    out.bandPixels = band;
    out.bandUm = band * pxUm;
    renderProfiles(heightCrop, cd58Crop, maskCrop, profile, xUm, yUm, ...
        lateralLabel, band, band * pxUm, cfg, heightMaps, acquisition, out.figure);
    fprintf('  %s\n', out.figure);
end

% ------------------------------------------------------------------ rendering

function renderProfiles(height, cd58, mask, profile, xUm, yUm, lateralLabel, ...
        band, bandUm, cfg, heightMaps, acquisition, outputFile)
    fs = cfg.fontSize;
    colourX = [0.30 0.80 1.00];      % cut along x
    colourY = [1.00 0.62 0.20];      % cut along y

    visibility = 'off';
    if cfg.showFigures; visibility = 'on'; end
    screen = get(groot, 'ScreenSize');
    wanted = [1350 820];
    figureSize = [min(wanted(1), screen(3) - 80), min(wanted(2), screen(4) - 120)];
    h = figure('Color', 'w', 'Visible', visibility, ...
        'Position', [40 40 figureSize]);

    % ---- the two maps, with the cuts drawn on them -----------------------
    axHeight = axes('Parent', h, 'Position', [0.055 0.60 0.26 0.34]);
    drawMapWithCuts(axHeight, xUm, yUm, height, ...
        perceptualColormapOrBuiltin(cfg.colormap), profile, colourX, colourY, ...
        lateralLabel, 'MIET height [nm]', fs, robustLimits(height, cfg.displayClip));
    axCd58 = axes('Parent', h, 'Position', [0.395 0.60 0.26 0.34]);
    drawMapWithCuts(axCd58, xUm, yUm, cd58, ...
        perceptualColormapOrBuiltin(cfg.cd58Colormap), profile, colourX, ...
        colourY, lateralLabel, 'CD58 [photons]', fs, robustLimits(cd58, cfg.displayClip));

    % ---- the profiles ----------------------------------------------------
    ax = axes('Parent', h, 'Position', [0.075 0.09 0.86 0.40]);
    yyaxis(ax, 'left');
    plot(ax, xUm, profile.heightAlongX, '-', 'Color', colourX, 'LineWidth', 2);
    hold(ax, 'on');
    plot(ax, yUm, profile.heightAlongY, '-', 'Color', colourY, 'LineWidth', 2);
    ylabel(ax, 'MIET height [nm]', 'FontSize', fs + 1);
    set(ax, 'YColor', [0.15 0.15 0.15]);

    yyaxis(ax, 'right');
    plot(ax, xUm, profile.cd58AlongX, '--', 'Color', colourX, 'LineWidth', 1.6);
    plot(ax, yUm, profile.cd58AlongY, '--', 'Color', colourY, 'LineWidth', 1.6);
    ylabel(ax, 'CD58-Atto488 [photons]', 'FontSize', fs + 1);
    set(ax, 'YColor', [0.35 0.35 0.35]);

    xlabel(ax, sprintf('distance [%s]', lateralLabel), 'FontSize', fs + 1);
    set(ax, 'FontSize', fs, 'Box', 'off', 'LineWidth', 1.0);
    grid(ax, 'on');
    legend(ax, { sprintf('height, cut along x (y = %.2f)', profile.cutYUm), ...
                 sprintf('height, cut along y (x = %.2f)', profile.cutXUm), ...
                 'CD58, cut along x', 'CD58, cut along y' }, ...
        'Location', 'northoutside', 'Orientation', 'horizontal', ...
        'FontSize', fs - 1, 'Box', 'off');

    calib = heightMaps.calibration;
    annotation(h, 'textbox', [0.02 0.955 0.96 0.04], ...
        'String', sprintf(['%s   cross-sections   |   %s   |   band %d px ' ...
            '(%.2f um) x %d px along   |   tau_0 %.2f ns, qy %.3f'], regexprep(acquisition, '^_', ''), ...
            prettyBinning(cfg.binning), band, bandUm, ...
            profile.smoothAlongPixels, ...
            calibField(calib, 'tauFreeNs'), calibField(calib, 'quantumYield')), ...
        'EdgeColor', 'none', 'FontSize', fs + 2, 'FontWeight', 'bold', ...
        'Interpreter', 'none', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');

    set(h, 'PaperPositionMode', 'auto');
    print(h, outputFile, '-dpng', '-r150');
    if ~cfg.showFigures; close(h); end
end

function drawMapWithCuts(ax, xUm, yUm, data, map, profile, colourX, colourY, ...
        lateralLabel, barLabel, fs, limits)
    handle = imagesc(ax, [min(xUm) max(xUm)], [min(yUm) max(yUm)], data);
    set(handle, 'AlphaData', isfinite(data));
    colormap(ax, map);
    if all(isfinite(limits)) && limits(2) > limits(1)
        caxis(ax, limits); %#ok<CAXIS>
    end
    axis(ax, 'image');
    set(ax, 'YDir', 'reverse', 'FontSize', fs - 1, 'Box', 'off', ...
        'Color', [0.92 0.92 0.92]);
    hold(ax, 'on');
    plot(ax, [min(xUm) max(xUm)], profile.cutYUm * [1 1], '-', ...
        'Color', colourX, 'LineWidth', 2);
    plot(ax, profile.cutXUm * [1 1], [min(yUm) max(yUm)], '-', ...
        'Color', colourY, 'LineWidth', 2);
    xlabel(ax, sprintf('x [%s]', lateralLabel), 'FontSize', fs);
    ylabel(ax, sprintf('y [%s]', lateralLabel), 'FontSize', fs);
    bar = colorbar(ax, 'FontSize', fs - 2);
    bar.Label.String = barLabel;
    bar.Label.FontSize = fs - 1;
end

% ------------------------------------------------------------------- helpers

function values = bandMedian(block, dimension)
% Median across the band, ignoring the NaNs that mark pixels with no height.
    values = median(block, dimension, 'omitnan');
    values = values(:).';
end

function values = smoothAlong(values, window)
% Moving median along the cut. movmedian keeps edges and NaN gaps intact, so a
% stretch of the cut that lies outside the cell stays a gap rather than being
% bridged by neighbouring values.
    if window <= 1; return; end
    values = movmedian(values, window, 'omitnan', 'Endpoints', 'shrink');
end

function limits = robustLimits(data, clip)
    values = double(data(:));
    values = values(isfinite(values));
    if isempty(values); limits = [0 1]; return; end
    limits = [quantile(values, clip(1)), quantile(values, clip(2))];
    if ~(limits(2) > limits(1)); limits = [min(values), max(values) + eps]; end
end

function index = clampIndex(value, limit)
    index = max(1, min(limit, round(value)));
end

function cd58 = loadCd58(acquisitionDir, imageSize)
    ptu = fullfile(acquisitionDir, 'RawImage.ptu');
    if ~isfile(ptu)
        warning('immune_cell_MIET_profile_figure:NoPtu', ...
            'No RawImage.ptu in %s; the CD58 profile will be blank.', acquisitionDir);
        cd58 = nan(imageSize);
        return;
    end
    try
        result = immune_cell_MIET_cd58_image(ptu, struct( ...
            'imageSize', imageSize, 'excitationNm', 485));
        cd58 = double(result.intensity);
    catch err
        warning('immune_cell_MIET_profile_figure:Cd58Failed', ...
            'CD58 unavailable (%s).', err.message);
        cd58 = nan(imageSize);
    end
end

function map = perceptualColormapOrBuiltin(name)
    switch lower(char(name))
        case {'viridis', 'magma', 'cividis', 'coolwarm', 'gray'}
            map = perceptualColormap(name, 256);
        otherwise
            map = feval(char(name), 256);
    end
end

function [rows, cols] = croppedRange(mask, margin)
    rowsAny = find(any(mask, 2));
    colsAny = find(any(mask, 1));
    if isempty(rowsAny); rowsAny = 1:size(mask, 1); end
    if isempty(colsAny); colsAny = 1:size(mask, 2); end
    rows = max(1, rowsAny(1) - margin):min(size(mask, 1), rowsAny(end) + margin);
    cols = max(1, colsAny(1) - margin):min(size(mask, 2), colsAny(end) + margin);
end

function folder = acquisitionFolder(heightMapsMat)
    folder = fileparts(heightMapsMat);
    for step = 1:3
        [parent, name] = fileparts(folder);
        if ~isempty(regexp(name, '^_\d{8}-\d{6}$', 'once')); return; end
        if isempty(parent) || strcmp(parent, folder); break; end
        folder = parent;
    end
end

function pxUm = pixelSizeFromHeader(acquisitionDir)
    pxUm = [];
    ptu = fullfile(acquisitionDir, 'RawImage.ptu');
    if ~isfile(ptu); return; end
    try
        head = PTU_Read_Head(ptu);
    catch
        return;
    end
    if ~isfield(head, 'ImgHdr_PixResol'); return; end
    value = double(head.ImgHdr_PixResol(1));
    if isfinite(value) && value > 0; pxUm = value; end
end

function value = calibField(calib, name)
    value = NaN;
    if isstruct(calib) && isfield(calib, 'params') && isfield(calib.params, name)
        value = double(calib.params.(name));
    elseif isstruct(calib) && isfield(calib, name)
        value = double(calib.(name));
    end
end

function text = prettyBinning(binning)
    switch binning
        case 'native';     text = 'native pixels';
        case 'sliding2x2'; text = '2x2 sliding TCSPC';
        case 'sliding4x4'; text = '4x4 sliding TCSPC';
        otherwise;         text = binning;
    end
end

function cfg = withDefaults(cfg, defaults)
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(cfg, names{index}) || isempty(cfg.(names{index}))
            cfg.(names{index}) = defaults.(names{index});
        end
    end
end
