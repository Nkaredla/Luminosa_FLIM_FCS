function out = immune_cell_MIET_topography_figure(heightMapsMat, cfg)
%IMMUNE_CELL_MIET_TOPOGRAPHY_FIGURE Rendered 3D membrane topography per cell.
%
% out = immune_cell_MIET_topography_figure(heightMapsMat)
% out = immune_cell_MIET_topography_figure(heightMapsMat, cfg)
%
% Renders one cell's MIET height map as a lit 3D surface, floating above a
% second layer that shows the CD58-Atto488 intensity image of the bilayer
% underneath it. One PNG per colormap.
%
% This follows makenima ("make nice images") from the MIET-GUI, which pairs a
% lit height surface with a shadow underneath, but replaces the shadow with a
% real measurement: the 485-excited CD58 image. The two layers are then
% physically what they look like - the ligand-presenting bilayer, and the cell
% membrane standing off it by the measured distance.
%
% TWO COLORMAPS IN ONE AXES
%
% makenima packs a grey ramp and a copper ramp into a single 64-entry colormap
% and splits it by CLim, because one axes could only hold one colormap. That
% constraint is gone, but an axes still has ONE colormap, so the CD58 layer is
% converted to true RGB with its own map and drawn as a texture. The axes
% colormap is then free for the height surface, and neither layer's colour
% scale distorts the other's.
%
% HOLES ARE FILLED FOR RENDERING ONLY
%
% A MIET height map has no value where the lifetime fell outside the invertible
% window, and a surface with holes cannot be lit convincingly. Those pixels are
% filled by diffusion from their neighbours and the whole surface is smoothed,
% purely so the rendering reads as a surface. The filled fraction is printed on
% the figure, the fill NEVER leaves this function, and nothing quantitative is
% taken from the rendered array.
%
% cfg fields
%   binning         'sliding4x4'
%   colormaps       {'turbo','parula','copper'}
%   cd58Colormap    'bone' for the 485 window
%   membraneColormap 'hot' for the 640 window (MemGlow)
%   outputDir       '' writes beside the height MAT
%   namePrefix      '' derived from the acquisition folder
%   smoothPixels    3.5 - Gaussian sigma for the render-only smoothing
%   renderDownsample 2 - block-average factor before rendering. A 4x4 sliding
%                   map is already correlated over 4 pixels, so rendering every
%                   pixel adds spikes without adding information
%   heightClip      [0.02 0.98] quantiles for the colour and z limits
%   zExaggeration   [] auto so the height range fills ~35% of the frame width
%   view            [-52 20]
%   layerZ          [-200 -700] - z of the two flat planes below the surface:
%                   the MIET height projection, then the CD58 image
%   showMembraneLayer false - add the 640 nm MemGlow image as a third plane
%   outlineColor    [1 1 1] - cell outline traced on the FLAT planes only
%   outlineWidth    1.6
%   clipLayersToFootprint true - drop the flat layers outside the cell outline
%                   so whatever is behind shows through
%   fontSize        14 - base font size; labels and titles scale from it
%   pixelSizeUm     [] read from the PTU header (ImgHdr_PixResol)
%   showFigures     false
%
% See also IMMUNE_CELL_MIET_CD58_IMAGE, MAKENIMA,
% RUN_BATCH_IMMUNE_CELL_MIET_TOPOGRAPHY.

    if nargin < 2 || isempty(cfg); cfg = struct(); end
    cfg = withDefaults(cfg, struct( ...
        'binning', 'sliding4x4', ...
        'colormaps', {{'turbo', 'parula', 'copper'}}, ...
        'cd58Colormap', 'bone', 'membraneColormap', 'hot', ...
        'outputDir', '', 'namePrefix', '', ...
        'smoothPixels', 3.5, 'renderDownsample', 2, ...
        'heightClip', [0.02 0.98], ...
        'zExaggeration', [], 'view', [-32 15], 'fontSize', 14, ...
        'layerZ', [-300 -600], 'showMembraneLayer', false, ...
        'outlineColor', [1 1 1], 'outlineWidth', 1.6, ...
        'clipLayersToFootprint', true, ...
        'pixelSizeUm', [], 'showFigures', true));

    heightMapsMat = char(heightMapsMat);
    if ~isfile(heightMapsMat)
        error('immune_cell_MIET_topography_figure:MissingMat', ...
            'Height MAT not found: %s', heightMapsMat);
    end
    loaded = load(heightMapsMat, 'heightMaps');
    heightMaps = loaded.heightMaps;
    if ~isfield(heightMaps.stages, cfg.binning)
        error('immune_cell_MIET_topography_figure:NoBinning', ...
            'Binning %s is not in %s', cfg.binning, heightMapsMat);
    end
    stage = heightMaps.stages.(cfg.binning);
    if ~stage.available
        error('immune_cell_MIET_topography_figure:StageUnavailable', ...
            'Binning %s is unavailable in %s', cfg.binning, heightMapsMat);
    end

    acquisitionDir = acquisitionFolder(heightMapsMat);
    [~, acquisition] = fileparts(acquisitionDir);
    if isempty(cfg.outputDir); cfg.outputDir = fileparts(heightMapsMat); end
    if ~isfolder(cfg.outputDir); mkdir(cfg.outputDir); end
    if isempty(cfg.namePrefix)
        cfg.namePrefix = sprintf('immune_cell_MIET_topography_%s_%s', ...
            regexprep(acquisition, '^_', ''), cfg.binning);
    end

    heightNm = double(stage.heightNm);
    heightNm(~stage.heightMask) = NaN;
    footprint = logical(stage.displayMask);
    if ~any(footprint(:)); footprint = isfinite(heightNm); end

    cd58 = loadPieImage(acquisitionDir, size(heightNm), 485);
    membrane = loadPieImage(acquisitionDir, size(heightNm), 640);
    if isempty(cfg.pixelSizeUm)
        cfg.pixelSizeUm = pixelSizeFromHeader(acquisitionDir);
    end

    [rows, cols] = croppedRange(footprint, 6);
    heightCrop = heightNm(rows, cols);
    maskCrop = footprint(rows, cols);
    cd58Crop = cd58(rows, cols);
    membraneCrop = membrane(rows, cols);

    valid = isfinite(heightCrop) & maskCrop;
    if nnz(valid) < 50
        error('immune_cell_MIET_topography_figure:TooFewPixels', ...
            'Only %d usable height pixels in %s', nnz(valid), heightMapsMat);
    end
    filledFraction = 1 - nnz(valid) / max(nnz(maskCrop), 1);
    rendered = fillAndSmooth(heightCrop, valid, maskCrop, cfg.smoothPixels);
    measured = double(valid);
    [rendered, maskCrop, cd58Crop, measured, pixelScale] = ...
        downsampleForRender(rendered, maskCrop, cd58Crop, measured, ...
        cfg.renderDownsample);
    membraneCrop = blockMean(membraneCrop, cfg.renderDownsample);
    cfg.pixelSizeUm = cfg.pixelSizeUm * pixelScale;
    if filledFraction > 0.5
        warning('immune_cell_MIET_topography_figure:MostlyFilled', ...
            ['%s: only %.0f%% of the footprint was inverted, so most of the ' ...
             'rendered surface is interpolated. The figure marks it, but read ' ...
             'the shape with care.'], acquisition, 100 * (1 - filledFraction));
    end

    limits = [quantile(heightCrop(valid), cfg.heightClip(1)), ...
              quantile(heightCrop(valid), cfg.heightClip(2))];
    if ~(limits(2) > limits(1)); limits = [min(heightCrop(valid)) max(heightCrop(valid))]; end

    out = struct();
    out.heightMapsMat = heightMapsMat;
    out.acquisition = acquisition;
    out.binning = cfg.binning;
    out.filledFraction = filledFraction;
    out.heightLimitsNm = limits;
    out.medianHeightNm = median(heightCrop(valid));
    out.figures = {};

    for index = 1:numel(cfg.colormaps)
        name = cfg.colormaps{index};
        file = fullfile(cfg.outputDir, sprintf('%s_%s.png', cfg.namePrefix, name));
        renderOne(rendered, maskCrop, cd58Crop, membraneCrop, measured, ...
            limits, name, cfg, heightMaps, stage, acquisition, ...
            filledFraction, file);
        out.figures{end + 1} = file; %#ok<AGROW>
        fprintf('  %s\n', file);
    end
end

% ------------------------------------------------------------------ rendering

function renderOne(height, mask, cd58, membrane, measured, limits, mapName, ...
        cfg, heightMaps, stage, acquisition, filledFraction, outputFile)
    height(~mask) = NaN;
    [nRow, nCol] = size(height);
    fs = cfg.fontSize;
    pxUm = cfg.pixelSizeUm;
    if isempty(pxUm) || ~isfinite(pxUm) || pxUm <= 0; pxUm = NaN; end
    if isfinite(pxUm)
        xUm = (0:nCol - 1) * pxUm;  yUm = (0:nRow - 1) * pxUm;
        lateralLabel = '\mum';
    else
        xUm = 0:nCol - 1;           yUm = 0:nRow - 1;
        lateralLabel = 'pixels';
    end
    [X, Y] = meshgrid(xUm, yUm);
    xSpan = max(max(xUm) - min(xUm), eps);
    ySpan = max(max(yUm) - min(yUm), eps);

    zRange = max(limits(2) - limits(1), eps);
    zTop = limits(2) + 0.08 * zRange;
    layerZ = sort(double(cfg.layerZ), 'descend');
    wantedLayers = 2 + double(cfg.showMembraneLayer);
    if numel(layerZ) < wantedLayers
        % Continue the spacing of the last given pair rather than guessing.
        step = abs(diff(layerZ(max(1, end - 1):end)));
        if isempty(step) || step == 0; step = 250; end
        while numel(layerZ) < wantedLayers
            layerZ(end + 1) = layerZ(end) - step; %#ok<AGROW>
        end
    end
    layerZ = layerZ(1:wantedLayers);
    zBottom = min(layerZ) - 0.10 * (zTop - min(layerZ));
    zStack = zTop - zBottom;
    % The whole stack - topography plus both flat layers - is scaled to a fixed
    % fraction of the lateral extent. Dropping a plane further down lengthens
    % the column, so the fraction has to be generous or the relief itself gets
    % squashed to make room for empty space between the planes.
    if isempty(cfg.zExaggeration)
        zAspect = zStack / (0.85 * max(xSpan, ySpan));
    else
        zAspect = 1 / cfg.zExaggeration;
    end
    if isfinite(pxUm)
        zExagFactor = 1000 / zAspect;
    else
        zExagFactor = NaN;
    end

    visibility = 'off';
    if cfg.showFigures; visibility = 'on'; end
    % MATLAB silently CLAMPS a figure larger than the screen and the layout is
    % then squeezed by a different factor in each direction - a 900x2700
    % request came back as 900x808. Ask for a size that fits.
    screen = get(groot, 'ScreenSize');
    wanted = [1250 1150];
    available = [screen(3) - 80, screen(4) - 120];
    figureSize = [min(wanted(1), available(1)), min(wanted(2), available(2))];
    h = figure('Color', 'k', 'Visible', visibility, ...
        'Position', [40 40 figureSize], 'InvertHardcopy', 'off');
    ax = axes('Parent', h, 'Position', [0.08 0.07 0.74 0.83]);
    hold(ax, 'on');

    % ---- flat layers, as textures at the requested z ---------------------
    % Order matters and is the caller's: directly under the surface sits the
    % SAME height data as a flat projection, so the 3D relief and its map can
    % be read against each other, and the 485 nm CD58 image goes underneath.
    heightRgb = heightProjectionRgb(height, mask, limits, mapName);
    layerRgb = { heightRgb, toRgb(cd58, cfg.cd58Colormap, [0.02 0.99]) };
    layerName = { 'MIET height, 2D projection', 'CD58-Atto488, 485 nm' };
    if cfg.showMembraneLayer
        layerRgb{end + 1} = toRgb(membrane, cfg.membraneColormap, [0.02 0.99]);
        layerName{end + 1} = 'MemGlow, 640 nm';
    end
    % The cell outline goes on the FLAT planes only. On the 3D relief it would
    % have to follow the surface in z, where it reads as a contour line of the
    % topography rather than as the footprint boundary.
    outline = maskOutline(xUm, yUm, mask);
    % Clipping is done with NaN in ZData, not with AlphaData: a surface with
    % NaN there simply has no face, which is true transparency and independent
    % of the renderer. Surface alpha silently dropped the whole object under
    % headless software OpenGL when it was tried on the relief.
    planeZ = ones(nRow, nCol);
    if cfg.clipLayersToFootprint
        planeZ(~mask) = NaN;
    end
    for index = 1:numel(layerRgb)
        surface(ax, X, Y, layerZ(index) * planeZ, layerRgb{index}, ...
            'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
            'FaceLighting', 'none', 'AmbientStrength', 1);
        for piece = 1:numel(outline)
            segment = outline{piece};
            plot3(ax, segment(1, :), segment(2, :), ...
                layerZ(index) * ones(1, size(segment, 2)), '-', ...
                'Color', cfg.outlineColor, 'LineWidth', cfg.outlineWidth);
        end
    end

    % ---- membrane topography --------------------------------------------
    faceRgb = desaturateFilled(height, measured, mask, limits, mapName);
    surf(ax, X, Y, height, faceRgb, 'EdgeColor', 'none', 'FaceColor', 'interp', ...
        'FaceLighting', 'gouraud', 'BackFaceLighting', 'unlit', ...
        'AmbientStrength', 0.45, 'DiffuseStrength', 0.75, ...
        'SpecularStrength', 0.28, 'SpecularExponent', 18);

    colormap(ax, resolveColormap(mapName));
    caxis(ax, limits); %#ok<CAXIS>
    view(ax, cfg.view);
    camlight(ax, -30, 55); camlight(ax, 35, 18);
    material(ax, 'dull');
    % Both aspect ratios must be pinned. DataAspectRatio alone keeps the DATA
    % isotropic but lets the plot box stretch to the axes rectangle, which is
    % what made a non-square crop look elongated.
    daspect(ax, [1 1 zAspect]);
    pbaspect(ax, [xSpan, ySpan, zStack / zAspect]);
    % Image convention: row 1 is the TOP of the map, so y must run downwards.
    set(ax, 'YDir', 'reverse', 'Color', 'k', 'XColor', [0.8 0.8 0.8], ...
        'YColor', [0.8 0.8 0.8], 'ZColor', [0.8 0.8 0.8], ...
        'GridColor', [0.45 0.45 0.45], 'Box', 'off', 'FontSize', fs, ...
        'LineWidth', 1.1);
    grid(ax, 'on');
    xlim(ax, [min(xUm) max(xUm)]); ylim(ax, [min(yUm) max(yUm)]);
    zlim(ax, [zBottom, zTop]);
    % Tick the height range only. The negative part of the axis exists purely
    % to position the display planes - labelling it would invite reading -700
    % as a measured depth, when the CD58 layer is a flat image parked there.
    set(ax, 'ZTick', heightTicks(zTop));
    xlabel(ax, sprintf('x [%s]', lateralLabel), 'FontSize', fs + 1);
    ylabel(ax, sprintf('y [%s]', lateralLabel), 'FontSize', fs + 1);
    zlabel(ax, 'z [nm]', 'FontSize', fs + 1);
    barMain = colorbar(ax, 'Color', [0.88 0.88 0.88], 'FontSize', fs, ...
        'Position', [0.90 0.30 0.022 0.40]);
    barMain.Label.String = 'MIET height [nm]';
    barMain.Label.Color = [0.92 0.92 0.92];
    barMain.Label.FontSize = fs + 1;

    calib = heightMaps.calibration;
    annotation(h, 'textbox', [0.02 0.925 0.96 0.06], 'Color', 'w', ...
        'String', sprintf('%s    MIET topography over CD58', ...
            regexprep(acquisition, '^_', '')), ...
        'EdgeColor', 'none', 'FontSize', fs + 4, 'FontWeight', 'bold', ...
        'Interpreter', 'none', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');
    annotation(h, 'textbox', [0.02 0.885 0.96 0.045], 'Color', [0.78 0.78 0.78], ...
        'String', sprintf(['%s  |  median %.1f nm  |  z x%.0f  |  %.0f%% ' ...
            'inverted, rest interpolated  |  tau_0 %.2f ns, qy %.3f'], ...
            prettyBinning(cfg.binning), median(height(isfinite(height))), ...
            zExagFactor, 100 * (1 - filledFraction), ...
            calibField(calib, 'tauFreeNs'), calibField(calib, 'quantumYield')), ...
        'EdgeColor', 'none', 'FontSize', fs, 'Interpreter', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    % A fixed-position legend for the three layers. Text placed in DATA
    % coordinates disappeared behind the planes as soon as the view rotated.
    legendLines = { 'MIET membrane topography   (3D relief)' };
    for index = 1:numel(layerName)
        legendLines{end + 1} = sprintf('%-28s z = %g nm', ...
            layerName{index}, layerZ(index)); %#ok<AGROW>
    end
%     annotation(h, 'textbox', [0.015 0.60 0.32 0.22], 'Color', [0.92 0.92 0.92], ...
%         'String', legendLines, 'EdgeColor', 'none', 'FontSize', fs, ...
%         'Interpreter', 'none', 'HorizontalAlignment', 'left', ...
%         'VerticalAlignment', 'top', 'FitBoxToText', 'off');

    % exportgraphics rescaled width and height by different factors on this
    % figure (1000x2100 came out 1056x4106), which squashed every panel. print
    % with PaperPositionMode 'auto' reproduces the on-screen geometry exactly.
    set(h, 'PaperPositionMode', 'auto');
    print(h, outputFile, '-dpng', '-r150');
    if ~cfg.showFigures; close(h); end
end

function ticks = heightTicks(zTop)
    span = max(zTop, eps);
    raw = span / 4;
    magnitude = 10^floor(log10(raw));
    normalised = raw / magnitude;
    if normalised < 1.5;      step = 1 * magnitude;
    elseif normalised < 3.5;  step = 2 * magnitude;
    elseif normalised < 7.5;  step = 5 * magnitude;
    else;                     step = 10 * magnitude;
    end
    ticks = 0:step:zTop;
end

function segments = maskOutline(xUm, yUm, mask)
% Footprint boundary as polylines in data coordinates. contourc at the 0.5
% level of the logical mask is base MATLAB - bwboundaries would pull in the
% Image Processing Toolbox for something this simple, and it returns pixel
% indices that would still have to be mapped back to microns.
    segments = {};
    if ~any(mask(:)) || any(size(mask) < 2); return; end
    matrix = contourc(xUm, yUm, double(mask), [0.5 0.5]);
    column = 1;
    while column < size(matrix, 2)
        count = matrix(2, column);
        if count < 2 || column + count > size(matrix, 2); break; end
        segments{end + 1} = matrix(:, column + 1:column + count); %#ok<AGROW>
        column = column + count + 1;
    end
end

function rgb = heightProjectionRgb(height, mask, limits, mapName)
% The flat height map uses the SAME colormap and colour limits as the relief
% above it, so the colourbar serves both. Pixels outside the cell are painted
% the figure background rather than the colormap's lowest colour, which would
% read as a real, very close contact.
    map = imresizeMap(resolveColormap(mapName), 256);
    scaled = (height - limits(1)) / max(limits(2) - limits(1), eps);
    scaled(~isfinite(scaled)) = 0;
    scaled = min(max(scaled, 0), 1);
    rgb = ind2rgb(uint8(round(scaled * 255)) + 1, map);
    outside = repmat(~mask | ~isfinite(height), 1, 1, 3);
    rgb(outside) = 0;
end

function rgb = desaturateFilled(height, measured, mask, limits, mapName)
% Colour by height through the chosen colormap, then blend interpolated pixels
% toward a neutral grey in proportion to how much of the block was filled.
    map = imresizeMap(resolveColormap(mapName), 256);
    scaled = (height - limits(1)) / max(limits(2) - limits(1), eps);
    scaled(~isfinite(scaled)) = 0;
    scaled = min(max(scaled, 0), 1);
    rgb = ind2rgb(uint8(round(scaled * 255)) + 1, map);
    weight = 0.35 + 0.65 * min(max(measured, 0), 1);
    weight(~mask) = 1;
    neutral = 0.38;
    rgb = rgb .* weight + neutral * (1 - weight);
end

function rgb = toRgb(image, mapName, clipQuantiles)
% Convert an intensity image to RGB through its own colormap, so it can be
% drawn as a texture in an axes whose colormap belongs to something else.
    image = double(image);
    finiteValues = image(isfinite(image));
    if isempty(finiteValues); finiteValues = 0; end
    lo = quantile(finiteValues, clipQuantiles(1));
    hi = quantile(finiteValues, clipQuantiles(2));
    if ~(hi > lo); lo = min(finiteValues); hi = max(finiteValues) + eps; end
    scaled = (image - lo) / (hi - lo);
    scaled(~isfinite(scaled)) = 0;
    scaled = min(max(scaled, 0), 1);
    map = resolveColormap(mapName);
    rgb = ind2rgb(uint8(round(scaled * 255)) + 1, imresizeMap(map, 256));
end

function map = imresizeMap(map, n)
    if size(map, 1) == n; return; end
    map = interp1(linspace(0, 1, size(map, 1)), map, linspace(0, 1, n), 'linear');
end

function map = resolveColormap(name)
    name = char(name);
    switch lower(name)
        case {'viridis', 'magma', 'cividis', 'coolwarm', 'gray'}
            map = perceptualColormap(name, 256);
        otherwise
            map = feval(name, 256);
    end
end

% ------------------------------------------------------------------- helpers

function filled = fillAndSmooth(height, valid, mask, sigma)
% Diffusion fill: repeatedly replace missing pixels by the mean of their
% neighbours. Cheap, edge-preserving enough for a rendering, and it cannot
% invent values outside the range of the surrounding data.
    filled = height;
    filled(~valid) = NaN;
    kernel = [0 1 0; 1 0 1; 0 1 0] / 4;
    for iteration = 1:200
        missing = ~isfinite(filled) & mask;
        if ~any(missing(:)); break; end
        padded = filled;
        padded(~isfinite(padded)) = 0;
        weight = double(isfinite(filled));
        numerator = conv2(padded, kernel, 'same');
        denominator = conv2(weight, kernel, 'same');
        candidate = numerator ./ max(denominator, eps);
        take = missing & denominator > 0;
        filled(take) = candidate(take);
    end
    remaining = ~isfinite(filled) & mask;
    if any(remaining(:))
        filled(remaining) = median(height(valid));
    end
    if sigma > 0
        filled(~isfinite(filled)) = median(height(valid));
        filled = imgaussfiltLocal(filled, sigma);
    end
    filled(~mask) = NaN;
end

function smoothed = imgaussfiltLocal(image, sigma)
    radius = max(1, ceil(3 * sigma));
    grid = -radius:radius;
    kernel = exp(-(grid.^2) / (2 * sigma^2));
    kernel = kernel / sum(kernel);
    smoothed = conv2(kernel', kernel, padReplicate(image, radius), 'valid');
end

function padded = padReplicate(image, radius)
    padded = image([ones(1, radius), 1:end, end * ones(1, radius)], ...
                   [ones(1, radius), 1:end, end * ones(1, radius)]);
end

function [height, mask, cd58, measured, factor] = ...
        downsampleForRender(height, mask, cd58, measured, factor)
    factor = max(1, round(factor));
    if factor == 1; return; end
    [nRow, nCol] = size(height);
    rows = 1:factor * floor(nRow / factor);
    cols = 1:factor * floor(nCol / factor);
    block = @(A, fun) fun(reshape(permute(reshape(A(rows, cols), factor, [], ...
        factor, numel(cols) / factor), [1 3 2 4]), factor * factor, []), 1);
    shape = [numel(rows) / factor, numel(cols) / factor];
    heightBlocks = reshape(block(height, @(B, d) mean(B, d, 'omitnan')), shape);
    maskBlocks = reshape(block(double(mask), @(B, d) mean(B, d)), shape) > 0.5;
    cd58Blocks = reshape(block(cd58, @(B, d) mean(B, d, 'omitnan')), shape);
    measuredBlocks = reshape(block(measured, @(B, d) mean(B, d)), shape);
    height = heightBlocks; mask = maskBlocks; cd58 = cd58Blocks;
    measured = measuredBlocks;
end

function [rows, cols] = croppedRange(mask, margin)
    rowsAny = find(any(mask, 2));
    colsAny = find(any(mask, 1));
    if isempty(rowsAny); rowsAny = 1:size(mask, 1); end
    if isempty(colsAny); colsAny = 1:size(mask, 2); end
    rows = max(1, rowsAny(1) - margin):min(size(mask, 1), rowsAny(end) + margin);
    cols = max(1, colsAny(1) - margin):min(size(mask, 2), colsAny(end) + margin);
end

function image = loadPieImage(acquisitionDir, imageSize, excitationNm)
    ptu = fullfile(acquisitionDir, 'RawImage.ptu');
    if ~isfile(ptu)
        warning('immune_cell_MIET_topography_figure:NoPtu', ...
            'No RawImage.ptu in %s; the %g nm layer will be blank.', ...
            acquisitionDir, excitationNm);
        image = zeros(imageSize);
        return;
    end
    try
        result = immune_cell_MIET_cd58_image(ptu, struct( ...
            'imageSize', imageSize, 'excitationNm', excitationNm));
        image = double(result.intensity);
    catch err
        warning('immune_cell_MIET_topography_figure:PieImageFailed', ...
            'The %g nm layer is unavailable (%s); drawing a blank panel.', ...
            excitationNm, err.message);
        image = zeros(imageSize);
    end
end

function out = blockMean(image, factor)
    factor = max(1, round(factor));
    if factor == 1; out = image; return; end
    [nRow, nCol] = size(image);
    rows = 1:factor * floor(nRow / factor);
    cols = 1:factor * floor(nCol / factor);
    trimmed = image(rows, cols);
    shape = [numel(rows) / factor, numel(cols) / factor];
    stacked = reshape(permute(reshape(trimmed, factor, shape(1), factor, ...
        shape(2)), [1 3 2 4]), factor * factor, []);
    out = reshape(mean(stacked, 1, 'omitnan'), shape);
end

function folder = acquisitionFolder(heightMapsMat)
% .../<acquisition>/<resultsFolder>/<cfgHash>/immune_cell_MIET_height_maps*.mat
    folder = fileparts(heightMapsMat);
    for step = 1:3
        [parent, name] = fileparts(folder);
        if ~isempty(regexp(name, '^_\d{8}-\d{6}$', 'once')); return; end
        if isempty(parent) || strcmp(parent, folder); break; end
        folder = parent;
    end
end

function pxUm = pixelSizeFromHeader(acquisitionDir)
% ImgHdr_PixResol is exactly what batch_immune_cell_MIET records as
% pixelSizeUm, and reading it here avoids the batch summary CSV, whose quoted
% free-text columns contain commas that defeat readtable's header detection.
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
    if isstruct(calib) && isfield(calib, 'params') && ...
            isfield(calib.params, name)
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

function value = pxUm0OrNaN(cfg)
    value = NaN;
    if isfield(cfg, 'pixelSizeUm') && ~isempty(cfg.pixelSizeUm)
        value = double(cfg.pixelSizeUm);
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
