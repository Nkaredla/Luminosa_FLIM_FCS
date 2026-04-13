function fig = luminosa_miet_nice_image(shadowSource, surfaceData, titleStr, varargin)
%LUMINOSA_MIET_NICE_IMAGE Render a shaded 3D MIET surface and side profile.

    opts = parseOptions(varargin{:});

    surfaceData = double(surfaceData);
    if isempty(surfaceData)
        error('Surface data is empty.');
    end

    dataMask = isfinite(surfaceData);
    if ~any(dataMask(:))
        error('Surface data contains no finite values.');
    end

    if nargin < 1 || isempty(shadowSource) || ~isequal(size(shadowSource), size(surfaceData))
        shadowSource = ones(size(surfaceData));
    else
        shadowSource = double(shadowSource);
    end
    shadowSource(~isfinite(shadowSource)) = 0;

    mask = dataMask;
    positiveMask = shadowSource > 0;
    if any(positiveMask(dataMask))
        mask = dataMask & positiveMask;
    end
    if ~any(mask(:))
        mask = dataMask;
    end

    if isempty(opts.ValueRange)
        vals = surfaceData(mask);
        valueRange = [min(vals) max(vals)];
    else
        valueRange = double(opts.ValueRange(:)).';
    end
    if numel(valueRange) ~= 2 || any(~isfinite(valueRange))
        vals = surfaceData(mask);
        valueRange = [min(vals) max(vals)];
    end
    if valueRange(2) <= valueRange(1)
        valueRange(2) = valueRange(1) + eps;
    end

    [cmap, shadowFraction] = luminosa_miet_colormap(opts.ColormapName, 256, 'nice');
    shadowData = buildShadowPlane(shadowSource, mask, valueRange, shadowFraction);
    surfacePlot = surfaceData;
    surfacePlot(~mask) = NaN;

    [xGrid, yGrid] = meshgrid(1:size(surfacePlot, 2), 1:size(surfacePlot, 1));
    zBase = ones(size(surfacePlot)) * valueRange(1);

    fig = figure('Name', titleStr, 'NumberTitle', 'off', 'Color', [0 0 0], 'Colormap', cmap);

    axBottom = axes('Parent', fig, ...
        'Position', [0.06522 0.06634 0.8934 0.1388], ...
        'Layer', 'top', ...
        'Color', [0 0 0], ...
        'XColor', [1 1 1], ...
        'YColor', [1 1 1], ...
        'ZColor', [1 1 1], ...
        'CLim', valueRange);
    view(axBottom, [0 0]);
    box(axBottom, 'on');
    hold(axBottom, 'on');
    surf(axBottom, xGrid, yGrid, surfacePlot, surfacePlot, ...
        'SpecularExponent', 15, ...
        'SpecularStrength', 0.6, ...
        'FaceLighting', 'phong', ...
        'LineStyle', 'none', ...
        'FaceColor', 'interp');
    axis(axBottom, 'tight');
    zlabel(axBottom, axisLabel(opts.Mode), 'Color', [1 1 1]);
    hold(axBottom, 'off');

    axTop = axes('Parent', fig, ...
        'Position', [0.08799 0.305 0.7588 0.6424], ...
        'PlotBoxAspectRatio', [3.5 3.5 1.156], ...
        'DataAspectRatio', [0.35 0.35 1], ...
        'Layer', 'top', ...
        'Color', [0 0 0], ...
        'XColor', [1 1 1], ...
        'YColor', [1 1 1], ...
        'ZColor', [1 1 1], ...
        'CLim', valueRange, ...
        'CameraViewAngle', 7.162);
    view(axTop, [-79 62]);
    hold(axTop, 'on');
    surf(axTop, xGrid, yGrid, zBase, shadowData, ...
        'LineStyle', 'none', ...
        'FaceLighting', 'none', ...
        'FaceColor', 'interp');
    light(axTop);
    surf(axTop, xGrid, yGrid, surfacePlot, surfacePlot, ...
        'SpecularExponent', 20, ...
        'SpecularStrength', 0.3, ...
        'DiffuseStrength', 0.3, ...
        'AmbientStrength', 0.7, ...
        'BackFaceLighting', 'unlit', ...
        'FaceLighting', 'phong', ...
        'LineStyle', 'none', ...
        'FaceColor', 'interp');
    axis(axTop, 'tight');
    title(axTop, titleStr, 'Color', [1 1 1], 'Interpreter', 'none');
    zlabel(axTop, axisLabel(opts.Mode), 'Color', [1 1 1]);
    hold(axTop, 'off');
end

function opts = parseOptions(varargin)
    opts = struct('Mode', 'height', 'ColormapName', 'MIET Bronze', 'ValueRange', []);
    if mod(numel(varargin), 2) ~= 0
        error('Options must be provided as name/value pairs.');
    end

    for idx = 1:2:numel(varargin)
        name = lower(char(string(varargin{idx})));
        value = varargin{idx + 1};
        switch name
            case 'mode'
                opts.Mode = char(string(value));
            case 'colormapname'
                opts.ColormapName = char(string(value));
            case 'valuerange'
                opts.ValueRange = value;
            otherwise
                error('Unknown option: %s', varargin{idx});
        end
    end
end

function shadowPlane = buildShadowPlane(shadowSource, mask, valueRange, shadowFraction)
    shadowPlane = nan(size(shadowSource));
    rangeWidth = valueRange(2) - valueRange(1);
    if rangeWidth <= 0
        rangeWidth = 1;
    end

    vals = double(shadowSource(mask));
    if isempty(vals)
        shadowNorm = zeros(size(shadowSource));
        shadowNorm(mask) = 1;
    else
        lo = prctile(vals, 5);
        hi = prctile(vals, 95);
        if ~isfinite(lo)
            lo = min(vals);
        end
        if ~isfinite(hi) || hi <= lo
            hi = max(vals);
        end
        if hi <= lo
            shadowNorm = zeros(size(shadowSource));
            shadowNorm(mask) = 1;
        else
            shadowNorm = (shadowSource - lo) / (hi - lo);
            shadowNorm = min(max(shadowNorm, 0), 1);
        end
    end

    shadowSpan = max(0.12, min(0.48, shadowFraction)) * rangeWidth;
    shadowPlane(mask) = valueRange(1) + shadowSpan * shadowNorm(mask);
end

function label = axisLabel(modeName)
    if strcmpi(modeName, 'lifetime')
        label = 'Lifetime (ns)';
    else
        label = 'Height (nm)';
    end
end
