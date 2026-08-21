function ax = plotDetectorHexMap(detXY, values, varargin)
%PLOTDETECTORHEXMAP Draw detector-channel values as hexagonal cells.
%
%   plotDetectorHexMap(detXY, values)
%   plotDetectorHexMap(..., 'CLim', [lo hi], 'CellRadius', r)
%
%   detXY is [nDetector x 2]. values is [nDetector x 1]. This is intended
%   for honeycomb detector layouts where a sparse imagesc grid is visually
%   misleading.

    p = inputParser;
    addParameter(p, 'Parent', []);
    addParameter(p, 'CLim', []);
    addParameter(p, 'CellRadius', []);
    addParameter(p, 'CellScale', 1.01);
    addParameter(p, 'EdgeColor', [0.55 0.55 0.55]);
    addParameter(p, 'LineWidth', 0.7);
    addParameter(p, 'BackgroundColor', [1 1 1]);
    parse(p, varargin{:});
    opts = p.Results;

    if isempty(opts.Parent)
        ax = gca;
    else
        ax = opts.Parent;
    end

    detXY = double(detXY);
    values = double(values(:));
    if size(detXY,2) ~= 2 || size(detXY,1) ~= numel(values)
        error('plotDetectorHexMap:BadInput', ...
            'detXY must be [nDetector x 2] and values must contain one value per detector.');
    end

    pitch = estimateNearestDetectorPitch(detXY);
    if isempty(opts.CellRadius)
        % Point-up hexagons touch for horizontal pitch p when R = p/sqrt(3).
        % A slight scale-up avoids antialiasing gaps in exported figures.
        cellRadius = opts.CellScale * pitch / sqrt(3);
    else
        cellRadius = opts.CellRadius;
    end

    cla(ax);
    hold(ax, 'on');
    set(ax, 'Color', opts.BackgroundColor);

    theta = (0:5) * pi/3 + pi/6;
    hx = cellRadius * cos(theta);
    hy = cellRadius * sin(theta);

    for k = 1:numel(values)
        if ~isfinite(values(k))
            continue;
        end

        patch(ax, ...
            detXY(k,1) + hx, ...
            detXY(k,2) + hy, ...
            values(k), ...
            'FaceColor', 'flat', ...
            'EdgeColor', opts.EdgeColor, ...
            'LineWidth', opts.LineWidth);
    end

    if ~isempty(opts.CLim)
        caxis(ax, opts.CLim);
    end

    axis(ax, 'equal');
    axis(ax, 'off');
    set(ax, 'YDir', 'normal');

    pad = max(cellRadius * 1.4, pitch * 0.25);
    xlim(ax, [min(detXY(:,1))-pad, max(detXY(:,1))+pad]);
    ylim(ax, [min(detXY(:,2))-pad, max(detXY(:,2))+pad]);
    hold(ax, 'off');
end

function pitch = estimateNearestDetectorPitch(detXY)
    n = size(detXY, 1);
    if n < 2
        pitch = 1;
        return;
    end

    nearest = nan(n, 1);
    for k = 1:n
        dx = detXY(:,1) - detXY(k,1);
        dy = detXY(:,2) - detXY(k,2);
        d = sqrt(dx.^2 + dy.^2);
        d = d(d > eps(max(d)));
        if ~isempty(d)
            nearest(k) = min(d);
        end
    end

    nearest = nearest(isfinite(nearest) & nearest > 0);
    if isempty(nearest)
        pitch = 1;
    else
        pitch = median(nearest);
    end
end
