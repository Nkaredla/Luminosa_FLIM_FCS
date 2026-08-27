function publicationMapPanel(ax, data, mapName, titleText, unitText, fontSize)
%PUBLICATIONMAPPANEL One map panel, styled for a figure rather than for debugging.
%
% publicationMapPanel(ax, data, mapName, titleText, unitText, fontSize)
%
% NaN is drawn transparent over a neutral ground, the axes carry no ticks or
% box, and the colour limits come from the 2nd-98th percentile of the FINITE
% data so a handful of extreme pixels cannot flatten everything else.
%
% caxis is set BEFORE the colorbar is created. Doing it the other way leaves the
% colorbar's limit-update listener pending, and it then fires inside
% exportgraphics mid-traversal and throws "Attempt to modify the tree during an
% update traversal".

    if nargin < 6 || isempty(fontSize); fontSize = 12; end

    handle = imagesc(ax, data);
    set(handle, 'AlphaData', isfinite(data));

    finite = data(isfinite(data));
    if numel(finite) > 10
        lo = quantileLocalBiexp(finite, 0.02);
        hi = quantileLocalBiexp(finite, 0.98);
        if hi > lo
            caxis(ax, [lo hi]);
        end
    end

    if ischar(mapName) || isstring(mapName)
        switch lower(char(mapName))
            case {'viridis', 'magma', 'cividis', 'coolwarm', 'gray'}
                colormap(ax, perceptualColormap(mapName));
            otherwise
                colormap(ax, char(mapName));
        end
    else
        colormap(ax, mapName);
    end

    % Square pixels AND a square-for-square-data box; axis('image') alone leaves
    % the axes box stretched by the tile and by the colorbar.
    [nRow, nCol] = size(data);
    axis(ax, 'image');
    set(ax, 'PlotBoxAspectRatio', [nCol nRow 1]);
    xlim(ax, [0.5, nCol + 0.5]);
    ylim(ax, [0.5, nRow + 0.5]);
    set(ax, 'XTick', [], 'YTick', [], 'Box', 'off', ...
        'Color', [0.97 0.97 0.97], 'FontSize', fontSize);

    c = colorbar(ax);
    set(c, 'FontSize', fontSize - 1, 'Box', 'off', 'TickDirection', 'out', ...
        'LineWidth', 0.8);
    if ~isempty(unitText)
        % Written as a tick-adjacent title rather than through c.Label, which
        % triggers the same listener problem noted above.
        set(get(c, 'Title'), 'String', unitText, 'FontSize', fontSize - 1);
    end

    title(ax, titleText, 'FontSize', fontSize + 1, 'FontWeight', 'normal');
end
