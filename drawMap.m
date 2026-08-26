function handle = drawMap(ax, data, mapName)
%DRAWMAP Show one map with NaN rendered transparent.
%
% handle = drawMap(ax, data, mapName)
%
% AlphaData on the finite entries is what keeps unfitted pixels from being
% painted as the colormap's lowest colour, which would otherwise read as a
% genuine short lifetime rather than as no measurement.

    handle = imagesc(ax, data);
    set(handle, 'AlphaData', isfinite(data));
    colormap(ax, mapName);

    % Square pixels AND a square-for-square-data axes box. axis('image') alone
    % only fixes the DATA aspect ratio; inside a tiledlayout the axes BOX is
    % still stretched to whatever the tile gives it, and adding a colorbar
    % narrows it further, so the map ends up visibly distorted. Pinning
    % PlotBoxAspectRatio to the data dimensions is what actually prevents that.
    [nRow, nCol] = size(data);
    axis(ax, 'image');
    set(ax, 'PlotBoxAspectRatio', [nCol nRow 1]);
    xlim(ax, [0.5, nCol + 0.5]);
    ylim(ax, [0.5, nRow + 0.5]);

    set(ax, 'XTick', [], 'YTick', [], 'Color', [0.15 0.15 0.15]);
end
