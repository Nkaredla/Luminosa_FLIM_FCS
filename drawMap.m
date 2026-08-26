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
    axis(ax, 'image');
    colormap(ax, mapName);
    set(ax, 'XTick', [], 'YTick', [], 'Color', [0.15 0.15 0.15]);
end
