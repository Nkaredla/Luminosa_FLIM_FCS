function showMap(map, titleText, unitText, colourMap)
%SHOWMAP Display one transposed FLIM map with transparent NaN pixels.

    imagesc(map', 'AlphaData', isfinite(map'));
    axis image off;
    colormap(gca, colourMap);
    title(titleText);
    set(gca, 'Color', [0.08 0.08 0.08]);
    colourBar = colorbar;
    if ~isempty(unitText)
        colourBar.Label.String = unitText;
    end
end
