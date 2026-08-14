function showBoundaryOverlay(intensity, mask, colour, maskIsBoundary)
%SHOWBOUNDARYOVERLAY Display log intensity with a coloured boundary overlay.

    imagesc(log1p(double(intensity')));
    axis image off;
    colormap(gca, gray(256));
    hold on;
    if ~any(mask(:))
        return;
    end
    if maskIsBoundary
        boundary = mask;
    else
        boundary = bwperim(mask, 8);
    end
    [row, column] = find(boundary');
    plot(column, row, '.', 'Color', colour, 'MarkerSize', 3);
end
