function [hBand, hLine] = plot_shaded_mean_curve(ax, x, yMean, yStd, color, label, faceAlpha, showLine)
% PLOT_SHADED_MEAN_CURVE
% Plot mean +/- std as a shaded band plus line.

    if nargin < 7 || isempty(faceAlpha)
        faceAlpha = 0.18;
    end
    if nargin < 8 || isempty(showLine)
        showLine = true;
    end

    x = double(x(:));
    yMean = double(yMean(:));
    yStd = double(yStd(:));
    lo = yMean - yStd;
    hi = yMean + yStd;

    hold(ax, 'on');
    hBand = fill(ax, [x; flipud(x)], [lo; flipud(hi)], color, ...
        'FaceAlpha', faceAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    if showLine
        hLine = plot(ax, x, yMean, 'Color', color, 'LineWidth', 1.4, 'DisplayName', label);
    else
        hLine = gobjects(0);
    end
end
