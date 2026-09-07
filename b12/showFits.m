function figureHandle = showFits(timeNs, names, curves, fittedCurves, ...
        exponentialCounts, lifetimes, messages, maximumPlots, visibility)
%SHOWFITS Display a selected, limited set of connected-region TCSPC fits.

    plotCount = min(numel(names), maximumPlots);
    columnCount = min(3, max(1, ceil(sqrt(plotCount))));
    rowCount = max(1, ceil(plotCount / columnCount));
    figureHandle = figure('Name', 'Selected connected-region B12 TCSPC fits', ...
        'Color', 'w', 'Visible', visibility);
    tiledlayout(figureHandle, rowCount, columnCount, ...
        'Padding', 'compact', 'TileSpacing', 'compact');

    for plotIndex = 1:plotCount
        nexttile;
        semilogy(timeNs, max(curves{plotIndex}, 0.5), 'k-', 'LineWidth', 1);
        hold on;
        if exponentialCounts(plotIndex) > 0
            fitToPlot = fittedCurves{plotIndex};
            finiteBins = isfinite(fitToPlot);
            fitToPlot(finiteBins) = max(fitToPlot(finiteBins), 0.5);
            semilogy(timeNs, fitToPlot, ...
                'r-', 'LineWidth', 1.4);
            validTau = lifetimes(plotIndex, isfinite(lifetimes(plotIndex, :)));
            tauText = strtrim(sprintf('%.2f ', validTau));
            title(sprintf('%s: %d-exp, tau [%s] ns', names{plotIndex}, ...
                exponentialCounts(plotIndex), tauText));
        else
            title(sprintf('%s: %s', names{plotIndex}, messages{plotIndex}), ...
                'Interpreter', 'none');
        end
        xlabel('Time (ns)');
        ylabel('Counts/bin');
        grid on;
    end
end
