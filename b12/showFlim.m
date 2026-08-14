function figureHandle = showFlim(b12Intensity, tauMean, tauStd, modelOrder, ...
        posteriorEvaluated, compartmentMap, visibility)
%SHOWFLIM Display B12 lifetime, uncertainty, and posterior coverage maps.

    figureHandle = figure('Name', 'B12 Bayesian FLIM', 'Color', 'w', 'Visible', visibility);
    tiledlayout(figureHandle, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    imagesc(log1p(double(b12Intensity')));
    axis image off; colorbar; title('B12 intensity (log counts)');

    nexttile;
    showMap(tauMean, 'Lifetime map (posterior or regional fallback)', 'ns', turbo(256));

    nexttile;
    showMap(tauStd, 'Pixel posterior lifetime uncertainty', 'ns', hot(256));

    nexttile;
    imagesc(modelOrder'); axis image off; caxis([0 2]); colorbar;
    colormap(gca, [0 0 0; 0.1 0.7 0.9; 0.95 0.35 0.1]);
    title('Selected model order');

    nexttile;
    imagesc(posteriorEvaluated'); axis image off; caxis([0 1]); colorbar;
    colormap(gca, [0.15 0.15 0.15; 0.2 0.9 0.35]);
    title('Green: pixel posterior evaluated');

    nexttile;
    imagesc(compartmentMap'); axis image off; caxis([0 3]); colorbar;
    colormap(gca, [0 0 0; 0.9 0.25 0.05; 0 0.75 0.95; 0.85 0.78 0.20]);
    title('Complete filled-cell assignment');
end
