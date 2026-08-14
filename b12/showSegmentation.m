function figureHandle = showSegmentation(mitoIntensity, b12Intensity, mitoMask, ...
        lysoMask, cellBoundary, compartmentMap, mitoRegionCount, ...
        lysoRegionCount, visibility)
%SHOWSEGMENTATION Display channel sources, boundaries, and compartment labels.

    figureHandle = figure('Name', 'B12 compartment segmentation', ...
        'Color', 'w', 'Visible', visibility);
    tiledlayout(figureHandle, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    showBoundaryOverlay(mitoIntensity, mitoMask, [1 0.25 0], false);
    title(sprintf('640 nm Det2: mitochondria (%d regions)', mitoRegionCount));

    nexttile;
    showBoundaryOverlay(b12Intensity, lysoMask, [0 0.8 1], false);
    title(sprintf('560 nm Det1 B12: lysosomes (%d spots)', lysoRegionCount));

    nexttile;
    imagesc(compartmentMap');
    axis image off;
    caxis([0 3]);
    colormap(gca, [0 0 0; 0.9 0.25 0.05; 0 0.75 0.95; 0.85 0.78 0.20]);
    title('1 mitochondria, 2 lysosome, 3 cytosol');

    nexttile;
    showBoundaryOverlay(b12Intensity, cellBoundary, [0.2 1 0.2], true);
    title('Filled cell boundary used for fitting');
end
