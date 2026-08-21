function fig = plotSummary(trajectories, frames, detector, ...
        localizationTable, ~, curvature, diffusion, opts)
    import membrane_tracking.curved_miet.*

    fig = figure('Color', 'w', 'Name', ...
        'Curved membrane MIET particle tracking', ...
        'Position', [80 80 1250 850]);
    layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax1 = nexttile(layout, 1);
    gridVector = linspace(-opts.membraneRadiusUm, ...
        opts.membraneRadiusUm, 101);
    [gridX, gridY] = meshgrid(gridVector, gridVector);
    gridR2 = gridX.^2 + gridY.^2;
    gridZ = opts.tipHeightUm + 0.5 * opts.curvaturePerUm * gridR2;
    gridZ(gridR2 > opts.membraneRadiusUm^2) = NaN;
    gridLifetime = opts.lifetimeAtSubstrateNs + ...
        opts.lifetimeSlopeNsPerUm * gridZ;
    surf(ax1, gridX, gridY, gridZ, gridLifetime, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.78);
    hold(ax1, 'on');
    truePath = squeeze(trajectories.positionsXYZUm(1,:,:)).';
    plot3(ax1, truePath(:,1), truePath(:,2), truePath(:,3), ...
        '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 1.2);
    if height(localizationTable) > 0 && curvature.fitSucceeded
        scatter3(ax1, localizationTable.xUm, localizationTable.yUm, ...
            localizationTable.fittedSurfaceHeightUm, 16, ...
            localizationTable.lifetimeNs, 'filled', ...
            'MarkerEdgeColor', [0.25 0.25 0.25]);
    end
    plot3(ax1, 0, 0, opts.tipHeightUm, 'o', ...
        'MarkerSize', 7, 'MarkerFaceColor', [0.85 0.20 0.15], ...
        'MarkerEdgeColor', 'w');
    hold(ax1, 'off');
    axis(ax1, 'equal');
    axis(ax1, 'tight');
    grid(ax1, 'on');
    view(ax1, 38, 28);
    xlabel(ax1, 'x [um]');
    ylabel(ax1, 'y [um]');
    zlabel(ax1, 'height [um]');
    title(ax1, 'Intrinsic diffusion on the MIET membrane');
    colorbar1 = colorbar(ax1, 'Location', 'eastoutside');
    ylabel(colorbar1, 'lifetime [ns]');

    ax2 = nexttile(layout, 2);
    [~, brightFrame] = max(sum(frames,1));
    plotDetectorHexMap(detector.detXY, frames(:,brightFrame), ...
        'Parent', ax2, 'CellRadius', detector.detectorHexRadiusUm);
    colorbar2 = colorbar(ax2, 'Location', 'eastoutside');
    ylabel(colorbar2, 'photons');
    title(ax2, sprintf('23-channel microimage, frame %d', brightFrame));

    ax3 = nexttile(layout, 3);
    hold(ax3, 'on');
    if height(localizationTable) > 0
        radialPosition = localizationTable.radialPositionUm;
        scatter(ax3, radialPosition, localizationTable.lifetimeNs, ...
            18, localizationTable.heightSigmaUm, 'filled', ...
            'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0.35);
        if curvature.fitSucceeded
            fitRadius = linspace(0, max(radialPosition), 250).';
            fitHeight = curvature.tipHeightUm + ...
                0.5 * curvature.curvaturePerUm * fitRadius.^2;
            fitLifetime = opts.lifetimeAtSubstrateNs + ...
                opts.lifetimeSlopeNsPerUm * fitHeight;
            plot(ax3, fitRadius, fitLifetime, 'k-', 'LineWidth', 2);
        end
    end
    hold(ax3, 'off');
    grid(ax3, 'on');
    xlabel(ax3, 'radial detector position [um]');
    ylabel(ax3, 'fitted lifetime [ns]');
    title(ax3, 'MIET lifetime versus molecular position');
    if height(localizationTable) > 0
        colorbar3 = colorbar(ax3, 'Location', 'eastoutside');
        ylabel(colorbar3, 'height CRB sigma [um]');
    end

    ax4 = nexttile(layout, 4);
    msd = diffusion.msdTable;
    hold(ax4, 'on');
    validRaw = isfinite(msd.surfaceMsdUm2);
    validCorrected = isfinite(msd.surfaceMsdCorrectedUm2);
    if any(validRaw)
        plot(ax4, msd.lagTimeS(validRaw), msd.surfaceMsdUm2(validRaw), ...
            'o-', 'Color', [0.20 0.48 0.75], 'LineWidth', 1.2, ...
            'MarkerFaceColor', [0.20 0.48 0.75]);
    end
    if any(validCorrected)
        plot(ax4, msd.lagTimeS(validCorrected), ...
            msd.surfaceMsdCorrectedUm2(validCorrected), 's-', ...
            'Color', [0.83 0.29 0.20], 'LineWidth', 1.4, ...
            'MarkerFaceColor', [0.83 0.29 0.20]);
    end
    truthLineTime = [0; max(msd.lagTimeS)];
    plot(ax4, truthLineTime, 4*opts.diffusionUm2PerS*truthLineTime, ...
        'k--', 'LineWidth', 1.3);
    hold(ax4, 'off');
    grid(ax4, 'on');
    xlabel(ax4, 'lag time [s]');
    ylabel(ax4, 'surface MSD [um^2]');
    title(ax4, 'Surface-metric diffusion estimate');
    legend(ax4, {'raw surface MSD', 'Fisher-noise corrected', ...
        'ground truth 4D tau'}, 'Location', 'northwest');

    colormap(fig, parula(256));
    if curvature.fitSucceeded
        summaryTitle = sprintf(['MIET curvature %.3g (true %.3g) 1/um; ' ...
            'D_{MLE} %.3g (true %.3g) um^2/s'], ...
            curvature.curvaturePerUm, opts.curvaturePerUm, ...
            diffusion.DFisherMLEUm2PerS, opts.diffusionUm2PerS);
    else
        summaryTitle = 'Curvature was not estimable from accepted photons';
    end
    sgtitle(layout, summaryTitle, 'FontWeight', 'bold');
end
