function fig = plotSummary(trajectories, frames, detector, localizationTable, ...
        trackTable, diffusion, experimentalBnp, opts)
    import membrane_tracking.focused_ism.*

    fig = figure('Color', 'w', 'Position', [80 80 1100 780]);
    tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', ...
        'TileSpacing', 'compact');

    ax1 = nexttile(tl);
    hold(ax1, 'on');
    boxSize = trajectories.boxSizeUm;
    rectangle(ax1, 'Position', [-boxSize(1)/2, -boxSize(2)/2, ...
        boxSize(1), boxSize(2)], 'EdgeColor', [0.2 0.2 0.2], ...
        'LineWidth', 1);
    rectangle(ax1, 'Position', [-opts.laserWaistUm, -opts.laserWaistUm, ...
        2*opts.laserWaistUm, 2*opts.laserWaistUm], ...
        'Curvature', [1 1], 'EdgeColor', [0.85 0.33 0.10], ...
        'LineStyle', '--', 'LineWidth', 1.2);
    cmapMol = lines(size(trajectories.positionsUm, 1));
    for m = 1:size(trajectories.positionsUm, 1)
        x = squeeze(trajectories.positionsUm(m,1,:));
        y = squeeze(trajectories.positionsUm(m,2,:));
        [xPlot, yPlot] = breakPeriodicSegments(x, y, opts);
        plot(ax1, xPlot, yPlot, '-', 'Color', cmapMol(m,:), ...
            'LineWidth', 0.8);
        plot(ax1, x(1), y(1), 'o', 'Color', cmapMol(m,:), ...
            'MarkerSize', 3);
    end
    axis(ax1, 'image');
    xlim(ax1, [-boxSize(1)/2 boxSize(1)/2]);
    ylim(ax1, [-boxSize(2)/2 boxSize(2)/2]);
    xlabel(ax1, 'x [um]');
    ylabel(ax1, 'y [um]');
    title(ax1, sprintf('True membrane diffusion, D = %.3g um^2/s', ...
        opts.diffusionUm2PerS));
    grid(ax1, 'on');

    ax2 = nexttile(tl);
    frameSums = sum(frames, 1);
    [~, exampleFrame] = max(frameSums);
    microimageCounts = double(frames(:,exampleFrame));
    colorLimit = [0 max(max(microimageCounts), 1)];
    plotDetectorHexMap(detector.detXY, microimageCounts, ...
        'Parent', ax2, 'CLim', colorLimit, ...
        'CellRadius', 1.01 * detector.detectorHexRadiusUm, ...
        'EdgeColor', [0.35 0.35 0.35], 'LineWidth', 0.8, ...
        'BackgroundColor', [0.08 0.08 0.08]);
    colormap(ax2, 'hot');
    cb = colorbar(ax2);
    cb.Label.String = 'photons / frame';
    hold(ax2, 'on');
    plot(ax2, 0, 0, '+', 'Color', [0.2 0.85 1], ...
        'MarkerSize', 8, 'LineWidth', 1.2);
    locFrame = localizationTable(localizationTable.frame == exampleFrame, :);
    if height(locFrame) > 0
        plot(ax2, locFrame.xDetectorImageUm, locFrame.yDetectorImageUm, ...
            'wo', 'MarkerSize', 9, 'LineWidth', 1.4);
        plot(ax2, locFrame.xDetectorImageUm, locFrame.yDetectorImageUm, ...
            'w+', 'MarkerSize', 7, 'LineWidth', 1.2);
    end
    title(ax2, sprintf('%s microimage, frame %d, %d photons, K = %d', ...
        detector.layout, exampleFrame, round(frameSums(exampleFrame)), ...
        height(locFrame)));

    ax3 = nexttile(tl);
    hold(ax3, 'on');
    trackIds = unique(trackTable.trackId);
    cmapTrack = lines(max(numel(trackIds), 1));
    for k = 1:numel(trackIds)
        tr = trackTable(trackTable.trackId == trackIds(k), :);
        if height(tr) < opts.minTrackLength
            continue;
        end
        c = cmapTrack(1 + mod(k-1, size(cmapTrack, 1)), :);
        [xPlot, yPlot] = breakPeriodicSegments(tr.xUm, tr.yUm, opts);
        plot(ax3, xPlot, yPlot, '-o', 'Color', c, ...
            'MarkerSize', 2.5, 'LineWidth', 0.8);
    end
    rejected = trackTable.associationProbability < ...
        opts.minAssociationProbabilityForDiffusion | ...
        trackTable.emitterCountPosteriorProbability < ...
        opts.minEmitterCountPosteriorProbabilityForDiffusion;
    if any(rejected)
        plot(ax3, trackTable.xUm(rejected), trackTable.yUm(rejected), ...
            'rx', 'MarkerSize', 5, 'LineWidth', 1);
    end
    axis(ax3, 'image');
    xlim(ax3, opts.maxLocalizationRadiusUm * [-1 1]);
    ylim(ax3, opts.maxLocalizationRadiusUm * [-1 1]);
    xlabel(ax3, 'localized x [um]');
    ylabel(ax3, 'localized y [um]');
    title(ax3, sprintf('Linked tracks (%d used, %.1f%% rejected)', ...
        diffusion.nTracksUsed, ...
        100 * diffusion.fractionRejectedForDiffusion));
    grid(ax3, 'on');

    ax4 = nexttile(tl);
    msdTable = diffusion.msdTable;
    valid = isfinite(msdTable.msdCorrectedUm2) & msdTable.nPairs > 0;
    hold(ax4, 'on');
    if any(valid)
        plot(ax4, msdTable.lagTimeS(valid), msdTable.msdUm2(valid), ...
            'o-', 'Color', [0.45 0.45 0.45], 'LineWidth', 1);
        plot(ax4, msdTable.lagTimeS(valid), ...
            msdTable.msdCorrectedUm2(valid), 'o-', 'Color', [0 0.45 0.74], ...
            'LineWidth', 1.2);
        plot(ax4, msdTable.lagTimeS(valid), ...
            4 * opts.diffusionUm2PerS * msdTable.lagTimeS(valid), ...
            '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.2);
        legendEntries = {'raw MSD', 'Fisher-corrected MSD', ...
            'ground-truth 4D\tau'};
        if all(isfinite(diffusion.msdSlopeCorrected))
            tauFit = linspace(0, max(msdTable.lagTimeS(valid)), 100).';
            yFit = diffusion.msdSlopeCorrected(1) + ...
                diffusion.msdSlopeCorrected(2) * tauFit;
            plot(ax4, tauFit, yFit, '-', 'Color', [0.85 0.33 0.10], ...
                'LineWidth', 1.2);
            legendEntries{end+1} = 'short-lag zero-intercept fit';
        end
        legend(ax4, legendEntries, 'Location', 'northwest');
    end
    xlabel(ax4, 'lag time [s]');
    ylabel(ax4, 'MSD [um^2]');
    if experimentalBnp.estimateAvailable
        title(ax4, sprintf( ...
            'D_{MSD} %.3g, D_{step} %.3g, D_{joint} %.3g, true %.3g um^2/s', ...
            diffusion.DMsdCorrectedUm2PerS, ...
            diffusion.DFisherMLEUm2PerS, experimentalBnp.posteriorMeanD, ...
            opts.diffusionUm2PerS));
    else
        title(ax4, sprintf('D_{MSD} %.3g, D_{step} %.3g, true %.3g um^2/s', ...
            diffusion.DMsdCorrectedUm2PerS, ...
            diffusion.DFisherMLEUm2PerS, opts.diffusionUm2PerS));
    end
    grid(ax4, 'on');

    title(tl, ['Focused-excitation ISM diffusion, joint detector-channel ' ...
        'localization, Fisher information, and SPT'], ...
        'FontWeight', 'bold');
end
