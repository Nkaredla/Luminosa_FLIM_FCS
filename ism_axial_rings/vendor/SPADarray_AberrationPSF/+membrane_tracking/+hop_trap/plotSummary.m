function fig = plotSummary(trajectories, mesh, frames, detector, ...
        localizationTable, trackTable, stepTable, msdTable, ...
        modelComparison, stickyBnp, classification, opts)
%PLOTSUMMARY Plot motion truth, ISM data, tracks, and model evidence.

    import membrane_tracking.hop_trap.*

    fig = figure('Color', 'w', 'Position', [60 60 1350 820]);
    layout = tiledlayout(fig, 2, 3, 'Padding', 'compact', ...
        'TileSpacing', 'compact');

    ax = nexttile(layout);
    hold(ax, 'on');
    plotPeriodicVoronoiMesh(ax, mesh, opts.boxSizeUm);
    color = lines(max(opts.nMolecules, 1));
    for molecule = 1:opts.nMolecules
        x = squeeze(trajectories.positionsUm(molecule,1,:));
        y = squeeze(trajectories.positionsUm(molecule,2,:));
        [xPlot, yPlot] = membrane_tracking.focused_ism.breakPeriodicSegments( ...
            x, y, opts);
        plot(ax, xPlot, yPlot, '-', 'Color', color(molecule,:), ...
            'LineWidth', 0.8);
        trapped = trajectories.trapState(molecule,:) == 2;
        if any(trapped)
            plot(ax, x(trapped), y(trapped), '.', ...
                'Color', [0.84 0.15 0.16], 'MarkerSize', 7);
        end
    end
    rectangle(ax, 'Position', [-opts.laserWaistUm, -opts.laserWaistUm, ...
        2*opts.laserWaistUm, 2*opts.laserWaistUm], ...
        'Curvature', [1 1], 'LineStyle', '--', ...
        'EdgeColor', [0.85 0.33 0.10], 'LineWidth', 1.1);
    axis(ax, 'image');
    xlim(ax, 0.5 * opts.boxSizeUm(1) * [-1 1]);
    ylim(ax, 0.5 * opts.boxSizeUm(2) * [-1 1]);
    xlabel(ax, 'x [um]');
    ylabel(ax, 'y [um]');
    title(ax, sprintf('Truth: %s, preset %s', ...
        opts.diffusionMode, opts.referencePresetResolved));
    grid(ax, 'on');

    ax = nexttile(layout);
    frameSums = sum(frames, 1);
    [~, exampleFrame] = max(frameSums);
    counts = double(frames(:,exampleFrame));
    plotDetectorHexMap(detector.detXY, counts, 'Parent', ax, ...
        'CLim', [0 max(max(counts), 1)], ...
        'CellRadius', 1.01 * detector.detectorHexRadiusUm, ...
        'EdgeColor', [0.35 0.35 0.35], 'LineWidth', 0.8, ...
        'BackgroundColor', [0.08 0.08 0.08]);
    colormap(ax, 'hot');
    colorbar(ax);
    hold(ax, 'on');
    loc = localizationTable(localizationTable.frame == exampleFrame, :);
    if height(loc) > 0
        plot(ax, loc.xDetectorImageUm, loc.yDetectorImageUm, ...
            'wo', 'MarkerSize', 9, 'LineWidth', 1.4);
    end
    title(ax, sprintf('Focused-ISM microimage, frame %d', exampleFrame));

    ax = nexttile(layout);
    hold(ax, 'on');
    trackIds = unique(trackTable.trackId).';
    colors = lines(max(numel(trackIds), 1));
    for k = 1:numel(trackIds)
        tr = sortrows(trackTable(trackTable.trackId == trackIds(k), :), ...
            'frame');
        if height(tr) < opts.minTrackLength
            continue;
        end
        [xPlot, yPlot] = membrane_tracking.focused_ism.breakPeriodicSegments( ...
            tr.xUm, tr.yUm, opts);
        plot(ax, xPlot, yPlot, '-o', 'Color', colors(k,:), ...
            'MarkerSize', 2, 'LineWidth', 0.8);
    end
    axis(ax, 'image');
    radius = opts.maxLocalizationRadiusUm;
    xlim(ax, radius * [-1 1]);
    ylim(ax, radius * [-1 1]);
    xlabel(ax, 'localized x [um]');
    ylabel(ax, 'localized y [um]');
    title(ax, sprintf('%d localizations, %d independent steps', ...
        height(localizationTable), height(stepTable)));
    grid(ax, 'on');

    ax = nexttile(layout);
    hold(ax, 'on');
    valid = isfinite(msdTable.correctedMsdUm2);
    errorbar(ax, msdTable.lagTimeS(valid), ...
        msdTable.correctedMsdUm2(valid), ...
        msdTable.standardErrorUm2(valid), 'o', ...
        'Color', [0.20 0.20 0.20], 'MarkerSize', 4);
    hop = modelComparison.hopMsd;
    if hop.estimateAvailable
        plot(ax, hop.lagTimeS, hop.freePredictionUm2, '--', ...
            'Color', [0 0.45 0.74], 'LineWidth', 1.2);
        plot(ax, hop.lagTimeS, hop.hopPredictionUm2, '-', ...
            'Color', [0.85 0.33 0.10], 'LineWidth', 1.4);
        legend(ax, {'corrected MSD', 'free model', 'hop model'}, ...
            'Location', 'northwest');
    end
    xlabel(ax, 'lag time [s]');
    ylabel(ax, 'MSD [um^2]');
    title(ax, sprintf('MSD confinement evidence: Delta qBIC = %.2f', ...
        modelComparison.deltaHopQuasiBic));
    grid(ax, 'on');

    ax = nexttile(layout);
    evidence = [modelComparison.deltaTrapBic, ...
        modelComparison.deltaHopQuasiBic];
    evidenceBars = bar(ax, evidence, 0.62, 'FaceColor', 'flat');
    evidenceBars.CData = [0.30 0.55 0.75; 0.85 0.42 0.15];
    hold(ax, 'on');
    yline(ax, opts.bicEvidenceThreshold, 'k--', 'threshold');
    yline(ax, 0, '-', 'Color', [0.4 0.4 0.4]);
    set(ax, 'XTick', [1 2], ...
        'XTickLabel', {'trap step BIC', 'hop MSD qBIC'});
    ylabel(ax, 'evidence difference');
    title(ax, sprintf('Classification: %s', classification.label));
    grid(ax, 'on');

    ax = nexttile(layout);
    hold(ax, 'on');
    stateHandles = gobjects(0);
    stateLabels = cell(0, 1);
    if modelComparison.trapHmm.estimateAvailable
        stateHandles(end+1) = plot(ax, stepTable.timeCenterS, ...
            modelComparison.trapHmm.posteriorSlowProbability, '-', ...
            'Color', [0.20 0.55 0.30], 'LineWidth', 1.0);
        stateLabels{end+1} = 'two-state HMM';
    end
    if stickyBnp.estimateAvailable
        stateHandles(end+1) = plot(ax, stepTable.timeCenterS, ...
            stickyBnp.posteriorSlowStateProbability, '-', ...
            'Color', [0.60 0.20 0.65], 'LineWidth', 1.0);
        stateLabels{end+1} = 'sticky BNP';
    end
    if ~isempty(stateHandles)
        legend(ax, stateHandles, stateLabels, 'Location', 'best');
    end
    xlabel(ax, 'time [s]');
    ylabel(ax, 'P(slow state)');
    ylim(ax, [0 1]);
    if stickyBnp.estimateAvailable
        title(ax, sprintf('BNP state count: mode %d, mean %.2f', ...
            stickyBnp.posteriorModeStateCount, ...
            stickyBnp.posteriorMeanStateCount));
    else
        title(ax, sprintf('Mobility states: %s', stickyBnp.status));
    end
    grid(ax, 'on');

    title(layout, ['Flat-membrane hop/trap diffusion with focused ISM, ' ...
        'BIC, and sticky BNP'], 'FontWeight', 'bold');
end
