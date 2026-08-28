function name = pie_full_period_figure(out, outputDir)
%PIE_FULL_PERIOD_FIGURE Full-period decay, fit, residuals and lifetime spectrum.
%
% Three panels, because three separate questions are being asked: does the
% two-pulse model describe the data, where is the misfit, and what lifetimes are
% present in each colour.

    dt = out.basis.dtNs;
    t = (0:out.basis.nBin - 1)' * dt;
    h = figure('Color', 'w', 'Visible', 'off', 'Position', [40 40 1250 1050]);
    layout = tiledlayout(h, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax = nexttile(layout);
    semilogy(ax, t, max(out.decay, 0.5), 'k.', 'MarkerSize', 4); hold(ax, 'on');
    semilogy(ax, t, max(out.model, 0.5), '-', 'Color', [0.85 0.33 0.10], ...
        'LineWidth', 1.4);
    if ~isempty(out.discrete)
        semilogy(ax, t, max(out.discrete.model, 0.5), '-', ...
            'Color', [0 0.45 0.74], 'LineWidth', 1.2);
    end
    for b = out.pulseBins
        xline(ax, (b - 1) * dt, ':', 'Color', [0.4 0.4 0.4]);
    end
    yline(ax, out.background, '--', 'Color', [0.5 0.5 0.5]);
    grid(ax, 'on'); ylabel(ax, 'photons per bin');
    legendEntries = {sprintf('data, %.3g photons', out.photons), ...
        sprintf('distribution, \chi^2_\nu = %.1f', out.reducedDeviance)};
    if ~isempty(out.discrete)
        legendEntries{end + 1} = sprintf('discrete, \chi^2_\nu = %.1f', ...
            out.discrete.reducedDeviance);
    end
    legend(ax, legendEntries, 'Location', 'northeast');
    title(ax, 'full PIE period, both excitation pulses fitted together');
    hold(ax, 'off');

    ax = nexttile(layout);
    plot(ax, t, out.residual, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.9);
    hold(ax, 'on');
    if ~isempty(out.discrete)
        plot(ax, t, out.discrete.residual, '-', 'Color', [0 0.45 0.74], ...
            'LineWidth', 0.9);
    end
    yline(ax, 0, 'k-'); yline(ax, 3, ':'); yline(ax, -3, ':');
    grid(ax, 'on'); ylabel(ax, 'Poisson residual (\sigma)');
    xlabel(ax, 'time in period (ns)');
    hold(ax, 'off');

    ax = nexttile(layout);
    colours = {[0.85 0.33 0.10], [0 0.45 0.74]};
    hold(ax, 'on');
    for p = 1:size(out.weightsByPulse, 1)
        stem(ax, out.tauNs, out.weightsByPulse(p, :), 'filled', ...
            'Color', colours{min(p, 2)}, 'MarkerSize', 3);
    end
    set(ax, 'XScale', 'log', 'YScale', 'log');
    xline(ax, dt, '--', 'Color', [0.5 0.5 0.5]);
    grid(ax, 'on');
    xlabel(ax, 'lifetime (ns)'); ylabel(ax, 'photons');
    legend(ax, [out.pulseNames, {'one bin width'}], 'Location', 'northwest');
    title(ax, 'lifetime spectrum by excitation pulse (NNLS concentrates weight - read the modes)');
    hold(ax, 'off');

    title(layout, sprintf(['PIE full-period analysis  -  %s  -  background ' ...
        '%.2f/bin (%.2f%% of photons)'], strrep(out.ptuFile, '\', '/'), ...
        out.background, 100 * out.backgroundPhotons / max(out.photons, 1)), ...
        'Interpreter', 'none');

    name = fullfile(outputDir, 'pie_full_period.png');
    drawnow;
    exportgraphics(h, name, 'Resolution', 150);
    close(h);
end
