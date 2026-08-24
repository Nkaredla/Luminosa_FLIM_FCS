function immune_cell_MIET_explorer_select(fig, rowCentre, colCentre)
%IMMUNE_CELL_MIET_EXPLORER_SELECT Pool a region's TCSPC, fit it, redraw.
%
% immune_cell_MIET_explorer_select(fig, rowCentre, colCentre)
%
% Called by the map's ButtonDownFcn and by the coordinate boxes. Pools the
% per-pixel TCSPC over a square region centred on the clicked pixel, fits two
% and three exponentials, and updates the decay and residual axes.
%
% The cube is read from the MAT file one region at a time with matfile rather
% than being loaded whole. The native cube is 602 x 602 x 156 uint16, about
% 113 MB, and this machine runs with C: chronically full - which surfaces as
% MATLAB "Out of memory" rather than as a disk error - so holding it resident
% alongside the maps is a real risk for no benefit when every selection touches
% at most a few hundred pixels.

    state = fig.UserData;
    if isempty(state) || ~isfield(state, 'cube'); return; end

    half = floor(state.regionSize.Value / 2);
    rows = max(1, rowCentre - half):min(state.imageSize(1), rowCentre + half);
    cols = max(1, colCentre - half):min(state.imageSize(2), colCentre + half);

    try
        block = state.cube.(state.cubeVariable)(rows, cols, :);
    catch readError
        state.status.Text = sprintf(['Could not read the TCSPC cube: %s. ' ...
            'The analysis MAT may have been written without ' ...
            'saveTcspcPix = true.'], readError.message);
        return;
    end
    decay = double(reshape(sum(sum(block, 1), 2), [], 1));
    photons = sum(decay);

    % ---- fits ----------------------------------------------------------
    fits = immune_cell_MIET_explorer_fit(decay, state.irf, state.dtNs, ...
        state.periodNs);

    % ---- decay panel ---------------------------------------------------
    ax = state.decayAxes;
    cla(ax);
    hold(ax, 'on');
    timeNs = (0:numel(decay) - 1)' * state.dtNs;
    semilogy(ax, timeNs, max(decay, 0.5), 'k.', 'MarkerSize', 7, ...
        'DisplayName', sprintf('data, %.0f photons', photons));
    colours = {[0.85 0.33 0.10], [0 0.45 0.74]};
    for k = 1:numel(fits)
        if ~fits(k).ok || isempty(fits(k).model); continue; end
        semilogy(ax, timeNs, max(fits(k).model, 0.5), '-', ...
            'Color', colours{k}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%d-exp, \\chi^2_\\nu = %.2f', ...
            fits(k).nExp, fits(k).chiSquare));
    end
    set(ax, 'YScale', 'log');
    grid(ax, 'on');
    xlabel(ax, 'time (ns)');
    ylabel(ax, 'photons per bin');
    title(ax, sprintf('rows %d-%d, cols %d-%d  (%d x %d pixels)', ...
        rows(1), rows(end), cols(1), cols(end), numel(rows), numel(cols)));
    legend(ax, 'Location', 'northeast');
    hold(ax, 'off');

    % ---- residual panel ------------------------------------------------
    ax = state.residualAxes;
    cla(ax);
    hold(ax, 'on');
    for k = 1:numel(fits)
        if ~fits(k).ok || isempty(fits(k).residual); continue; end
        plot(ax, timeNs, fits(k).residual, '-', 'Color', colours{k}, ...
            'LineWidth', 1, 'DisplayName', sprintf('%d-exp', fits(k).nExp));
    end
    yline(ax, 0, 'k-');
    yline(ax, 3, ':', 'Color', [0.6 0.6 0.6]);
    yline(ax, -3, ':', 'Color', [0.6 0.6 0.6]);
    grid(ax, 'on');
    xlabel(ax, 'time (ns)');
    ylabel(ax, 'Poisson residual');
    ylim(ax, [-8 8]);
    hold(ax, 'off');

    % ---- readout -------------------------------------------------------
    lines = {sprintf('region rows %d-%d, cols %d-%d, %.0f photons', ...
        rows(1), rows(end), cols(1), cols(end), photons)};
    for k = 1:numel(fits)
        if ~fits(k).ok
            lines{end + 1} = sprintf('%d-exp: FAILED (%s)', ...
                fits(k).nExp, fits(k).message);
            continue;
        end
        parts = arrayfun(@(t, a) sprintf('%.3f ns (%.0f%%)', t, 100 * a), ...
            fits(k).tauNs, fits(k).amplitude, 'UniformOutput', false);
        lines{end + 1} = sprintf(['%d-exp: %s   offset %.1f   shift %.3f ' ...
            'ns   chi2nu %.2f'], fits(k).nExp, strjoin(parts, ', '), ...
            fits(k).offset, fits(k).shiftNs, fits(k).chiSquare);
    end
    if numel(fits) == 2 && all([fits.ok])
        lines{end + 1} = sprintf(['max |residual|: 2-exp %.1f sigma, ' ...
            '3-exp %.1f sigma'], max(abs(fits(1).residual)), ...
            max(abs(fits(2).residual)));
        lines{end + 1} = ['A lower chi2 for 3-exp is not evidence for a ' ...
            'third component - look for STRUCTURE in the 2-exp residual.'];
    end
    state.readout.Value = lines;

    % ---- marker on the map ---------------------------------------------
    if isfield(state, 'marker') && isgraphics(state.marker)
        set(state.marker, 'XData', mean(cols), 'YData', mean(rows));
    end
    state.rowBox.Value = rowCentre;
    state.colBox.Value = colCentre;
    % uilabel exposes Text; only uitextarea/uieditfield have Value.
    state.status.Text = sprintf('selected row %d, col %d', ...
        rowCentre, colCentre);
    fig.UserData = state;
end
