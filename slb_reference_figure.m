function name = slb_reference_figure(out, outputDir)
%SLB_REFERENCE_FIGURE How the SLB anchor was arrived at, in one page.
%
% name = slb_reference_figure(out, outputDir)
%
% The anchor is the number the whole cell analysis is pinned to, so every step
% that produced it should be checkable here without reopening the MAT:
% what the decay looks like, whether the model describes it, WHICH pixels were
% used, and what sets the prior width.
%
% Three choices worth stating, because each replaces something that misled:
%
%   - the pixel-selection panel is drawn as an INTENSITY IMAGE with the selected
%     pixels outlined, not as a colour-mapped 0/1 array. A binary mask under a
%     continuous colormap reads as though the two values were far apart on a
%     scale, and the excluded region - the cell - ends up as the most visually
%     dominant object in a panel that is about the region OUTSIDE it.
%   - the mono-exponential fit is DRAWN, not just named. It is the model being
%     rejected, and the rejection is the point of the panel; naming it in a
%     legend without a curve leaves labels that do not correspond to lines.
%   - pixels pinned at the lifetime grid floor are counted and reported rather
%     than drawn as if they were measurements. They are censored, not fitted.

    pooled = out.pooled;
    if isfield(out, 'dtNs') && ~isempty(out.dtNs); dt = out.dtNs; else; dt = 0.16; end
    model = pooled.model(:);
    resid = pooled.residual(:);
    decay = model + resid .* sqrt(max(model, 1));
    t = (0:numel(model) - 1)' * dt;

    h = figure('Color', 'w', 'Visible', 'off', 'Position', [40 40 1500 1000]);
    layout = tiledlayout(h, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    fs = 11;

    % ---- 1. the decay, and the model that was rejected -------------------
    ax = nexttile(layout);
    hold(ax, 'on');
    plot(ax, t, max(decay, 0.5), '.', 'Color', [0.35 0.35 0.35], ...
        'MarkerSize', 5);
    monoDrawn = false;
    if isfield(pooled, 'monoTauNs') && ~isempty(pooled.monoTauNs)
        try
            basis = biexp_slb_basis(out.pooledIrf, dt, out.periodNs, numel(model));
            pm = biexp_slb_pattern_batch(basis, pooled.monoTauNs);
            design = [ones(numel(model), 1), pm];
            bm = poisson_nnls_whitened(design, decay);
            plot(ax, t, max(design * bm, 0.5), '-', 'Color', [0.85 0.33 0.10], ...
                'LineWidth', 1.6);
            monoDrawn = true;
        catch
            monoDrawn = false;
        end
    end
    plot(ax, t, max(model, 0.5), '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.8);
    set(ax, 'YScale', 'log', 'FontSize', fs, 'Box', 'off', 'TickDir', 'out', ...
        'LineWidth', 0.9);
    xlim(ax, [0 t(end)]);
    ylim(ax, [max(1, 0.4 * min(decay(decay > 0))), 2 * max(decay)]);
    xlabel(ax, 'time (ns)'); ylabel(ax, 'photons per bin');
    entries = {sprintf('data  (%.3g photons)', pooled.photons)};
    if monoDrawn
        entries{end + 1} = sprintf('1 exp: %.4f ns,  \\chi^2_\\nu = %.0f', ...
            pooled.monoTauNs, pooled.monoReducedDeviance);
    end
    entries{end + 1} = sprintf('2 exp: %.4f + %.4f ns,  \\chi^2_\\nu = %.0f', ...
        pooled.tau1Ns, pooled.tau2Ns, pooled.reducedDeviance);
    legend(ax, entries, 'Location', 'northeast', 'FontSize', fs - 1, ...
        'Box', 'off');
    title(ax, sprintf('anchor = %.4f ns   (short component, %.1f%% of photons)', ...
        out.anchorNs, 100 * pooled.photonFraction(1)), 'FontSize', fs + 1);
    hold(ax, 'off');

    % ---- 2. residuals ----------------------------------------------------
    ax = nexttile(layout);
    plot(ax, t, resid, '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.1);
    hold(ax, 'on');
    yline(ax, 0, '-', 'Color', [0.2 0.2 0.2]);
    set(ax, 'FontSize', fs, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.9);
    xlim(ax, [0 t(end)]);
    xlabel(ax, 'time (ns)'); ylabel(ax, 'Poisson residual (\sigma)');
    title(ax, sprintf(['residuals of the 2-exp fit   (max %.0f\\sigma - the ' ...
        'model is not exact)'], max(abs(resid))), 'FontSize', fs + 1);
    hold(ax, 'off');

    % ---- 3. which pixels, on the intensity image ------------------------
    ax = nexttile(layout);
    if isfield(out, 'slbOnlyMask') && ~isempty(out.slbOnlyMask)
        maskUsed = logical(out.slbOnlyMask);
        if isfield(out, 'intensityMap') && ~isempty(out.intensityMap)
            base = log10(double(out.intensityMap) + 1);
        else
            base = double(maskUsed);
        end
        lo = quantileLocalBiexp(base(isfinite(base)), 0.02);
        hi = quantileLocalBiexp(base(isfinite(base)), 0.995);
        imagesc(ax, base, [lo hi]);
        colormap(ax, perceptualColormap('gray'));
        axis(ax, 'image');
        set(ax, 'XTick', [], 'YTick', []);
        hold(ax, 'on');
        % Outline rather than fill, so the underlying image stays readable.
        contour(ax, double(maskUsed), [0.5 0.5], 'LineColor', [0.98 0.75 0.10], ...
            'LineWidth', 1.1);
        hold(ax, 'off');
        title(ax, sprintf(['SLB-only pixels outlined:  %d of %d  (%.0f%% of ' ...
            'the geometric region)'], nnz(maskUsed), out.geometricPixelCount, ...
            100 * nnz(maskUsed) / max(out.geometricPixelCount, 1)), ...
            'FontSize', fs + 1);
    else
        axis(ax, 'off');
        title(ax, 'pixel mask not stored', 'FontSize', fs + 1);
    end

    % ---- 4. what sets the prior width ------------------------------------
    ax = nexttile(layout);
    hold(ax, 'on');
    floorNs = out.opts.tauRangeNs(1);
    labels = {};
    if isfield(out, 'perPixel') && isfield(out.perPixel, 'tauNs') && ...
            ~isempty(out.perPixel.tauNs)
        v = out.perPixel.tauNs(:);
        censored1 = mean(v <= floorNs * 1.001);
        v = v(v > floorNs * 1.001);          % drop the censored pile-up
        histogram(ax, v, 90, 'Normalization', 'pdf', ...
            'FaceColor', [0.62 0.62 0.62], 'EdgeColor', 'none');
        labels{end + 1} = sprintf('1x1  (shot noise; %.0f%% censored at the floor)', ...
            100 * censored1);
    end
    if isfield(out, 'binned') && isfield(out.binned, 'tauNs') && ...
            ~isempty(out.binned.tauNs)
        v = out.binned.tauNs(:);
        censored4 = mean(v <= floorNs * 1.001);
        v = v(v > floorNs * 1.001);
        histogram(ax, v, 90, 'Normalization', 'pdf', ...
            'FaceColor', [0.00 0.45 0.74], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        labels{end + 1} = sprintf('%dx%d  (real variation; %.0f%% censored)', ...
            out.opts.binForVariation, out.opts.binForVariation, 100 * censored4);
    end
    xline(ax, out.anchorNs, '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 1.8);
    xline(ax, out.anchorNs - out.suggestedPrior.slbSigmaNs, ':', ...
        'Color', [0.85 0.10 0.10], 'LineWidth', 1.2);
    xline(ax, out.anchorNs + out.suggestedPrior.slbSigmaNs, ':', ...
        'Color', [0.85 0.10 0.10], 'LineWidth', 1.2);
    labels = [labels, {'anchor', '\pm\sigma'}];
    set(ax, 'FontSize', fs, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.9);
    xlabel(ax, 'per-pixel lifetime (ns)'); ylabel(ax, 'density');
    legend(ax, labels, 'Location', 'northeast', 'FontSize', fs - 1, 'Box', 'off');
    title(ax, sprintf('prior width \\sigma = %.4f ns', ...
        out.suggestedPrior.slbSigmaNs), 'FontSize', fs + 1);
    hold(ax, 'off');

    [~, acq] = fileparts(fileparts(fileparts(out.analysisMat)));
    title(layout, sprintf('bare-SLB reference   -   %s', acq), ...
        'FontSize', fs + 3, 'FontWeight', 'bold');

    name = fullfile(outputDir, 'slb_reference.png');
    drawnow;
    exportgraphics(h, name, 'Resolution', 300);
    close(h);
end
