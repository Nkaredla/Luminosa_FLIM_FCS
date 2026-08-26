function name = immune_cell_MIET_biexp_component_figure(out, outputDir)
%IMMUNE_CELL_MIET_BIEXP_COMPONENT_FIGURE FLIM images of both lifetime components.
%
% name = immune_cell_MIET_biexp_component_figure(out, outputDir)
%
% Writes biexp_slb_component_lifetimes.png: the SLB (tau1) and long (tau2)
% lifetime maps side by side, with the photon share and the photon-weighted mean
% lifetime beneath them. The overview figure deliberately shows the fit-quality
% diagnostics; this one is the physical picture, and it exists because tau1 was
% never plotted anywhere despite being saved in every result.
%
% ON THE COLOUR SCALES
%
% tau1 is scaled to the tau1 GRID, not to a quantile of the data, because tau1
% can only take grid values and the grid spans just +/-2 sigma of the prior. A
% quantile scale would stretch that tiny range across the whole colormap and
% manufacture the appearance of structure where the fit has at most five
% distinct values to choose from. The tick labels name the actual nodes, so how
% coarse the quantisation is stays visible rather than hidden.
%
% tau2 is scaled to a robust quantile range, but pixels where tau2 was NOT
% refined off its grid node are hatched out in the companion panel so that
% quantisation is visible there too.

    maps = out.maps;
    h = figure('Color', 'w', 'Visible', 'off', 'Position', [40 40 1400 1050]);
    layout = tiledlayout(h, 2, 2, 'Padding', 'compact', ...
        'TileSpacing', 'compact');

    grid1 = sort(out.tau1Grid(:));
    finite2 = maps.tau2Ns(isfinite(maps.tau2Ns));

    % ---- tau1: the SLB component ---------------------------------------
    ax = nexttile(layout);
    drawMap(ax, maps.tau1Ns, 'parula');
    if numel(grid1) > 1 && grid1(end) > grid1(1)
        caxis(ax, [grid1(1) grid1(end)]);
        c = colorbar(ax);
        % Row vectors, and set Limits before Ticks: setting column-shaped ticks
        % on a colorbar whose limits have not settled is what raised
        % ColorBar/MarkDirty out of exportgraphics.
        c.Limits = [grid1(1) grid1(end)];
        c.Ticks = reshape(grid1, 1, []);
        c.TickLabels = reshape(arrayfun(@(v) sprintf('%.3f', v), ...
            grid1, 'UniformOutput', false), 1, []);
        c.Label.String = 'ns  (grid nodes)';
    else
        colorbar(ax);
    end
    atEdge = 0;
    if isfield(maps, 'tau1AtGridEdge')
        valid = isfinite(maps.tau1Ns);
        atEdge = 100 * mean(double(maps.tau1AtGridEdge(valid)));
    end
    title(ax, sprintf(['\\tau_1  SLB component  (prior %.4f \\pm %.3f ns, ' ...
        '%.1f%% at grid edge)'], out.opts.slbTauNs, out.opts.slbSigmaNs, ...
        atEdge));

    % ---- tau2: the long component --------------------------------------
    ax = nexttile(layout);
    drawMap(ax, maps.tau2Ns, 'parula');
    if numel(finite2) > 10
        lo = quantileLocalBiexp(finite2, 0.02);
        hi = quantileLocalBiexp(finite2, 0.98);
        if hi > lo; caxis(ax, [lo hi]); end
    end
    c = colorbar(ax); c.Label.String = 'ns';
    refinedPct = 100;
    if isfield(maps, 'tau2Refined')
        valid = isfinite(maps.tau2Ns);
        refinedPct = 100 * mean(double(maps.tau2Refined(valid)));
    end
    title(ax, sprintf(['\\tau_2  long component  (median %.3f ns, ' ...
        '%.0f%% refined off-grid)'], median(finite2), refinedPct));

    % ---- photon share of tau2 ------------------------------------------
    ax = nexttile(layout);
    drawMap(ax, maps.photonFraction2, 'hot');
    caxis(ax, [0 1]);
    c = colorbar(ax); c.Label.String = 'fraction';
    share = maps.photonFraction2(isfinite(maps.photonFraction2));
    title(ax, sprintf('photon share of \\tau_2  (median %.3f)', ...
        median(share)));

    % ---- photon-weighted mean lifetime ---------------------------------
    ax = nexttile(layout);
    drawMap(ax, maps.tauMeanNs, 'parula');
    mean2 = maps.tauMeanNs(isfinite(maps.tauMeanNs));
    if numel(mean2) > 10
        lo = quantileLocalBiexp(mean2, 0.02);
        hi = quantileLocalBiexp(mean2, 0.98);
        if hi > lo; caxis(ax, [lo hi]); end
    end
    c = colorbar(ax); c.Label.String = 'ns';
    title(ax, sprintf(['photon-weighted \\langle\\tau\\rangle  ' ...
        '(median %.3f ns)'], median(mean2)));

    binText = '1x1';
    if out.opts.binSize > 1
        binText = sprintf('%dx%d sliding, step 1', out.opts.binSize, ...
            out.opts.binSize);
    end
    title(layout, sprintf(['component lifetimes, soft-SLB biexponential  -  ' ...
        '%s  -  %d pixel(s), median %.0f photons'], binText, ...
        numel(out.pixelIndex), median(out.amplitude1 + out.amplitude2)));

    name = fullfile(outputDir, 'biexp_slb_component_lifetimes.png');
    exportgraphics(h, name, 'Resolution', 160);
    close(h);
end
