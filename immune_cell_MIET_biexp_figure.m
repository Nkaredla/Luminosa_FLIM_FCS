function name = immune_cell_MIET_biexp_figure(out, outputDir)
%IMMUNE_CELL_MIET_BIEXP_FIGURE Six-panel overview of the soft-SLB biexponential.
%
% The panels are chosen so the two things that could invalidate the fit are
% visible without opening the MAT: the prior-pull map shows whether the soft
% constraint is fighting the data, and the reduced-deviance map shows where the
% two-component model fails. A figure of only lifetimes and fractions would look
% fine in both cases.

    maps = out.maps;
    h = figure('Color', 'w', 'Visible', 'off', 'Position', [40 40 1500 900]);
    layout = tiledlayout(h, 2, 3, 'Padding', 'compact', ...
        'TileSpacing', 'compact');

    panels = { ...
        maps.tau2Ns,             'free component \tau_2 (ns)',      'parula', []; ...
        maps.photonFraction2,    'photon share of \tau_2',           'parula', [0 1]; ...
        maps.tauMeanNs,          'photon-weighted mean \tau (ns)',   'parula', []; ...
        maps.slbPullSigma,       'SLB prior pull (\sigma)',          'coolwarm', [-2 2]; ...
        maps.reducedDeviance,    'reduced deviance',                 'hot',    []; ...
        log10(maps.intensity + 1), 'log_{10} intensity',             'gray',   []};

    for k = 1:size(panels, 1)
        ax = nexttile(layout);
        data = panels{k, 1};
        handle = imagesc(ax, data);
        set(handle, 'AlphaData', isfinite(data));
        axis(ax, 'image');
        title(ax, panels{k, 2});
        if strcmp(panels{k, 3}, 'coolwarm')
            colormap(ax, buildDivergingMap());
        else
            colormap(ax, panels{k, 3});
        end
        % ORDER MATTERS. caxis FIRST, colorbar second. Creating the colorbar
        % and then changing the axes colour limits leaves the colorbar's
        % limit-update listener pending; it then fires inside exportgraphics,
        % mid-traversal, and throws
        %   "Attempt to modify the tree during an update traversal"
        % out of ColorBar/attachAxesListeners. Setting the limits before the
        % colorbar exists means there is nothing left to notify.
        limits = panels{k, 4};
        finite = data(isfinite(data));
        if ~isempty(limits)
            caxis(ax, limits);
        elseif numel(finite) > 10
            lo = quantileLocalBiexp(finite, 0.02);
            hi = quantileLocalBiexp(finite, 0.98);
            if hi > lo; caxis(ax, [lo hi]); end
        end
        colorbar(ax);
        set(ax, 'XTick', [], 'YTick', []);
    end

    title(layout, sprintf(['biexponential with soft-fixed SLB: prior ' ...
        '%.4f +/- %.3f ns, bin %dx%d, %d pixels'], out.opts.slbTauNs, ...
        out.opts.slbSigmaNs, out.opts.binSize, out.opts.binSize, ...
        numel(out.pixelIndex)));

    name = fullfile(outputDir, 'biexp_slb_maps.png');
    % Flush every pending graphics update before exporting, so exportgraphics
    % never starts a traversal with work still queued.
    drawnow;
    exportgraphics(h, name, 'Resolution', 160);
    close(h);
end
