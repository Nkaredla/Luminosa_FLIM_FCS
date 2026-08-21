function out = estimate_slb_lifetime_global(source, opts)
%ESTIMATE_SLB_LIFETIME_GLOBAL Shortest lifetime from all photons in one image.
%
% out = estimate_slb_lifetime_global(analysisMatFile)
% out = estimate_slb_lifetime_global(tcspcCube, opts)
%
% Sums every photon in one acquisition into a single decay, fits it with a
% range of component counts, and reports the shortest lifetime as a first
% approximation of the SLB lifetime - to be used for segmentation and as the
% fixed reference for downstream analysis.
%
% WHY DO IT THIS WAY
%
% The pipeline currently fits the SLB lifetime inside a segmented
% outside-cell region. But segmentation itself is driven by lifetime
% contrast, so segmentation depends on the lifetime and the lifetime depends
% on the segmentation. Estimating it from all photons removes segmentation
% from that loop, and uses the full photon budget (order 1e8) instead of a
% sub-region.
%
% It also matters empirically: refitted per acquisition from a segmented
% region, this lifetime came out between 0.108 and 0.399 ns across six
% nominally identical acquisitions, twice below the 0.16 ns bin width and so
% not determinable from the data at all. Fitted globally on one acquisition,
% two disjoint regions agreed to 4% (0.334 and 0.348 ns).
%
% THE CHECK THAT MAKES IT TRUSTWORTHY
%
% "The shortest lifetime" is not well defined on its own - it depends on how
% many components you fit. So the fit is repeated over a range of component
% counts and the shortest lifetime is reported for each. If it is stable
% across counts, it is a real component and can be trusted as a reference.
% If it drifts with the number of components, it is absorbing whatever the
% model cannot otherwise describe and should NOT be used as a fixed
% reference. The stability spread is reported explicitly, and so is the
% amplitude fraction it carries: a shortest component holding a negligible
% fraction is a fitting artefact, not the SLB.
%
% opts fields (all optional)
%   componentCounts   component counts to try, default 1:4
%   tauRangeNs        seed span for the log-spaced initial lifetimes,
%                     default [0.2 5]
%   stabilityTolNs    spread across counts (>=2) below which the shortest
%                     lifetime is called stable, default 0.05
%   minAmplitudeFrac  amplitude fraction the shortest component must carry to
%                     be credible, default 0.05
%   dtNs, pulsePeriodNs, irf
%                     required when source is a raw cube; read from the
%                     result struct when source is an analysis MAT
%   outputDir         default alongside the MAT, or pwd for a raw cube
%   makeFigure        default true

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    defaults = struct( ...
        'componentCounts', 1:4, ...
        'tauRangeNs', [0.2 5], ...
        'stabilityTolNs', 0.05, ...
        'minAmplitudeFrac', 0.05, ...
        'dtNs', [], 'pulsePeriodNs', [], 'irf', [], ...
        'outputDir', '', 'makeFigure', true);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if exist('Fluofit', 'file') ~= 2
        error('estimate_slb_lifetime_global:MissingDependency', ...
            'Fluofit.m must be on the MATLAB path.');
    end

    tauSlbPipeline = NaN;
    if ischar(source) || isstring(source)
        matPath = char(source);
        if ~isfile(matPath)
            error('estimate_slb_lifetime_global:MissingFile', ...
                'File not found: %s', matPath);
        end
        info = whos('-file', matPath);
        variables = {info.name};
        if ~ismember('tcspc_pix', variables)
            error('estimate_slb_lifetime_global:NoCube', ...
                ['%s has no tcspc_pix. Re-run the pipeline with ' ...
                 'cfg.saveTcspcPix = true.'], matPath);
        end
        store = matfile(matPath);
        cube = store.tcspc_pix;
        result = store.result;
        opts.irf = double(result.irf.curve(:));
        if isfield(result.bayesian, 'compact')
            if isfield(result.bayesian.compact, 'dtNs')
                opts.dtNs = double(result.bayesian.compact.dtNs);
            end
            if isfield(result.bayesian.compact, 'pulsePeriodNs')
                opts.pulsePeriodNs = ...
                    double(result.bayesian.compact.pulsePeriodNs);
            end
        end
        if isempty(opts.dtNs) || ~isfinite(opts.dtNs)
            opts.dtNs = double(result.config.tcspcBinNs);
        end
        tauSlbPipeline = double(result.slbReference.fixedLifetimeNs);
        if isempty(opts.outputDir); opts.outputDir = fileparts(matPath); end
        [~, label] = fileparts(matPath);
    else
        cube = source;
        if isempty(opts.dtNs) || isempty(opts.pulsePeriodNs) || isempty(opts.irf)
            error('estimate_slb_lifetime_global:MissingInstrument', ...
                ['A raw cube needs opts.dtNs, opts.pulsePeriodNs and ' ...
                 'opts.irf.']);
        end
        if isempty(opts.outputDir); opts.outputDir = pwd; end
        label = 'cube';
    end
    if ~isfolder(opts.outputDir); mkdir(opts.outputDir); end

    % Every photon in the acquisition, with no mask and no segmentation.
    decay = double(squeeze(sum(sum(double(cube), 1), 2)));
    decay = decay(:);
    photonTotal = sum(decay);
    nBins = numel(decay);
    timeNs = (0:nBins-1)' * opts.dtNs;

    fprintf('\nestimate_slb_lifetime_global\n');
    fprintf('  %s\n', label);
    fprintf(['  whole-image decay: %d bins, dt %.4g ns, period %.4g ns, ' ...
        '%.4g photons\n'], nBins, opts.dtNs, opts.pulsePeriodNs, photonTotal);
    if isfinite(tauSlbPipeline)
        fprintf('  pipeline fitted tauSLB (segmented region): %.4g ns\n', ...
            tauSlbPipeline);
    end
    fprintf('\n  nExp   shortest tau   its amplitude frac   all tau (ns)          chi2\n');
    fprintf(  '  ----   ------------   ------------------   --------------------  ------\n');

    counts = opts.componentCounts(:).';
    rows = struct([]);
    fits = cell(1, numel(counts));
    for k = 1:numel(counts)
        nExp = counts(k);
        if nExp == 1
            seed = sqrt(prod(opts.tauRangeNs));
        else
            seed = logspace(log10(opts.tauRangeNs(1)), ...
                log10(opts.tauRangeNs(2)), nExp);
        end
        fit = struct('nExp', nExp, 'ok', false);
        try
            [~, offset, A, tau, ~, dtau, ~, ~, ~, chi] = ...
                Fluofit(opts.irf, decay, opts.pulsePeriodNs, opts.dtNs, ...
                seed, [], 0);
            tau = double(tau(:).');
            A = double(A(:).');
            % A can carry a leading background/offset term depending on the
            % Fluofit variant; keep only as many amplitudes as lifetimes.
            if numel(A) > numel(tau)
                A = A(end - numel(tau) + 1:end);
            end
            [tauSorted, order] = sort(tau, 'ascend');
            ampSorted = nan(size(tauSorted));
            if numel(A) == numel(tau); ampSorted = A(order); end
            ampFraction = ampSorted / max(sum(ampSorted(ampSorted > 0)), eps);
            fit = struct('nExp', nExp, 'ok', true, 'tau', tauSorted, ...
                'amplitude', ampSorted, 'amplitudeFraction', ampFraction, ...
                'chi', double(chi(1)), 'offset', double(offset(1)), ...
                'dtau', double(dtau(:).'));
            fprintf('  %4d   %12.4f   %18.3f   %-20s  %6.4g\n', nExp, ...
                tauSorted(1), ampFraction(1), ...
                mat2str(round(tauSorted, 3)), fit.chi);
        catch fitError
            fprintf('  %4d   fit failed: %s\n', nExp, fitError.message);
            fit.message = fitError.message;
        end
        fits{k} = fit;
        entry = struct('nExp', nExp, 'ok', fit.ok, ...
            'shortestTauNs', NaN, 'shortestAmplitudeFraction', NaN, ...
            'chi', NaN);
        if fit.ok
            entry.shortestTauNs = fit.tau(1);
            entry.shortestAmplitudeFraction = fit.amplitudeFraction(1);
            entry.chi = fit.chi;
        end
        if isempty(rows); rows = entry; else; rows(end+1) = entry; end %#ok<AGROW>
    end

    summary = struct2table(rows);
    out = struct('summary', summary, 'fits', {fits}, 'decay', decay, ...
        'timeNs', timeNs, 'opts', opts, 'photonTotal', photonTotal, ...
        'tauSlbPipeline', tauSlbPipeline);

    % Non-converged fits must be excluded before judging stability, or one
    % broken fit destroys the verdict. Two objective symptoms, both seen on
    % real data: chi2 worse than a fit with FEWER components (nested models
    % cannot fit worse when converged), and a component carrying essentially
    % zero amplitude (the optimiser split one species and starved a half).
    converged = summary.ok;
    bestChiSoFar = Inf;
    for k = 1:height(summary)
        if ~summary.ok(k); continue; end
        degenerate = any(fits{k}.amplitudeFraction < 1e-3);
        worseThanSimpler = summary.chi(k) > bestChiSoFar;
        if degenerate || worseThanSimpler
            converged(k) = false;
            fprintf('  nExp %d rejected (chi2 %.4g): worseThanSimpler=%d, zeroAmplitudeComponent=%d\n', ...
                summary.nExp(k), summary.chi(k), worseThanSimpler, degenerate);
        else
            bestChiSoFar = min(bestChiSoFar, summary.chi(k));
        end
    end
    out.converged = converged;

    % Counts of 1 are excluded regardless: a single exponential averages
    % everything and is not expected to agree with the rest.
    multi = summary.nExp >= 2 & converged;
    shortest = summary.shortestTauNs(multi);
    fprintf('\n');
    if numel(shortest) < 2
        fprintf('  not enough successful multi-component fits to judge stability\n');
        out.stable = false;
        out.recommendedTauSlbNs = NaN;
    else
        spread = max(shortest) - min(shortest);
        out.shortestSpreadNs = spread;
        out.stable = spread <= opts.stabilityTolNs;
        fprintf(['  shortest lifetime across nExp >= 2: %s ns, spread ' ...
            '%.4f ns (tol %.3g)\n'], mat2str(round(shortest(:).', 4)), ...
            spread, opts.stabilityTolNs);
        credible = summary.shortestAmplitudeFraction(multi) >= ...
            opts.minAmplitudeFrac;
        if out.stable && all(credible)
            out.recommendedTauSlbNs = median(shortest);
            fprintf(['  STABLE. Recommended tauSLB = %.4f ns, usable as the ' ...
                'segmentation seed and\n  as the fixed reference.\n'], ...
                out.recommendedTauSlbNs);
        else
            out.recommendedTauSlbNs = NaN;
            if ~out.stable
                fprintf(['  NOT STABLE. The shortest lifetime moves with ' ...
                    'the number of components, so it is\n  absorbing ' ...
                    'unmodelled structure rather than measuring one ' ...
                    'species. Do not fix it.\n']);
            end
            if ~all(credible)
                fprintf(['  At least one fit puts a negligible amplitude ' ...
                    'fraction (< %.3g) on the shortest\n  component, which ' ...
                    'means it is a fitting artefact rather than the SLB.\n'], ...
                    opts.minAmplitudeFrac);
            end
        end
        if isfinite(tauSlbPipeline)
            fprintf('  pipeline value was %.4f ns (difference %+.4f ns)\n', ...
                tauSlbPipeline, tauSlbPipeline - median(shortest));
        end
    end

    if isfinite(out.recommendedTauSlbNs)
        fprintf(['\n  Use as:  cfg.slbTauSeedNs = %.4f;  and floor the fit ' ...
            'with\n           cfg.lifetimeBoundsNs = [%.3f 10];\n'], ...
            out.recommendedTauSlbNs, 2 * opts.dtNs);
    end

    if opts.makeFigure
        out.figure = plotGlobalDecayFits(out, label);
        fprintf('  wrote %s\n', out.figure);
    end
end

function name = plotGlobalDecayFits(out, label)
    opts = out.opts;
    decay = out.decay;
    timeNs = out.timeNs;
    fits = out.fits;

    h = figure('Color', 'w', 'Visible', 'off', 'Position', [60 60 1100 620]);
    layout = tiledlayout(h, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Decay with each fit overlaid, spanning both columns of the first two rows.
    ax = nexttile(layout, 1, [2 2]);
    hold(ax, 'on');
    stairs(ax, timeNs, max(decay, 0.5), 'Color', [0.45 0.45 0.45], ...
        'LineWidth', 0.7, 'DisplayName', 'all photons');
    colours = lines(numel(fits));
    for k = 1:numel(fits)
        fit = fits{k};
        if ~fit.ok; continue; end
        model = rebuildModel(fit, opts.irf, timeNs, opts.pulsePeriodNs);
        if isempty(model); continue; end
        model = model * (sum(decay) / max(sum(model), eps));
        plot(ax, timeNs, max(model, 0.5), '-', 'LineWidth', 1.3, ...
            'Color', colours(k, :), ...
            'DisplayName', sprintf('%d exp, shortest %.3f ns', ...
                fit.nExp, fit.tau(1)));
    end
    hold(ax, 'off');
    set(ax, 'YScale', 'log');
    upper = max(4, 2 * max(decay));
    ylim(ax, [0.5 upper]);
    set(ax, 'YTick', 10 .^ (0:ceil(log10(upper))));
    xlim(ax, [0 timeNs(end)]);
    grid(ax, 'on'); ax.GridAlpha = 0.12;
    ylabel(ax, 'photons / bin');
    xlabel(ax, 'time [ns]');
    legend(ax, 'Location', 'northeast', 'FontSize', 7, 'Box', 'off');
    title(ax, sprintf('Whole-image decay, %.3g photons', out.photonTotal), ...
        'FontSize', 9);

    % Shortest lifetime versus component count: the stability check.
    ax = nexttile(layout, 5);
    ok = out.summary.ok;
    plot(ax, out.summary.nExp(ok), out.summary.shortestTauNs(ok), ...
        'o-', 'LineWidth', 1.4, 'Color', [0.00 0.45 0.74]);
    if isfinite(out.tauSlbPipeline)
        yline(ax, out.tauSlbPipeline, 'r--', 'pipeline');
    end
    grid(ax, 'on');
    xlabel(ax, 'number of components');
    ylabel(ax, 'shortest tau [ns]');
    title(ax, 'Stability of the shortest lifetime', 'FontSize', 9);

    % Amplitude fraction carried by the shortest component.
    ax = nexttile(layout, 6);
    plot(ax, out.summary.nExp(ok), ...
        out.summary.shortestAmplitudeFraction(ok), ...
        's-', 'LineWidth', 1.4, 'Color', [0.85 0.33 0.10]);
    yline(ax, opts.minAmplitudeFrac, 'k--', 'credibility floor');
    grid(ax, 'on');
    ylim(ax, [0 1]);
    xlabel(ax, 'number of components');
    ylabel(ax, 'amplitude fraction');
    title(ax, 'Weight on the shortest component', 'FontSize', 9);

    title(layout, 'SLB lifetime from all photons in one acquisition', ...
        'FontWeight', 'bold');
    subtitle(layout, ['no segmentation used; a shortest lifetime that ' ...
        'moves with the component count is not a species'], 'FontSize', 8);

    name = fullfile(opts.outputDir, sprintf('%s_global_slb_lifetime.png', ...
        label));
    exportgraphics(h, name, 'Resolution', 200);
    exportgraphics(h, strrep(name, '.png', '.pdf'), 'ContentType', 'vector');
    close(h);
end

function model = rebuildModel(fit, irf, timeNs, periodNs)
    model = zeros(numel(timeNs), 1);
    irf = double(irf(:));
    for k = 1:numel(fit.tau)
        weight = 1;
        if numel(fit.amplitude) == numel(fit.tau) && ...
                isfinite(fit.amplitude(k)) && fit.amplitude(k) > 0
            weight = fit.amplitude(k);
        end
        model = model + weight * ...
            periodicDecay(irf, timeNs, periodNs, fit.tau(k));
    end
    if all(model == 0); model = []; end
end

function pattern = periodicDecay(irf, timeNs, periodNs, tauNs)
    nBins = numel(timeNs);
    decay = zeros(nBins, 1);
    for repeat = 0:3
        decay = decay + exp(-(timeNs + repeat * periodNs) / tauNs);
    end
    full = conv(irf, decay);
    pattern = zeros(nBins, 1);
    for start = 1:nBins:numel(full)
        stop = min(start + nBins - 1, numel(full));
        span = stop - start + 1;
        pattern(1:span) = pattern(1:span) + full(start:stop);
    end
    total = sum(pattern);
    if total > 0; pattern = pattern / total; end
end
