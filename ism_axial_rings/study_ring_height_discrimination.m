function out = study_ring_height_discrimination(opts)
%STUDY_RING_HEIGHT_DISCRIMINATION Is the long-lifetime pool ABOVE the membrane?
%
% out = study_ring_height_discrimination()
% out = study_ring_height_discrimination(opts)
%
% THE QUESTION THIS ANSWERS, WHICH THE SUMMED DECAY CANNOT
%
% In the immune-cell MIET data the third (longest) lifetime component is
% interpreted as internalised or endocytosed membrane dye - dye that has left
% the contact plane. But a third lifetime is equally consistent with a second
% chemically distinct population sitting AT the membrane. The detector-summed
% TCSPC decay cannot distinguish those two cases even in principle: both are a
% third exponential with some amplitude, and amplitude is confounded with
% abundance.
%
% The detector rings can, because a molecule at height z spreads its photons
% across rings with a weight pattern w(:, z) that changes with z while the
% decay shape does not.
%
% THE EXPERIMENT IS BUILT SO THE SUMMED DECAY HAS EXACTLY ZERO INFORMATION
%
% Two truths are simulated with the SAME three lifetimes:
%
%   displaced   third pool at heightUm = displacedHeightUm (default 0.40)
%   coplanar    third pool at the membrane height (default 0.02)
%
% and the third pool's molecule count is chosen SEPARATELY in each case so that
% the third pool contributes the same fraction of DETECTED photons in both.
% Because the detector-summed decay depends only on lifetimes and detected
% photon fractions, the two truths then have an identical expected summed
% decay - the study verifies this numerically and prints the residual. Any
% ability to tell them apart therefore comes from the ring dimension and from
% nothing else.
%
% This is a deliberately harder and more honest test than asking whether a
% third component exists at all. It also sidesteps a limitation established by
% diagnose_ring_height_identifiability: the normalised ring pattern saturates
% above roughly 0.3 um, so height is NOT identifiable as a continuous quantity
% at realistic photon counts. A binary displaced-versus-coplanar decision is
% still available, and that is what is measured here.
%
% THE TEST STATISTIC
%
% Two height-constrained fits of the same three-component model:
%
%   H0  all three components share ONE fitted height  (heightGroups = [1 1 1])
%   H1  each component has its own fitted height       (heightGroups = [1 2 3])
%
% dD = D(H0) - D(H1) is large when some component is axially displaced from the
% others. Neither fit needs to know the true heights, so the statistic is
% usable on real data. The threshold is calibrated on the coplanar truth at the
% (1-alpha) quantile, then the detection rate is measured on the displaced
% truth - identical false-positive rate by construction, so only power differs.
%
% For reference the same calibrate-then-measure procedure is applied to the
% summed decay using the only statistic available to it, the fitted amplitude
% fraction of the longest component. Its power should come out at alpha,
% because by construction it is looking at noise.
%
% opts fields (all optional)
%   tauNs                three lifetimes, default [0.34 1.90 3.30]
%   membraneHeightUm     default 0.02
%   displacedHeightUm    default 0.40
%   thirdPoolFractions   third pool's share of detected photons,
%                        default [0.05 0.10 0.20]
%   photonTotals         default [3e3 1e4 3e4]
%   repeats              default 60
%   alpha                default 0.05
%   restarts             fit restarts, default 3
%   dtNs, periodNs, irfFwhmNs, nBins, outputDir, makeFigure, seed

    thisDir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(opts); opts = struct(); end
    defaults = struct( ...
        'ringWeight', [], ...
        'tauNs', [0.34 1.90 3.30], ...
        'membraneHeightUm', 0.02, 'displacedHeightUm', 0.40, ...
        'slbMolecules', 1.00, 'membraneMolecules', 0.60, ...
        'thirdPoolFractions', [0.05 0.10 0.20], ...
        'photonTotals', [3e3 1e4 3e4], ...
        'repeats', 60, 'alpha', 0.05, 'restarts', 3, ...
        'dtNs', 50/156, 'periodNs', 50, 'irfFwhmNs', 0.35, 'nBins', 156, ...
        'outputDir', thisDir, 'makeFigure', true, 'seed', 23);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if ~isfolder(opts.outputDir); mkdir(opts.outputDir); end
    rng(opts.seed);

    if isempty(opts.ringWeight)
        fprintf('computing the ring weight table...\n');
        ring = simulate_ring_weight_vs_height(struct('makeFigure', false, ...
            'heightsUm', 0:0.02:1.0, 'nx', 41, 'fovUm', 2.6));
        weights = struct('table', ring.ringWeight ./ max(ring.ringWeight(:)), ...
            'heights', ring.heightsUm);
    else
        weights = struct('table', opts.ringWeight.ringWeight ./ ...
            max(opts.ringWeight.ringWeight(:)), ...
            'heights', opts.ringWeight.heightsUm);
    end
    nRing = size(weights.table, 1);

    timeNs = (0:opts.nBins-1)' * opts.dtNs;
    irf = exp(-0.5 * ((timeNs - 12 * opts.dtNs) / ...
        (opts.irfFwhmNs / (2 * sqrt(2 * log(2))))) .^ 2);
    irf = irf / sum(irf);

    wSlb = sum(ring_interpolate_weights(weights, 0));
    wMem = sum(ring_interpolate_weights(weights, opts.membraneHeightUm));
    wDisp = sum(ring_interpolate_weights(weights, opts.displacedHeightUm));
    baseSignal = opts.slbMolecules * wSlb + opts.membraneMolecules * wMem;

    fprintf('\nstudy_ring_height_discrimination\n');
    fprintf('  %d rings, %d bins, dt %.4g ns, period %.3g ns\n', ...
        nRing, opts.nBins, opts.dtNs, opts.periodNs);
    fprintf('  lifetimes %s ns; SLB at z=0, membrane at z=%.3f um\n', ...
        mat2str(opts.tauNs), opts.membraneHeightUm);
    fprintf(['  third pool: displaced z=%.3f um (collection %.4f) vs ' ...
        'coplanar z=%.3f um (collection %.4f)\n'], ...
        opts.displacedHeightUm, wDisp, opts.membraneHeightUm, wMem);
    fprintf('  collection ratio %.2fx, so equal photon share needs %.2fx the molecules\n', ...
        wMem / wDisp, wMem / wDisp);
    fprintf('  %d repeats, %d fit restarts, false-positive rate fixed at %.0f%%\n\n', ...
        opts.repeats, opts.restarts, 100 * opts.alpha);

    nFrac = numel(opts.thirdPoolFractions);
    nTotal = numel(opts.photonTotals);
    rows = struct([]);

    for iF = 1:nFrac
        fraction = opts.thirdPoolFractions(iF);
        % Detected photons from the third pool must equal fraction of the
        % total in BOTH truths, so its molecule count differs between them.
        thirdSignal = fraction / (1 - fraction) * baseSignal;
        displaced = makeTruth(opts, thirdSignal / wDisp, opts.displacedHeightUm);
        coplanar  = makeTruth(opts, thirdSignal / wMem,  opts.membraneHeightUm);

        % Verify the summed decays really are indistinguishable.
        eDisp = expectedCounts(displaced, weights, irf, timeNs, opts.periodNs);
        eCop  = expectedCounts(coplanar,  weights, irf, timeNs, opts.periodNs);
        sDisp = sum(eDisp, 1) / sum(eDisp(:));
        sCop  = sum(eCop, 1)  / sum(eCop(:));
        summedGap = max(abs(sDisp - sCop)) / max(sDisp);
        ringGap = norm(normalise(sum(eDisp, 2)) - normalise(sum(eCop, 2)));
        fprintf(['  third pool = %4.1f%% of photons | summed-decay ' ...
            'difference %.2e (relative) | ring-pattern difference %.4f\n'], ...
            100 * fraction, summedGap, ringGap);

        for iN = 1:nTotal
            photonTotal = opts.photonTotals(iN);
            tic;
            nullStat = nan(opts.repeats, 1);
            altStat = nan(opts.repeats, 1);
            nullSummed = nan(opts.repeats, 1);
            altSummed = nan(opts.repeats, 1);
            altHeightSpread = nan(opts.repeats, 1);
            for rep = 1:opts.repeats
                cNull = ringPoisson(eCop * (photonTotal / sum(eCop(:))));
                cAlt = ringPoisson(eDisp * (photonTotal / sum(eDisp(:))));
                [nullStat(rep), ~] = displacementStatistic(cNull, irf, ...
                    timeNs, opts, weights);
                [altStat(rep), altHeightSpread(rep)] = ...
                    displacementStatistic(cAlt, irf, timeNs, opts, weights);
                nullSummed(rep) = summedStatistic(cNull, irf, timeNs, opts);
                altSummed(rep) = summedStatistic(cAlt, irf, timeNs, opts);
            end
            ringThreshold = quantileLocal(nullStat, 1 - opts.alpha);
            summedThreshold = quantileLocal(nullSummed, 1 - opts.alpha);
            entry = struct('thirdPoolFraction', fraction, ...
                'photonTotal', photonTotal, ...
                'ringThreshold', ringThreshold, ...
                'ringPower', mean(altStat > ringThreshold), ...
                'summedThreshold', summedThreshold, ...
                'summedPower', mean(altSummed > summedThreshold), ...
                'medianNullStat', median(nullStat, 'omitnan'), ...
                'medianAltStat', median(altStat, 'omitnan'), ...
                'medianAltHeightSpreadUm', median(altHeightSpread, 'omitnan'), ...
                'summedDecayGap', summedGap, 'ringPatternGap', ringGap, ...
                'seconds', toc);
            if isempty(rows); rows = entry; else; rows(end+1) = entry; end %#ok<AGROW>
            fprintf(['      N=%6d  ring power %.2f (dD %.2f -> %.2f) | ' ...
                'summed power %.2f | fitted height spread %.3f um | %4.0f s\n'], ...
                photonTotal, entry.ringPower, entry.medianNullStat, ...
                entry.medianAltStat, entry.summedPower, ...
                entry.medianAltHeightSpreadUm, entry.seconds);
        end
    end

    summary = struct2table(rows);
    out = struct('summary', summary, 'opts', opts, 'weights', weights, ...
        'irf', irf, 'timeNs', timeNs);
    csvFile = fullfile(opts.outputDir, 'ring_height_discrimination.csv');
    writetable(summary, csvFile);
    fprintf('\n  wrote %s\n', csvFile);
    if opts.makeFigure
        out.figure = plotDiscrimination(out);
        fprintf('  wrote %s\n', out.figure);
    end

    fprintf(['\n  Reading this: the summed-decay column is a NULL CONTROL. ' ...
        'The two truths have\n  identical expected summed decays (see the ' ...
        'relative difference above), so the\n  summed power must land at ' ...
        'the %.0f%% false-positive rate. Any ring power above\n  that is ' ...
        'information the detector array provides and the summed decay ' ...
        'cannot.\n'], 100 * opts.alpha);
end

% ------------------------------------------------------------------ statistic

function [deltaD, heightSpread] = displacementStatistic(counts, irf, ...
        timeNs, opts, weights)
%DISPLACEMENTSTATISTIC Deviance gain from letting components differ in height.
% H0: one shared height for all components. H1: a height per component.
% Neither uses knowledge of the true heights, so this runs unchanged on data.
    nComponent = numel(opts.tauNs);
    common = struct('restarts', opts.restarts);
    tied = common; tied.heightGroups = ones(1, nComponent);
    free = common; free.heightGroups = 1:nComponent;
    fitTied = ring_fit_height_constrained(counts, irf, timeNs, ...
        opts.periodNs, nComponent, weights, tied);
    fitFree = ring_fit_height_constrained(counts, irf, timeNs, ...
        opts.periodNs, nComponent, weights, free);
    deltaD = fitTied.deviance - fitFree.deviance;
    heightSpread = max(fitFree.height) - min(fitFree.height);
end

function statistic = summedStatistic(counts, irf, timeNs, opts)
%SUMMEDSTATISTIC The best the summed decay can do: the long component's share.
% Included as a control. Because the two truths were constructed to have the
% same expected summed decay, this statistic is measuring noise, and its
% "power" should equal alpha. If it does not, the construction is broken.
    nComponent = numel(opts.tauNs);
    fit = ring_fit_free_amplitudes(sum(counts, 1), irf, timeNs, ...
        opts.periodNs, nComponent, struct('restarts', opts.restarts));
    total = sum(fit.amplitude);
    if total <= 0; statistic = 0; return; end
    statistic = fit.amplitude(end) / total;
end

% ----------------------------------------------------------------- helpers

function truth = makeTruth(opts, thirdMolecules, thirdHeight)
    truth = struct( ...
        'tauNs', {opts.tauNs(1), opts.tauNs(2), opts.tauNs(3)}, ...
        'heightUm', {0, opts.membraneHeightUm, thirdHeight}, ...
        'molecules', {opts.slbMolecules, opts.membraneMolecules, ...
                      thirdMolecules});
end

function expected = expectedCounts(truth, weights, irf, timeNs, periodNs)
    nRing = size(weights.table, 1);
    expected = zeros(nRing, numel(timeNs));
    for j = 1:numel(truth)
        w = ring_interpolate_weights(weights, truth(j).heightUm);
        pattern = ring_periodic_decay(irf, timeNs, periodNs, truth(j).tauNs);
        expected = expected + truth(j).molecules * (w(:) * pattern(:)');
    end
end

function v = normalise(v)
    total = sum(v);
    if total > 0; v = v / total; end
end

function counts = ringPoisson(lambda)
    % Exact inversion for small means, normal approximation only above 100
    % where its relative error is negligible. The statistic is a deviance
    % computed under an exact-Poisson assumption, so the sampler should not be
    % the coarser of the two.
    counts = zeros(size(lambda));
    large = lambda > 100;
    if any(large(:))
        draw = lambda(large) + sqrt(lambda(large)) .* randn(nnz(large), 1);
        counts(large) = max(0, round(draw));
    end
    small = find(~large & lambda > 0);
    for k = 1:numel(small)
        limit = exp(-lambda(small(k)));
        product = rand();
        n = 0;
        while product > limit
            product = product * rand();
            n = n + 1;
        end
        counts(small(k)) = n;
    end
end

function q = quantileLocal(values, level)
    values = sort(values(~isnan(values)));
    if isempty(values); q = NaN; return; end
    if numel(values) == 1; q = values; return; end
    position = min(max(level * numel(values) + 0.5, 1), numel(values));
    low = floor(position); high = ceil(position);
    if low == high
        q = values(low);
    else
        q = values(low) + (position - low) * (values(high) - values(low));
    end
end

% ------------------------------------------------------------------ figure

function name = plotDiscrimination(out)
    summary = out.summary;
    fractions = unique(summary.thirdPoolFraction, 'stable');
    h = figure('Color', 'w', 'Visible', 'off', 'Position', [60 60 1100 430]);
    layout = tiledlayout(h, 1, 2, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    colours = lines(numel(fractions));

    ax = nexttile(layout);
    hold(ax, 'on');
    for k = 1:numel(fractions)
        keep = summary.thirdPoolFraction == fractions(k);
        plot(ax, summary.photonTotal(keep), summary.ringPower(keep), '-o', ...
            'Color', colours(k, :), 'LineWidth', 1.7, ...
            'MarkerFaceColor', colours(k, :), ...
            'DisplayName', sprintf('rings, %.0f%% of photons', ...
            100 * fractions(k)));
        plot(ax, summary.photonTotal(keep), summary.summedPower(keep), '--s', ...
            'Color', colours(k, :), 'LineWidth', 1.1, ...
            'DisplayName', sprintf('summed, %.0f%%', 100 * fractions(k)));
    end
    yline(ax, out.opts.alpha, ':k', 'chance');
    set(ax, 'XScale', 'log'); grid(ax, 'on'); ylim(ax, [0 1.05]);
    xlabel(ax, 'total detected photons');
    ylabel(ax, sprintf('detection rate at %.0f%% FP', 100 * out.opts.alpha));
    title(ax, sprintf('displaced (%.2f um) vs coplanar third pool', ...
        out.opts.displacedHeightUm));
    legend(ax, 'Location', 'eastoutside');

    ax = nexttile(layout);
    hold(ax, 'on');
    for k = 1:numel(fractions)
        keep = summary.thirdPoolFraction == fractions(k);
        plot(ax, summary.photonTotal(keep), ...
            summary.medianAltHeightSpreadUm(keep), '-o', ...
            'Color', colours(k, :), 'LineWidth', 1.7, ...
            'MarkerFaceColor', colours(k, :), ...
            'DisplayName', sprintf('%.0f%% of photons', 100 * fractions(k)));
    end
    yline(ax, out.opts.displacedHeightUm - out.opts.membraneHeightUm, '--k', ...
        'true separation');
    set(ax, 'XScale', 'log'); grid(ax, 'on');
    xlabel(ax, 'total detected photons');
    ylabel(ax, 'fitted height spread (\mum)');
    title(ax, 'spread of the three fitted heights');
    legend(ax, 'Location', 'best');

    name = fullfile(out.opts.outputDir, 'ring_height_discrimination.png');
    exportgraphics(h, name, 'Resolution', 200);
    close(h);
end
