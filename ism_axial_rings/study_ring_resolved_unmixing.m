function out = study_ring_resolved_unmixing(opts)
%STUDY_RING_RESOLVED_UNMIXING Ring-resolved global fitting for lifetime splitting.
%
% out = study_ring_resolved_unmixing()
% out = study_ring_resolved_unmixing(opts)
%
% Does fitting the ring-resolved TCSPC jointly make the detection of a third
% lifetime pool (SLB / membrane / internalised) more reliable than fitting the
% detector-summed decay at the same total photon count?
%
% WHY THE RING DIMENSION COULD HELP
%
% Each detector ring collects a molecule at height z with a different weight
% w(ring, z), computed by simulate_ring_weight_vs_height. A population at a
% given height therefore imprints a specific, predictable amplitude pattern
% across rings, so height becomes an observable that is independent of
% lifetime. The M-open pathology measured earlier in this project - the false
% three-component rate RISING with photon count, 9% to 54% - happens because
% the likelihood alone cannot distinguish a real third exponential from two
% on-grid exponentials approximating one off-grid one. A physical constraint
% that the artefact cannot satisfy is the way out.
%
% THREE MODELS ARE COMPARED, AND THE DIFFERENCE BETWEEN THEM IS THE POINT
%
%   summed        rings summed first, then a free multiexponential fit.
%                 This is what the current pipeline does.
%                 params: nComp lifetimes + nComp amplitudes
%
%   ring-free     shared lifetimes across rings, amplitudes free and
%                 non-negative per ring. Uses the ring dimension as DATA but
%                 imposes no physics on it.
%                 params: nComp lifetimes + nComp*nRing amplitudes
%
%   ring-w(z)     shared lifetimes, and each component's ring pattern forced
%                 to a_j * w(:, z_j) with z_j fitted. The ring pattern is now
%                 a one-parameter physical family instead of nRing free
%                 numbers.
%                 params: nComp lifetimes + nComp heights + nComp amplitudes
%
% The middle model is included deliberately, because it isolates the two
% effects. Going from summed to ring-free adds information but also adds
% parameters; going from ring-free to ring-w(z) removes parameters by imposing
% the optics. If only the constrained model helps, the gain comes from the
% physics and not merely from resolving the detector.
%
% Splitting photons across rings does not create photons, it partitions them,
% so every ring decay is noisier than the summed one. All three models are
% therefore given the SAME total photons.
%
% HOW ROBUSTNESS IS MEASURED - CALIBRATED, NOT ASSUMED
%
% For every model and photon count, both a 2-component and a 3-component fit
% are run and the Poisson deviance drop dD = D2 - D3 is recorded. dD is then
% thresholded to decide "three pools present".
%
% The threshold is CALIBRATED on the two-pool truth: it is set to the
% (1-alpha) quantile of the null dD distribution, so every model is held to
% the SAME false-positive rate by construction. What is then compared is the
% detection rate on the three-pool truth at that fixed false-positive rate.
% This is the only comparison that is fair across models with different
% parameter counts, and it sidesteps the question of what an asymptotic
% chi-square threshold would be for a boundary-constrained problem - a
% question with no clean answer here.
%
% An earlier version of this study scored components by how well their ring
% pattern matched w(:, z), with a tolerance derived from the spacing of the
% height grid. That was wrong: the grid step is 20 nm, adjacent patterns are
% nearly identical, so the tolerance collapsed towards zero and every noisy
% pattern was flagged. Any threshold on a noisy statistic has to be calibrated
% against the noise, never against the model grid.
%
% THE FITTED HEIGHT IS NOT A MEASUREMENT - READ THIS BEFORE USING z3
%
% diagnose_ring_height_identifiability establishes, without invoking the
% optimiser at all, that the fitted height is not identifiable at these photon
% counts:
%
%   - the NORMALISED ring pattern saturates above about 0.3 um. Its distance
%     from the pattern at 0.40 um is 0.039 at z = 0.30 but only 0.037 at
%     z = 1.00, so heights above 0.3 um are nearly indistinguishable from one
%     another in shape. Absolute collection efficiency still falls steeply with
%     height, but that is degenerate with the free amplitude.
%   - scanning the height on a one-dimensional deviance profile, with the
%     lifetimes and other heights pinned at truth and no optimiser involved,
%     gives an interquartile range of about 0.6 um at N = 1e4 and no real
%     improvement at 3e4.
%   - the noiseless profile does put its minimum exactly at the true 0.40 um,
%     so the forward model and the design matrix are correct. The problem is
%     information, not implementation.
%
% So the z3 column below is reported with its spread and with the fraction of
% fits sitting at the edge of the height table, and it must NOT be read as a
% height measurement. What the ring dimension does support is the strictly
% weaker, and still useful, binary question - is some component axially
% displaced from the others - measured in study_ring_height_discrimination.
%
% opts fields (all optional)
%   ringWeight       struct with fields ringWeight [nRing x nz] and heightsUm;
%                    if omitted, simulate_ring_weight_vs_height is called
%   truth            struct array with fields tauNs, heightUm, molecules
%                    Default: SLB 0.34 ns at 0, membrane 1.9 ns at 0.02 um,
%                    internalised 3.3 ns at 0.4 um
%   nullTruth        same, without the third pool
%   photonTotals     total detected photons per local region,
%                    default [1e3 3e3 1e4 3e4]
%   repeats          noise realisations per condition, default 60
%   alpha            calibrated false-positive rate, default 0.05
%   restarts         fit restarts per model, default 3. Single-start fits were
%                    measured to miss the better optimum in 8 of 15 cases for
%                    the constrained model, which would penalise it unfairly.
%   evalBudget       give every model the same ABSOLUTE optimiser budget
%                    instead of the same budget per parameter. Set it to check
%                    that a conclusion is not an artefact of search effort.
%   dtNs, periodNs, irfFwhmNs, nBins  instrument, default 50/156, 50, 0.35, 156
%   outputDir, makeFigure, seed

    thisDir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(opts); opts = struct(); end
    defaults = struct( ...
        'ringWeight', [], ...
        'truth', [], 'nullTruth', [], ...
        'photonTotals', [1e3 3e3 1e4 3e4], ...
        'repeats', 60, 'alpha', 0.05, 'restarts', 3, 'evalBudget', [], ...
        'dtNs', 50/156, 'periodNs', 50, 'irfFwhmNs', 0.35, 'nBins', 156, ...
        'outputDir', thisDir, 'makeFigure', true, 'seed', 17);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if isempty(opts.truth)
        opts.truth = struct( ...
            'tauNs',     {0.34, 1.90, 3.30}, ...
            'heightUm',  {0.00, 0.02, 0.40}, ...
            'molecules', {1.00, 0.60, 0.50});
    end
    if isempty(opts.nullTruth)
        opts.nullTruth = opts.truth(1:2);
    end
    if ~isfolder(opts.outputDir); mkdir(opts.outputDir); end
    rng(opts.seed);

    % ---- ring weight table -------------------------------------------------
    if isempty(opts.ringWeight)
        fprintf('computing the ring weight table...\n');
        ring = simulate_ring_weight_vs_height(struct('makeFigure', false, ...
            'heightsUm', 0:0.02:1.0, 'nx', 41, 'fovUm', 2.6));
        weightTable = ring.ringWeight;
        weightHeights = ring.heightsUm;
    else
        weightTable = opts.ringWeight.ringWeight;
        weightHeights = opts.ringWeight.heightsUm;
    end
    weightTable = weightTable ./ max(weightTable(:));
    nRing = size(weightTable, 1);

    timeNs = (0:opts.nBins-1)' * opts.dtNs;
    irf = exp(-0.5 * ((timeNs - 12 * opts.dtNs) / ...
        (opts.irfFwhmNs / (2 * sqrt(2 * log(2))))) .^ 2);
    irf = irf / sum(irf);

    weights = struct('table', weightTable, 'heights', weightHeights);
    models = {'summed', 'ring-free', 'ring-w(z)'};
    nModel = numel(models);

    fprintf('\nstudy_ring_resolved_unmixing\n');
    fprintf('  %d rings, %d bins, dt %.4g ns, period %.3g ns, IRF %.2f ns\n', ...
        nRing, opts.nBins, opts.dtNs, opts.periodNs, opts.irfFwhmNs);
    printTruth('  three-pool truth', opts.truth, weights);
    printTruth('  two-pool truth  ', opts.nullTruth, weights);
    fprintf('  %d repeats per condition, false-positive rate fixed at %.0f%%\n', ...
        opts.repeats, 100 * opts.alpha);
    fprintf(['  models: summed (rings collapsed) | ring-free (shared tau, ' ...
        'free per-ring amps)\n          | ring-w(z) (shared tau, ring ' ...
        'pattern forced to w(:,z_j))\n\n']);

    nTotal = numel(opts.photonTotals);
    % deltaD(model, photonIndex, repeat, scenario) with scenario 1 = null.
    deltaD = nan(nModel, nTotal, opts.repeats, 2);
    tauErr3 = nan(nModel, nTotal, opts.repeats, 2);
    heightHat = nan(nModel, nTotal, opts.repeats, 3);   % ring-w(z) only
    atBound = false(nModel, nTotal, opts.repeats);
    evals = nan(nModel, nTotal, opts.repeats, 2);
    fitOpts = struct('restarts', opts.restarts, 'evalBudget', opts.evalBudget);
    scenarios = {opts.nullTruth, opts.truth};
    scenarioName = {'two-pool (null)', 'three-pool'};

    for sIndex = 1:2
        truth = scenarios{sIndex};
        trueTau = [truth.tauNs];
        for iN = 1:nTotal
            photonTotal = opts.photonTotals(iN);
            tic;
            for rep = 1:opts.repeats
                counts = simulateRingCounts(truth, weights, irf, timeNs, ...
                    opts.periodNs, photonTotal);
                summed = sum(counts, 1);
                for m = 1:nModel
                    switch models{m}
                        case 'summed'
                            f2 = ring_fit_free_amplitudes(summed, irf, ...
                                timeNs, opts.periodNs, 2, fitOpts);
                            f3 = ring_fit_free_amplitudes(summed, irf, ...
                                timeNs, opts.periodNs, 3, fitOpts);
                        case 'ring-free'
                            f2 = ring_fit_free_amplitudes(counts, irf, ...
                                timeNs, opts.periodNs, 2, fitOpts);
                            f3 = ring_fit_free_amplitudes(counts, irf, ...
                                timeNs, opts.periodNs, 3, fitOpts);
                        otherwise
                            f2 = ring_fit_height_constrained(counts, irf, ...
                                timeNs, opts.periodNs, 2, weights, fitOpts);
                            f3 = ring_fit_height_constrained(counts, irf, ...
                                timeNs, opts.periodNs, 3, weights, fitOpts);
                            heightHat(m, iN, rep, :) = f3.height(:)';
                            atBound(m, iN, rep) = any(f3.atBound);
                    end
                    deltaD(m, iN, rep, sIndex) = f2.deviance - f3.deviance;
                    tauErr3(m, iN, rep, sIndex) = tauMatchError(f3.tau, trueTau);
                    evals(m, iN, rep, sIndex) = f2.funcCount + f3.funcCount;
                end
            end
            fprintf('  %-16s N=%6d  done in %5.1f s\n', ...
                scenarioName{sIndex}, photonTotal, toc);
        end
    end

    % ---- calibrate on the null, then measure power -------------------------
    rows = struct([]);
    fprintf(['\n  model         N   dD_thresh   power   median|dTau| ns' ...
        '     z3 (um, IQR)   evals\n']);
    fprintf(['  --------------------------------------------------------' ...
        '-------------------------\n']);
    for m = 1:nModel
        for iN = 1:nTotal
            nullD = squeeze(deltaD(m, iN, :, 1));
            altD  = squeeze(deltaD(m, iN, :, 2));
            threshold = quantileLocal(nullD, 1 - opts.alpha);
            power = mean(altD > threshold);
            entry = struct('model', models{m}, ...
                'photonTotal', opts.photonTotals(iN), ...
                'deltaDThreshold', threshold, ...
                'powerAtFixedFP', power, ...
                'achievedFP', mean(nullD > threshold), ...
                'medianNullDeltaD', median(nullD, 'omitnan'), ...
                'medianAltDeltaD', median(altD, 'omitnan'), ...
                'tauErrNs', median(squeeze(tauErr3(m, iN, :, 2)), 'omitnan'), ...
                'height3Um', NaN);
            entry.height3IqrUm = NaN;
            entry.atBoundFraction = NaN;
            if strcmp(models{m}, 'ring-w(z)')
                zAll = squeeze(heightHat(m, iN, :, :));
                highest = max(zAll, [], 2);
                entry.height3Um = median(highest, 'omitnan');
                entry.height3IqrUm = quantileLocal(highest, 0.75) - ...
                    quantileLocal(highest, 0.25);
                entry.atBoundFraction = mean(squeeze(atBound(m, iN, :)));
            end
            entry.medianEvals = median(squeeze(evals(m, iN, :, 2)), 'omitnan');
            if isempty(rows); rows = entry; else; rows(end+1) = entry; end %#ok<AGROW>
            if isnan(entry.height3Um)
                heightText = '       -';
            else
                heightText = sprintf('%.3f (%.3f)', entry.height3Um, ...
                    entry.height3IqrUm);
            end
            fprintf('  %-10s %6d  %9.2f  %6.2f  %14.3f  %15s  %6d\n', ...
                models{m}, opts.photonTotals(iN), threshold, power, ...
                entry.tauErrNs, heightText, round(entry.medianEvals));
        end
    end

    summary = struct2table(rows);
    out = struct('summary', summary, 'opts', opts, 'models', {models}, ...
        'weightTable', weightTable, 'weightHeights', weightHeights, ...
        'deltaD', deltaD, 'tauErr3', tauErr3, 'heightHat', heightHat, ...
        'atBound', atBound, 'evals', evals, 'irf', irf, 'timeNs', timeNs);
    csvFile = fullfile(opts.outputDir, 'ring_resolved_unmixing.csv');
    writetable(summary, csvFile);
    fprintf('\n  wrote %s\n', csvFile);

    if opts.makeFigure
        out.figure = plotUnmixing(out);
        fprintf('  wrote %s\n', out.figure);
    end

    fprintf(['\n  Reading this: every model is held to a %.0f%% false-' ...
        'positive rate by\n  construction, so the only figure of merit is ' ...
        'POWER - the rate at which a\n  genuine third pool is detected. ' ...
        'The ring models earn their keep only if\n  their power exceeds ' ...
        'the summed model at the same photon count.\n'], 100 * opts.alpha);
    fprintf(['  The z3 column is NOT a height measurement - its ' ...
        'interquartile range is\n  comparable to its own value, because the ' ...
        'normalised ring pattern saturates\n  above about 0.3 um. See ' ...
        'diagnose_ring_height_identifiability, and use\n  ' ...
        'study_ring_height_discrimination for the binary displaced-versus-' ...
        'coplanar\n  question the rings genuinely can answer. True third-' ...
        'pool height %.2f um.\n'], max([opts.truth.heightUm]));
end

% --------------------------------------------------------------- simulation

function counts = simulateRingCounts(truth, weights, irf, timeNs, ...
        periodNs, photonTotal)
    nRing = size(weights.table, 1);
    expected = zeros(nRing, numel(timeNs));
    for j = 1:numel(truth)
        w = ring_interpolate_weights(weights, truth(j).heightUm);
        pattern = ring_periodic_decay(irf, timeNs, periodNs, truth(j).tauNs);
        expected = expected + truth(j).molecules * (w(:) * pattern(:)');
    end
    total = sum(expected(:));
    if total <= 0
        counts = zeros(nRing, numel(timeNs));
        return;
    end
    expected = expected * (photonTotal / total);
    counts = ring_poisson_sample(expected);
end

function err = tauMatchError(tauFit, tauTrue)
    % Greedy nearest matching of fitted to true lifetimes, then the median
    % absolute error over the true components.
    tauFit = sort(tauFit(:)); tauTrue = sort(tauTrue(:));
    used = false(numel(tauFit), 1);
    errors = nan(numel(tauTrue), 1);
    for j = 1:numel(tauTrue)
        candidate = abs(tauFit - tauTrue(j));
        candidate(used) = Inf;
        [errors(j), pick] = min(candidate);
        if isfinite(errors(j)); used(pick) = true; end
    end
    err = median(errors, 'omitnan');
end

function q = quantileLocal(values, level)
    values = sort(values(~isnan(values)));
    if isempty(values); q = NaN; return; end
    if numel(values) == 1; q = values; return; end
    position = level * numel(values) + 0.5;
    position = min(max(position, 1), numel(values));
    low = floor(position); high = ceil(position);
    if low == high
        q = values(low);
    else
        q = values(low) + (position - low) * (values(high) - values(low));
    end
end

function printTruth(label, truth, weights)
    fprintf('%s:', label);
    for j = 1:numel(truth)
        w = ring_interpolate_weights(weights, truth(j).heightUm);
        fprintf(' [tau %.2f ns, z %.2f um, mol %.2f, detected %.3f]', ...
            truth(j).tauNs, truth(j).heightUm, truth(j).molecules, ...
            truth(j).molecules * sum(w));
    end
    fprintf('\n');
end

% ------------------------------------------------------------------ figures

function name = plotUnmixing(out)
    summary = out.summary;
    models = out.models;
    h = figure('Color', 'w', 'Visible', 'off', 'Position', [60 60 1280 420]);
    layout = tiledlayout(h, 1, 3, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    colours = lines(numel(models));

    ax = nexttile(layout);
    hold(ax, 'on');
    for m = 1:numel(models)
        keep = strcmp(summary.model, models{m});
        plot(ax, summary.photonTotal(keep), summary.powerAtFixedFP(keep), ...
            '-o', 'Color', colours(m, :), 'LineWidth', 1.6, ...
            'MarkerFaceColor', colours(m, :), 'DisplayName', models{m});
    end
    set(ax, 'XScale', 'log'); grid(ax, 'on'); ylim(ax, [0 1.05]);
    xlabel(ax, 'total detected photons');
    ylabel(ax, sprintf('detection rate at %.0f%% FP', 100 * out.opts.alpha));
    title(ax, 'power to detect the third pool');
    legend(ax, 'Location', 'southeast');

    ax = nexttile(layout);
    hold(ax, 'on');
    for m = 1:numel(models)
        keep = strcmp(summary.model, models{m});
        plot(ax, summary.photonTotal(keep), summary.tauErrNs(keep), ...
            '-o', 'Color', colours(m, :), 'LineWidth', 1.6, ...
            'MarkerFaceColor', colours(m, :), 'DisplayName', models{m});
    end
    set(ax, 'XScale', 'log', 'YScale', 'log'); grid(ax, 'on');
    xlabel(ax, 'total detected photons');
    ylabel(ax, 'median |\Delta\tau| (ns)');
    title(ax, 'lifetime accuracy, 3-component fit');

    ax = nexttile(layout);
    hold(ax, 'on');
    keep = strcmp(summary.model, 'ring-w(z)');
    plot(ax, summary.photonTotal(keep), summary.height3Um(keep), '-o', ...
        'Color', colours(3, :), 'LineWidth', 1.6, ...
        'MarkerFaceColor', colours(3, :), 'DisplayName', 'recovered');
    trueZ = max([out.opts.truth.heightUm]);
    yline(ax, trueZ, '--k', sprintf('truth %.2f um', trueZ));
    set(ax, 'XScale', 'log'); grid(ax, 'on');
    xlabel(ax, 'total detected photons');
    ylabel(ax, 'highest fitted z (\mum)');
    title(ax, 'height of the third pool, ring-w(z)');

    name = fullfile(out.opts.outputDir, 'ring_resolved_unmixing.png');
    exportgraphics(h, name, 'Resolution', 200);
    close(h);
end
