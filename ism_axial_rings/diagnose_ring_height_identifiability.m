function out = diagnose_ring_height_identifiability(opts)
%DIAGNOSE_RING_HEIGHT_IDENTIFIABILITY Is the fitted height real or an artefact?
%
% out = diagnose_ring_height_identifiability()
%
% study_ring_resolved_unmixing reports a recovered third-pool height of
% 0.488 um at every photon count from 1e3 to 3e4, when the truth is 0.40 um.
% Two things make that suspicious:
%
%   1. the bias is essentially independent of photon count (0.4876 to 0.4891
%      across a 30-fold range), which is the signature of model or optimiser
%      mis-specification, not of statistical error
%   2. the default height seed for a three-component fit is
%      linspace(zLow+0.01, 0.5, 3) = [0.01 0.255 0.5], whose third element is
%      0.5 - uncomfortably close to the reported 0.488
%
% So the reported height could be (a) a real measurement with a systematic
% bias, (b) the seed, barely moved because the objective is flat in z, or
% (c) a genuine flat direction, meaning height simply is not identifiable up
% there. These have completely different consequences and must be separated.
%
% THE TESTS, ORDERED SO THE DECISIVE ONE USES NO OPTIMISER
%
% 1. PATTERN SEPARATION. Print the normalised ring pattern w(:,z)/sum(w(:,z))
%    at a range of heights and its distance from the pattern at the true
%    0.40 um. If the pattern barely changes above 0.3 um then height is not
%    identifiable there and no fitter could recover it - a property of the
%    optics, not of the code.
%
% 2. NOISELESS PROFILE. Fix the lifetimes and the other two heights at truth,
%    scan the third height across the whole table, and find the minimum of the
%    deviance against NOISELESS expected counts. fminsearch is not involved.
%    If this minimum is at 0.40, the model is correctly specified and the bias
%    comes from noise or from the optimiser. If it is at 0.488, the forward
%    model or the interpolation is wrong. If the profile is flat, see test 1.
%
% 3. NOISY PROFILE. The same one-dimensional scan on Poisson data, repeated, to
%    get the distribution of the profile minimum. Still no optimiser. This is
%    the honest measure of how well height can be estimated at a given photon
%    count.
%
% 4. SEED SENSITIVITY. Full fits from deliberately different height seeds,
%    including one seeded at the truth. If the recovered height tracks the
%    seed, the reported value is the seed.
%
% 5. MULTISTART. Full fits with eight random restarts, keeping the best
%    deviance. If the answer moves toward 0.40, the single-start fit was stuck.

    thisDir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(opts); opts = struct(); end
    defaults = struct( ...
        'truth', [], 'photonTotals', [1e4 3e4], 'repeats', 15, ...
        'dtNs', 50/156, 'periodNs', 50, 'irfFwhmNs', 0.35, 'nBins', 156, ...
        'restarts', 8, 'outputDir', thisDir, 'seed', 5);
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
    rng(opts.seed);

    fprintf('computing the ring weight table...\n');
    ring = simulate_ring_weight_vs_height(struct('makeFigure', false, ...
        'heightsUm', 0:0.02:1.0, 'nx', 41, 'fovUm', 2.6));
    weights = struct('table', ring.ringWeight ./ max(ring.ringWeight(:)), ...
        'heights', ring.heightsUm);
    nRing = size(weights.table, 1);

    timeNs = (0:opts.nBins-1)' * opts.dtNs;
    irf = exp(-0.5 * ((timeNs - 12 * opts.dtNs) / ...
        (opts.irfFwhmNs / (2 * sqrt(2 * log(2))))) .^ 2);
    irf = irf / sum(irf);

    trueTau = [opts.truth.tauNs];
    trueZ = [opts.truth.heightUm];
    targetZ = trueZ(end);

    fprintf('\ndiagnose_ring_height_identifiability\n');
    fprintf('  %d rings, height grid %.2f to %.2f um in %.3f um steps\n', ...
        nRing, weights.heights(1), weights.heights(end), ...
        weights.heights(2) - weights.heights(1));
    fprintf('  truth tau %s ns at z %s um\n', mat2str(trueTau), mat2str(trueZ));

    % ---- 1. pattern separation --------------------------------------------
    fprintf('\n  [1] ring pattern vs height (normalised), and distance from z=%.2f\n', ...
        targetZ);
    fprintf('      z(um)   sum(w)   %s   |pattern - pattern(%.2f)|\n', ...
        sprintf('r%-6d', 1:nRing), targetZ);
    probe = [0 0.02 0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.0];
    reference = normalisePattern(ring_interpolate_weights(weights, targetZ));
    patternDistance = nan(numel(probe), 1);
    for k = 1:numel(probe)
        w = ring_interpolate_weights(weights, probe(k));
        p = normalisePattern(w);
        patternDistance(k) = norm(p - reference);
        fprintf('      %5.2f   %6.3f   %s   %.5f\n', probe(k), sum(w), ...
            sprintf('%-7.4f', p), patternDistance(k));
    end
    % How far must z move from the truth before the pattern changes by as much
    % as photon noise would? Poisson noise on a fraction f from N photons has
    % standard deviation roughly sqrt(f/N), so this sets the resolution scale.
    gridDistance = arrayfun(@(z) norm(normalisePattern( ...
        ring_interpolate_weights(weights, z)) - reference), weights.heights);
    fprintf('      pattern distance from z=%.2f: max over the whole grid %.4f\n', ...
        targetZ, max(gridDistance));
    fprintf('      distance to z=0.49 (the reported value) %.5f\n', ...
        norm(normalisePattern(ring_interpolate_weights(weights, 0.49)) - reference));

    % ---- 2. noiseless profile ---------------------------------------------
    expected = expectedCounts(opts.truth, weights, irf, timeNs, opts.periodNs);
    noiseless = expected * (1e6 / sum(expected(:)));
    [profileZ, profileD] = heightProfile(noiseless, trueTau, trueZ, weights, ...
        irf, timeNs, opts.periodNs);
    [~, iMin] = min(profileD);
    fprintf(['\n  [2] noiseless profile over the third height, no optimiser:' ...
        '\n      minimum at z = %.3f um (truth %.3f), deviance %.4g\n'], ...
        profileZ(iMin), targetZ, profileD(iMin));
    fprintf('      deviance at truth %.4g, at 0.49 um %.4g, span over grid %.4g\n', ...
        interp1(profileZ, profileD, targetZ), ...
        interp1(profileZ, profileD, 0.49), max(profileD) - min(profileD));

    % ---- 3. noisy profile --------------------------------------------------
    fprintf('\n  [3] noisy profile minimum, no optimiser (median and IQR over %d reps)\n', ...
        opts.repeats);
    noisyZ = nan(numel(opts.photonTotals), opts.repeats);
    for iN = 1:numel(opts.photonTotals)
        scaled = expected * (opts.photonTotals(iN) / sum(expected(:)));
        for rep = 1:opts.repeats
            counts = ringPoisson(scaled);
            [zGrid, dGrid] = heightProfile(counts, trueTau, trueZ, weights, ...
                irf, timeNs, opts.periodNs);
            [~, j] = min(dGrid);
            noisyZ(iN, rep) = zGrid(j);
        end
        values = noisyZ(iN, :);
        fprintf('      N=%6d  median z3 %.3f um, IQR [%.3f %.3f], at grid top %d/%d\n', ...
            opts.photonTotals(iN), median(values), ...
            prctileLocal(values, 25), prctileLocal(values, 75), ...
            nnz(values >= weights.heights(end) - 1e-9), opts.repeats);
    end

    % ---- 4. seed sensitivity ----------------------------------------------
    fprintf('\n  [4] seed sensitivity of the full three-component fit, N=%d\n', ...
        opts.photonTotals(1));
    seedSets = {[0.01 0.06 0.12], [0.01 0.255 0.50], [0.02 0.30 0.90], trueZ};
    seedLabel = {'low', 'study default', 'high', 'at the truth'};
    scaled = expected * (opts.photonTotals(1) / sum(expected(:)));
    seedResult = nan(numel(seedSets), opts.repeats);
    seedDeviance = nan(numel(seedSets), opts.repeats);
    % The SAME realisations are reused across seeds. Drawing fresh data per seed
    % would confound the seed effect with sampling noise and would make the
    % deviance column incomparable between rows.
    pairedCounts = cell(opts.repeats, 1);
    for rep = 1:opts.repeats
        pairedCounts{rep} = ringPoisson(scaled);
    end
    for s = 1:numel(seedSets)
        for rep = 1:opts.repeats
            counts = pairedCounts{rep};
            fit = ring_fit_height_constrained(counts, irf, timeNs, ...
                opts.periodNs, 3, weights, ...
                struct('seedHeightUm', seedSets{s}));
            seedResult(s, rep) = max(fit.height);
            seedDeviance(s, rep) = fit.deviance;
        end
        fprintf('      seed %-14s z3 seed %.3f -> recovered median %.3f um, deviance %.1f\n', ...
            seedLabel{s}, max(seedSets{s}), median(seedResult(s, :)), ...
            median(seedDeviance(s, :)));
    end

    % ---- 5. multistart ----------------------------------------------------
    fprintf('\n  [5] multistart, %d restarts, N=%d\n', opts.restarts, ...
        opts.photonTotals(1));
    singleZ = nan(opts.repeats, 1); multiZ = nan(opts.repeats, 1);
    singleD = nan(opts.repeats, 1); multiD = nan(opts.repeats, 1);
    for rep = 1:opts.repeats
        counts = pairedCounts{rep};
        one = ring_fit_height_constrained(counts, irf, timeNs, ...
            opts.periodNs, 3, weights, struct('restarts', 1));
        many = ring_fit_height_constrained(counts, irf, timeNs, ...
            opts.periodNs, 3, weights, struct('restarts', opts.restarts));
        singleZ(rep) = max(one.height); multiZ(rep) = max(many.height);
        singleD(rep) = one.deviance;    multiD(rep) = many.deviance;
    end
    fprintf('      single start  median z3 %.3f um, median deviance %.2f\n', ...
        median(singleZ), median(singleD));
    fprintf('      %2d restarts   median z3 %.3f um, median deviance %.2f\n', ...
        opts.restarts, median(multiZ), median(multiD));
    fprintf('      restarts found a strictly better optimum in %d/%d reps\n', ...
        nnz(multiD < singleD - 1e-6), opts.repeats);

    out = struct('weights', weights, 'irf', irf, 'timeNs', timeNs, ...
        'opts', opts, 'probeHeights', probe, ...
        'patternDistance', patternDistance, ...
        'profileZ', profileZ, 'profileDeviance', profileD, ...
        'noisyProfileZ', noisyZ, 'seedSets', {seedSets}, ...
        'seedRecoveredZ', seedResult, 'seedDeviance', seedDeviance, ...
        'singleStartZ', singleZ, 'multiStartZ', multiZ, ...
        'singleStartDeviance', singleD, 'multiStartDeviance', multiD);
end

% ---------------------------------------------------------------- helpers

function p = normalisePattern(w)
    total = sum(w);
    if total <= 0; p = zeros(size(w))'; return; end
    p = (w / total)';
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

function [zGrid, deviance] = heightProfile(counts, tau, heightTruth, ...
        weights, irf, timeNs, periodNs)
%HEIGHTPROFILE Deviance as the LAST height is scanned, everything else fixed.
% No optimiser is involved, so the result reflects the information in the data
% and nothing about fminsearch.
    zGrid = weights.heights(:)';
    deviance = nan(size(zGrid));
    height = heightTruth;
    for k = 1:numel(zGrid)
        height(end) = zGrid(k);
        deviance(k) = ring_constrained_deviance(tau, height, counts, irf, ...
            timeNs, periodNs, weights);
    end
end

function counts = ringPoisson(lambda)
    % Exact Poisson draws by inversion for small means and a corrected normal
    % approximation for large ones. The normal branch is only used where the
    % relative error is negligible, because the statistic being studied is a
    % deviance computed under an exact-Poisson assumption.
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

function value = prctileLocal(values, level)
    values = sort(values(~isnan(values)));
    if isempty(values); value = NaN; return; end
    if numel(values) == 1; value = values; return; end
    position = min(max(level / 100 * numel(values) + 0.5, 1), numel(values));
    low = floor(position); high = ceil(position);
    if low == high
        value = values(low);
    else
        value = values(low) + (position - low) * (values(high) - values(low));
    end
end
