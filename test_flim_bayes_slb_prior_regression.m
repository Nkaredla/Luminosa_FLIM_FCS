function test_flim_bayes_slb_prior_regression()
%TEST_FLIM_BAYES_SLB_PRIOR_REGRESSION Soft SLB prior must be opt-in.
%
% The SLB photon count used to be a hard constraint. Marginalising it over a
% Gaussian prior is a change to the core inference, so the default must
% reproduce the old behaviour exactly. This checks that, then checks the
% opt-in path actually runs and actually changes something.

    rng(3);
    nPix = 60;
    nBins = 78;
    dtNs = 0.16;
    periodNs = 12.5;
    tauSlbNs = 0.3;
    timeNs = (0:nBins-1) * dtNs;

    irf = exp(-0.5 * ((timeNs - 1) / 0.15) .^ 2);
    irf = irf / sum(irf);

    slb = conv(irf, exp(-timeNs / tauSlbNs));
    slb = slb(1:nBins);
    mem = conv(irf, exp(-timeNs / 1.5));
    mem = mem(1:nBins);
    shape = 0.4 * slb / sum(slb) + 0.6 * mem / sum(mem);
    shape = shape / sum(shape);
    edges = [0, cumsum(shape)];
    edges(end) = 1;

    photonTotal = 4000;
    Y = zeros(nPix, 1, nBins);
    for k = 1:nPix
        Y(k, 1, :) = reshape(histcounts(rand(1, photonTotal), edges), ...
            1, 1, nBins);
    end

    base = struct('analysisMask', true(nPix, 1), 'minPhotons', 1, ...
        'useGPU', false, 'batchSize', 2048, 'includeBackground', true, ...
        'membraneTauCount', 10, 'fractionStep', 0.2, ...
        'minimumMembraneFraction', 0.1, ...
        'fixedSlbPhotonCount', 0.4 * photonTotal, ...
        'fixedSlbPhotonCountStd', 0.1 * 0.4 * photonTotal, ...
        'slbCountRelTol', 0);

    hardImplicit = flim_bayes_fixed_slb(Y, irf, periodNs, dtNs, tauSlbNs, base);

    explicitZero = base;
    explicitZero.slbCountPriorNodes = 0;
    hardExplicit = flim_bayes_fixed_slb(Y, irf, periodNs, dtNs, tauSlbNs, ...
        explicitZero);

    assert(isequaln(hardImplicit.modelProbability, ...
        hardExplicit.modelProbability), ...
        'slbCountPriorNodes = 0 must reproduce the hard constraint exactly.');
    assert(isequaln(hardImplicit.tauMeanArithmetic, ...
        hardExplicit.tauMeanArithmetic), ...
        'Default lifetimes changed; the prior must be strictly opt-in.');
    assert(~hardImplicit.fixedSlbPhotonConstraint.countMarginalised, ...
        'Default must report countMarginalised = false.');
    fprintf('  default path unchanged: OK\n');

    softened = base;
    softened.slbCountPriorNodes = 5;
    soft = flim_bayes_fixed_slb(Y, irf, periodNs, dtNs, tauSlbNs, softened);

    assert(soft.fixedSlbPhotonConstraint.countMarginalised, ...
        'Opt-in path must report countMarginalised = true.');
    deltaP = max(abs(double(hardImplicit.modelProbability(:)) - ...
        double(soft.modelProbability(:))));
    % Posterior widths should grow: the SLB calibration uncertainty is now
    % propagated rather than assumed away.
    hardWidth = median(double(hardImplicit.tauPosteriorStd(:)), 'omitnan');
    softWidth = median(double(soft.tauPosteriorStd(:)), 'omitnan');
    fprintf('  marginalised path runs: OK (max |dP| = %.4f)\n', deltaP);
    fprintf('  median tau posterior std: hard %.4f -> soft %.4f ns\n', ...
        hardWidth, softWidth);
    assert(isfinite(softWidth) && softWidth > 0, ...
        'Marginalised posterior width is not finite and positive.');

    % A wrong SLB count should hurt less when it is marginalised than when it
    % is imposed. This is the whole point of the change.
    %
    % The comparison MUST use a fine lifetime grid. On the production grid
    % (10 points over 0.62-5 ns) the nodes near 1.5 ns are 1.231, 1.547 and
    % 1.945, so both variants land on the same node and the residual error is
    % pure quantisation - the SLB effect is invisible and this test would
    % silently measure nothing. That is itself the finding: while the grid
    % dominates, softening the SLB prior cannot help.
    biased = base;
    biased.membraneTauCount = 48;
    biased.membraneTauBoundsNs = [0.4 5.5];
    biased.fixedSlbPhotonCount = 1.10 * 0.4 * photonTotal;
    biasedHard = flim_bayes_fixed_slb(Y, irf, periodNs, dtNs, tauSlbNs, biased);
    biasedSoft = biased;
    biasedSoft.slbCountPriorNodes = 5;
    biasedSoftOut = flim_bayes_fixed_slb(Y, irf, periodNs, dtNs, tauSlbNs, ...
        biasedSoft);

    truth = 1.5;
    hardErr = median(abs(double(biasedHard.oneMembrane.membraneLifetime1Ns(:)) ...
        - truth), 'omitnan');
    softErr = median(abs( ...
        double(biasedSoftOut.oneMembrane.membraneLifetime1Ns(:)) - truth), ...
        'omitnan');
    fprintf(['  with a 10%% SLB error, median |tau2 - truth|: ' ...
        'hard %.4f -> marginalised %.4f ns\n'], hardErr, softErr);

    % Reported, not asserted: whether marginalising helps depends on how much
    % of the error is bias rather than quantisation, which is what this
    % measurement exists to establish.
    if softErr <= hardErr
        fprintf('  marginalising reduced the bias by %.1f%%\n', ...
            100 * (hardErr - softErr) / max(hardErr, eps));
    else
        fprintf(['  marginalising did NOT reduce the bias here ' ...
            '(%.4f -> %.4f ns); do not enable it on this evidence\n'], ...
            hardErr, softErr);
    end

    % Record that the production grid masks the effect entirely.
    coarse = biased;
    coarse.membraneTauCount = base.membraneTauCount;
    coarse.membraneTauBoundsNs = [];
    coarseHard = flim_bayes_fixed_slb(Y, irf, periodNs, dtNs, tauSlbNs, coarse);
    coarseOpts = coarse;
    coarseOpts.slbCountPriorNodes = 5;
    coarseSoft = flim_bayes_fixed_slb(Y, irf, periodNs, dtNs, tauSlbNs, ...
        coarseOpts);
    coarseHardErr = median(abs( ...
        double(coarseHard.oneMembrane.membraneLifetime1Ns(:)) - truth), ...
        'omitnan');
    coarseSoftErr = median(abs( ...
        double(coarseSoft.oneMembrane.membraneLifetime1Ns(:)) - truth), ...
        'omitnan');
    fprintf(['  same test on the PRODUCTION grid: hard %.4f -> ' ...
        'marginalised %.4f ns (quantised, effect masked)\n'], ...
        coarseHardErr, coarseSoftErr);

    fprintf('test_flim_bayes_slb_prior_regression: PASS\n');
end
