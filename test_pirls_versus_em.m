function results = test_pirls_versus_em(mapsMat, analysisMat, sampleSize)
%TEST_PIRLS_VERSUS_EM Head-to-head amplitude solve: PIRLSnonneg vs poisson_nnls_em.
%
% results = test_pirls_versus_em()
% results = test_pirls_versus_em(mapsMat, analysisMat, sampleSize)
%
% Both solvers are handed the SAME decay and the SAME (tau1, tau2), so the only
% difference is how the three non-negative amplitudes [B; a1; a2] are found. They
% are minimising the same thing - the Poisson deviance - so the one reaching the
% lower value is simply better, and no judgement call is involved.
%
% WHY THIS COMPARISON EXISTS
%
% PIRLSnonneg solves lsqnonneg(x'*w*x, x'*w*y). Passing the GRAM MATRIX to
% lsqnonneg as though it were a design matrix minimises the norm of the
% NORMAL-EQUATION RESIDUAL, not the weighted residual. Unconstrained the two
% coincide, so this is invisible in easy cases; with a bound active they do not.
% Here the background is driven to zero on roughly half the pixels, so the
% bound-active case is the common one rather than an edge case.
%
% A USEFUL SELF-CHECK IS BUILT IN: poisson_nnls_em converges to the global
% optimum of a concave problem, so PIRLSnonneg should NEVER beat it by more than
% numerical noise. If this test reports PIRLS winning anywhere by a real margin,
% the EM solver is not converging and that is the finding, not the other way
% round.
%
% PIRLSnonneg is NOT modified by this or any other file in this change - other
% pipelines depend on it. This only measures.

    if nargin < 1 || isempty(mapsMat)
        mapsMat = ['D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_20260813_1' ...
            '\_20260813-160712\immune_cell_MIET\biexp_slb_v2\biexp_slb_maps.mat'];
    end
    if nargin < 2 || isempty(analysisMat)
        analysisMat = ['D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_' ...
            '20260813_1\_20260813-160712\immune_cell_MIET\' ...
            'immune_cell_MIET_640nm_red_analysis.mat'];
    end
    if nargin < 3 || isempty(sampleSize); sampleSize = 1500; end

    L = load(mapsMat, 'out');
    out = L.out;
    m = matfile(analysisMat);
    [nRow, nCol, ~] = size(m, 'tcspc_pix');

    nFit = numel(out.pixelIndex);
    % Deterministic, evenly spread sample - no RNG, so this is reproducible.
    pick = round(linspace(1, nFit, min(sampleSize, nFit)));
    pick = unique(pick);
    nPick = numel(pick);

    devEM = nan(nPick, 1);
    devPI = nan(nPick, 1);
    betaEM = nan(3, nPick);
    betaPI = nan(3, nPick);
    photons = nan(nPick, 1);
    pirlsFailed = 0;

    fprintf('\ntest_pirls_versus_em\n  %d pixel(s) sampled from %d fitted\n', ...
        nPick, nFit);
    started = tic;
    for k = 1:nPick
        row = pick(k);
        [r, c] = ind2sub([nRow nCol], out.pixelIndex(row));
        y = double(squeeze(m.tcspc_pix(r, c, :)));
        photons(k) = sum(y);
        p1 = biexp_slb_pattern(out.irf, out.dtNs, out.periodNs, ...
            out.tau1Ns(row), out.nBin);
        p2 = biexp_slb_pattern(out.irf, out.dtNs, out.periodNs, ...
            out.tau2Ns(row), out.nBin);
        design = [ones(out.nBin, 1), p1, p2];

        [b, d] = poisson_nnls_em(design, y, struct('maxIter', 2000, ...
            'tol', 1e-14, 'checkEvery', 25));
        betaEM(:, k) = b;
        devEM(k) = d;

        try
            bp = PIRLSnonneg(design, y);
            bp = max(bp(:), 0);
            betaPI(:, k) = bp;
            devPI(k) = poisson_nnls_em_deviance(y, ...
                max(design * bp, 1e-12));
        catch
            pirlsFailed = pirlsFailed + 1;
        end
        if mod(k, 250) == 0
            fprintf('    %d/%d (%.0f s)\n', k, nPick, toc(started));
        end
    end

    ok = isfinite(devEM) & isfinite(devPI);
    gap = devPI(ok) - devEM(ok);          % positive means EM is better
    atBound = betaEM(1, ok)' <= 1e-9;

    results = struct();
    results.sampled = nPick;
    results.pirlsFailed = pirlsFailed;
    results.emBetterFraction = mean(gap > 1e-6);
    results.pirlsBetterFraction = mean(gap < -1e-6);
    results.gapMedian = median(gap);
    results.gapP90 = prctile(gap, 90);
    results.gapMax = max(gap);
    results.gapWorstForPirls = min(gap);
    results.gapMedianAtBound = median(gap(atBound));
    results.gapMedianOffBound = median(gap(~atBound));
    results.boundActiveFraction = mean(atBound);
    results.a1RatioMedian = median(betaPI(2, ok)' ./ ...
        max(betaEM(2, ok)', eps));
    results.a2RatioMedian = median(betaPI(3, ok)' ./ ...
        max(betaEM(3, ok)', eps));
    results.photonsMedian = median(photons);

    fprintf('\n  RESULT (positive gap = EM reaches the lower deviance)\n');
    fprintf('    background at the bound      %.1f%% of sampled pixels\n', ...
        100 * results.boundActiveFraction);
    fprintf('    EM strictly better           %.1f%%\n', ...
        100 * results.emBetterFraction);
    fprintf('    PIRLS strictly better        %.2f%%  (should be ~0)\n', ...
        100 * results.pirlsBetterFraction);
    fprintf('    deviance gap  median %.3f | p90 %.3f | max %.3f\n', ...
        results.gapMedian, results.gapP90, results.gapMax);
    fprintf('      where bound ACTIVE   median %.3f\n', ...
        results.gapMedianAtBound);
    fprintf('      where bound INACTIVE median %.3f\n', ...
        results.gapMedianOffBound);
    fprintf('    amplitude ratio PIRLS/EM     a1 %.4f | a2 %.4f\n', ...
        results.a1RatioMedian, results.a2RatioMedian);
    if results.pirlsBetterFraction > 0.01
        fprintf(['    WARNING: PIRLS beats EM on more than 1%% of pixels, ' ...
            'so the EM solver is\n    not converging - that is the finding ' ...
            'here, not a PIRLS defect.\n']);
    end
    fprintf('    worst case for PIRLS (most negative gap) %.3e\n', ...
        results.gapWorstForPirls);
end
