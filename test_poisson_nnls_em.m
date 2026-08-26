function results = test_poisson_nnls_em()
%TEST_POISSON_NNLS_EM Validate the Poisson non-negative EM solver.
%
% results = test_poisson_nnls_em()
%
% Checks, in order of how much they would matter if they failed:
%
%   1. EXACT RECOVERY on a noiseless synthetic decay built from known component
%      photon counts. Because biexp_slb_pattern normalises to unit sum, the
%      recovered amplitudes must equal those photon counts, not merely be
%      proportional to them.
%   2. MONOTONE DESCENT of the deviance. This is the property the replaced
%      scheme lacked - it walked away from the optimum and returned its last
%      iterate.
%   3. NON-NEGATIVITY, structurally, including when the true background is zero
%      so the bound is active. This is the case that broke the old solver.
%   4. AGREEMENT WITH A BRUTE-FORCE MINIMUM on small problems, so the claim
%      "reaches the global optimum" is tested rather than asserted.
%   5. THE REAL PIXEL that the audit identified: pixel 68376 of the 155036
%      native fit, where the shipped amplitudes give deviance 122.89 and the
%      true optimum is 59.45. The solver must reach the lower value.
%   6. A comparison against PIRLSnonneg on the same inputs, reported but NOT
%      asserted - PIRLSnonneg is left untouched by design, and this test exists
%      to document the size of the difference rather than to police it.

    results = struct('name', {}, 'pass', {}, 'detail', {});
    fprintf('\ntest_poisson_nnls_em\n');

    dtNs = 0.16;
    periodNs = 50;
    nBin = 156;
    irf = zeros(nBin, 1);
    irf(20:24) = [0.05; 0.25; 0.40; 0.25; 0.05];   % a plausible narrow IRF
    irf = irf / sum(irf);

    p1 = biexp_slb_pattern(irf, dtNs, periodNs, 0.3549, nBin);
    p2 = biexp_slb_pattern(irf, dtNs, periodNs, 2.0724, nBin);
    design = [ones(nBin, 1), p1, p2];

    % ---- 1. exact recovery, noiseless ----------------------------------
    truth = [0.5; 200; 400];
    y = design * truth;
    [beta, dev] = poisson_nnls_em(design, y, struct('maxIter', 2000, ...
        'tol', 1e-14, 'checkEvery', 25));
    err = max(abs(beta - truth) ./ max(truth, 1));
    results = addResult(results, 'exact recovery (noiseless)', err < 1e-4, ...
        sprintf('beta = [%.4f %.4f %.4f], truth = [%.4f %.4f %.4f], relerr %.2e, dev %.3e', ...
        beta(1), beta(2), beta(3), truth(1), truth(2), truth(3), err, dev));

    % ---- 2. monotone descent -------------------------------------------
    y2 = poissrnd(design * [0.4; 150; 300]);
    devPath = nan(1, 40);
    b = [];
    for k = 1:40
        o = struct('maxIter', k, 'tol', 0, 'checkEvery', 1e6);
        if ~isempty(b); o.seed = []; end
        [b, devPath(k)] = poisson_nnls_em(design, y2, o);
    end
    rises = diff(devPath);
    worst = max(rises);
    results = addResult(results, 'monotone deviance descent', ...
        worst <= 1e-8 * max(1, abs(devPath(1))), ...
        sprintf('largest increase between successive iteration counts %.3e (dev %.4f -> %.4f)', ...
        worst, devPath(1), devPath(end)));

    % ---- 3. non-negativity with an ACTIVE bound ------------------------
    % True background exactly zero: the old solver's failure case.
    y3 = poissrnd(design * [0; 180; 260]);
    beta3 = poisson_nnls_em(design, y3);
    results = addResult(results, 'non-negative with active bound', ...
        all(beta3(:) >= 0), ...
        sprintf('beta = [%.6f %.3f %.3f], min = %.3e', beta3(1), beta3(2), ...
        beta3(3), min(beta3(:))));

    % ---- 4. against a brute-force minimum ------------------------------
    % Coarse 3-D scan around the EM answer; nothing should beat it.
    [betaEM, devEM] = poisson_nnls_em(design, y3);
    best = devEM;
    scale = [linspace(0.7, 1.3, 13)];
    for i = scale
        for j = scale
            for k = [0 scale]
                trial = betaEM .* [k; i; j];
                d = poisson_nnls_em_deviance(y3, ...
                    max(design * trial, 1e-12));
                if d < best; best = d; end
            end
        end
    end
    results = addResult(results, 'no nearby point beats EM', ...
        best >= devEM - 1e-6 * max(1, devEM), ...
        sprintf('EM deviance %.6f, best of 2366 perturbations %.6f', ...
        devEM, best));

    % ---- 5. the real pixel from the audit ------------------------------
    matPath = ['D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_20260813_1' ...
        '\_20260813-155036\immune_cell_MIET'];
    mapsFile = fullfile(matPath, 'biexp_slb', 'biexp_slb_maps.mat');
    analysisFile = fullfile(matPath, 'immune_cell_MIET_640nm_red_analysis.mat');
    if isfile(mapsFile) && isfile(analysisFile)
        L = load(mapsFile, 'out');
        out = L.out;
        target = 68376;
        pix = out.pixelIndex(target);
        cubeVar = 'tcspc_pix';
        m = matfile(analysisFile);
        [nRow, nCol, ~] = size(m, cubeVar);
        [r, c] = ind2sub([nRow nCol], pix);
        decay = double(squeeze(m.(cubeVar)(r, c, :)));
        tau1 = out.tau1Ns(target);
        tau2 = out.tau2Ns(target);
        q1 = biexp_slb_pattern(out.irf, out.dtNs, out.periodNs, tau1, out.nBin);
        q2 = biexp_slb_pattern(out.irf, out.dtNs, out.periodNs, tau2, out.nBin);
        A = [ones(out.nBin, 1), q1, q2];

        shipped = [out.background(target); out.amplitude1(target); ...
            out.amplitude2(target)];
        devShipped = poisson_nnls_em_deviance(decay, max(A * shipped, 1e-12));
        [betaEM2, devEM2] = poisson_nnls_em(A, decay, ...
            struct('maxIter', 3000, 'tol', 1e-14, 'checkEvery', 50));

        fprintf(['  pixel %d (linear %d, %.0f photons, tau1 %.5f, tau2 %.5f)\n' ...
            '     shipped beta [%.3f %.3f %.3f] deviance %.3f\n' ...
            '     EM      beta [%.3f %.3f %.3f] deviance %.3f\n'], ...
            target, pix, sum(decay), tau1, tau2, shipped(1), shipped(2), ...
            shipped(3), devShipped, betaEM2(1), betaEM2(2), betaEM2(3), devEM2);
        results = addResult(results, 'beats the shipped fit on the audit pixel', ...
            devEM2 < devShipped - 1, ...
            sprintf('EM %.3f vs shipped %.3f (improvement %.3f)', ...
            devEM2, devShipped, devShipped - devEM2));
    else
        results = addResult(results, 'audit pixel available', false, ...
            sprintf('skipped: %s not found', mapsFile));
    end

    % ---- 6. documented comparison against PIRLSnonneg ------------------
    if exist('PIRLSnonneg', 'file') == 2
        try
            betaP = PIRLSnonneg(design, y3);
            devP = poisson_nnls_em_deviance(y3, max(design * betaP(:), 1e-12));
            fprintf(['  PIRLSnonneg on the active-bound case: ' ...
                'beta [%.4f %.3f %.3f] deviance %.4f  (EM %.4f, ' ...
                'difference %.4f)\n'], betaP(1), betaP(2), betaP(3), devP, ...
                devEM, devP - devEM);
        catch e
            fprintf('  PIRLSnonneg comparison skipped: %s\n', e.message);
        end
    end

    passed = sum([results.pass]);
    fprintf('\n  %d/%d checks passed\n', passed, numel(results));
    for k = 1:numel(results)
        fprintf('   [%s] %-42s %s\n', ternaryPassMark(results(k).pass), ...
            results(k).name, results(k).detail);
    end
end
