function out = immune_cell_MIET_biexp_run(cube, mask, pixelIndex, intensity, ...
        irf, dtNs, periodNs, nBin, tau1Grid, tau2Grid, opts)
%IMMUNE_CELL_MIET_BIEXP_RUN Fit every selected pixel; return maps.
%
% Four stages. There are four rather than two because of a defect found by audit
% in the previous version, which biased the reported tau2 upward by 5-7.4%.
%
% WHAT WENT WRONG BEFORE
%
% The old code ranked the (tau1, tau2) grid by
%
%     sum((Y - model).^2 ./ max(model, 1))
%
% with amplitudes from an unweighted least-squares solve clamped at zero, and
% then reported tau2Grid(winner) directly - the later refit re-solved only the
% three amplitudes at the already-chosen node. So THE RANKING ERROR WAS THE
% REPORTED LIFETIME. At 350-700 photons over 156 bins the model sits below one
% count in roughly 80% of bins, and max(model, 1) discards the correct 1/model
% Poisson weight in exactly those bins. Against injected truth the bias was
% one-directional: truth 1.9302 -> 2.0724 (+7.4%), truth 1.9716 -> 2.0724
% (+5.1%), collapsing only at ~100x the photon count. Re-ranking the same grid
% by the Poisson deviance moves 84.3% of pixels and the median tau2 from 2.0724
% to 1.9302 ns.
%
% That the bias shrinks with photon count is why this had to be fixed before any
% binned comparison: spatial binning raises the photon count, so the artefact
% would have masqueraded as a real tau2 trend with bin size.
%
% THE STAGES
%
%   1. SHORTLIST using the cheap score above - kept only to rank candidates,
%      never to report one. A cheap score is legitimate for that, since it only
%      has to keep the true optimum inside a short list; how often it fails is
%      measured and reported as shortlistEdgeFraction.
%   2. EXACT RE-RANK. For each grid pair, gather the pixels that shortlisted it,
%      solve the amplitudes by Poisson maximum likelihood (poisson_nnls_em) and
%      score with the exact Poisson deviance plus the tau1 prior penalty.
%      Looping over PAIRS rather than pixels is what makes this affordable: each
%      pair's design is shared across its member pixels, so the work scales with
%      the shortlist length instead of the whole grid.
%   3. PARABOLIC REFINEMENT of tau2 between nodes, in log tau2, from the
%      deviances at the winner and its two neighbours. This removes the 7.37%
%      node quantisation that made six separate acquisitions all report
%      2.07237064230134. Guarded for convexity and kept inside the bracket.
%   4. FINAL SOLVE at the refined tau2 - the source of the reported amplitudes,
%      background and deviance.
%
% Amplitudes come from poisson_nnls_em throughout rather than PIRLSnonneg,
% because that function passes a Gram matrix to lsqnonneg, which minimises the
% norm of the normal-equation residual instead of the weighted residual and so
% returns the wrong answer whenever a bound is active - here the background
% clamps to zero on about 54% of pixels. PIRLSnonneg is deliberately left
% unmodified, since other pipelines depend on it.

    nPixel = numel(pixelIndex);
    n1 = numel(tau1Grid);
    n2 = numel(tau2Grid);
    shortlist = max(1, min(opts.shortlist, n1 * n2));

    % Patterns and their PRE-NORMALISATION sums. The sums are needed for the
    % species fraction: for a unit-sum basis the pre-exponential weight is
    % a_i / S_i, not a_i / tau_i.
    pattern1 = zeros(nBin, n1);
    sum1 = zeros(n1, 1);
    for a = 1:n1
        [pattern1(:, a), sum1(a)] = ...
            biexp_slb_pattern(irf, dtNs, periodNs, tau1Grid(a), nBin);
    end
    pattern2 = zeros(nBin, n2);
    for b = 1:n2
        pattern2(:, b) = biexp_slb_pattern(irf, dtNs, periodNs, ...
            tau2Grid(b), nBin);
    end
    penalty1 = ((tau1Grid(:) - opts.slbTauNs) / max(opts.slbSigmaNs, eps)) .^ 2;

    % Pair bookkeeping: linear index q <-> (a, b).
    [aGrid, bGrid] = ndgrid(1:n1, 1:n2);
    aOf = aGrid(:); bOf = bGrid(:);
    pairOk = tau2Grid(bOf)' > tau1Grid(aOf)';
    pairOk = pairOk(:);
    nPair = numel(aOf);
    logTau2 = log(tau2Grid(:));
    logStep = mean(diff(logTau2));

    tau1Hat = nan(nPixel, 1);
    tau2Hat = nan(nPixel, 1);
    tau1Index = nan(nPixel, 1);
    amp1 = nan(nPixel, 1);
    amp2 = nan(nPixel, 1);
    background = nan(nPixel, 1);
    deviance = nan(nPixel, 1);
    speciesRaw1 = nan(nPixel, 1);
    speciesRaw2 = nan(nPixel, 1);
    refined = false(nPixel, 1);
    atEdgeOfShortlist = false(nPixel, 1);

    ones1 = ones(nBin, 1);
    cubeFlat = reshape(cube, [], nBin);
    blockCount = ceil(nPixel / opts.blockSize);
    fprintf(['  fitting %d pixel(s) in %d block(s); shortlist %d of %d ' ...
        'valid pair(s)\n'], nPixel, blockCount, shortlist, nnz(pairOk));
    started = tic;

    for blockIndex = 1:blockCount
        lo = (blockIndex - 1) * opts.blockSize + 1;
        hi = min(nPixel, blockIndex * opts.blockSize);
        idx = pixelIndex(lo:hi);
        nThis = numel(idx);
        Y = double(cubeFlat(idx, :))';

        % ---- stage 1: shortlist -----------------------------------------
        score = inf(nPair, nThis);
        for q = 1:nPair
            if ~pairOk(q); continue; end
            design = [ones1, pattern1(:, aOf(q)), pattern2(:, bOf(q))];
            beta = max((design' * design) \ (design' * Y), 0);
            model = design * beta;
            model(model < 1e-9) = 1e-9;
            score(q, :) = sum((Y - model) .^ 2 ./ max(model, 1), 1) + ...
                penalty1(aOf(q));
        end
        [~, ranked] = sort(score, 1);
        top = ranked(1:shortlist, :);

        % ---- stage 2: exact Poisson re-rank over the shortlist ----------
        bestDev = inf(1, nThis);
        bestQ = ones(1, nThis);
        bestRank = ones(1, nThis);
        for q = 1:nPair
            if ~pairOk(q); continue; end
            member = any(top == q, 1);
            if ~any(member); continue; end
            cols = find(member);
            design = [ones1, pattern1(:, aOf(q)), pattern2(:, bOf(q))];
            [~, dev] = poisson_nnls_em(design, Y(:, cols));
            total = dev + penalty1(aOf(q));
            better = total < bestDev(cols);
            if any(better)
                hit = cols(better);
                bestDev(hit) = total(better);
                bestQ(hit) = q;
                [~, rankOfQ] = max(top(:, hit) == q, [], 1);
                bestRank(hit) = rankOfQ;
            end
        end
        atEdgeOfShortlist(lo:hi) = bestRank(:) == shortlist;

        aBest = aOf(bestQ);
        bBest = bOf(bestQ);
        tau2Fit = tau2Grid(bBest)';
        tau2Fit = tau2Fit(:)';

        % ---- stage 3: parabolic refinement of tau2 ----------------------
        if opts.refineTau2
            devLeft = inf(1, nThis);
            devRight = inf(1, nThis);
            for q = 1:nPair
                if ~pairOk(q); continue; end
                sameTau1 = aBest(:)' == aOf(q);
                isLeftNode = sameTau1 & (bBest(:)' == bOf(q) + 1);
                isRightNode = sameTau1 & (bBest(:)' == bOf(q) - 1);
                cols = find(isLeftNode | isRightNode);
                if isempty(cols); continue; end
                design = [ones1, pattern1(:, aOf(q)), pattern2(:, bOf(q))];
                [~, dev] = poisson_nnls_em(design, Y(:, cols));
                dev = dev + penalty1(aOf(q));
                pickLeft = isLeftNode(cols);
                devLeft(cols(pickLeft)) = dev(pickLeft);
                devRight(cols(~pickLeft)) = dev(~pickLeft);
            end
            y0 = devLeft; y1 = bestDev; y2 = devRight;
            curvature = y0 - 2 * y1 + y2;
            usable = isfinite(y0) & isfinite(y2) & curvature > 0;
            shift = zeros(1, nThis);
            shift(usable) = 0.5 * logStep * ...
                (y0(usable) - y2(usable)) ./ curvature(usable);
            % never leave the bracket the three points define
            shift = max(min(shift, logStep), -logStep);
            shift(~usable) = 0;
            tau2Fit = exp(logTau2(bBest)' + shift);
            tau2Fit = tau2Fit(:)';
            refined(lo:hi) = usable(:);
        end

        % ---- stage 4: final solve at the refined tau2 -------------------
        for p = 1:nThis
            a = aBest(p);
            [q2, s2] = biexp_slb_pattern(irf, dtNs, periodNs, ...
                tau2Fit(p), nBin);
            design = [ones1, pattern1(:, a), q2];
            [beta, dev] = poisson_nnls_em(design, Y(:, p), ...
                struct('maxIter', 400, 'tol', 1e-12, 'checkEvery', 20));
            row = lo + p - 1;
            tau1Hat(row) = tau1Grid(a);
            tau1Index(row) = a;
            tau2Hat(row) = tau2Fit(p);
            background(row) = beta(1);
            amp1(row) = beta(2);
            amp2(row) = beta(3);
            deviance(row) = dev;
            speciesRaw1(row) = beta(2) / max(sum1(a), eps);
            speciesRaw2(row) = beta(3) / max(s2, eps);
        end
        fprintf('    block %d/%d done (%.0f s elapsed)\n', blockIndex, ...
            blockCount, toc(started));
    end

    % ---- assemble maps ---------------------------------------------------
    imageSize = size(mask);
    blank = @() nan(imageSize);
    maps = struct();
    maps.tau1Ns = blank(); maps.tau1Ns(pixelIndex) = tau1Hat;
    maps.tau2Ns = blank(); maps.tau2Ns(pixelIndex) = tau2Hat;
    % True where tau1 landed on the first or last available node: the fit wanted
    % to keep going and the grid stopped it. A two-sided sigma threshold cannot
    % see this, because the positivity filter drops only LOW nodes and so can
    % leave a pull of -2 unreachable while +2 stays available.
    maps.tau1AtGridEdge = false(imageSize);
    maps.tau1AtGridEdge(pixelIndex) = tau1Index == 1 | tau1Index == n1;
    maps.slbPullSigma = blank();
    maps.slbPullSigma(pixelIndex) = ...
        (tau1Hat - opts.slbTauNs) / max(opts.slbSigmaNs, eps);

    % PHOTON fraction. biexp_slb_pattern normalises every column to unit sum, so
    % a fitted amplitude already IS that component's photon count - verified to
    % 2e-16 by recovering 200 and 400 from a decay synthesised with those
    % counts. It must therefore not be weighted by tau again.
    photonTotal = amp1 + amp2;
    maps.photonFraction1 = blank();
    maps.photonFraction1(pixelIndex) = amp1 ./ max(photonTotal, eps);
    maps.photonFraction2 = blank();
    maps.photonFraction2(pixelIndex) = amp2 ./ max(photonTotal, eps);

    % SPECIES (pre-exponential) fraction: a_i / S_i with S_i the
    % PRE-NORMALISATION pattern sum. Dividing by tau instead is only
    % asymptotically right - at tau = 0.058 ns, tau/dt is 0.34 times the true S,
    % so that shortcut was wrong by a factor of about 2.4 at short lifetimes.
    speciesTotal = speciesRaw1 + speciesRaw2;
    maps.speciesFraction1 = blank();
    maps.speciesFraction1(pixelIndex) = speciesRaw1 ./ max(speciesTotal, eps);
    maps.speciesFraction2 = blank();
    maps.speciesFraction2(pixelIndex) = speciesRaw2 ./ max(speciesTotal, eps);

    % Photon-weighted (intensity-weighted) mean lifetime, sum(a*tau)/sum(a).
    maps.tauMeanNs = blank();
    maps.tauMeanNs(pixelIndex) = ...
        (amp1 .* tau1Hat + amp2 .* tau2Hat) ./ max(photonTotal, eps);

    maps.background = blank(); maps.background(pixelIndex) = background;
    maps.backgroundFraction = blank();
    maps.backgroundFraction(pixelIndex) = ...
        background * nBin ./ max(intensity(pixelIndex), 1);
    maps.intensity = intensity;
    maps.reducedDeviance = blank();
    maps.reducedDeviance(pixelIndex) = deviance / max(nBin - 5, 1);
    maps.tau2Refined = false(imageSize);
    maps.tau2Refined(pixelIndex) = refined;
    maps.fittedMask = mask;

    out = struct('maps', maps, 'pixelIndex', pixelIndex, ...
        'tau1Ns', tau1Hat, 'tau2Ns', tau2Hat, 'amplitude1', amp1, ...
        'amplitude2', amp2, 'background', background, ...
        'deviance', deviance, 'tau1Index', tau1Index, ...
        'speciesRaw1', speciesRaw1, 'speciesRaw2', speciesRaw2, ...
        'dtNs', dtNs, 'periodNs', periodNs, 'nBin', nBin, 'irf', irf, ...
        'tau1Grid', tau1Grid, 'tau2Grid', tau2Grid, ...
        'shortlist', shortlist, ...
        'shortlistEdgeFraction', mean(atEdgeOfShortlist), ...
        'tau2RefinedFraction', mean(refined), ...
        'elapsedSeconds', toc(started));
    fprintf('  fitting done in %.0f s\n', out.elapsedSeconds);
    fprintf(['    tau2 refined off-grid for %.1f%% of pixel(s); the winning ' ...
        'pair sat at the\n    shortlist edge for %.2f%%\n'], ...
        100 * out.tau2RefinedFraction, 100 * out.shortlistEdgeFraction);
    if out.shortlistEdgeFraction > 0.02
        fprintf(['    WARNING: the shortlist may be too short, so the true ' ...
            'optimum could lie\n    outside it. Raise opts.shortlist and ' ...
            'compare.\n']);
    end
end
