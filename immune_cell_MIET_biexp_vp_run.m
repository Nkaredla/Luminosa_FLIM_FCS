function out = immune_cell_MIET_biexp_vp_run(cube, mask, pixelIndex, ...
        intensity, irf, dtNs, periodNs, nBin, opts)
%IMMUNE_CELL_MIET_BIEXP_VP_RUN Variable-projection fit of every selected pixel.
%
% Replaces the grid search. The outer BFGS searches [log tau1, log tau2] and the
% three non-negative amplitudes are profiled out by Poisson maximum likelihood at
% every trial point - the separable-nonlinear-least-squares structure the
% single-decay fitter always used, now applied per pixel and vectorised.
%
% WHY THIS REPLACED THE GRID
%
% The grid needed a ranking metric to choose among nodes, and the cheap metric it
% used was biased: it weighted residuals by max(model,1), which at 350-700
% photons over 156 bins discards the correct 1/model Poisson weight in the ~80%
% of bins holding under one count. Because tau2 was reported straight from the
% winning node, that ranking error WAS the answer, biasing tau2 upward by
% 5-7.4% in a way that shrank with photon count - which would have masqueraded
% as a real trend in any comparison across bin sizes.
%
% Variable projection has no ranking metric, no nodes and no shortlist, so none
% of those failure modes exist. Measured against the grid on the same pixels and
% the same penalised deviance, it reaches a strictly lower objective on 100% of
% pixels (median 0.145 lower) and the grid wins on none - and it does so at
% 2.8 ms/pixel against the grid's 6.2, because the shared circulant lets all
% pixels' patterns come from one matrix product and the unit-sum normalisation
% makes the inner EM pure elementwise work.
%
% It also returns CONTINUOUS tau1. The grid could only ever report one of five
% values spanning +/-2 sigma of the prior, which is most of why the tau1 maps
% looked like noise.

    nPixel = numel(pixelIndex);
    basis = biexp_slb_basis(irf, dtNs, periodNs, nBin);
    cubeFlat = reshape(cube, [], nBin);

    tau1Hat = nan(nPixel, 1);
    tau2Hat = nan(nPixel, 1);
    amp1 = nan(nPixel, 1);
    amp2 = nan(nPixel, 1);
    background = nan(nPixel, 1);
    deviance = nan(nPixel, 1);
    speciesRaw1 = nan(nPixel, 1);
    speciesRaw2 = nan(nPixel, 1);
    converged = false(nPixel, 1);
    gradNorm = nan(nPixel, 1);
    evaluations = nan(nPixel, 1);
    % Adequacy of the two-component model, per pixel. Reduced deviance alone
    % cannot answer this at a few hundred photons: with most bins holding 0-3
    % counts the deviance sits well below 1 whether the model is right or not,
    % which is why every fit so far has reported ~0.5 and told us nothing. A
    % STRUCTURED misfit instead shows up as correlation between neighbouring
    % residuals - a model that is the wrong shape misses several adjacent bins
    % the same way. Lag-1 residual autocorrelation is ~0 for a correct model at
    % any count level, so it is comparable across binnings.
    residAcf1 = nan(nPixel, 1);
    residMaxAbs = nan(nPixel, 1);

    blockCount = ceil(nPixel / opts.blockSize);
    fprintf(['  fitting %d pixel(s) in %d block(s) by variable projection ' ...
        '(BFGS + Poisson EM)\n'], nPixel, blockCount);
    started = tic;

    for blockIndex = 1:blockCount
        lo = (blockIndex - 1) * opts.blockSize + 1;
        hi = min(nPixel, blockIndex * opts.blockSize);
        rows = lo:hi;
        Y = double(cubeFlat(pixelIndex(rows), :))';

        vp = biexp_slb_bfgs_batch(Y, basis, opts);

        tau1Hat(rows) = vp.tau1Ns(:);
        tau2Hat(rows) = vp.tau2Ns(:);
        background(rows) = vp.beta(1, :)';
        amp1(rows) = vp.beta(2, :)';
        amp2(rows) = vp.beta(3, :)';
        deviance(rows) = vp.deviance(:);
        % Species (pre-exponential) weight is amplitude divided by the
        % PRE-normalisation pattern sum, not by tau: tau/dt is only the
        % asymptotic form and is 0.34x the true sum at tau = 0.058 ns.
        speciesRaw1(rows) = vp.beta(2, :)' ./ max(vp.patternSum1(:), eps);
        speciesRaw2(rows) = vp.beta(3, :)' ./ max(vp.patternSum2(:), eps);
        converged(rows) = vp.converged(:);
        gradNorm(rows) = vp.gradInfNorm(:);
        evaluations(rows) = vp.evaluations(:);

        % Rebuild the fitted model for this block and score the residuals.
        P1 = biexp_slb_pattern_batch(basis, vp.tau1Ns);
        P2 = biexp_slb_pattern_batch(basis, vp.tau2Ns);
        model = max(vp.beta(1, :) + vp.beta(2, :) .* P1 + ...
            vp.beta(3, :) .* P2, 1e-12);
        R = (Y - model) ./ sqrt(model);
        R = R - mean(R, 1);
        num = sum(R(1:end - 1, :) .* R(2:end, :), 1);
        den = sum(R .^ 2, 1);
        residAcf1(rows) = (num ./ max(den, eps))';
        residMaxAbs(rows) = max(abs(R), [], 1)';

        fprintf('    block %d/%d done (%.0f s elapsed, %.1f%% converged)\n', ...
            blockIndex, blockCount, toc(started), 100 * mean(vp.converged));
    end

    % ---- assemble maps ---------------------------------------------------
    imageSize = size(mask);
    blank = @() nan(imageSize);
    maps = struct();
    maps.tau1Ns = blank(); maps.tau1Ns(pixelIndex) = tau1Hat;
    maps.tau2Ns = blank(); maps.tau2Ns(pixelIndex) = tau2Hat;
    maps.slbPullSigma = blank();
    maps.slbPullSigma(pixelIndex) = ...
        (tau1Hat - opts.slbTauNs) / max(opts.slbSigmaNs, eps);
    % No grid, so no grid edge. Kept as an all-false map so every consumer of
    % these results keeps working unchanged.
    maps.tau1AtGridEdge = false(imageSize);
    maps.converged = false(imageSize);
    maps.converged(pixelIndex) = converged;

    photonTotal = amp1 + amp2;
    maps.photonFraction1 = blank();
    maps.photonFraction1(pixelIndex) = amp1 ./ max(photonTotal, eps);
    maps.photonFraction2 = blank();
    maps.photonFraction2(pixelIndex) = amp2 ./ max(photonTotal, eps);
    speciesTotal = speciesRaw1 + speciesRaw2;
    maps.speciesFraction1 = blank();
    maps.speciesFraction1(pixelIndex) = speciesRaw1 ./ max(speciesTotal, eps);
    maps.speciesFraction2 = blank();
    maps.speciesFraction2(pixelIndex) = speciesRaw2 ./ max(speciesTotal, eps);
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
    % tau2 is continuous here, so it is "refined" everywhere by construction.
    maps.tau2Refined = false(imageSize);
    maps.tau2Refined(pixelIndex) = true;
    maps.residualAcf1 = blank();
    maps.residualAcf1(pixelIndex) = residAcf1;
    maps.residualMaxAbs = blank();
    maps.residualMaxAbs(pixelIndex) = residMaxAbs;
    maps.fittedMask = mask;

    out = struct('maps', maps, 'pixelIndex', pixelIndex, ...
        'tau1Ns', tau1Hat, 'tau2Ns', tau2Hat, 'amplitude1', amp1, ...
        'amplitude2', amp2, 'background', background, ...
        'deviance', deviance, 'speciesRaw1', speciesRaw1, ...
        'speciesRaw2', speciesRaw2, 'converged', converged, ...
        'gradInfNorm', gradNorm, 'evaluations', evaluations, ...
        'residualAcf1', residAcf1, 'residualMaxAbs', residMaxAbs, ...
        'dtNs', dtNs, 'periodNs', periodNs, 'nBin', nBin, 'irf', irf, ...
        'method', 'variable-projection', ...
        'tau1Grid', [], 'tau2Grid', [], ...
        'convergedFraction', mean(converged), ...
        'tau2RefinedFraction', 1, ...
        'elapsedSeconds', toc(started));
    fprintf('  fitting done in %.0f s (%.2f ms/pixel)\n', ...
        out.elapsedSeconds, 1000 * out.elapsedSeconds / max(nPixel, 1));
    fprintf('    converged %.2f%% | objective evaluations median %.0f\n', ...
        100 * out.convergedFraction, median(evaluations));
    good = isfinite(residAcf1);
    fprintf(['    BIEXP ADEQUACY: lag-1 residual autocorrelation median ' ...
        '%+.4f (p90 %+.4f);\n                    max |residual| median ' ...
        '%.2f sigma, p90 %.2f\n'], median(residAcf1(good)), ...
        prctile(residAcf1(good), 90), median(residMaxAbs(good)), ...
        prctile(residMaxAbs(good), 90));
    if median(residAcf1(good)) > 0.05
        fprintf(['                    NOTE: positive median correlation ' ...
            'means the residuals are\n                    structured, so ' ...
            'two components are NOT describing the shape.\n']);
    end
    if out.convergedFraction < 0.99
        fprintf(['    WARNING: %.2f%% of pixels did not meet the gradient ' ...
            'tolerance.\n'], 100 * (1 - out.convergedFraction));
    end
end
