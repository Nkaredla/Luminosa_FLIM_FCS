function out = immune_cell_MIET_biexp_run(cube, mask, pixelIndex, intensity, ...
        irf, dtNs, periodNs, nBin, tau1Grid, tau2Grid, opts)
%IMMUNE_CELL_MIET_BIEXP_RUN Fit every selected pixel; return maps.
%
% Two stages, for the reason given in immune_cell_MIET_biexp_slb: a vectorised
% grid ranks the lifetime pairs cheaply, then PIRLSnonneg produces the reported
% amplitudes and background once per pixel at its winning pair.
%
% ON THE RANKING METRIC
%
% The grid ranks by chi-square plus the Gaussian prior penalty, not by the
% Poisson deviance, purely because ranking 200 pairs times 100k pixels with the
% deviance would mean about 3e9 logarithms. Chi-square and deviance are both
% twice a negative log-likelihood up to terms that do not depend on the model, so
% they order candidate pairs almost identically, and the pair that wins is then
% refitted with the exact Poisson-weighted solver. The tau2 grid spacing is
% roughly 5%, far finer than the lifetime uncertainty at these photon counts, so
% grid resolution is not the limiting error.
%
% Amplitudes from the grid stage are clamped at zero before the residual is
% formed, so a pair is never ranked using a model with negative photon counts.

    nPixel = numel(pixelIndex);
    n1 = numel(tau1Grid);
    n2 = numel(tau2Grid);

    % Patterns depend only on the lifetime, so build them once.
    pattern1 = zeros(nBin, n1);
    for a = 1:n1
        pattern1(:, a) = biexp_slb_pattern(irf, dtNs, periodNs, tau1Grid(a), nBin);
    end
    pattern2 = zeros(nBin, n2);
    for b = 1:n2
        pattern2(:, b) = biexp_slb_pattern(irf, dtNs, periodNs, tau2Grid(b), nBin);
    end
    penalty1 = ((tau1Grid - opts.slbTauNs) / max(opts.slbSigmaNs, eps)) .^ 2;

    tau1Hat = nan(nPixel, 1);
    tau2Hat = nan(nPixel, 1);
    amp1 = nan(nPixel, 1);
    amp2 = nan(nPixel, 1);
    background = nan(nPixel, 1);
    deviance = nan(nPixel, 1);

    ones156 = ones(nBin, 1);
    cubeFlat = reshape(cube, [], nBin);
    blockCount = ceil(nPixel / opts.blockSize);
    fprintf('  fitting %d pixel(s) in %d block(s)\n', nPixel, blockCount);
    started = tic;

    for blockIndex = 1:blockCount
        lo = (blockIndex - 1) * opts.blockSize + 1;
        hi = min(nPixel, blockIndex * opts.blockSize);
        idx = pixelIndex(lo:hi);
        nThis = numel(idx);

        % One row-index into the flattened cube, rather than pulling 156
        % full 602x602 slices per block.
        Y = double(cubeFlat(idx, :))';

        bestScore = inf(1, nThis);
        best1 = ones(1, nThis);
        best2 = ones(1, nThis);
        for a = 1:n1
            for b = 1:n2
                if tau2Grid(b) <= tau1Grid(a); continue; end
                design = [ones156, pattern1(:, a), pattern2(:, b)];
                % One small solve for the whole block: the design is shared.
                beta = max((design' * design) \ (design' * Y), 0);
                model = design * beta;
                model(model < 1e-9) = 1e-9;
                score = sum((Y - model) .^ 2 ./ max(model, 1), 1) + penalty1(a);
                better = score < bestScore;
                if any(better)
                    bestScore(better) = score(better);
                    best1(better) = a;
                    best2(better) = b;
                end
            end
        end

        % Poisson-weighted, non-negative refit at each pixel's winning pair.
        for p = 1:nThis
            design = [ones156, pattern1(:, best1(p)), pattern2(:, best2(p))];
            y = Y(:, p);
            try
                beta = PIRLSnonneg(design, y);
            catch
                beta = max(design \ y, 0);
            end
            model = max(design * beta, 1e-12);
            row = lo + p - 1;
            tau1Hat(row) = tau1Grid(best1(p));
            tau2Hat(row) = tau2Grid(best2(p));
            background(row) = beta(1);
            amp1(row) = beta(2);
            amp2(row) = beta(3);
            deviance(row) = biexp_slb_deviance(y, model);
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
    maps.slbPullSigma = blank();
    maps.slbPullSigma(pixelIndex) = ...
        (tau1Hat - opts.slbTauNs) / max(opts.slbSigmaNs, eps);
    total = amp1 + amp2;
    maps.speciesFraction1 = blank();
    maps.speciesFraction1(pixelIndex) = amp1 ./ max(total, eps);
    maps.speciesFraction2 = blank();
    maps.speciesFraction2(pixelIndex) = amp2 ./ max(total, eps);
    % Photon share weights amplitude by lifetime: what each population actually
    % contributes to the measured brightness.
    photon1 = amp1 .* tau1Hat;
    photon2 = amp2 .* tau2Hat;
    photonTotal = photon1 + photon2;
    maps.photonFraction1 = blank();
    maps.photonFraction1(pixelIndex) = photon1 ./ max(photonTotal, eps);
    maps.photonFraction2 = blank();
    maps.photonFraction2(pixelIndex) = photon2 ./ max(photonTotal, eps);
    maps.tauMeanNs = blank();
    maps.tauMeanNs(pixelIndex) = ...
        (photon1 .* tau1Hat + photon2 .* tau2Hat) ./ max(photonTotal, eps);
    maps.background = blank(); maps.background(pixelIndex) = background;
    maps.backgroundFraction = blank();
    maps.backgroundFraction(pixelIndex) = ...
        background * nBin ./ max(intensity(pixelIndex), 1);
    maps.intensity = intensity;
    maps.reducedDeviance = blank();
    maps.reducedDeviance(pixelIndex) = deviance / max(nBin - 5, 1);
    maps.fittedMask = mask;

    out = struct('maps', maps, 'pixelIndex', pixelIndex, ...
        'tau1Ns', tau1Hat, 'tau2Ns', tau2Hat, 'amplitude1', amp1, ...
        'amplitude2', amp2, 'background', background, ...
        'deviance', deviance, 'dtNs', dtNs, 'periodNs', periodNs, ...
        'nBin', nBin, 'irf', irf, 'tau1Grid', tau1Grid, ...
        'tau2Grid', tau2Grid, 'elapsedSeconds', toc(started));
    fprintf('  fitting done in %.0f s\n', out.elapsedSeconds);
end
