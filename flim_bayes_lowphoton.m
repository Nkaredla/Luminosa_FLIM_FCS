function out = flim_bayes_lowphoton(tcspc_pix, irf, pulsePeriodNs, dtNs, tau0, opts)
% FLIM_BAYES_LOWPHOTON
% Low-photon Bayesian FLIM map using an IRF-convolved mixture model.
%
% This adapts the low_photon Bayesian framework from
% bryankaye1/bayesian-analysis-of-fluorescent-lifetime-data to the
% Luminosa GUI workflow:
%   1) use the GUI-computed IRF on the current TCSPC grid
%   2) fit a global one- or two-state model to seed lifetimes
%   3) use PIRLS pattern matching for seed/background estimation
%   4) compute a per-pixel posterior mean lifetime on a low-photon grid
%
% INPUTS
%   tcspc_pix     : [nx x ny x t] or [nx x ny x t x nCh] TCSPC cube
%   irf           : [t x 1] IRF already sampled on the same time axis
%   pulsePeriodNs : hardware repetition period in ns (stored for reference)
%   dtNs          : TCSPC bin width in ns
%   tau0          : 1- or 2-state seed lifetimes in ns. If more are given,
%                   the shortest/longest pair is used.
%   opts fields
%       .useGPU             : true/false, default auto
%       .batchSize          : pixel batch size, default 2048
%       .includeBackground  : true/false, default true
%       .optimizeTau        : true/false, default true
%       .signalGrid         : signal-fraction grid, default linspace(0.05,1,25)
%       .fractionGrid       : species-1 fraction grid, default linspace(0,1,41)
%       .singleExpTauGrid   : explicit tau grid for 1-state mode
%       .shiftBounds        : integer IRF shift search bounds in bins, default [-5 5]
%
% OUTPUT
%   out struct with posterior mean lifetime map and diagnostics.

    if nargin < 6 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'useGPU') || isempty(opts.useGPU)
        opts.useGPU = (gpuDeviceCount > 0);
    end
    if ~isfield(opts, 'batchSize') || isempty(opts.batchSize)
        opts.batchSize = 2048;
    end
    if ~isfield(opts, 'includeBackground') || isempty(opts.includeBackground)
        opts.includeBackground = true;
    end
    if ~isfield(opts, 'optimizeTau') || isempty(opts.optimizeTau)
        opts.optimizeTau = true;
    end
    if ~isfield(opts, 'signalGrid') || isempty(opts.signalGrid)
        opts.signalGrid = linspace(0.0, 1.0, 26);
    end
    if ~isfield(opts, 'fractionGrid') || isempty(opts.fractionGrid)
        opts.fractionGrid = linspace(0.0, 1.0, 41);
    end
    if ~isfield(opts, 'singleExpTauGrid')
        opts.singleExpTauGrid = [];
    end
    if ~isfield(opts, 'shiftBounds') || isempty(opts.shiftBounds)
        opts.shiftBounds = [-5 5];
    end

    Ypix = single(tcspc_pix);
    if ndims(Ypix) == 4
        Ypix = sum(Ypix, 4);
    end
    if ndims(Ypix) ~= 3
        error('tcspc_pix must be a 3D or 4D TCSPC cube.');
    end

    [nx, ny, nt] = size(Ypix);
    modelPeriodNs = nt * dtNs;

    irf = max(double(irf(:)), 0);
    if numel(irf) ~= nt
        error('IRF length (%d) must match the TCSPC time axis (%d).', numel(irf), nt);
    end
    if any(irf > 0)
        irf = irf ./ max(sum(irf), eps);
    else
        irf = zeros(nt, 1);
        irf(1) = 1;
    end

    tau0Requested = sort(double(tau0(:)).', 'ascend');
    if isempty(tau0Requested)
        error('tau0 must contain at least one lifetime guess.');
    end
    if numel(tau0Requested) > 2
        tau0Use = [tau0Requested(1), tau0Requested(end)];
    else
        tau0Use = tau0Requested;
    end

    globalDecay = squeeze(sum(sum(Ypix, 1), 2));
    globalDecay = double(globalDecay(:));

    seedFit = fit_global_bayes_seed(globalDecay, irf, modelPeriodNs, dtNs, tau0Use, opts);

    [AmpSeed, ~, pmSeedInfo] = PatternMatchIm_matlab(Ypix, seedFit.modelMatrix, 'PIRLS', opts.useGPU, opts.batchSize);
    if opts.includeBackground
        ampSeedSum = sum(AmpSeed, 3);
        bgFracSeed = AmpSeed(:, :, 1) ./ max(ampSeedSum, eps('single'));
    else
        bgFracSeed = zeros(nx, ny, 'single');
    end

    if numel(seedFit.tauFit) == 1
        posterior = bayes_single_state_map(Ypix, seedFit, tau0Requested, opts);
    else
        posterior = bayes_two_state_map(Ypix, seedFit, opts);
    end

    out = struct();
    out.globalDecay = globalDecay;
    out.irf = irf;
    out.irfShifted = seedFit.irfShifted;
    out.pulsePeriodNs = pulsePeriodNs;
    out.modelPeriodNs = modelPeriodNs;
    out.dtNs = dtNs;
    out.tau0Requested = tau0Requested(:).';
    out.tau0Used = tau0Use(:).';
    out.includeBackground = logical(opts.includeBackground);
    out.optimizeTau = logical(opts.optimizeTau);

    out.globalFit = seedFit;
    out.seedAmp = AmpSeed;
    out.seedBackgroundFraction = bgFracSeed;
    out.seedPatternMatchInfo = pmSeedInfo;

    out.tauMeanArithmetic = posterior.tauMeanArithmetic;
    out.tauPosteriorStd = posterior.tauPosteriorStd;
    out.tauMap = posterior.tauMap;
    out.intensity = sum(Ypix, 3);
    out.signalFractionMean = posterior.signalFractionMean;
    out.stateFractionMean = posterior.stateFractionMean;
    out.posteriorInfo = posterior.posteriorInfo;
end

function seedFit = fit_global_bayes_seed(globalDecay, irf, modelPeriodNs, dtNs, tau0, opts)
    tau0 = double(tau0(:));
    seedMat = build_tau_seed_matrix(tau0);
    shiftMin = floor(opts.shiftBounds(1));
    shiftMax = ceil(opts.shiftBounds(end));
    shiftGrid = shiftMin:shiftMax;
    if isempty(shiftGrid)
        shiftGrid = 0;
    end

    best.err = inf;
    best.tauFit = tau0(:).';
    best.coeff = [];
    best.fitCounts = [];
    best.modelMatrix = [];
    best.stats = struct();
    best.shiftBins = 0;
    best.irfShifted = irf(:);

    for shiftBins = shiftGrid
        irfShifted = local_circshift(irf(:), shiftBins);
        for seedIdx = 1:size(seedMat, 2)
            tauSeed = seedMat(:, seedIdx);
            if opts.optimizeTau
                p0 = tauSeed;
                xmin = max(0.03, p0 / 10);
                xmax = max(p0 * 10, p0 + 0.05);
                tol = 1e-5;
                steps = max(250, 180 * numel(p0));
                [tauCand, ~] = Simplex(@global_bayes_err, p0, xmin, xmax, tol, steps, ...
                    globalDecay, modelPeriodNs, dtNs, irfShifted, opts.includeBackground);
            else
                tauCand = tauSeed;
            end

            tauCand = sort(double(tauCand(:)), 'ascend');
            [err, coeff, fitCounts, modelMatrix] = global_bayes_err( ...
                tauCand, globalDecay, modelPeriodNs, dtNs, irfShifted, opts.includeBackground);
            stats = calc_bayes_fit_stats(globalDecay, fitCounts, numel(tauCand), opts.includeBackground);
            if err < best.err
                best.err = err;
                best.tauFit = tauCand(:).';
                best.coeff = coeff(:).';
                best.fitCounts = fitCounts(:);
                best.modelMatrix = modelMatrix;
                best.stats = stats;
                best.shiftBins = shiftBins;
                best.irfShifted = irfShifted(:);
            end
        end
    end

    seedFit = best;
    seedFit.dtNs = dtNs;
    seedFit.modelPeriodNs = modelPeriodNs;
end

function [err, coeff, fitCounts, modelMatrix] = global_bayes_err(tau, counts, modelPeriodNs, dtNs, irf, includeBackground)
    tau = max(double(tau(:)), 1e-6);
    counts = max(double(counts(:)), 0);
    irf = max(double(irf(:)), 0);
    irf = irf ./ max(sum(irf), eps);

    modelMatrix = build_bayes_model_matrix(irf, modelPeriodNs, dtNs, tau, includeBackground);
    coeff = solve_bayes_coefficients_pirls(modelMatrix, counts);
    fitCounts = modelMatrix * coeff;
    err = sum((counts - fitCounts).^2 ./ max(fitCounts, 1));
end

function M = build_bayes_model_matrix(irf, modelPeriodNs, dtNs, tau, includeBackground)
    irf = max(double(irf(:)), 0);
    irf = irf ./ max(sum(irf), eps);
    tau = max(double(tau(:)), 1e-6);

    n = numel(irf);
    tNs = ((0:n-1)') * dtNs;
    periodNs = max(modelPeriodNs, n * dtNs);
    M = zeros(n, numel(tau) + double(includeBackground));

    if includeBackground
        M(:, 1) = 1;
    end

    for k = 1:numel(tau)
        decay = exp(-tNs ./ tau(k)) ./ max(1 - exp(-periodNs / tau(k)), eps);
        convSig = Convol(irf, decay);
        M(:, double(includeBackground) + k) = max(double(convSig(:)), 0);
    end

    colScale = sum(M, 1);
    colScale(~isfinite(colScale) | colScale <= 0) = 1;
    M = M ./ colScale;
end

function coeff = solve_bayes_coefficients_pirls(M, counts)
    counts = max(double(counts(:)), 0);
    M = max(double(M), 0);
    try
        coeff = double(PIRLSnonneg(M, counts, 10));
        coeff = coeff(:);
        return;
    catch
    end

    try
        [beta, ~] = PIRLSnonneg_batch_gpu_matlab(M, counts, 10, 25, 1, false);
        coeff = double(beta(:, 1));
        coeff = coeff(:);
        return;
    catch
    end

    coeff = lsqnonneg(M, counts);
    coeff = coeff(:);
end

function seedMat = build_tau_seed_matrix(tau0)
    tau0 = double(tau0(:));
    scales = [1.0, 0.75, 1.35];
    seedMat = zeros(numel(tau0), numel(scales));
    for idx = 1:numel(scales)
        seedMat(:, idx) = max(tau0 * scales(idx), 0.03);
    end
end

function stats = calc_bayes_fit_stats(counts, fitCounts, nExp, includeBackground)
    counts = max(double(counts(:)), 0);
    fitCounts = max(double(fitCounts(:)), eps);
    n = numel(counts);
    k = nExp + double(includeBackground) + 1;
    chi2 = sum((counts - fitCounts).^2 ./ max(fitCounts, 1));
    ndf = max(n - k, 1);
    stats = struct();
    stats.chi2 = chi2;
    stats.chi2red = chi2 / ndf;
    stats.n = n;
    stats.k = k;
    stats.includeBackground = logical(includeBackground);
end

function posterior = bayes_two_state_map(Ypix, seedFit, opts)
    M = max(double(seedFit.modelMatrix), 0);
    tauFit = double(seedFit.tauFit(:)).';
    if opts.includeBackground
        bg = M(:, 1);
        p1 = M(:, 2);
        p2 = M(:, 3);
    else
        bg = [];
        p1 = M(:, 1);
        p2 = M(:, 2);
    end

    signalGrid = double(opts.signalGrid(:));
    fractionGrid = double(opts.fractionGrid(:));
    if ~opts.includeBackground
        signalGrid = 1;
    end

    [Pgrid, tauPerGrid, signalPerGrid, fracPerGrid] = build_two_state_grid( ...
        bg, p1, p2, tauFit(1), tauFit(2), signalGrid, fractionGrid, opts.includeBackground);

    [tauMean, tauStd, tauMap, signalMean, fracMean, postInfo] = evaluate_posterior_grid( ...
        Ypix, Pgrid, tauPerGrid, signalPerGrid, fracPerGrid, opts.useGPU, opts.batchSize);

    posterior = struct();
    posterior.tauMeanArithmetic = tauMean;
    posterior.tauPosteriorStd = tauStd;
    posterior.tauMap = tauMap;
    posterior.signalFractionMean = signalMean;
    posterior.stateFractionMean = fracMean;
    posterior.posteriorInfo = postInfo;
    posterior.posteriorInfo.method = 'Bayesian low-photon 2-state posterior';
    posterior.posteriorInfo.signalGrid = signalGrid(:).';
    posterior.posteriorInfo.fractionGrid = fractionGrid(:).';
end

function posterior = bayes_single_state_map(Ypix, seedFit, tau0Requested, opts)
    tauFit = double(seedFit.tauFit(1));
    if isempty(opts.singleExpTauGrid)
        tauMin = max(0.05, min([0.5 * tauFit, 0.5 * min(tau0Requested)]));
        tauMax = max([1.5 * tauFit, 1.5 * max(tau0Requested), tauMin + 0.05]);
        tauGrid = linspace(tauMin, tauMax, 41);
    else
        tauGrid = double(opts.singleExpTauGrid(:)).';
    end

    if opts.includeBackground
        bg = double(seedFit.modelMatrix(:, 1));
    else
        bg = [];
    end
    signalGrid = double(opts.signalGrid(:));
    if ~opts.includeBackground
        signalGrid = 1;
    end

    [Pgrid, tauPerGrid, signalPerGrid] = build_single_state_grid( ...
        seedFit.irfShifted, bg, seedFit.modelPeriodNs, seedFit.dtNs, tauGrid, signalGrid, opts.includeBackground);

    [tauMean, tauStd, tauMap, signalMean, ~, postInfo] = evaluate_posterior_grid( ...
        Ypix, Pgrid, tauPerGrid, signalPerGrid, zeros(size(signalPerGrid)), opts.useGPU, opts.batchSize);

    posterior = struct();
    posterior.tauMeanArithmetic = tauMean;
    posterior.tauPosteriorStd = tauStd;
    posterior.tauMap = tauMap;
    posterior.signalFractionMean = signalMean;
    posterior.stateFractionMean = [];
    posterior.posteriorInfo = postInfo;
    posterior.posteriorInfo.method = 'Bayesian low-photon 1-state posterior';
    posterior.posteriorInfo.signalGrid = signalGrid(:).';
    posterior.posteriorInfo.tauGrid = tauGrid(:).';
end

function [Pgrid, tauPerGrid, signalPerGrid, fracPerGrid] = build_two_state_grid(bg, p1, p2, tau1, tau2, signalGrid, fractionGrid, includeBackground)
    signalGrid = double(signalGrid(:));
    fractionGrid = double(fractionGrid(:));
    sigPatterns = p2 + (p1 - p2) * reshape(fractionGrid.', 1, []);

    if includeBackground
        Pgrid = zeros(numel(p1), numel(signalGrid) * numel(fractionGrid));
        tauPerGrid = zeros(1, size(Pgrid, 2));
        signalPerGrid = zeros(1, size(Pgrid, 2));
        fracPerGrid = zeros(1, size(Pgrid, 2));
        col = 1;
        for sIdx = 1:numel(signalGrid)
            s = signalGrid(sIdx);
            cols = col:(col + numel(fractionGrid) - 1);
            Pgrid(:, cols) = s * sigPatterns + (1 - s) * bg;
            tauPerGrid(cols) = fractionGrid * tau1 + (1 - fractionGrid) * tau2;
            signalPerGrid(cols) = s;
            fracPerGrid(cols) = fractionGrid;
            col = col + numel(fractionGrid);
        end
    else
        Pgrid = sigPatterns;
        tauPerGrid = fractionGrid * tau1 + (1 - fractionGrid) * tau2;
        signalPerGrid = ones(size(tauPerGrid));
        fracPerGrid = fractionGrid;
    end

    colScale = sum(Pgrid, 1);
    colScale(colScale <= 0) = 1;
    Pgrid = Pgrid ./ colScale;
end

function [Pgrid, tauPerGrid, signalPerGrid] = build_single_state_grid(irfShifted, bg, modelPeriodNs, dtNs, tauGrid, signalGrid, includeBackground)
    tauGrid = double(tauGrid(:));
    signalGrid = double(signalGrid(:));
    n = numel(irfShifted);

    sigPatterns = zeros(n, numel(tauGrid));
    for idx = 1:numel(tauGrid)
        Mtmp = build_bayes_model_matrix(irfShifted, modelPeriodNs, dtNs, tauGrid(idx), false);
        sigPatterns(:, idx) = Mtmp(:, 1);
    end

    if includeBackground
        Pgrid = zeros(n, numel(signalGrid) * numel(tauGrid));
        tauPerGrid = zeros(1, size(Pgrid, 2));
        signalPerGrid = zeros(1, size(Pgrid, 2));
        col = 1;
        for sIdx = 1:numel(signalGrid)
            s = signalGrid(sIdx);
            cols = col:(col + numel(tauGrid) - 1);
            Pgrid(:, cols) = s * sigPatterns + (1 - s) * bg;
            tauPerGrid(cols) = tauGrid;
            signalPerGrid(cols) = s;
            col = col + numel(tauGrid);
        end
    else
        Pgrid = sigPatterns;
        tauPerGrid = tauGrid(:).';
        signalPerGrid = ones(size(tauPerGrid));
    end

    colScale = sum(Pgrid, 1);
    colScale(colScale <= 0) = 1;
    Pgrid = Pgrid ./ colScale;
end

function [tauMeanMap, tauStdMap, tauMapMap, signalMeanMap, fracMeanMap, info] = ...
        evaluate_posterior_grid(Ypix, Pgrid, tauPerGrid, signalPerGrid, fracPerGrid, useGPU, batchSize)
    [nx, ny, nt] = size(Ypix);
    Y = reshape(permute(single(Ypix), [3 1 2]), nt, []);
    nPix = size(Y, 2);
    validPix = sum(Y, 1) > 0;

    tauMean = nan(1, nPix, 'single');
    tauStd = nan(1, nPix, 'single');
    tauMap = nan(1, nPix, 'single');
    signalMean = nan(1, nPix, 'single');
    fracMean = nan(1, nPix, 'single');

    tauPerGrid = double(tauPerGrid(:));
    tauSqPerGrid = tauPerGrid .^ 2;
    signalPerGrid = double(signalPerGrid(:));
    fracPerGrid = double(fracPerGrid(:));

    logP = max(double(Pgrid), realmin('double'));
    logP = log(logP);

    gpuUsed = false;
    if useGPU
        try
            logPg = gpuArray(single(logP));
            gpuUsed = true;
        catch
            gpuUsed = false;
        end
    end

    validIdx = find(validPix);
    for i0 = 1:batchSize:numel(validIdx)
        idx = validIdx(i0:min(i0 + batchSize - 1, numel(validIdx)));
        Ychunk = Y(:, idx);

        if gpuUsed
            Yg = gpuArray(Ychunk);
            ll = gather(logPg.' * Yg);
        else
            ll = logP.' * double(Ychunk);
        end

        ll = double(ll);
        ll = ll - max(ll, [], 1);
        w = exp(ll);
        wsum = sum(w, 1);
        wsum(wsum <= 0) = 1;
        w = w ./ wsum;

        tauMean(idx) = single(tauPerGrid.' * w);
        tauVar = tauSqPerGrid.' * w - double(tauMean(idx)) .^ 2;
        tauStd(idx) = single(sqrt(max(tauVar, 0)));
        signalMean(idx) = single(signalPerGrid.' * w);
        if ~isempty(fracPerGrid)
            fracMean(idx) = single(fracPerGrid.' * w);
        end

        [~, mapIdx] = max(w, [], 1);
        tauMap(idx) = single(tauPerGrid(mapIdx));
    end

    tauMeanMap = reshape(tauMean, nx, ny);
    tauStdMap = reshape(tauStd, nx, ny);
    tauMapMap = reshape(tauMap, nx, ny);
    signalMeanMap = reshape(signalMean, nx, ny);
    if isempty(fracPerGrid) || all(fracPerGrid == 0)
        fracMeanMap = [];
    else
        fracMeanMap = reshape(fracMean, nx, ny);
    end

    info = struct();
    info.usedGPU = gpuUsed;
    info.nGrid = size(Pgrid, 2);
    info.batchSize = batchSize;
end

function y = local_circshift(x, shiftBins)
    x = x(:);
    n = numel(x);
    if n == 0
        y = x;
        return;
    end
    shiftBins = round(shiftBins);
    idx = mod((1:n) + n - shiftBins - 1, n) + 1;
    y = x(idx);
end
