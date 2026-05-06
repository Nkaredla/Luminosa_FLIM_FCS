function results = simulate_flim_bayes_monoexp_error(opts)
% SIMULATE_FLIM_BAYES_MONOEXP_ERROR
% Monte Carlo benchmark for FLIM_BAYES_LOWPHOTON on monoexponential TCSPC.
%
% The function generates IRF-convolved monoexponential decays, draws noisy
% TCSPC histograms with exactly N photons per replicate, fits them with
% FLIM_BAYES_LOWPHOTON in 1-state mode, and summarizes the fitting error
% across a grid of true lifetimes and photon counts.
%
% Example
%   results = simulate_flim_bayes_monoexp_error();
%
%   opts = struct();
%   opts.trueTauNs = 2.8;
%   opts.tau0Ns = 2.0;
%   opts.nReplicates = 500;
%   opts.useGPU = false;
%   results = simulate_flim_bayes_monoexp_error(opts);
%
% INPUT
%   opts fields (all optional)
%       .trueTauNs         scalar or vector of true monoexponential lifetimes [ns]
%       .tau0Ns            scalar or vector lifetime seed(s) passed to flim_bayes [ns]
%       .photonCounts      vector of photons per decay
%       .nReplicates       number of Monte Carlo repeats per photon count
%       .nTimeBins         number of TCSPC bins
%       .dtNs              TCSPC bin width [ns]
%       .pulsePeriodNs     laser repetition period [ns]
%       .backgroundFraction fraction of uniform background photons [0..1]
%       .irf               explicit IRF sampled on the TCSPC grid
%       .irfCenterBin      Gaussian IRF center bin if .irf is omitted
%       .irfSigmaBins      Gaussian IRF sigma in bins if .irf is omitted
%       .includeBackground pass background component to flim_bayes
%       .optimizeTau       pass through to flim_bayes (default true)
%       .signalGrid        pass through to flim_bayes
%       .singleExpTauGrid  pass through to flim_bayes
%       .shiftBounds       pass through to flim_bayes
%       .batchSize         pass through to flim_bayes
%       .useGPU            pass through to flim_bayes
%       .showPlot          draw summary figure (default true)
%       .rngSeed           RNG seed for reproducibility (default 1)
%       .calibrationMatFile MIET calibration MAT file, default calibrationCurve.mat
%       .curvePlotTauIdx   selected tau indices for curve subplots
%       .curvePlotPhotonIdx selected photon-count indices for curve subplots
%       .curvePlotMaxTau   default representative tau count if indices omitted
%       .curvePlotMaxPhoton default representative photon count if indices omitted
%
% OUTPUT
%   results struct with fields:
%       .summaryTable      table of bias / spread / RMSE by (tau, photons)
%       .replicates        [nTau x nPhot] replicate results
%       .summaryMaps       2D maps over (tau, photons)
%       .uncertaintyMapNs  empirical std(error) map over (tau, photons)
%       .trueDecayProb     [t x nTau] noiseless normalized decay probabilities
%       .irf               normalized IRF used for simulation
%       .calibration       loaded calibration struct
%       .optsUsed          resolved option struct

    if nargin < 1 || isempty(opts)
        opts = struct();
    end
    opts = fill_default_options(opts);

    rng(opts.rngSeed, 'twister');

    nt = opts.nTimeBins;
    dtNs = opts.dtNs;
    pulsePeriodNs = opts.pulsePeriodNs;
    tauTrueNsList = opts.trueTauNs(:).';
    tau0NsList = opts.tau0Ns(:).';

    irf = build_irf(opts);
    calib = load_lifetime_height_calibration(opts.calibrationMatFile);
    bgPattern = ones(nt, 1) ./ nt;

    fitOpts = struct();
    fitOpts.useGPU = logical(opts.useGPU);
    fitOpts.batchSize = opts.batchSize;
    fitOpts.includeBackground = logical(opts.includeBackground);
    fitOpts.optimizeTau = logical(opts.optimizeTau);
    fitOpts.signalGrid = opts.signalGrid;
    fitOpts.fractionGrid = 0;
    fitOpts.singleExpTauGrid = opts.singleExpTauGrid;
    fitOpts.shiftBounds = opts.shiftBounds;

    nTau = numel(tauTrueNsList);
    nPhot = numel(opts.photonCounts);
    replicates = repmat(struct( ...
        'trueTauNs', [], ...
        'tau0Ns', [], ...
        'photonCount', [], ...
        'estimatedTauNs', [], ...
        'estimatedHeightNm', [], ...
        'posteriorStdNs', [], ...
        'errorNs', [], ...
        'errorHeightNm', [], ...
        'meanDataCurve', [], ...
        'stdDataCurve', [], ...
        'meanFitCurve', [], ...
        'stdFitCurve', [], ...
        'summary', [], ...
        'fitOutput', []), nTau, nPhot);

    trueDecayProb = zeros(nt, nTau);
    trueHeightNm = nan(1, nTau);
    meanTauMap = nan(nTau, nPhot);
    biasMap = nan(nTau, nPhot);
    stdErrorMap = nan(nTau, nPhot);
    rmseMap = nan(nTau, nPhot);
    maeMap = nan(nTau, nPhot);
    medianAbsErrorMap = nan(nTau, nPhot);
    meanPosteriorStdMap = nan(nTau, nPhot);
    q05Map = nan(nTau, nPhot);
    q50Map = nan(nTau, nPhot);
    q95Map = nan(nTau, nPhot);
    meanHeightMap = nan(nTau, nPhot);
    biasHeightMap = nan(nTau, nPhot);
    stdErrorHeightMap = nan(nTau, nPhot);
    rmseHeightMap = nan(nTau, nPhot);

    summaryRows = cell(nTau * nPhot, 1);
    rowIdx = 1;
    for tauIdx = 1:nTau
        tauTrueNs = tauTrueNsList(tauIdx);
        tau0Ns = tau0NsList(tauIdx);
        decayProb = build_monoexp_probability(irf, tauTrueNs, pulsePeriodNs, dtNs, opts.backgroundFraction);
        trueDecayProb(:, tauIdx) = decayProb(:);
        trueHeightHere = lifetime_to_height_from_calibration(tauTrueNs, calib);
        trueHeightNm(tauIdx) = trueHeightHere;

        for photIdx = 1:nPhot
            photonCount = opts.photonCounts(photIdx);
            counts = simulate_exact_photon_histograms(decayProb, photonCount, opts.nReplicates);
            tcspc_pix = permute(reshape(counts, nt, opts.nReplicates, 1), [2 3 1]);

            out = flim_bayes_lowphoton(tcspc_pix, irf, pulsePeriodNs, dtNs, tau0Ns, fitOpts);
            estTauNs = double(out.tauMeanArithmetic(:));
            postStdNs = double(out.tauPosteriorStd(:));
            errNs = estTauNs - tauTrueNs;
            estHeightNm = lifetime_to_height_from_calibration(estTauNs, calib);

            totalCounts = sum(double(counts), 1);
            signalFrac = double(out.signalFractionMean(:));
            fitCounts = reconstruct_monoexp_fit_counts( ...
                estTauNs, signalFrac, totalCounts, irf, pulsePeriodNs, dtNs, bgPattern);

            row = summarize_errors(photonCount, tauTrueNs, tau0Ns, estTauNs, postStdNs);
            row = append_height_summary(row, trueHeightHere, estHeightNm);

            replicates(tauIdx, photIdx).trueTauNs = tauTrueNs;
            replicates(tauIdx, photIdx).tau0Ns = tau0Ns;
            replicates(tauIdx, photIdx).photonCount = photonCount;
            replicates(tauIdx, photIdx).estimatedTauNs = estTauNs;
            replicates(tauIdx, photIdx).estimatedHeightNm = estHeightNm;
            replicates(tauIdx, photIdx).posteriorStdNs = postStdNs;
            replicates(tauIdx, photIdx).errorNs = errNs;
            replicates(tauIdx, photIdx).errorHeightNm = estHeightNm - trueHeightHere;
            replicates(tauIdx, photIdx).meanDataCurve = mean(double(counts), 2);
            replicates(tauIdx, photIdx).stdDataCurve = std(double(counts), 0, 2);
            replicates(tauIdx, photIdx).meanFitCurve = mean(fitCounts, 2);
            replicates(tauIdx, photIdx).stdFitCurve = std(fitCounts, 0, 2);
            replicates(tauIdx, photIdx).summary = row;
            replicates(tauIdx, photIdx).fitOutput = out;
            summaryRows{rowIdx} = row;

            meanTauMap(tauIdx, photIdx) = row.MeanEstimatedTauNs;
            biasMap(tauIdx, photIdx) = row.BiasNs;
            stdErrorMap(tauIdx, photIdx) = row.StdErrorNs;
            rmseMap(tauIdx, photIdx) = row.RMSE_Ns;
            maeMap(tauIdx, photIdx) = row.MAE_Ns;
            medianAbsErrorMap(tauIdx, photIdx) = row.MedianAbsErrorNs;
            meanPosteriorStdMap(tauIdx, photIdx) = row.MeanPosteriorStdNs;
            q05Map(tauIdx, photIdx) = row.Q05EstimatedTauNs;
            q50Map(tauIdx, photIdx) = row.Q50EstimatedTauNs;
            q95Map(tauIdx, photIdx) = row.Q95EstimatedTauNs;
            meanHeightMap(tauIdx, photIdx) = row.MeanEstimatedHeightNm;
            biasHeightMap(tauIdx, photIdx) = row.BiasHeightNm;
            stdErrorHeightMap(tauIdx, photIdx) = row.StdErrorHeightNm;
            rmseHeightMap(tauIdx, photIdx) = row.RMSE_HeightNm;

            rowIdx = rowIdx + 1;
        end
    end

    summaryTable = vertcat(summaryRows{:});

    results = struct();
    results.summaryTable = summaryTable;
    results.replicates = replicates;
    results.trueTauNs = tauTrueNsList;
    results.trueHeightNm = trueHeightNm(:);
    results.tau0Ns = tau0NsList;
    results.tAxisNs = ((0:nt-1)') * dtNs;
    results.trueDecayProb = trueDecayProb;
    results.irf = irf(:);
    results.calibration = calib;
    results.summaryMaps = struct( ...
        'MeanEstimatedTauNs', meanTauMap, ...
        'BiasNs', biasMap, ...
        'StdErrorNs', stdErrorMap, ...
        'RMSE_Ns', rmseMap, ...
        'MAE_Ns', maeMap, ...
        'MedianAbsErrorNs', medianAbsErrorMap, ...
        'MeanPosteriorStdNs', meanPosteriorStdMap, ...
        'Q05EstimatedTauNs', q05Map, ...
        'Q50EstimatedTauNs', q50Map, ...
        'Q95EstimatedTauNs', q95Map, ...
        'MeanEstimatedHeightNm', meanHeightMap, ...
        'BiasHeightNm', biasHeightMap, ...
        'StdErrorHeightNm', stdErrorHeightMap, ...
        'RMSE_HeightNm', rmseHeightMap, ...
        'tauAxisNs', tauTrueNsList(:), ...
        'heightAxisNm', trueHeightNm(:), ...
        'photonAxis', opts.photonCounts(:).');
    results.uncertaintyMapNs = stdErrorMap;
    results.posteriorUncertaintyMapNs = meanPosteriorStdMap;
    results.rmseMapNs = rmseMap;
    results.biasMapNs = biasMap;
    results.lifetimeDeviationMapNs = biasMap;
    results.heightBiasMapNm = biasHeightMap;
    results.heightUncertaintyMapNm = stdErrorHeightMap;
    results.heightRmseMapNm = rmseHeightMap;
    results.optsUsed = opts;

    disp(summaryTable);

    if opts.showPlot
        plot_summary(results);
    end
end

function opts = fill_default_options(opts)
    if ~isfield(opts, 'trueTauNs') || isempty(opts.trueTauNs)
        opts.trueTauNs = 0.5:0.5:5.0;
    end
    if ~isfield(opts, 'tau0Ns') || isempty(opts.tau0Ns)
        opts.tau0Ns = [];
    end
    if ~isfield(opts, 'photonCounts') || isempty(opts.photonCounts)
        opts.photonCounts = [50 100 500 1000];
    end
    if ~isfield(opts, 'nReplicates') || isempty(opts.nReplicates)
        opts.nReplicates = 250;
    end
    if ~isfield(opts, 'nTimeBins') || isempty(opts.nTimeBins)
        opts.nTimeBins = 256;
    end
    if ~isfield(opts, 'dtNs') || isempty(opts.dtNs)
        opts.dtNs = 0.05;
    end
    if ~isfield(opts, 'pulsePeriodNs') || isempty(opts.pulsePeriodNs)
        opts.pulsePeriodNs = opts.nTimeBins * opts.dtNs;
    end
    if ~isfield(opts, 'backgroundFraction') || isempty(opts.backgroundFraction)
        opts.backgroundFraction = 0;
    end
    if ~isfield(opts, 'irf')
        opts.irf = [];
    end
    if ~isfield(opts, 'irfCenterBin') || isempty(opts.irfCenterBin)
        opts.irfCenterBin = 8;
    end
    if ~isfield(opts, 'irfSigmaBins') || isempty(opts.irfSigmaBins)
        opts.irfSigmaBins = 1.5;
    end
    if ~isfield(opts, 'includeBackground') || isempty(opts.includeBackground)
        opts.includeBackground = (opts.backgroundFraction > 0);
    end
    if ~isfield(opts, 'optimizeTau') || isempty(opts.optimizeTau)
        opts.optimizeTau = true;
    end
    if ~isfield(opts, 'signalGrid') || isempty(opts.signalGrid)
        if opts.includeBackground
            opts.signalGrid = linspace(0.0, 1.0, 26);
        else
            opts.signalGrid = 1;
        end
    end
    if ~isfield(opts, 'singleExpTauGrid') || isempty(opts.singleExpTauGrid)
        opts.singleExpTauGrid = [];
    end
    if ~isfield(opts, 'shiftBounds') || isempty(opts.shiftBounds)
        opts.shiftBounds = [-5 5];
    end
    if ~isfield(opts, 'batchSize') || isempty(opts.batchSize)
        opts.batchSize = 2048;
    end
    if ~isfield(opts, 'useGPU') || isempty(opts.useGPU)
        opts.useGPU = false;
    end
    if ~isfield(opts, 'showPlot') || isempty(opts.showPlot)
        opts.showPlot = true;
    end
    if ~isfield(opts, 'rngSeed') || isempty(opts.rngSeed)
        opts.rngSeed = 1;
    end
    if ~isfield(opts, 'calibrationMatFile') || isempty(opts.calibrationMatFile)
        opts.calibrationMatFile = 'calibrationCurve.mat';
    end
    if ~isfield(opts, 'curvePlotTauIdx')
        opts.curvePlotTauIdx = [];
    end
    if ~isfield(opts, 'curvePlotPhotonIdx')
        opts.curvePlotPhotonIdx = [];
    end
    if ~isfield(opts, 'curvePlotMaxTau') || isempty(opts.curvePlotMaxTau)
        opts.curvePlotMaxTau = 3;
    end
    if ~isfield(opts, 'curvePlotMaxPhoton') || isempty(opts.curvePlotMaxPhoton)
        opts.curvePlotMaxPhoton = 4;
    end

    opts.trueTauNs = double(opts.trueTauNs(:)).';
    if ~isempty(opts.tau0Ns)
        opts.tau0Ns = double(opts.tau0Ns(:)).';
    else
        opts.tau0Ns = [];
    end
    opts.photonCounts = max(1, round(double(opts.photonCounts(:)).'));
    opts.nReplicates = max(1, round(double(opts.nReplicates(1))));
    opts.nTimeBins = max(8, round(double(opts.nTimeBins(1))));
    opts.dtNs = double(opts.dtNs(1));
    opts.pulsePeriodNs = double(opts.pulsePeriodNs(1));
    opts.backgroundFraction = min(max(double(opts.backgroundFraction(1)), 0), 1);
    opts.irfCenterBin = double(opts.irfCenterBin(1));
    opts.irfSigmaBins = max(double(opts.irfSigmaBins(1)), eps);
    opts.signalGrid = double(opts.signalGrid(:)).';
    if ~isempty(opts.singleExpTauGrid)
        opts.singleExpTauGrid = double(opts.singleExpTauGrid(:)).';
    else
        opts.singleExpTauGrid = [];
    end
    opts.shiftBounds = double(opts.shiftBounds(:)).';
    opts.batchSize = max(1, round(double(opts.batchSize(1))));
    opts.rngSeed = round(double(opts.rngSeed(1)));
    opts.calibrationMatFile = char(opts.calibrationMatFile);
    opts.curvePlotTauIdx = round(double(opts.curvePlotTauIdx(:)).');
    opts.curvePlotPhotonIdx = round(double(opts.curvePlotPhotonIdx(:)).');
    opts.curvePlotMaxTau = max(1, round(double(opts.curvePlotMaxTau(1))));
    opts.curvePlotMaxPhoton = max(1, round(double(opts.curvePlotMaxPhoton(1))));

    opts.tau0Ns = normalize_tau0_seeds(opts.tau0Ns, opts.trueTauNs);
    if isempty(opts.singleExpTauGrid)
        tauMin = max(0.05, 0.4 * min([opts.trueTauNs opts.tau0Ns]));
        tauMax = max([1.8 * max([opts.trueTauNs opts.tau0Ns]), tauMin + 0.05]);
        opts.singleExpTauGrid = linspace(tauMin, tauMax, 61);
    end
end

function tau0Ns = normalize_tau0_seeds(tau0Ns, trueTauNs)
    if isempty(tau0Ns)
        tau0Ns = trueTauNs;
        return;
    end
    if isscalar(tau0Ns)
        tau0Ns = repmat(double(tau0Ns), size(trueTauNs));
        return;
    end
    if numel(tau0Ns) ~= numel(trueTauNs)
        error('opts.tau0Ns must be scalar or have the same length as opts.trueTauNs.');
    end
    tau0Ns = double(tau0Ns(:)).';
end

function irf = build_irf(opts)
    if ~isempty(opts.irf)
        irf = max(double(opts.irf(:)), 0);
        if numel(irf) ~= opts.nTimeBins
            error('opts.irf must have exactly %d bins.', opts.nTimeBins);
        end
    else
        bins = (1:opts.nTimeBins)';
        irf = exp(-0.5 .* ((bins - opts.irfCenterBin) ./ opts.irfSigmaBins) .^ 2);
    end

    irfSum = sum(irf);
    if ~(isfinite(irfSum) && irfSum > 0)
        error('IRF must contain at least one positive finite value.');
    end
    irf = irf ./ irfSum;
end

function prob = build_monoexp_probability(irf, tauNs, pulsePeriodNs, dtNs, backgroundFraction)
    nt = numel(irf);
    tNs = ((0:nt-1)') * dtNs;
    tauNs = max(double(tauNs), eps);
    pulsePeriodNs = max(double(pulsePeriodNs), nt * dtNs);

    decay = exp(-tNs ./ tauNs) ./ max(1 - exp(-pulsePeriodNs ./ tauNs), eps);
    signal = Convol(irf(:), decay(:));
    signal = max(double(signal(:)), 0);
    signal = signal ./ max(sum(signal), eps);

    if backgroundFraction > 0
        background = ones(nt, 1) ./ nt;
        prob = (1 - backgroundFraction) .* signal + backgroundFraction .* background;
    else
        prob = signal;
    end

    prob = max(prob, 0);
    prob = prob ./ max(sum(prob), eps);
end

function counts = simulate_exact_photon_histograms(prob, photonCount, nReplicates)
    prob = double(prob(:));
    prob = prob ./ max(sum(prob), eps);
    cdf = cumsum(prob);
    cdf(end) = 1;
    edges = [0; cdf];
    edges(end) = 1 + eps(1);

    nBins = numel(prob);
    counts = zeros(nBins, nReplicates, 'single');
    for rep = 1:nReplicates
        u = rand(photonCount, 1);
        counts(:, rep) = single(histcounts(u, edges).');
    end
end

function fitCounts = reconstruct_monoexp_fit_counts(estTauNs, signalFrac, totalCounts, irf, pulsePeriodNs, dtNs, bgPattern)
    nRep = numel(estTauNs);
    nt = numel(irf);
    fitCounts = zeros(nt, nRep);
    signalFrac = reshape(double(signalFrac), [], 1);
    totalCounts = reshape(double(totalCounts), [], 1);
    for rep = 1:nRep
        sigPattern = build_monoexp_probability(irf, estTauNs(rep), pulsePeriodNs, dtNs, 0);
        pattern = signalFrac(rep) .* sigPattern + (1 - signalFrac(rep)) .* bgPattern;
        pattern = pattern ./ max(sum(pattern), eps);
        fitCounts(:, rep) = totalCounts(rep) .* pattern(:);
    end
end

function row = summarize_errors(photonCount, tauTrueNs, tau0Ns, estTauNs, postStdNs)
    valid = isfinite(estTauNs);
    estTauNs = estTauNs(valid);
    postStdNs = postStdNs(valid);
    errNs = estTauNs - tauTrueNs;

    meanTau = mean(estTauNs);
    biasNs = mean(errNs);
    stdErrNs = std(errNs);
    rmseNs = sqrt(mean(errNs .^ 2));
    maeNs = mean(abs(errNs));
    medAbsErrNs = median(abs(errNs));
    meanPostStdNs = mean(postStdNs);
    q = local_quantiles(estTauNs, [0.05 0.50 0.95]);

    row = table( ...
        photonCount, ...
        numel(estTauNs), ...
        tauTrueNs, ...
        tau0Ns, ...
        meanTau, ...
        biasNs, ...
        stdErrNs, ...
        rmseNs, ...
        maeNs, ...
        medAbsErrNs, ...
        meanPostStdNs, ...
        q(1), ...
        q(2), ...
        q(3), ...
        'VariableNames', { ...
        'PhotonCount', ...
        'NReplicates', ...
        'TrueTauNs', ...
        'Tau0Ns', ...
        'MeanEstimatedTauNs', ...
        'BiasNs', ...
        'StdErrorNs', ...
        'RMSE_Ns', ...
        'MAE_Ns', ...
        'MedianAbsErrorNs', ...
        'MeanPosteriorStdNs', ...
        'Q05EstimatedTauNs', ...
        'Q50EstimatedTauNs', ...
        'Q95EstimatedTauNs'});
end

function row = append_height_summary(row, trueHeightNm, estHeightNm)
    valid = isfinite(estHeightNm);
    estHeightNm = double(estHeightNm(valid));
    errHeightNm = estHeightNm - trueHeightNm;

    row.TrueHeightNm = trueHeightNm;
    row.MeanEstimatedHeightNm = mean(estHeightNm);
    row.BiasHeightNm = mean(errHeightNm);
    row.StdErrorHeightNm = std(errHeightNm);
    row.RMSE_HeightNm = sqrt(mean(errHeightNm .^ 2));
end

function q = local_quantiles(x, probs)
    x = sort(double(x(:)));
    probs = double(probs(:));
    if isempty(x)
        q = nan(size(probs));
        return;
    end
    if numel(x) == 1
        q = repmat(x, size(probs));
        return;
    end

    pos = 1 + (numel(x) - 1) .* min(max(probs, 0), 1);
    lo = floor(pos);
    hi = ceil(pos);
    w = pos - lo;
    q = (1 - w) .* x(lo) + w .* x(hi);
end

function plot_summary(results)
    if numel(results.trueTauNs) == 1
        plot_scalar_summary(results);
    else
        plot_grid_summary(results);
    end
    plot_lifetime_deviation_summary(results);
    plot_height_summary(results);
    plot_condition_curves(results);
end

function plot_scalar_summary(results)
    summary = results.summaryTable;
    photonCounts = summary.PhotonCount;
    tauTrueNs = results.trueTauNs(1);

    figure('Color', 'w', 'Name', 'FLIM_bayes monoexponential Monte Carlo');

    subplot(2, 1, 1);
    hold on;
    colors = lines(numel(results.replicates));
    for idx = 1:numel(results.replicates)
        photonCount = results.replicates(idx).photonCount;
        estTauNs = results.replicates(idx).estimatedTauNs;
        xj = photonCount .* 10 .^ (0.045 .* (rand(size(estTauNs)) - 0.5));
        plot(xj, estTauNs, '.', 'Color', colors(idx, :), 'MarkerSize', 8);
        plot(photonCount, summary.MeanEstimatedTauNs(idx), 'o', ...
            'Color', colors(idx, :), 'MarkerFaceColor', colors(idx, :), 'MarkerSize', 6);
        plot([photonCount photonCount], ...
            summary.MeanEstimatedTauNs(idx) + [-1 1] .* summary.StdErrorNs(idx), ...
            '-', 'Color', colors(idx, :), 'LineWidth', 1.5);
    end
    plot([min(photonCounts) max(photonCounts)], [tauTrueNs tauTrueNs], 'k--', 'LineWidth', 1.2);
    set(gca, 'XScale', 'log');
    xlabel('Photons per decay');
    ylabel('Estimated lifetime (ns)');
    title(sprintf('FLIM\\_bayes monoexponential fits (true \\tau = %.3f ns)', tauTrueNs));
    grid on;
    box on;

    subplot(2, 1, 2);
    hold on;
    loglog(photonCounts, summary.RMSE_Ns, '-o', 'LineWidth', 1.5, 'DisplayName', 'RMSE');
    loglog(photonCounts, summary.StdErrorNs, '-s', 'LineWidth', 1.5, 'DisplayName', 'Std(error)');
    loglog(photonCounts, abs(summary.BiasNs), '-d', 'LineWidth', 1.5, 'DisplayName', '|Bias|');
    loglog(photonCounts, summary.MeanPosteriorStdNs, '-^', 'LineWidth', 1.5, 'DisplayName', 'Mean posterior std');
    xlabel('Photons per decay');
    ylabel('Lifetime error (ns)');
    title('Error summary versus photon count');
    grid on;
    box on;
    legend('Location', 'southwest');
end

function plot_grid_summary(results)
    maps = results.summaryMaps;
    photonCounts = maps.photonAxis;
    tauAxisNs = maps.tauAxisNs(:);
    [X, Y] = meshgrid(photonCounts, tauAxisNs);

    figure('Color', 'w', 'Name', 'FLIM_bayes monoexponential uncertainty maps');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    render_uncertainty_map(X, Y, maps.StdErrorNs, photonCounts, 'Empirical std(error)', 'Std(error) [ns]');

    nexttile;
    render_uncertainty_map(X, Y, maps.MeanPosteriorStdNs, photonCounts, 'Mean posterior std', 'Posterior std [ns]');

    nexttile;
    render_uncertainty_map(X, Y, maps.RMSE_Ns, photonCounts, 'RMSE', 'RMSE [ns]');

    nexttile;
    render_uncertainty_map(X, Y, maps.BiasNs, photonCounts, 'Bias', 'Bias [ns]');
end

function render_uncertainty_map(X, Y, Z, photonCounts, titleStr, colorbarStr)
    surface(X, Y, zeros(size(Z)), Z, 'EdgeColor', 'none');
    view(2);
    set(gca, 'XScale', 'log', 'XTick', photonCounts);
    xlabel('Photons per decay');
    ylabel('True lifetime (ns)');
    title(titleStr);
    grid on;
    box on;
    cb = colorbar;
    cb.Label.String = colorbarStr;
end

function plot_condition_curves(results)
    tauIdx = choose_plot_indices(numel(results.trueTauNs), results.optsUsed.curvePlotTauIdx, results.optsUsed.curvePlotMaxTau);
    photIdx = choose_plot_indices(numel(results.optsUsed.photonCounts), results.optsUsed.curvePlotPhotonIdx, results.optsUsed.curvePlotMaxPhoton);
    nTau = numel(tauIdx);
    nPhot = numel(photIdx);
    tAxisNs = results.tAxisNs(:);

    figure('Color', 'w', 'Name', 'FLIM_bayes monoexponential mean data and fit curves', ...
        'Units', 'normalized', 'Position', [0.02 0.04 0.96 0.90]);
    tiledlayout(nTau, nPhot, 'TileSpacing', 'compact', 'Padding', 'compact');

    for iTau = 1:nTau
        for iPhot = 1:nPhot
            nexttile;
            tauHere = tauIdx(iTau);
            photHere = photIdx(iPhot);
            entry = results.replicates(tauHere, photHere);
            plot_shaded_mean_curve(gca, tAxisNs, entry.meanDataCurve, entry.stdDataCurve, [0.20 0.45 0.85], 'Data', 0.14, false);
            plot_shaded_mean_curve(gca, tAxisNs, entry.meanFitCurve, entry.stdFitCurve, [0.85 0.33 0.10], 'Fit', 0.28, false);
            plot(gca, tAxisNs, entry.meanDataCurve, 'Color', [0.20 0.45 0.85], 'LineWidth', 1.4, 'DisplayName', 'Data');
            plot(gca, tAxisNs, entry.meanFitCurve, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6, 'DisplayName', 'Fit');
            title(sprintf('\\tau=%.2f ns | h=%.1f nm | N=%d', ...
                results.trueTauNs(tauHere), results.trueHeightNm(tauHere), entry.photonCount), 'FontSize', 9);
            grid on;
            box on;
            if iTau == nTau
                xlabel('Time (ns)');
            end
            if iPhot == 1
                ylabel('Counts');
            end
            if iTau == 1 && iPhot == 1
                legend('Location', 'northeast');
            end
        end
    end
end

function plot_height_summary(results)
    maps = results.summaryMaps;
    photonCounts = maps.photonAxis;
    heightAxisNm = maps.heightAxisNm(:);
    [X, Y] = meshgrid(photonCounts, heightAxisNm);

    figure('Color', 'w', 'Name', 'FLIM_bayes monoexponential height summary');
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    render_uncertainty_map(X, Y, maps.BiasHeightNm, photonCounts, 'Height bias: estimated - true', 'Bias (nm)');
    ylabel(gca, 'True height (nm)');

    nexttile;
    render_uncertainty_map(X, Y, maps.StdErrorHeightNm, photonCounts, 'Height std(error)', 'Std(error) (nm)');
    ylabel(gca, 'True height (nm)');

    nexttile;
    render_uncertainty_map(X, Y, maps.MeanEstimatedHeightNm, photonCounts, 'Mean estimated height', 'Height (nm)');
    ylabel(gca, 'True height (nm)');
end

function plot_lifetime_deviation_summary(results)
    maps = results.summaryMaps;
    photonCounts = maps.photonAxis;
    tauAxisNs = maps.tauAxisNs(:);
    [X, Y] = meshgrid(photonCounts, tauAxisNs);

    figure('Color', 'w', 'Name', 'FLIM_bayes monoexponential lifetime deviation');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    render_uncertainty_map(X, Y, maps.BiasNs, photonCounts, ...
        'Lifetime deviation: estimated - true', 'Deviation (ns)');

    nexttile;
    hold on;
    colors = lines(numel(photonCounts));
    for photIdx = 1:numel(photonCounts)
        plot(tauAxisNs, maps.MeanEstimatedTauNs(:, photIdx), '-', ...
            'Color', colors(photIdx, :), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%d photons', photonCounts(photIdx)));
    end
    plot(tauAxisNs, tauAxisNs, 'k--', 'LineWidth', 1.2, 'DisplayName', 'True = fitted');
    xlabel('True lifetime (ns)');
    ylabel('Mean fitted lifetime (ns)');
    title('Fitted lifetime vs true lifetime');
    grid on;
    box on;
    legend('Location', 'northwest');
end

function idx = choose_plot_indices(n, requestedIdx, maxCount)
    if isempty(requestedIdx)
        if n <= maxCount
            idx = 1:n;
        else
            idx = unique(round(linspace(1, n, maxCount)));
        end
    else
        idx = unique(requestedIdx(requestedIdx >= 1 & requestedIdx <= n));
        if isempty(idx)
            idx = 1:min(n, maxCount);
        end
    end
end
