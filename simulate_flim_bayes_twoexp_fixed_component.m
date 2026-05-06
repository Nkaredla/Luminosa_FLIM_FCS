function results = simulate_flim_bayes_twoexp_fixed_component(opts)
% SIMULATE_FLIM_BAYES_TWOEXP_FIXED_COMPONENT
% Monte Carlo benchmark for separating two lifetime components when one
% lifetime is fixed and known.
%
% This companion simulator uses the same IRF-convolved Bayesian grid
% posterior idea as FLIM_BAYES_LOWPHOTON, but keeps one component fixed at
% a known lifetime (default 5.5 ns). The second component is swept over a
% user-defined range, and the two signal components contribute equal photon
% numbers in every replicate.
%
% Example
%   results = simulate_flim_bayes_twoexp_fixed_component();
%
%   opts = struct();
%   opts.trueTauUnknownNs = 0.5:0.25:5.0;
%   opts.totalPhotonCounts = [100 200 500 1000];
%   opts.nReplicates = 500;
%   results = simulate_flim_bayes_twoexp_fixed_component(opts);
%
% INPUT
%   opts fields (all optional)
%       .trueTauUnknownNs    vector of true unknown-component lifetimes [ns]
%       .knownTauNs          fixed known lifetime [ns], default 5.5
%       .totalPhotonCounts   total photons per decay (both components summed)
%       .nReplicates         Monte Carlo repeats per condition
%       .tauGridNs           candidate grid for the unknown lifetime [ns]
%       .fractionGrid        candidate grid for unknown-component fraction
%       .backgroundFraction  fraction of total photons assigned to uniform bg
%       .includeBackground   include bg in the fit model, default auto
%       .signalGrid          signal-fraction grid if background is modeled
%       .nTimeBins           number of TCSPC bins
%       .dtNs                TCSPC bin width [ns]
%       .pulsePeriodNs       laser repetition period [ns]
%       .irf                 explicit IRF sampled on the TCSPC grid
%       .irfCenterBin        Gaussian IRF center if .irf is omitted
%       .irfSigmaBins        Gaussian IRF sigma in bins if .irf is omitted
%       .batchSize           posterior evaluation batch size
%       .useGPU              use gpuArray for posterior matrix multiply
%       .showPlot            draw summary heatmaps, default true
%       .rngSeed             RNG seed, default 1
%       .calibrationMatFile  MIET calibration MAT file, default calibrationCurve.mat
%       .curvePlotTauIdx     selected unknown-lifetime indices for curve subplots
%       .curvePlotPhotonIdx  selected photon-count indices for curve subplots
%       .curvePlotMaxTau     default representative tau count if indices omitted
%       .curvePlotMaxPhoton  default representative photon count if indices omitted
%
% OUTPUT
%   results struct with fields:
%       .summaryTable               row-wise summary by (tau, photons)
%       .replicates                 [nTau x nPhot] replicate-level results
%       .summaryMaps                2D maps over (tau, photons)
%       .tauUncertaintyMapNs        empirical std(error) for unknown tau
%       .fractionUncertaintyMap     empirical std(error) for unknown fraction
%       .trueUnknownTauAxisNs       lifetime axis
%       .photonAxis                 photon axis
%       .knownTauNs                 fixed known lifetime
%       .trueFractionUnknown        true signal fraction of unknown component
%       .calibration                loaded calibration struct

    if nargin < 1 || isempty(opts)
        opts = struct();
    end
    opts = fill_default_options(opts);

    rng(opts.rngSeed, 'twister');

    nt = opts.nTimeBins;
    dtNs = opts.dtNs;
    pulsePeriodNs = opts.pulsePeriodNs;
    tauTrueList = opts.trueTauUnknownNs(:).';
    photonCounts = opts.totalPhotonCounts(:).';
    trueFracUnknown = 0.5;

    irf = build_irf(opts);
    calib = load_lifetime_height_calibration(opts.calibrationMatFile);
    bgPattern = ones(nt, 1) ./ nt;
    knownPattern = build_signal_pattern(irf, opts.knownTauNs, pulsePeriodNs, dtNs);

    [Pgrid, tauPerGrid, fracPerGrid, signalPerGrid] = build_fixed_component_grid( ...
        irf, bgPattern, pulsePeriodNs, dtNs, opts.knownTauNs, opts.tauGridNs, ...
        opts.fractionGrid, opts.signalGrid, opts.includeBackground);

    nTau = numel(tauTrueList);
    nPhot = numel(photonCounts);

    replicates = repmat(struct( ...
        'trueTauUnknownNs', [], ...
        'knownTauNs', [], ...
        'totalPhotonCount', [], ...
        'trueFractionUnknown', [], ...
        'estimatedTauUnknownNs', [], ...
        'estimatedHeightUnknownNm', [], ...
        'posteriorStdTauUnknownNs', [], ...
        'estimatedFractionUnknown', [], ...
        'posteriorStdFractionUnknown', [], ...
        'estimatedSignalFraction', [], ...
        'errorTauUnknownNs', [], ...
        'errorHeightUnknownNm', [], ...
        'errorFractionUnknown', [], ...
        'meanDataCurve', [], ...
        'stdDataCurve', [], ...
        'meanFitCurve', [], ...
        'stdFitCurve', [], ...
        'summary', [], ...
        'posteriorInfo', []), nTau, nPhot);

    tauMeanMap = nan(nTau, nPhot);
    tauBiasMap = nan(nTau, nPhot);
    tauStdErrorMap = nan(nTau, nPhot);
    tauRmseMap = nan(nTau, nPhot);
    tauMaeMap = nan(nTau, nPhot);
    tauPostStdMap = nan(nTau, nPhot);

    fracMeanMap = nan(nTau, nPhot);
    fracBiasMap = nan(nTau, nPhot);
    fracStdErrorMap = nan(nTau, nPhot);
    fracRmseMap = nan(nTau, nPhot);
    fracMaeMap = nan(nTau, nPhot);
    fracPostStdMap = nan(nTau, nPhot);

    signalMeanMap = nan(nTau, nPhot);
    trueHeightUnknownNm = nan(1, nTau);
    meanHeightMap = nan(nTau, nPhot);
    biasHeightMap = nan(nTau, nPhot);
    stdErrorHeightMap = nan(nTau, nPhot);
    rmseHeightMap = nan(nTau, nPhot);

    summaryRows = cell(nTau * nPhot, 1);
    rowIdx = 1;

    for tauIdx = 1:nTau
        tauTrue = tauTrueList(tauIdx);
        unknownPattern = build_signal_pattern(irf, tauTrue, pulsePeriodNs, dtNs);
        trueHeightHere = lifetime_to_height_from_calibration(tauTrue, calib);
        trueHeightUnknownNm(tauIdx) = trueHeightHere;

        for photIdx = 1:nPhot
            totalPhotons = photonCounts(photIdx);
            counts = simulate_equal_component_histograms( ...
                unknownPattern, knownPattern, bgPattern, totalPhotons, ...
                opts.backgroundFraction, opts.nReplicates);
            tcspc_pix = permute(reshape(counts, nt, opts.nReplicates, 1), [2 3 1]);

            posterior = evaluate_fixed_component_posterior( ...
                tcspc_pix, Pgrid, tauPerGrid, fracPerGrid, signalPerGrid, ...
                opts.useGPU, opts.batchSize);

            estTau = double(posterior.tauUnknownMean(:));
            estFrac = double(posterior.fractionUnknownMean(:));
            estSignal = double(posterior.signalFractionMean(:));
            postTauStd = double(posterior.tauUnknownStd(:));
            postFracStd = double(posterior.fractionUnknownStd(:));
            estHeightNm = lifetime_to_height_from_calibration(estTau, calib);

            totalCounts = sum(double(counts), 1);
            fitCounts = reconstruct_twoexp_fit_counts( ...
                estTau, estFrac, estSignal, totalCounts, irf, pulsePeriodNs, dtNs, ...
                opts.knownTauNs, bgPattern);

            row = summarize_twoexp_errors( ...
                totalPhotons, tauTrue, opts.knownTauNs, trueFracUnknown, ...
                estTau, postTauStd, estFrac, postFracStd, estSignal);
            row = append_height_summary(row, trueHeightHere, estHeightNm);

            replicates(tauIdx, photIdx).trueTauUnknownNs = tauTrue;
            replicates(tauIdx, photIdx).knownTauNs = opts.knownTauNs;
            replicates(tauIdx, photIdx).totalPhotonCount = totalPhotons;
            replicates(tauIdx, photIdx).trueFractionUnknown = trueFracUnknown;
            replicates(tauIdx, photIdx).estimatedTauUnknownNs = estTau;
            replicates(tauIdx, photIdx).estimatedHeightUnknownNm = estHeightNm;
            replicates(tauIdx, photIdx).posteriorStdTauUnknownNs = postTauStd;
            replicates(tauIdx, photIdx).estimatedFractionUnknown = estFrac;
            replicates(tauIdx, photIdx).posteriorStdFractionUnknown = postFracStd;
            replicates(tauIdx, photIdx).estimatedSignalFraction = estSignal;
            replicates(tauIdx, photIdx).errorTauUnknownNs = estTau - tauTrue;
            replicates(tauIdx, photIdx).errorHeightUnknownNm = estHeightNm - trueHeightHere;
            replicates(tauIdx, photIdx).errorFractionUnknown = estFrac - trueFracUnknown;
            replicates(tauIdx, photIdx).meanDataCurve = mean(double(counts), 2);
            replicates(tauIdx, photIdx).stdDataCurve = std(double(counts), 0, 2);
            replicates(tauIdx, photIdx).meanFitCurve = mean(fitCounts, 2);
            replicates(tauIdx, photIdx).stdFitCurve = std(fitCounts, 0, 2);
            replicates(tauIdx, photIdx).summary = row;
            replicates(tauIdx, photIdx).posteriorInfo = posterior.posteriorInfo;

            summaryRows{rowIdx} = row;
            tauMeanMap(tauIdx, photIdx) = row.MeanEstimatedTauUnknownNs;
            tauBiasMap(tauIdx, photIdx) = row.BiasTauUnknownNs;
            tauStdErrorMap(tauIdx, photIdx) = row.StdErrorTauUnknownNs;
            tauRmseMap(tauIdx, photIdx) = row.RMSE_TauUnknownNs;
            tauMaeMap(tauIdx, photIdx) = row.MAE_TauUnknownNs;
            tauPostStdMap(tauIdx, photIdx) = row.MeanPosteriorStdTauUnknownNs;

            fracMeanMap(tauIdx, photIdx) = row.MeanEstimatedFractionUnknown;
            fracBiasMap(tauIdx, photIdx) = row.BiasFractionUnknown;
            fracStdErrorMap(tauIdx, photIdx) = row.StdErrorFractionUnknown;
            fracRmseMap(tauIdx, photIdx) = row.RMSE_FractionUnknown;
            fracMaeMap(tauIdx, photIdx) = row.MAE_FractionUnknown;
            fracPostStdMap(tauIdx, photIdx) = row.MeanPosteriorStdFractionUnknown;

            signalMeanMap(tauIdx, photIdx) = row.MeanEstimatedSignalFraction;
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
    results.summaryMaps = struct( ...
        'MeanEstimatedTauUnknownNs', tauMeanMap, ...
        'BiasTauUnknownNs', tauBiasMap, ...
        'StdErrorTauUnknownNs', tauStdErrorMap, ...
        'RMSE_TauUnknownNs', tauRmseMap, ...
        'MAE_TauUnknownNs', tauMaeMap, ...
        'MeanPosteriorStdTauUnknownNs', tauPostStdMap, ...
        'MeanEstimatedFractionUnknown', fracMeanMap, ...
        'BiasFractionUnknown', fracBiasMap, ...
        'StdErrorFractionUnknown', fracStdErrorMap, ...
        'RMSE_FractionUnknown', fracRmseMap, ...
        'MAE_FractionUnknown', fracMaeMap, ...
        'MeanPosteriorStdFractionUnknown', fracPostStdMap, ...
        'MeanEstimatedSignalFraction', signalMeanMap, ...
        'MeanEstimatedHeightNm', meanHeightMap, ...
        'BiasHeightNm', biasHeightMap, ...
        'StdErrorHeightNm', stdErrorHeightMap, ...
        'RMSE_HeightNm', rmseHeightMap, ...
        'tauAxisNs', tauTrueList(:), ...
        'heightAxisNm', trueHeightUnknownNm(:), ...
        'photonAxis', photonCounts(:).');
    results.tauUncertaintyMapNs = tauStdErrorMap;
    results.fractionUncertaintyMap = fracStdErrorMap;
    results.tauRmseMapNs = tauRmseMap;
    results.tauBiasMapNs = tauBiasMap;
    results.lifetimeDeviationMapNs = tauBiasMap;
    results.fractionRmseMap = fracRmseMap;
    results.trueUnknownTauAxisNs = tauTrueList(:);
    results.trueUnknownHeightNm = trueHeightUnknownNm(:);
    results.photonAxis = photonCounts(:).';
    results.knownTauNs = opts.knownTauNs;
    results.trueFractionUnknown = trueFracUnknown;
    results.irf = irf(:);
    results.knownPattern = knownPattern(:);
    results.calibration = calib;
    results.heightBiasMapNm = biasHeightMap;
    results.heightUncertaintyMapNm = stdErrorHeightMap;
    results.heightRmseMapNm = rmseHeightMap;
    results.posteriorGrid = struct( ...
        'tauGridNs', opts.tauGridNs(:), ...
        'fractionGrid', opts.fractionGrid(:), ...
        'signalGrid', opts.signalGrid(:), ...
        'nGrid', size(Pgrid, 2));
    results.optsUsed = opts;

    disp(summaryTable);

    if opts.showPlot
        plot_twoexp_summary(results);
    end
end

function opts = fill_default_options(opts)
    if ~isfield(opts, 'trueTauUnknownNs') || isempty(opts.trueTauUnknownNs)
        opts.trueTauUnknownNs = 0.5:0.25:5.0;
    end
    if ~isfield(opts, 'knownTauNs') || isempty(opts.knownTauNs)
        opts.knownTauNs = 5.5;
    end
    if ~isfield(opts, 'totalPhotonCounts') || isempty(opts.totalPhotonCounts)
        opts.totalPhotonCounts = [26 50 76 100 200 500 750 1000];
    end
    if ~isfield(opts, 'nReplicates') || isempty(opts.nReplicates)
        opts.nReplicates = 250;
    end
    if ~isfield(opts, 'backgroundFraction') || isempty(opts.backgroundFraction)
        opts.backgroundFraction = 0;
    end
    if ~isfield(opts, 'includeBackground') || isempty(opts.includeBackground)
        opts.includeBackground = (opts.backgroundFraction > 0);
    end
    if ~isfield(opts, 'fractionGrid') || isempty(opts.fractionGrid)
        opts.fractionGrid = linspace(0.0, 1.0, 51);
    end
    if ~isfield(opts, 'signalGrid') || isempty(opts.signalGrid)
        if opts.includeBackground
            opts.signalGrid = linspace(0.0, 1.0, 26);
        else
            opts.signalGrid = 1;
        end
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
    if ~isfield(opts, 'irf')
        opts.irf = [];
    end
    if ~isfield(opts, 'irfCenterBin') || isempty(opts.irfCenterBin)
        opts.irfCenterBin = 8;
    end
    if ~isfield(opts, 'irfSigmaBins') || isempty(opts.irfSigmaBins)
        opts.irfSigmaBins = 1.5;
    end
    if ~isfield(opts, 'tauGridNs') || isempty(opts.tauGridNs)
        tauMaxDefault = max(0.25, opts.knownTauNs - 0.05);
        opts.tauGridNs = unique([0.2:0.05:tauMaxDefault, opts.trueTauUnknownNs(:).']);
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

    opts.trueTauUnknownNs = double(opts.trueTauUnknownNs(:)).';
    opts.knownTauNs = double(opts.knownTauNs(1));
    opts.totalPhotonCounts = max(2, round(double(opts.totalPhotonCounts(:)).'));
    opts.nReplicates = max(1, round(double(opts.nReplicates(1))));
    opts.backgroundFraction = min(max(double(opts.backgroundFraction(1)), 0), 1);
    opts.includeBackground = logical(opts.includeBackground);
    opts.fractionGrid = min(max(double(opts.fractionGrid(:)).', 0), 1);
    opts.signalGrid = min(max(double(opts.signalGrid(:)).', 0), 1);
    opts.nTimeBins = max(8, round(double(opts.nTimeBins(1))));
    opts.dtNs = double(opts.dtNs(1));
    opts.pulsePeriodNs = double(opts.pulsePeriodNs(1));
    opts.irfCenterBin = double(opts.irfCenterBin(1));
    opts.irfSigmaBins = max(double(opts.irfSigmaBins(1)), eps);
    opts.tauGridNs = unique(max(double(opts.tauGridNs(:)).', 0.03));
    opts.batchSize = max(1, round(double(opts.batchSize(1))));
    opts.useGPU = logical(opts.useGPU);
    opts.showPlot = logical(opts.showPlot);
    opts.rngSeed = round(double(opts.rngSeed(1)));
    opts.calibrationMatFile = char(opts.calibrationMatFile);
    opts.curvePlotTauIdx = round(double(opts.curvePlotTauIdx(:)).');
    opts.curvePlotPhotonIdx = round(double(opts.curvePlotPhotonIdx(:)).');
    opts.curvePlotMaxTau = max(1, round(double(opts.curvePlotMaxTau(1))));
    opts.curvePlotMaxPhoton = max(1, round(double(opts.curvePlotMaxPhoton(1))));

    if any(opts.trueTauUnknownNs <= 0)
        error('opts.trueTauUnknownNs must be positive.');
    end
    if opts.knownTauNs <= 0
        error('opts.knownTauNs must be positive.');
    end
    if any(opts.totalPhotonCounts <= 0)
        error('opts.totalPhotonCounts must be positive.');
    end
    nSignalTest = opts.totalPhotonCounts - round(opts.backgroundFraction * opts.totalPhotonCounts);
    if any(mod(nSignalTest, 2) ~= 0)
        error(['Total signal photons must be even to keep both components equal. ' ...
               'Choose even totalPhotonCounts or a backgroundFraction that preserves an even signal count.']);
    end
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
        error('IRF must contain positive finite values.');
    end
    irf = irf ./ irfSum;
end

function pattern = build_signal_pattern(irf, tauNs, pulsePeriodNs, dtNs)
    nt = numel(irf);
    tNs = ((0:nt-1)') * dtNs;
    tauNs = max(double(tauNs), eps);
    pulsePeriodNs = max(double(pulsePeriodNs), nt * dtNs);

    decay = exp(-tNs ./ tauNs) ./ max(1 - exp(-pulsePeriodNs ./ tauNs), eps);
    pattern = Convol(irf(:), decay(:));
    pattern = max(double(pattern(:)), 0);
    pattern = pattern ./ max(sum(pattern), eps);
end

function [Pgrid, tauPerGrid, fracPerGrid, signalPerGrid] = build_fixed_component_grid( ...
        irf, bgPattern, pulsePeriodNs, dtNs, knownTauNs, tauGridNs, fractionGrid, signalGrid, includeBackground)

    knownPattern = build_signal_pattern(irf, knownTauNs, pulsePeriodNs, dtNs);
    tauGridNs = double(tauGridNs(:));
    fractionGrid = double(fractionGrid(:));
    signalGrid = double(signalGrid(:));

    nTau = numel(tauGridNs);
    nFrac = numel(fractionGrid);
    nSig = numel(signalGrid);
    nt = numel(irf);

    if includeBackground
        nGrid = nTau * nFrac * nSig;
    else
        nGrid = nTau * nFrac;
        signalGrid = 1;
        nSig = 1;
    end

    Pgrid = zeros(nt, nGrid);
    tauPerGrid = zeros(nGrid, 1);
    fracPerGrid = zeros(nGrid, 1);
    signalPerGrid = zeros(nGrid, 1);

    col = 1;
    for tauIdx = 1:nTau
        unknownPattern = build_signal_pattern(irf, tauGridNs(tauIdx), pulsePeriodNs, dtNs);
        for fracIdx = 1:nFrac
            frac = fractionGrid(fracIdx);
            mixSignal = frac .* unknownPattern + (1 - frac) .* knownPattern;
            for sigIdx = 1:nSig
                sigFrac = signalGrid(sigIdx);
                if includeBackground
                    pattern = sigFrac .* mixSignal + (1 - sigFrac) .* bgPattern;
                else
                    pattern = mixSignal;
                end
                pattern = max(pattern(:), 0);
                pattern = pattern ./ max(sum(pattern), eps);
                Pgrid(:, col) = pattern;
                tauPerGrid(col) = tauGridNs(tauIdx);
                fracPerGrid(col) = frac;
                signalPerGrid(col) = sigFrac;
                col = col + 1;
            end
        end
    end
end

function counts = simulate_equal_component_histograms(unknownPattern, knownPattern, bgPattern, totalPhotons, backgroundFraction, nReplicates)
    totalPhotons = round(double(totalPhotons));
    nBg = round(backgroundFraction * totalPhotons);
    nSignal = totalPhotons - nBg;
    if mod(nSignal, 2) ~= 0
        error('Signal photon count must be even to keep the two components equal.');
    end

    nPerComp = nSignal / 2;
    nt = numel(unknownPattern);
    counts = zeros(nt, nReplicates, 'single');
    for rep = 1:nReplicates
        c1 = draw_histogram_counts(unknownPattern, nPerComp);
        c2 = draw_histogram_counts(knownPattern, nPerComp);
        cbg = draw_histogram_counts(bgPattern, nBg);
        counts(:, rep) = single(c1 + c2 + cbg);
    end
end

function counts = draw_histogram_counts(prob, nPhotons)
    prob = double(prob(:));
    prob = prob ./ max(sum(prob), eps);
    if nPhotons <= 0
        counts = zeros(size(prob));
        return;
    end
    cdf = cumsum(prob);
    cdf(end) = 1;
    edges = [0; cdf];
    edges(end) = 1 + eps(1);
    u = rand(nPhotons, 1);
    counts = histcounts(u, edges).';
end

function fitCounts = reconstruct_twoexp_fit_counts(estTau, estFrac, estSignal, totalCounts, irf, pulsePeriodNs, dtNs, knownTauNs, bgPattern)
    nRep = numel(estTau);
    nt = numel(irf);
    fitCounts = zeros(nt, nRep);
    estFrac = reshape(double(estFrac), [], 1);
    estSignal = reshape(double(estSignal), [], 1);
    totalCounts = reshape(double(totalCounts), [], 1);
    knownPattern = build_signal_pattern(irf, knownTauNs, pulsePeriodNs, dtNs);
    for rep = 1:nRep
        unknownPattern = build_signal_pattern(irf, estTau(rep), pulsePeriodNs, dtNs);
        mixSignal = estFrac(rep) .* unknownPattern + (1 - estFrac(rep)) .* knownPattern;
        pattern = estSignal(rep) .* mixSignal + (1 - estSignal(rep)) .* bgPattern;
        pattern = pattern ./ max(sum(pattern), eps);
        fitCounts(:, rep) = totalCounts(rep) .* pattern(:);
    end
end

function posterior = evaluate_fixed_component_posterior(Ypix, Pgrid, tauPerGrid, fracPerGrid, signalPerGrid, useGPU, batchSize)
    [nx, ny, nt] = size(Ypix);
    Y = reshape(permute(single(Ypix), [3 1 2]), nt, []);
    nPix = size(Y, 2);
    validPix = sum(Y, 1) > 0;

    tauMean = nan(1, nPix, 'single');
    tauStd = nan(1, nPix, 'single');
    tauMap = nan(1, nPix, 'single');
    fracMean = nan(1, nPix, 'single');
    fracStd = nan(1, nPix, 'single');
    fracMap = nan(1, nPix, 'single');
    signalMean = nan(1, nPix, 'single');

    tauPerGrid = double(tauPerGrid(:));
    tauSqPerGrid = tauPerGrid .^ 2;
    fracPerGrid = double(fracPerGrid(:));
    fracSqPerGrid = fracPerGrid .^ 2;
    signalPerGrid = double(signalPerGrid(:));

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

        fracMean(idx) = single(fracPerGrid.' * w);
        fracVar = fracSqPerGrid.' * w - double(fracMean(idx)) .^ 2;
        fracStd(idx) = single(sqrt(max(fracVar, 0)));

        signalMean(idx) = single(signalPerGrid.' * w);

        [~, mapIdx] = max(w, [], 1);
        tauMap(idx) = single(tauPerGrid(mapIdx));
        fracMap(idx) = single(fracPerGrid(mapIdx));
    end

    posterior = struct();
    posterior.tauUnknownMean = reshape(tauMean, nx, ny);
    posterior.tauUnknownStd = reshape(tauStd, nx, ny);
    posterior.tauUnknownMap = reshape(tauMap, nx, ny);
    posterior.fractionUnknownMean = reshape(fracMean, nx, ny);
    posterior.fractionUnknownStd = reshape(fracStd, nx, ny);
    posterior.fractionUnknownMap = reshape(fracMap, nx, ny);
    posterior.signalFractionMean = reshape(signalMean, nx, ny);
    posterior.posteriorInfo = struct( ...
        'usedGPU', gpuUsed, ...
        'nGrid', size(Pgrid, 2), ...
        'batchSize', batchSize, ...
        'method', 'Bayesian 2-state posterior with one fixed lifetime');
end

function row = summarize_twoexp_errors(totalPhotons, tauTrue, knownTau, fracTrue, estTau, postTauStd, estFrac, postFracStd, estSignal)
    valid = isfinite(estTau) & isfinite(estFrac);
    estTau = estTau(valid);
    postTauStd = postTauStd(valid);
    estFrac = estFrac(valid);
    postFracStd = postFracStd(valid);
    estSignal = estSignal(valid);

    tauErr = estTau - tauTrue;
    fracErr = estFrac - fracTrue;

    row = table( ...
        totalPhotons, ...
        numel(estTau), ...
        tauTrue, ...
        knownTau, ...
        fracTrue, ...
        mean(estTau), ...
        mean(tauErr), ...
        std(tauErr), ...
        sqrt(mean(tauErr .^ 2)), ...
        mean(abs(tauErr)), ...
        mean(postTauStd), ...
        mean(estFrac), ...
        mean(fracErr), ...
        std(fracErr), ...
        sqrt(mean(fracErr .^ 2)), ...
        mean(abs(fracErr)), ...
        mean(postFracStd), ...
        mean(estSignal), ...
        'VariableNames', { ...
        'TotalPhotonCount', ...
        'NReplicates', ...
        'TrueTauUnknownNs', ...
        'KnownTauNs', ...
        'TrueFractionUnknown', ...
        'MeanEstimatedTauUnknownNs', ...
        'BiasTauUnknownNs', ...
        'StdErrorTauUnknownNs', ...
        'RMSE_TauUnknownNs', ...
        'MAE_TauUnknownNs', ...
        'MeanPosteriorStdTauUnknownNs', ...
        'MeanEstimatedFractionUnknown', ...
        'BiasFractionUnknown', ...
        'StdErrorFractionUnknown', ...
        'RMSE_FractionUnknown', ...
        'MAE_FractionUnknown', ...
        'MeanPosteriorStdFractionUnknown', ...
        'MeanEstimatedSignalFraction'});
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

function plot_twoexp_summary(results)
    maps = results.summaryMaps;
    photonCounts = maps.photonAxis;
    tauAxisNs = maps.tauAxisNs(:);
    [X, Y] = meshgrid(photonCounts, tauAxisNs);

    figure('Color', 'w', 'Name', 'FLIM_bayes two-component fixed-lifetime separation');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    render_map(X, Y, maps.StdErrorTauUnknownNs, photonCounts, ...
        'Unknown lifetime std(error)', 'Std(error) [ns]');

    nexttile;
    render_map(X, Y, maps.RMSE_TauUnknownNs, photonCounts, ...
        'Unknown lifetime RMSE', 'RMSE [ns]');

    nexttile;
    render_map(X, Y, maps.StdErrorFractionUnknown, photonCounts, ...
        'Unknown fraction std(error)', 'Std(error)');

    nexttile;
    render_map(X, Y, maps.RMSE_FractionUnknown, photonCounts, ...
        'Unknown fraction RMSE', 'RMSE');
    plot_lifetime_deviation_summary(results);
    plot_height_summary(results);
    plot_condition_curves(results);
end

function render_map(X, Y, Z, photonCounts, titleStr, colorbarStr)
    surface(X, Y, zeros(size(Z)), Z, 'EdgeColor', 'none');
    view(2);
    set(gca, 'XScale', 'log', 'XTick', photonCounts);
    xlabel('Total photons per decay');
    ylabel('True unknown lifetime (ns)');
    title(titleStr);
    grid on;
    box on;
    cb = colorbar;
    cb.Label.String = colorbarStr;
end

function plot_condition_curves(results)
    tauIdx = choose_plot_indices(numel(results.trueUnknownTauAxisNs), results.optsUsed.curvePlotTauIdx, results.optsUsed.curvePlotMaxTau);
    photIdx = choose_plot_indices(numel(results.photonAxis), results.optsUsed.curvePlotPhotonIdx, results.optsUsed.curvePlotMaxPhoton);
    nTau = numel(tauIdx);
    nPhot = numel(photIdx);
    tAxisNs = ((0:numel(results.irf)-1)') * results.optsUsed.dtNs;

    figure('Color', 'w', 'Name', 'Fixed-component two-state mean data and fit curves', ...
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
            title(sprintf('\\tau_u=%.2f ns | h=%.1f nm | N=%d', ...
                results.trueUnknownTauAxisNs(tauHere), results.trueUnknownHeightNm(tauHere), entry.totalPhotonCount), ...
                'FontSize', 9);
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

    figure('Color', 'w', 'Name', 'Fixed-component two-state height summary');
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    render_map(X, Y, maps.BiasHeightNm, photonCounts, 'Height bias: estimated - true', 'Bias (nm)');
    ylabel(gca, 'True height (nm)');

    nexttile;
    render_map(X, Y, maps.StdErrorHeightNm, photonCounts, 'Height std(error)', 'Std(error) (nm)');
    ylabel(gca, 'True height (nm)');

    nexttile;
    render_map(X, Y, maps.MeanEstimatedHeightNm, photonCounts, 'Mean estimated height', 'Height (nm)');
    ylabel(gca, 'True height (nm)');
end

function plot_lifetime_deviation_summary(results)
    maps = results.summaryMaps;
    photonCounts = maps.photonAxis;
    tauAxisNs = maps.tauAxisNs(:);
    [X, Y] = meshgrid(photonCounts, tauAxisNs);

    figure('Color', 'w', 'Name', 'Fixed-component two-state lifetime deviation');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    render_map(X, Y, maps.BiasTauUnknownNs, photonCounts, ...
        'Unknown lifetime deviation: estimated - true', 'Deviation (ns)');

    nexttile;
    hold on;
    colors = lines(numel(photonCounts));
    for photIdx = 1:numel(photonCounts)
        plot(tauAxisNs, maps.MeanEstimatedTauUnknownNs(:, photIdx), '-', ...
            'Color', colors(photIdx, :), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%d photons', photonCounts(photIdx)));
    end
    plot(tauAxisNs, tauAxisNs, 'k--', 'LineWidth', 1.2, 'DisplayName', 'True = fitted');
    xlabel('True unknown lifetime (ns)');
    ylabel('Mean fitted unknown lifetime (ns)');
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
