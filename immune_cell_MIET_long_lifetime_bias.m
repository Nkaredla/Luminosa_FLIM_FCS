function analysis = immune_cell_MIET_long_lifetime_bias(referenceResultMat, cfg)
%IMMUNE_CELL_MIET_LONG_LIFETIME_BIAS Monte-Carlo bias of the MIET tau2 estimate.
%
% analysis = immune_cell_MIET_long_lifetime_bias(referenceResultMat, cfg)
%
% The simulation uses the measured IRF, TCSPC timing and fixed SLB lifetime
% from an immune_cell_MIET result, generates fixed-SLB + one-long-lifetime
% photon histograms, and sends them through flim_bayes_fixed_slb with the same
% one-membrane-state, spatially regularized-amplitude inference used by
% run_batch_immune_cell_MIET.
%
% Three different biases are deliberately kept separate:
%   posteriorMeanBias  E[tau2 | M2] - true tau2 for every replicate
%   mapBias            lifetime-grid MAP - true tau2
%   displayBias        posterior-mean bias only among replicates passing the
%                      production P(M2) and expected-photon display gates
%
% The last quantity includes selection bias and is the closest analogue of the
% lifetimes visible in the exported cell maps.
%
% Important cfg fields (all optional)
%   longLifetimeNs          true long lifetimes (default below)
%   slbPhotonCounts         exact detected short-component photons
%   longPhotonCounts        exact detected long-component photons
%   replicates              Monte-Carlo replicates per condition (100)
%   slbPriorScale           calibration mean / true SLB count (1)
%   slbPriorRelativeStd     Gaussian SLB-amplitude prior width (0.20)
%   slbCountPriorNodes      amplitude-prior quadrature nodes (5)
%   posteriorThreshold      displayed-component threshold (0.8)
%   minExpectedPhotons      displayed long-component photons (10)
%   randomSeed              reproducible random seed (13)
%   outputDir               default <reference result>/long_lifetime_bias
%   writeFigures            export the four diagnostic heatmaps (true)
%   showFigures             leave figures open (false)
%   bayesOverrides          fields merged into the reference Bayes options
%
% With no input, the function uses the reference acquisition discussed during
% development. Pass another result MAT to use another measured IRF/tau_SLB.

    if nargin < 1 || isempty(referenceResultMat)
        referenceResultMat = fullfile( ...
            'D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_20260813_1', ...
            '_20260813-155036', ...
            'immune_cell_MIET_slb_regularized_one_membrane_grid48', ...
            'cfg_589c3e5b', 'immune_cell_MIET_640nm_red_analysis.mat');
    end
    if nargin < 2 || isempty(cfg)
        cfg = struct();
    end
    referenceResultMat = char(referenceResultMat);
    if ~isfile(referenceResultMat)
        error('immune_cell_MIET_long_lifetime_bias:ReferenceResult', ...
            'Reference result MAT does not exist: %s', referenceResultMat);
    end
    cfg = fillDefaults(cfg, referenceResultMat);
    if ~isfolder(cfg.outputDir)
        mkdir(cfg.outputDir);
    end
    fprintf('immune_cell_MIET_long_lifetime_bias: results go to %s\n', ...
        cfg.outputDir);

    loaded = load(referenceResultMat, 'result');
    if ~isfield(loaded, 'result') || ~isstruct(loaded.result)
        error('immune_cell_MIET_long_lifetime_bias:MissingResult', ...
            'The reference MAT must contain a result struct.');
    end
    reference = loaded.result;
    instrument = referenceInstrument(reference);
    fprintf(['immune_cell_MIET_long_lifetime_bias: tau_SLB %.4g ns; dt %.4g ns; ' ...
        '%d TCSPC bins; %d conditions x %d replicates.\n'], ...
        instrument.tauSlbNs, instrument.dtNs, numel(instrument.irf), ...
        conditionCount(cfg), cfg.replicates);

    rng(cfg.randomSeed, 'twister');
    [condition, replicate] = simulateHistograms(instrument, cfg);
    bayesOpts = productionBayesOptions(reference, instrument, cfg, replicate);
    fprintf(['immune_cell_MIET_long_lifetime_bias: fitting %d synthetic pixels ' ...
        'with the production one-membrane Bayesian estimator.\n'], ...
        size(replicate.tcspc, 1));
    bayes = flim_bayes_fixed_slb(replicate.tcspc, instrument.irf, ...
        instrument.pulsePeriodNs, instrument.dtNs, instrument.tauSlbNs, ...
        bayesOpts);
    lifetimeGrid = double(bayes.membraneTauGridNs(:));
    fprintf(['immune_cell_MIET_long_lifetime_bias: lifetime grid %.3g--%.3g ns; ' ...
        'its weak-data prior mean is %.3g ns.\n'], min(lifetimeGrid), ...
        max(lifetimeGrid), mean(lifetimeGrid));

    [replicateTable, summaryTable] = measureBias( ...
        bayes, condition, replicate, instrument, cfg);

    csvFile = fullfile(cfg.outputDir, 'long_lifetime_bias_summary.csv');
    replicateCsvFile = fullfile(cfg.outputDir, ...
        'long_lifetime_bias_replicates.csv');
    matFile = fullfile(cfg.outputDir, 'long_lifetime_bias_analysis.mat');
    writetable(summaryTable, csvFile);
    writetable(replicateTable, replicateCsvFile);

    % Save all numerical results before invoking MATLAB graphics. This keeps
    % the bias analysis usable on headless systems with no display driver.
    figureFiles = {};
    analysis = struct();
    analysis.method = ['Monte-Carlo fixed-SLB + one-membrane TCSPC simulation ' ...
        'using the production Bayesian estimator'];
    analysis.referenceResultMat = referenceResultMat;
    analysis.instrument = instrument;
    analysis.config = cfg;
    analysis.estimator = struct('lifetimeGridNs', lifetimeGrid(:).', ...
        'weakDataPriorMeanLifetimeNs', mean(lifetimeGrid), ...
        'weakDataPriorMedianLifetimeNs', median(lifetimeGrid), ...
        'reportedEstimator', 'posterior mean conditional on M2', ...
        'displayPosteriorThreshold', cfg.posteriorThreshold, ...
        'displayMinimumExpectedLongPhotons', cfg.minExpectedPhotons);
    analysis.summary = summaryTable;
    analysis.replicates = replicateTable;
    analysis.outputFiles = struct('mat', matFile, 'summaryCsv', csvFile, ...
        'replicateCsv', replicateCsvFile, 'figures', {figureFiles});
    save(matFile, 'analysis', '-v7.3');

    if cfg.writeFigures
        try
            figureFiles = writeBiasFigures(summaryTable, instrument, cfg);
            analysis.outputFiles.figures = figureFiles;
            save(matFile, 'analysis', '-v7.3');
        catch figureError
            warning('immune_cell_MIET_long_lifetime_bias:FigureExport', ...
                ['Numerical results were saved, but figure export failed: %s. ' ...
                 'Set cfg.writeFigures=false on headless systems.'], ...
                figureError.message);
        end
    end

    fprintf(['immune_cell_MIET_long_lifetime_bias: complete. Summary CSV:\n' ...
        '  %s\n'], csvFile);
    printLargestBiases(summaryTable);
end

function cfg = fillDefaults(cfg, referenceResultMat)
    defaults = struct( ...
        'longLifetimeNs', [0.65 0.7 0.75 0.8 0.9 1.0 1.2 1.5 2.0 2.5 3.0 4.0], ...
        'slbPhotonCounts', [100 250 1000 4000 8000], ...
        'longPhotonCounts', [25 50 100 250 500 1000 2500 5000], ...
        'replicates', 100, ...
        'slbPriorScale', 1, ...
        'slbPriorRelativeStd', 0.20, ...
        'slbCountPriorNodes', 5, ...
        'slbCountRelTol', 0.0025, ...
        'posteriorThreshold', 0.8, ...
        'minExpectedPhotons', 10, ...
        'randomSeed', 13, ...
        'outputDir', fullfile(fileparts(referenceResultMat), ...
            'long_lifetime_bias'), ...
        'writeFigures', true, ...
        'showFigures', false, ...
        'bayesOverrides', struct());
    cfg = mergeStruct(defaults, cfg);
    cfg.longLifetimeNs = positiveVector(cfg.longLifetimeNs, ...
        'cfg.longLifetimeNs', false);
    cfg.slbPhotonCounts = positiveVector(cfg.slbPhotonCounts, ...
        'cfg.slbPhotonCounts', true);
    cfg.longPhotonCounts = positiveVector(cfg.longPhotonCounts, ...
        'cfg.longPhotonCounts', true);
    cfg.slbPriorScale = positiveVector(cfg.slbPriorScale, ...
        'cfg.slbPriorScale', false);
    validateattributes(cfg.replicates, {'numeric'}, ...
        {'scalar','integer','positive','finite'});
    validateattributes(cfg.slbPriorRelativeStd, {'numeric'}, ...
        {'scalar','real','finite','positive','<',1});
    validateattributes(cfg.slbCountPriorNodes, {'numeric'}, ...
        {'scalar','integer','positive','<=',15});
    if mod(cfg.slbCountPriorNodes, 2) == 0
        error('immune_cell_MIET_long_lifetime_bias:PriorNodes', ...
            'cfg.slbCountPriorNodes must be odd.');
    end
    validateattributes(cfg.slbCountRelTol, {'numeric'}, ...
        {'scalar','real','finite','nonnegative','<',1});
    validateattributes(cfg.posteriorThreshold, {'numeric'}, ...
        {'scalar','real','finite','>=',0,'<=',1});
    validateattributes(cfg.minExpectedPhotons, {'numeric'}, ...
        {'scalar','real','finite','nonnegative'});
    validateattributes(cfg.randomSeed, {'numeric'}, ...
        {'scalar','integer','nonnegative','finite'});
    validateattributes(cfg.showFigures, {'numeric','logical'}, {'scalar'});
    validateattributes(cfg.writeFigures, {'numeric','logical'}, {'scalar'});
    cfg.showFigures = logical(cfg.showFigures);
    cfg.writeFigures = logical(cfg.writeFigures);
    if ~(isstruct(cfg.bayesOverrides) && isscalar(cfg.bayesOverrides))
        error('immune_cell_MIET_long_lifetime_bias:BayesOverrides', ...
            'cfg.bayesOverrides must be a scalar struct.');
    end
    cfg.outputDir = char(cfg.outputDir);
end

function value = positiveVector(value, name, integerRequired)
    attributes = {'real','finite','vector','nonempty','positive'};
    if integerRequired
        attributes{end + 1} = 'integer';
    end
    validateattributes(value, {'numeric'}, attributes, mfilename, name);
    value = unique(double(value(:)).', 'stable');
end

function instrument = referenceInstrument(result)
    requiredTop = {'irf','slbReference','channel','bayesian'};
    for k = 1:numel(requiredTop)
        if ~isfield(result, requiredTop{k})
            error('immune_cell_MIET_long_lifetime_bias:ReferenceSchema', ...
                'Reference result is missing result.%s.', requiredTop{k});
        end
    end
    instrument = struct();
    instrument.irf = double(result.irf.curve(:));
    instrument.tauSlbNs = double(result.slbReference.fixedLifetimeNs);
    instrument.dtNs = double(result.channel.dtNs);
    compact = result.bayesian.compact;
    if isfield(compact, 'pulsePeriodNs')
        instrument.pulsePeriodNs = double(compact.pulsePeriodNs);
    else
        instrument.pulsePeriodNs = numel(instrument.irf) * instrument.dtNs;
    end
    if isfield(compact, 'irfShiftBins')
        instrument.irfShiftBins = double(compact.irfShiftBins);
    else
        instrument.irfShiftBins = 0;
    end
    if isfield(compact, 'convolutionMethod')
        instrument.convolutionMethod = char(compact.convolutionMethod);
    else
        instrument.convolutionMethod = 'linear';
    end
    if any(~isfinite(instrument.irf)) || sum(instrument.irf) <= 0
        error('immune_cell_MIET_long_lifetime_bias:ReferenceIrf', ...
            'Reference IRF is invalid.');
    end
    instrument.irf = instrument.irf / sum(instrument.irf);
end

function count = conditionCount(cfg)
    count = numel(cfg.longLifetimeNs) * numel(cfg.slbPhotonCounts) * ...
        numel(cfg.longPhotonCounts) * numel(cfg.slbPriorScale);
end

function [condition, replicate] = simulateHistograms(instrument, cfg)
    nConditions = conditionCount(cfg);
    nReplicateRows = nConditions * cfg.replicates;
    nBins = numel(instrument.irf);
    tcspc = zeros(nReplicateRows, 1, nBins, 'uint32');
    conditionId = zeros(nReplicateRows, 1, 'uint32');
    tauTrue = zeros(nReplicateRows, 1);
    slbPhotons = zeros(nReplicateRows, 1);
    longPhotons = zeros(nReplicateRows, 1);
    priorScale = zeros(nReplicateRows, 1);
    priorMean = zeros(nReplicateRows, 1);
    priorStd = zeros(nReplicateRows, 1);

    shiftedIrf = circshift(instrument.irf, round(instrument.irfShiftBins));
    slbPattern = decayPattern(shiftedIrf, instrument.pulsePeriodNs, ...
        instrument.dtNs, instrument.tauSlbNs, ...
        instrument.convolutionMethod);
    longPatterns = zeros(nBins, numel(cfg.longLifetimeNs));
    for tauIndex = 1:numel(cfg.longLifetimeNs)
        longPatterns(:, tauIndex) = decayPattern(shiftedIrf, ...
            instrument.pulsePeriodNs, instrument.dtNs, ...
            cfg.longLifetimeNs(tauIndex), instrument.convolutionMethod);
    end

    condition = struct('id', zeros(nConditions, 1, 'uint32'), ...
        'tauTrueNs', zeros(nConditions, 1), ...
        'slbPhotons', zeros(nConditions, 1), ...
        'longPhotons', zeros(nConditions, 1), ...
        'slbPriorScale', zeros(nConditions, 1), ...
        'crlbTauStdKnownFractionNs', zeros(nConditions, 1), ...
        'crlbTauStdUnknownFractionNs', zeros(nConditions, 1));
    row = 0;
    id = 0;
    for scaleIndex = 1:numel(cfg.slbPriorScale)
        scale = cfg.slbPriorScale(scaleIndex);
        for slbIndex = 1:numel(cfg.slbPhotonCounts)
            nSlb = cfg.slbPhotonCounts(slbIndex);
            for longIndex = 1:numel(cfg.longPhotonCounts)
                nLong = cfg.longPhotonCounts(longIndex);
                for tauIndex = 1:numel(cfg.longLifetimeNs)
                    id = id + 1;
                    tau = cfg.longLifetimeNs(tauIndex);
                    condition.id(id) = uint32(id);
                    condition.tauTrueNs(id) = tau;
                    condition.slbPhotons(id) = nSlb;
                    condition.longPhotons(id) = nLong;
                    condition.slbPriorScale(id) = scale;
                    [condition.crlbTauStdKnownFractionNs(id), ...
                        condition.crlbTauStdUnknownFractionNs(id)] = ...
                        lifetimeCrlb(slbPattern, longPatterns(:, tauIndex), ...
                        shiftedIrf, instrument, tau, nSlb, nLong);
                    for repetition = 1:cfg.replicates
                        row = row + 1;
                        histogram = sampleHistogram(slbPattern, nSlb) + ...
                            sampleHistogram(longPatterns(:, tauIndex), nLong);
                        tcspc(row, 1, :) = uint32(histogram);
                        conditionId(row) = uint32(id);
                        tauTrue(row) = tau;
                        slbPhotons(row) = nSlb;
                        longPhotons(row) = nLong;
                        priorScale(row) = scale;
                        priorMean(row) = scale * nSlb;
                        priorStd(row) = cfg.slbPriorRelativeStd * priorMean(row);
                    end
                end
            end
        end
    end
    replicate = struct('tcspc', tcspc, 'conditionId', conditionId, ...
        'tauTrueNs', tauTrue, 'slbPhotons', slbPhotons, ...
        'longPhotons', longPhotons, 'slbPriorScale', priorScale, ...
        'slbPriorMeanPhotons', priorMean, 'slbPriorStdPhotons', priorStd);
end

function opts = productionBayesOptions(reference, instrument, cfg, replicate)
    if isfield(reference, 'config') && isfield(reference.config, 'bayes') && ...
            isstruct(reference.config.bayes)
        opts = reference.config.bayes;
    else
        opts = struct();
    end
    opts = mergeStruct(opts, cfg.bayesOverrides);
    opts.analysisMask = true(size(replicate.tcspc, 1), 1);
    opts.minPhotons = 1;
    opts.useGPU = false;
    opts.batchSize = 2048;
    opts.includeBackground = true;
    opts.maxMembraneStates = 1;
    opts.fixedSlbPhotonCount = replicate.slbPriorMeanPhotons;
    opts.fixedSlbPhotonCountStd = replicate.slbPriorStdPhotons;
    opts.slbCountPriorNodes = cfg.slbCountPriorNodes;
    opts.slbCountRelTol = cfg.slbCountRelTol;
    opts.irfShiftBins = instrument.irfShiftBins;
    opts.convolutionMethod = instrument.convolutionMethod;
    if ~isfield(opts, 'membraneTauCount') || isempty(opts.membraneTauCount)
        opts.membraneTauCount = 48;
    end
    if ~isfield(opts, 'fractionStep') || isempty(opts.fractionStep)
        opts.fractionStep = 0.2;
    end
    if ~isfield(opts, 'minimumMembraneFraction') || ...
            isempty(opts.minimumMembraneFraction)
        opts.minimumMembraneFraction = 0.1;
    end
end

function [replicateTable, summaryTable] = measureBias( ...
        bayes, condition, replicate, instrument, cfg)
    posteriorMean = double(bayes.oneMembrane.membraneLifetime1Ns(:));
    lifetimeMap = double(bayes.oneMembrane.membraneLifetime1MAPNs(:));
    posteriorStd = double(bayes.oneMembrane.membraneLifetime1StdNs(:));
    probabilityM2 = double(bayes.probabilityBiexponential(:));
    totalPhotons = double(bayes.intensity(:));
    expectedLongPhotons = totalPhotons .* ...
        double(bayes.membrane1PhotonFraction(:));
    displaySelected = probabilityM2 >= cfg.posteriorThreshold & ...
        expectedLongPhotons >= cfg.minExpectedPhotons;
    mapSelectedM2 = double(bayes.completeExponentialCountMAP(:)) == 2;

    replicateTable = table(double(replicate.conditionId), ...
        replicate.tauTrueNs, replicate.slbPhotons, replicate.longPhotons, ...
        replicate.slbPriorScale, replicate.slbPriorMeanPhotons, ...
        replicate.slbPriorStdPhotons, posteriorMean, lifetimeMap, ...
        posteriorStd, probabilityM2, expectedLongPhotons, mapSelectedM2, ...
        displaySelected, ...
        'VariableNames', {'conditionId','trueLongLifetimeNs','slbPhotons', ...
        'longPhotons','slbPriorScale','slbPriorMeanPhotons', ...
        'slbPriorStdPhotons','posteriorMeanLifetimeNs','mapLifetimeNs', ...
        'posteriorStdLifetimeNs','probabilityTwoComponent', ...
        'posteriorExpectedLongPhotons','mapSelectedTwoComponent', ...
        'displaySelected'});

    nConditions = numel(condition.id);
    meanEstimate = nan(nConditions, 1);
    medianEstimate = nan(nConditions, 1);
    meanBias = nan(nConditions, 1);
    relativeBias = nan(nConditions, 1);
    rmse = nan(nConditions, 1);
    meanMapEstimate = nan(nConditions, 1);
    mapBias = nan(nConditions, 1);
    mapRelativeBias = nan(nConditions, 1);
    meanPosteriorStd = nan(nConditions, 1);
    meanProbabilityM2 = nan(nConditions, 1);
    meanExpectedLongPhotons = nan(nConditions, 1);
    medianExpectedLongPhotons = nan(nConditions, 1);
    expectedLongPhotonRelativeError = nan(nConditions, 1);
    mapSelectionRate = nan(nConditions, 1);
    displaySelectionRate = nan(nConditions, 1);
    displayMeanEstimate = nan(nConditions, 1);
    displayBias = nan(nConditions, 1);
    displayRelativeBias = nan(nConditions, 1);
    overestimateFraction = nan(nConditions, 1);
    displayOverestimateFraction = nan(nConditions, 1);
    posteriorEmpiricalStd = nan(nConditions, 1);
    posteriorBiasMonteCarloSE = nan(nConditions, 1);
    posteriorBiasCi95Low = nan(nConditions, 1);
    posteriorBiasCi95High = nan(nConditions, 1);
    displaySelectedCount = zeros(nConditions, 1);
    displayEmpiricalStd = nan(nConditions, 1);
    displayBiasMonteCarloSE = nan(nConditions, 1);
    displayBiasCi95Low = nan(nConditions, 1);
    displayBiasCi95High = nan(nConditions, 1);
    for id = 1:nConditions
        use = double(replicate.conditionId) == id;
        truth = condition.tauTrueNs(id);
        estimates = posteriorMean(use);
        mapEstimates = lifetimeMap(use);
        meanEstimate(id) = finiteMean(estimates);
        medianEstimate(id) = finiteMedian(estimates);
        meanBias(id) = meanEstimate(id) - truth;
        relativeBias(id) = 100 * meanBias(id) / truth;
        rmse(id) = sqrt(finiteMean((estimates - truth) .^ 2));
        finiteAll = isfinite(estimates);
        posteriorEmpiricalStd(id) = finiteStd(estimates);
        if any(finiteAll)
            posteriorBiasMonteCarloSE(id) = posteriorEmpiricalStd(id) / ...
                sqrt(nnz(finiteAll));
            posteriorBiasCi95Low(id) = meanBias(id) - ...
                1.96 * posteriorBiasMonteCarloSE(id);
            posteriorBiasCi95High(id) = meanBias(id) + ...
                1.96 * posteriorBiasMonteCarloSE(id);
        end
        meanMapEstimate(id) = finiteMean(mapEstimates);
        mapBias(id) = meanMapEstimate(id) - truth;
        mapRelativeBias(id) = 100 * mapBias(id) / truth;
        meanPosteriorStd(id) = finiteMean(posteriorStd(use));
        meanProbabilityM2(id) = finiteMean(probabilityM2(use));
        meanExpectedLongPhotons(id) = finiteMean(expectedLongPhotons(use));
        medianExpectedLongPhotons(id) = finiteMedian(expectedLongPhotons(use));
        expectedLongPhotonRelativeError(id) = 100 * ...
            (meanExpectedLongPhotons(id) - condition.longPhotons(id)) / ...
            condition.longPhotons(id);
        mapSelectionRate(id) = mean(mapSelectedM2(use));
        selected = use & displaySelected;
        selectedEstimates = posteriorMean(selected);
        selectedEstimates = selectedEstimates(isfinite(selectedEstimates));
        displaySelectedCount(id) = numel(selectedEstimates);
        displaySelectionRate(id) = nnz(selected) / nnz(use);
        displayMeanEstimate(id) = finiteMean(selectedEstimates);
        displayBias(id) = displayMeanEstimate(id) - truth;
        displayRelativeBias(id) = 100 * displayBias(id) / truth;
        displayEmpiricalStd(id) = finiteStd(selectedEstimates);
        if displaySelectedCount(id) > 0
            displayBiasMonteCarloSE(id) = displayEmpiricalStd(id) / ...
                sqrt(displaySelectedCount(id));
            displayBiasCi95Low(id) = displayBias(id) - ...
                1.96 * displayBiasMonteCarloSE(id);
            displayBiasCi95High(id) = displayBias(id) + ...
                1.96 * displayBiasMonteCarloSE(id);
        end
        if any(finiteAll)
            overestimateFraction(id) = mean(estimates(finiteAll) > truth);
        end
        if ~isempty(selectedEstimates)
            displayOverestimateFraction(id) = ...
                mean(selectedEstimates > truth);
        end
    end

    summaryTable = table(double(condition.id), condition.tauTrueNs, ...
        condition.slbPhotons, condition.longPhotons, ...
        condition.longPhotons ./ (condition.slbPhotons + condition.longPhotons), ...
        condition.slbPriorScale, repmat(cfg.replicates, nConditions, 1), ...
        meanEstimate, medianEstimate, meanBias, relativeBias, rmse, ...
        posteriorEmpiricalStd, posteriorBiasMonteCarloSE, ...
        posteriorBiasCi95Low, posteriorBiasCi95High, ...
        meanMapEstimate, mapBias, mapRelativeBias, meanPosteriorStd, ...
        meanProbabilityM2, meanExpectedLongPhotons, ...
        medianExpectedLongPhotons, expectedLongPhotonRelativeError, ...
        mapSelectionRate, displaySelectionRate, ...
        displayMeanEstimate, displayBias, displayRelativeBias, ...
        displaySelectedCount, displayEmpiricalStd, ...
        displayBiasMonteCarloSE, displayBiasCi95Low, displayBiasCi95High, ...
        overestimateFraction, displayOverestimateFraction, ...
        condition.crlbTauStdKnownFractionNs, ...
        condition.crlbTauStdUnknownFractionNs, ...
        'VariableNames', {'conditionId','trueLongLifetimeNs','slbPhotons', ...
        'longPhotons','trueLongPhotonFraction','slbPriorScale','replicates', ...
        'posteriorMeanEstimateNs','posteriorMedianEstimateNs', ...
        'posteriorMeanBiasNs','posteriorRelativeBiasPercent','posteriorRmseNs', ...
        'posteriorEmpiricalStdNs','posteriorBiasMonteCarloSENs', ...
        'posteriorBiasCi95LowNs','posteriorBiasCi95HighNs', ...
        'meanMapEstimateNs','mapBiasNs','mapRelativeBiasPercent', ...
        'meanPosteriorStdNs','meanProbabilityTwoComponent', ...
        'meanPosteriorExpectedLongPhotons', ...
        'medianPosteriorExpectedLongPhotons', ...
        'expectedLongPhotonRelativeErrorPercent', ...
        'mapTwoComponentRate','displaySelectionRate','displayMeanEstimateNs', ...
        'displayBiasNs','displayRelativeBiasPercent','displaySelectedCount', ...
        'displayEmpiricalStdNs','displayBiasMonteCarloSENs', ...
        'displayBiasCi95LowNs', ...
        'displayBiasCi95HighNs', ...
        'overestimateFraction','displayOverestimateFraction', ...
        'crlbTauStdKnownFractionNs', ...
        'crlbTauStdUnknownFractionNs'});
    summaryTable.fixedSlbLifetimeNs = ...
        repmat(instrument.tauSlbNs, nConditions, 1);
end

function files = writeBiasFigures(summary, instrument, cfg)
    files = { ...
        fullfile(cfg.outputDir, 'posterior_mean_relative_bias.png'), ...
        fullfile(cfg.outputDir, 'display_selected_relative_bias.png'), ...
        fullfile(cfg.outputDir, 'map_relative_bias.png'), ...
        fullfile(cfg.outputDir, 'display_selection_rate.png')};
    makeHeatmapFigure(summary, cfg, 'posteriorRelativeBiasPercent', ...
        'Posterior-mean lifetime bias [%]', true, instrument, files{1});
    makeHeatmapFigure(summary, cfg, 'displayRelativeBiasPercent', ...
        'Displayed-pixel lifetime bias [%]', true, instrument, files{2});
    makeHeatmapFigure(summary, cfg, 'mapRelativeBiasPercent', ...
        'Lifetime-grid MAP bias [%]', true, instrument, files{3});
    makeHeatmapFigure(summary, cfg, 'displaySelectionRate', ...
        'Fraction passing display gates', false, instrument, files{4});
end

function makeHeatmapFigure(summary, cfg, fieldName, colourLabel, ...
        symmetricColour, instrument, outputFile)
    [~, scaleIndex] = min(abs(cfg.slbPriorScale - 1));
    displayedScale = cfg.slbPriorScale(scaleIndex);
    useScale = summary.slbPriorScale == displayedScale;
    valuesForLimits = double(summary.(fieldName)(useScale));
    valuesForLimits = valuesForLimits(isfinite(valuesForLimits));
    visibility = 'off';
    if cfg.showFigures
        visibility = 'on';
    end
    h = figure('Color', 'w', 'Visible', visibility, ...
        'Position', [60 60 420 * numel(cfg.slbPhotonCounts) 520]);
    layout = tiledlayout(h, 1, numel(cfg.slbPhotonCounts), ...
        'Padding', 'compact', 'TileSpacing', 'compact');
    for slbIndex = 1:numel(cfg.slbPhotonCounts)
        nSlb = cfg.slbPhotonCounts(slbIndex);
        matrix = nan(numel(cfg.longPhotonCounts), numel(cfg.longLifetimeNs));
        for longIndex = 1:numel(cfg.longPhotonCounts)
            for tauIndex = 1:numel(cfg.longLifetimeNs)
                match = useScale & summary.slbPhotons == nSlb & ...
                    summary.longPhotons == cfg.longPhotonCounts(longIndex) & ...
                    summary.trueLongLifetimeNs == cfg.longLifetimeNs(tauIndex);
                if nnz(match) == 1
                    matrix(longIndex, tauIndex) = summary.(fieldName)(match);
                end
            end
        end
        ax = nexttile(layout);
        object = imagesc(ax, matrix);
        set(object, 'AlphaData', isfinite(matrix));
        set(ax, 'YDir', 'normal', 'Color', [0.85 0.85 0.85], ...
            'XTick', 1:numel(cfg.longLifetimeNs), ...
            'XTickLabel', compose('%g', cfg.longLifetimeNs), ...
            'YTick', 1:numel(cfg.longPhotonCounts), ...
            'YTickLabel', compose('%g', cfg.longPhotonCounts));
        xlabel(ax, 'True long lifetime [ns]');
        ylabel(ax, 'Long-component photons');
        title(ax, sprintf('SLB photons = %g', nSlb));
        colourBar = colorbar(ax);
        colourBar.Label.String = colourLabel;
        if symmetricColour
            limit = finitePercentile(abs(valuesForLimits), 0.98);
            limit = max(limit, 1);
            clim(ax, [-limit limit]);
            colormap(ax, blueWhiteRed(256));
        else
            clim(ax, [0 1]);
            colormap(ax, parula(256));
        end
    end
    title(layout, sprintf(['%s | tau_{SLB}=%.3g ns | SLB prior scale %.3g | ' ...
        '%d replicates'], colourLabel, instrument.tauSlbNs, displayedScale, ...
        cfg.replicates), 'FontWeight', 'bold');
    exportgraphics(h, outputFile, 'Resolution', 250);
    if ~cfg.showFigures
        close(h);
    end
end

function [knownFractionStd, unknownFractionStd] = lifetimeCrlb( ...
        slbPattern, longPattern, shiftedIrf, instrument, tau, nSlb, nLong)
    total = nSlb + nLong;
    fraction = nLong / total;
    probability = (1 - fraction) * slbPattern + fraction * longPattern;
    probability = max(probability, realmin('double'));
    step = max(1e-4, 1e-3 * tau);
    lowerPattern = decayPattern(shiftedIrf, instrument.pulsePeriodNs, ...
        instrument.dtNs, max(tau - step, eps), instrument.convolutionMethod);
    upperPattern = decayPattern(shiftedIrf, instrument.pulsePeriodNs, ...
        instrument.dtNs, tau + step, instrument.convolutionMethod);
    derivativeTau = fraction * (upperPattern - lowerPattern) / (2 * step);
    informationTau = total * sum((derivativeTau .^ 2) ./ probability);
    knownFractionStd = 1 / sqrt(max(informationTau, eps));

    derivativeFraction = longPattern - slbPattern;
    information = total * [ ...
        sum((derivativeTau .^ 2) ./ probability), ...
        sum((derivativeTau .* derivativeFraction) ./ probability); ...
        sum((derivativeTau .* derivativeFraction) ./ probability), ...
        sum((derivativeFraction .^ 2) ./ probability)];
    if rcond(information) < 1e-12
        unknownFractionStd = Inf;
    else
        covariance = inv(information);
        unknownFractionStd = sqrt(max(covariance(1, 1), 0));
    end
end

function histogram = sampleHistogram(pattern, photonCount)
    probability = max(double(pattern(:)), 1e-12);
    probability = probability / sum(probability);
    edges = [0; cumsum(probability)];
    edges(end) = 1;
    histogram = histcounts(rand(photonCount, 1), edges).';
end

function pattern = decayPattern(irf, periodNs, dtNs, tauNs, method)
    timeNs = (0:numel(irf)-1).' * dtNs;
    decay = exp(-timeNs / tauNs) / max(1 - exp(-periodNs / tauNs), eps);
    if strcmp(method, 'gui') && exist('Convol', 'file') == 2
        convolved = double(Convol(irf, decay));
        convolved = convolved(:);
    else
        convolved = conv(irf, decay, 'full');
    end
    convolved = convolved(1:numel(irf));
    pattern = max(real(convolved), 0);
    pattern = pattern / sum(pattern);
end

function value = finiteMean(data)
    data = double(data(isfinite(data)));
    if isempty(data)
        value = NaN;
    else
        value = mean(data);
    end
end

function value = finiteMedian(data)
    data = double(data(isfinite(data)));
    if isempty(data)
        value = NaN;
    else
        value = median(data);
    end
end

function value = finiteStd(data)
    data = double(data(isfinite(data)));
    if isempty(data)
        value = NaN;
    else
        value = std(data);
    end
end

function value = finitePercentile(data, probability)
    data = sort(double(data(isfinite(data))));
    if isempty(data)
        value = NaN;
        return;
    end
    position = 1 + (numel(data) - 1) * min(max(probability, 0), 1);
    low = floor(position);
    high = ceil(position);
    fraction = position - low;
    value = data(low) * (1 - fraction) + data(high) * fraction;
end

function map = blueWhiteRed(count)
    x = linspace(-1, 1, count).';
    red = interp1([-1 0 1], [0.10 1.00 0.80], x);
    green = interp1([-1 0 1], [0.30 1.00 0.10], x);
    blue = interp1([-1 0 1], [0.80 1.00 0.10], x);
    map = min(max([red green blue], 0), 1);
end

function printLargestBiases(summary)
    usable = isfinite(summary.displayRelativeBiasPercent) & ...
        summary.displaySelectionRate > 0;
    if any(usable)
        candidates = find(usable);
        [~, local] = max(abs(summary.displayRelativeBiasPercent(usable)));
        row = candidates(local);
        fprintf(['  Largest absolute displayed-pixel bias: %+.1f%% at tau %.3g ns, ' ...
            '%g SLB + %g long photons (selection rate %.1f%%).\n'], ...
            summary.displayRelativeBiasPercent(row), ...
            summary.trueLongLifetimeNs(row), summary.slbPhotons(row), ...
            summary.longPhotons(row), 100 * summary.displaySelectionRate(row));
    end
end

function out = mergeStruct(base, override)
    out = base;
    names = fieldnames(override);
    for index = 1:numel(names)
        out.(names{index}) = override.(names{index});
    end
end
