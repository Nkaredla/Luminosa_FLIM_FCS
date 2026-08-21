function [locs, modelOrder] = fitFocusedISMFrame(data, detector, opts)
    import membrane_tracking.focused_ism.*

    % K=0 is a shared-background-only microimage. Positive-K models add
    % independent monomer positions and photon amplitudes.
    backgroundOnly = max(mean(data), opts.minExpectedCount);
    mu0 = backgroundOnly * ones(detector.nChannels, 1);
    nll0 = poissonImageNll(data, mu0, opts);
    % A spatial Poisson process gives a Poisson prior for the number of
    % molecules in the effective Gaussian excitation area.
    priorMean = max(opts.emitterCountPriorConcentrationPerUm2 * ...
        detector.effectiveExcitationAreaUm2, 0.05);

    signalEstimate = max(sum(data - detector.backgroundMeanCounts), 0);
    photonSupportedK = max(1, floor(signalEstimate / ...
        max(opts.minPhotonsPerLocalization, 1)) + 1);
    maxK = min(opts.maxEmittersPerFitRegion, photonSupportedK);
    models = cell(maxK + 1, 1);
    scores = inf(maxK + 1, 1);
    kValues = (0:maxK).';
    models{1} = struct('nEmitters', 0, 'background', backgroundOnly, ...
        'mu', mu0, 'negLogLikelihood', nll0);
    scores(1) = bayesianModelOrderScore(nll0, 1, numel(data), ...
        0, priorMean);

    previousModel = [];
    seeds = zeros(0, 2);
    consecutiveFailures = 0;
    for k = 1:maxK
        if isempty(previousModel)
            residualMean = mu0;
        else
            residualMean = previousModel.mu;
        end
        % Build K+1 from the strongest spatial pattern left in the residual.
        % Warm-starting from K reduces label switching and optimizer cost.
        newSeed = focusedISMResidualSeed(data, residualMean, seeds, ...
            detector, opts);
        if isempty(newSeed)
            break;
        end
        seeds(end+1,:) = newSeed; %#ok<AGROW>

        [model, ok] = fitFocusedISMModel(data, seeds, detector, ...
            opts, previousModel);
        if ~ok
            % Failure at one order does not prove that every higher order is
            % absent. Retain the distinct seed and allow a bounded number of
            % higher-order attempts before stopping the ladder.
            consecutiveFailures = consecutiveFailures + 1;
            if consecutiveFailures > opts.maxConsecutiveFitFailures
                break;
            end
            continue;
        end
        consecutiveFailures = 0;
        models{k+1} = model;
        scores(k+1) = bayesianModelOrderScore( ...
            model.negLogLikelihood, 3*k + 1, numel(data), k, priorMean);
        previousModel = model;
    end

    finiteModels = isfinite(scores);
    finiteScores = scores(finiteModels);
    finiteK = kValues(finiteModels);
    finiteCells = models(finiteModels);
    if isempty(finiteScores)
        locs = struct([]);
        modelOrder = emptyModelOrderResult();
        return;
    end

    % The score is -2 log approximate evidence, so exp(-score/2) produces
    % normalized approximate posterior model probabilities.
    logWeights = -0.5 * (finiteScores - min(finiteScores));
    weights = exp(logWeights);
    weights = weights / sum(weights);
    [~, posteriorOrder] = sort(weights, 'descend');

    % Keep the MAP emitter count even when its joint position Fisher matrix
    % is singular. Falling back to a lower K turns close emitters into a
    % photon-weighted centroid and biases the apparent motion. Singular MAP
    % positions remain linking waypoints but are explicitly excluded from D.
    bestIndex = posteriorOrder(1);
    selectedModel = finiteCells{bestIndex};
    selectedProbability = weights(bestIndex);
    selectedK = selectedModel.nEmitters;
    selectedFisher = struct();
    positionsIdentifiable = true;
    if selectedK > 0
        selectedFisher = focusedISMFisher(selectedModel, detector, opts);
        positionsIdentifiable = selectedFisher.isValid;
    end

    modelOrder = struct('nEmitters', selectedK, ...
        'posteriorProbability', selectedProbability, ...
        'entropy', -sum(weights .* log(max(weights, realmin))), ...
        'candidateEmitterCounts', finiteK, 'posteriorWeights', weights, ...
        'positionsIdentifiable', positionsIdentifiable);
    if selectedK == 0
        locs = struct([]);
        return;
    end

    locs = repmat(struct(), selectedK, 1);
    for e = 1:selectedK
        if positionsIdentifiable
            covarianceXY = selectedFisher.covarianceXY(:,:,e);
        else
            % This covariance is deliberately conservative and is used only
            % by the linker. positionIdentifiable prevents it entering MSD
            % or step-likelihood inference.
            photonLimitedVariance = (opts.psfSigmaUm^2 + ...
                detector.detectorPitchUm^2 / 12) / ...
                max(selectedModel.amplitude(e), 1);
            covarianceXY = 4 * photonLimitedVariance * eye(2);
        end
        locs(e).xUm = selectedModel.xUm(e);
        locs(e).yUm = selectedModel.yUm(e);
        locs(e).signalPhotons = selectedModel.amplitude(e);
        locs(e).backgroundCountsPerDetector = selectedModel.background;
        locs(e).crbVarXUm2 = covarianceXY(1,1);
        locs(e).crbVarYUm2 = covarianceXY(2,2);
        locs(e).crbCovXYUm2 = covarianceXY(1,2);
        locs(e).crbSigmaUm = sqrt(max(0.5 * trace(covarianceXY), 0));
        locs(e).negLogLikelihood = selectedModel.negLogLikelihood;
        locs(e).positionIdentifiable = positionsIdentifiable;
    end
end
