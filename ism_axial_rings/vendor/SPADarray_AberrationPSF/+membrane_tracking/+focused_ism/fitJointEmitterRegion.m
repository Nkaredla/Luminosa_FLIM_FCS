function [locs, modelOrder] = fitJointEmitterRegion( ...
        image, candidates, detector, opts)
    import membrane_tracking.focused_ism.*

    radius = opts.fitWindowRadiusPx;
    rows = max(1, min(candidates(:,1)) - radius): ...
        min(detector.ny, max(candidates(:,1)) + radius);
    cols = max(1, min(candidates(:,2)) - radius): ...
        min(detector.nx, max(candidates(:,2)) + radius);
    data = double(image(rows, cols));

    [~, candidateOrder] = sort(candidates(:,3), 'descend');
    candidates = candidates(candidateOrder, :);
    seedX = detector.xCenters(candidates(:,2));
    seedY = detector.yCenters(candidates(:,1));
    seedXY = [seedX(:), seedY(:)];

    regionAreaUm2 = numel(rows) * numel(cols) * detector.pixelSizeUm^2;
    priorMean = max(opts.emitterCountPriorConcentrationPerUm2 * ...
        regionAreaUm2, 0.25);
    maxK = opts.maxEmittersPerFitRegion;

    models = cell(maxK + 1, 1);
    scores = inf(maxK + 1, 1);
    kValues = (0:maxK).';

    backgroundOnly = max(mean(data(:)), opts.minExpectedCount);
    mu0 = backgroundOnly * ones(size(data));
    nll0 = poissonImageNll(data, mu0, opts);
    models{1} = struct('nEmitters', 0, 'background', backgroundOnly, ...
        'mu', mu0, 'negLogLikelihood', nll0);
    scores(1) = bayesianModelOrderScore(nll0, 1, numel(data), ...
        0, priorMean);

    previousModel = [];
    activeSeeds = zeros(0, 2);
    for k = 1:maxK
        if k <= size(seedXY, 1)
            newSeed = seedXY(k,:);
        elseif ~isempty(previousModel)
            newSeed = residualEmitterSeed(data, previousModel.mu, ...
                rows, cols, activeSeeds, detector, opts);
            if isempty(newSeed)
                break;
            end
        else
            break;
        end
        activeSeeds(end+1, :) = newSeed; %#ok<AGROW>

        [model, ok] = fitJointGaussianModel(data, rows, cols, ...
            activeSeeds, detector, opts, previousModel);
        if ~ok
            break;
        end
        modelIndex = k + 1;
        models{modelIndex} = model;
        parameterCount = 3 * k + 1;
        scores(modelIndex) = bayesianModelOrderScore( ...
            model.negLogLikelihood, parameterCount, numel(data), ...
            k, priorMean);
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

    logWeights = -0.5 * (finiteScores - min(finiteScores));
    weights = exp(logWeights);
    weights = weights / sum(weights);
    [~, posteriorOrder] = sort(weights, 'descend');

    selectedModel = [];
    selectedProbability = NaN;
    selectedK = 0;
    selectedFisher = struct();
    for index = posteriorOrder(:).'
        candidateModel = finiteCells{index};
        if candidateModel.nEmitters == 0
            selectedModel = candidateModel;
            selectedProbability = weights(index);
            selectedK = 0;
            break;
        end
        fisher = multiEmitterLocalizationFisher(candidateModel, ...
            rows, cols, detector, opts);
        if fisher.isValid
            selectedModel = candidateModel;
            selectedProbability = weights(index);
            selectedK = candidateModel.nEmitters;
            selectedFisher = fisher;
            break;
        end
    end

    entropyValue = -sum(weights .* log(max(weights, realmin)));
    modelOrder = struct('nEmitters', selectedK, ...
        'posteriorProbability', selectedProbability, ...
        'entropy', entropyValue, 'candidateEmitterCounts', finiteK, ...
        'posteriorWeights', weights);

    if isempty(selectedModel) || selectedK == 0
        locs = struct([]);
        return;
    end

    locs = repmat(struct(), selectedK, 1);
    for e = 1:selectedK
        covarianceXY = selectedFisher.covarianceXY(:,:,e);
        probabilityInWindow = sum(selectedModel.probability(:,:,e), 'all');
        signalPhotons = selectedModel.amplitude(e) * probabilityInWindow;
        locs(e).xUm = selectedModel.xUm(e);
        locs(e).yUm = selectedModel.yUm(e);
        locs(e).xPixel = (selectedModel.xUm(e) - detector.xEdges(1)) / ...
            detector.pixelSizeUm + 0.5;
        locs(e).yPixel = (selectedModel.yUm(e) - detector.yEdges(1)) / ...
            detector.pixelSizeUm + 0.5;
        locs(e).signalPhotons = signalPhotons;
        locs(e).windowCounts = sum(data(:));
        locs(e).backgroundCountsPerPixel = selectedModel.background;
        locs(e).crbVarXUm2 = covarianceXY(1,1);
        locs(e).crbVarYUm2 = covarianceXY(2,2);
        locs(e).crbCovXYUm2 = covarianceXY(1,2);
        locs(e).crbSigmaUm = sqrt(max(0.5 * trace(covarianceXY), 0));
        locs(e).negLogLikelihood = selectedModel.negLogLikelihood;
    end
end
