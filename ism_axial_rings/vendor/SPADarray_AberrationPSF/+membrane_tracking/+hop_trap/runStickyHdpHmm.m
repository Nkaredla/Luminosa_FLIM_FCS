function result = runStickyHdpHmm(stepTable, freeModel, opts)
%RUNSTICKYHDPHMM Truncated sticky HDP-HMM posterior mobility-state sampler.

    import membrane_tracking.hop_trap.*

    result = emptyStickyHdpHmmResult(stepTable);
    if ~opts.runStickyBnp
        result.status = 'disabled';
        return;
    end
    if ~freeModel.estimateAvailable || ...
            height(stepTable) < opts.minimumInferenceSteps
        return;
    end

    n = height(stepTable);
    kMaximum = opts.stickyBnpMaximumStates;
    grid = logspace(-8, log10(max(100, ...
        20 * opts.freeDiffusionUm2PerS)), ...
        opts.stickyBnpDiffusionGridSize);
    logEmissionGrid = zeros(n, numel(grid));
    for g = 1:numel(grid)
        logEmissionGrid(:,g) = logGaussianStepDensity(stepTable, grid(g));
    end
    priorCenter = log(max(freeModel.diffusionUm2PerS, 1e-8));
    logGridPrior = -0.5 * ((log(grid) - priorCenter) / ...
        opts.stickyBnpLogDiffusionPriorSigma).^2;

    proxy = max((stepTable.dxUm.^2 + stepTable.dyUm.^2 - ...
        stepTable.noiseVarXUm2 - stepTable.noiseVarYUm2) ./ ...
        max(4 * stepTable.stepDtS, realmin), 1e-8);
    initialD = logspace(log10(max(min(proxy), 1e-8)), ...
        log10(max(max(proxy), freeModel.diffusionUm2PerS)), kMaximum);
    diffusionIndex = zeros(1, kMaximum);
    for state = 1:kMaximum
        [~, diffusionIndex(state)] = min(abs(log(grid) - ...
            log(initialD(state))));
    end
    [~, state] = min(abs(log(proxy) - log(grid(diffusionIndex))), [], 2);

    beta = ones(1, kMaximum) / kMaximum;
    transition = zeros(kMaximum);
    for row = 1:kMaximum
        concentration = opts.stickyBnpAlpha * beta;
        concentration(row) = concentration(row) + opts.stickyBnpKappa;
        transition(row,:) = concentration / sum(concentration);
    end
    nSaved = opts.stickyBnpSamples;
    stateCountTrace = zeros(nSaved, 1);
    diffusionTrace = nan(nSaved, kMaximum);
    occupancyTrace = zeros(nSaved, kMaximum);
    rankCounts = zeros(n, kMaximum);
    slowCounts = zeros(n, 1);
    totalSweeps = opts.stickyBnpBurnIn + ...
        opts.stickyBnpSamples * opts.stickyBnpThinning;
    saved = 0;

    for sweep = 1:totalSweeps
        transitionCounts = zeros(kMaximum);
        initialCounts = zeros(1, kMaximum);
        sequences = unique(stepTable.sequenceId, 'stable').';
        for sequence = sequences
            rows = find(stepTable.sequenceId == sequence);
            initialCounts(state(rows(1))) = ...
                initialCounts(state(rows(1))) + 1;
            for t = 1:numel(rows)-1
                i = state(rows(t));
                j = state(rows(t+1));
                transitionCounts(i,j) = transitionCounts(i,j) + 1;
            end
        end

        tableCounts = sampleHdpTableCounts(transitionCounts, beta, ...
            opts.stickyBnpAlpha, opts.stickyBnpKappa);
        beta = sampleDirichlet(opts.stickyBnpGamma / kMaximum + ...
            sum(tableCounts, 1));
        initialProbability = sampleDirichlet(beta + initialCounts);
        for row = 1:kMaximum
            concentration = opts.stickyBnpAlpha * beta + ...
                transitionCounts(row,:);
            concentration(row) = concentration(row) + opts.stickyBnpKappa;
            transition(row,:) = sampleDirichlet(concentration);
        end

        for k = 1:kMaximum
            assigned = state == k;
            logPosterior = logGridPrior;
            if any(assigned)
                logPosterior = logPosterior + ...
                    sum(logEmissionGrid(assigned,:), 1);
            end
            diffusionIndex(k) = sampleDiscreteLogWeights(logPosterior);
        end
        logEmission = logEmissionGrid(:,diffusionIndex);
        state = sampleHmmStatePaths(logEmission, stepTable.sequenceId, ...
            initialProbability, transition);

        if sweep > opts.stickyBnpBurnIn && ...
                mod(sweep - opts.stickyBnpBurnIn, ...
                opts.stickyBnpThinning) == 0
            saved = saved + 1;
            occupancy = accumarray(state, 1, [kMaximum 1]).' / n;
            active = find(occupancy >= ...
                opts.stickyBnpMinimumOccupancyFraction);
            if isempty(active)
                [~, active] = max(occupancy);
            end
            [~, order] = sort(grid(diffusionIndex(active)));
            active = active(order);
            stateCountTrace(saved) = numel(active);
            diffusionTrace(saved,1:numel(active)) = ...
                grid(diffusionIndex(active));
            occupancyTrace(saved,1:numel(active)) = occupancy(active);
            for rank = 1:numel(active)
                rankCounts(:,rank) = rankCounts(:,rank) + ...
                    double(state == active(rank));
            end
            if numel(active) >= 2
                slowCounts = slowCounts + double(state == active(1));
            end
        end
    end

    stateCountMode = mode(stateCountTrace);
    result.status = 'ok';
    result.estimateAvailable = true;
    result.method = 'finite sticky HDP-HMM approximation';
    result.maximumStates = kMaximum;
    result.nPosteriorSamples = nSaved;
    result.stateCountTrace = stateCountTrace;
    result.posteriorModeStateCount = stateCountMode;
    result.posteriorMeanStateCount = mean(stateCountTrace);
    result.diffusionTraceUm2PerS = diffusionTrace;
    result.occupancyTrace = occupancyTrace;
    result.posteriorRankProbability = rankCounts / nSaved;
    result.posteriorSlowStateProbability = slowCounts / nSaved;
    result.posteriorMedianDiffusionUm2PerS = ...
        median(diffusionTrace, 1, 'omitnan');
    result.posteriorMedianOccupancy = ...
        median(occupancyTrace, 1, 'omitnan');
    result.stateCountProbability = accumarray(stateCountTrace, 1, ...
        [kMaximum 1]).' / nSaved;
end
