function diffusion = estimateSurfaceDiffusion(trackTable, trajectories, ...
        surfaceModel, modes, opts)
    import membrane_tracking.fluctuating_miet.*

    if height(trackTable) == 0 || ~surfaceModel.fitSucceeded
        diffusion = emptyDiffusionResult(opts);
        diffusion.truthValidation = truthDiffusionDiagnostics( ...
            trajectories, trackTable, modes, opts);
        return;
    end
    curvature = surfaceModel.curvaturePerUm;
    usable = rowsUsableForDiffusion(trackTable, opts);
    trackIds = unique(trackTable.trackId).';
    usedTrackIds = zeros(0,1);
    for trackId = trackIds
        if sum(usable & trackTable.trackId == trackId) >= opts.minTrackLength
            usedTrackIds(end+1,1) = trackId; %#ok<AGROW>
        end
    end
    if isempty(usedTrackIds)
        diffusion = emptyDiffusionResult(opts);
        diffusion.nTracksTotal = numel(trackIds);
        diffusion.truthValidation = truthDiffusionDiagnostics( ...
            trajectories, trackTable, modes, opts);
        return;
    end

    lagFrames = (1:opts.maxMsdLag).';
    lagTimeS = lagFrames * opts.dtS;
    nPairs = zeros(opts.maxMsdLag,1);
    msd = nan(opts.maxMsdLag,1);
    noise = nan(opts.maxMsdLag,1);
    msdCorrected = nan(opts.maxMsdLag,1);
    for lag = 1:opts.maxMsdLag
        squared = zeros(0,1);
        noiseTrace = zeros(0,1);
        for k = 1:numel(usedTrackIds)
            rows = find(trackTable.trackId == usedTrackIds(k) & usable);
            [~, order] = sort(trackTable.frame(rows));
            rows = rows(order);
            for i = 1:numel(rows)
                target = trackTable.frame(rows(i)) + lag;
                j = find(trackTable.frame(rows) == target, 1, 'first');
                if isempty(j)
                    continue;
                end
                [step, stepNoise, ok] = tangentStep(trackTable, rows(i), ...
                    rows(j), surfaceModel, modes, opts);
                if ok
                    squared(end+1,1) = sum(step.^2); %#ok<AGROW>
                    noiseTrace(end+1,1) = trace(stepNoise); %#ok<AGROW>
                end
            end
        end
        nPairs(lag) = numel(squared);
        if nPairs(lag) >= opts.minMsdPairs
            msd(lag) = mean(squared);
            noise(lag) = mean(noiseTrace);
            msdCorrected(lag) = msd(lag) - noise(lag);
        end
    end
    fitMask = isfinite(msdCorrected) & lagFrames <= opts.maxMsdFitLag;
    Dmsd = fitMsdThroughOrigin(lagTimeS(fitMask), msdCorrected(fitMask), ...
        nPairs(fitMask));

    [stepXY, stepNoise, stepDt, driftTerm, gateRadialLimit] = ...
        collectIndependentSteps(trackTable, usable, usedTrackIds, ...
        surfaceModel, modes, opts);
    if ~opts.gateCensoringCorrection
        gateRadialLimit = zeros(0,opts.gateAcceptanceQuadraturePoints);
    end
    fisherEstimate = surfaceDiffusionMLE(stepXY, stepNoise, stepDt, ...
        driftTerm, gateRadialLimit, Dmsd, opts.dtS);

    diffusion = struct();
    diffusion.geometry = [opts.diffusionGeometry ...
        ' tangent metric from the estimated membrane'];
    diffusion.trueDiffusionUm2PerS = opts.diffusionUm2PerS;
    diffusion.curvatureUsedPerUm = curvature;
    diffusion.DMsdCorrectedUm2PerS = Dmsd;
    diffusion.DFisherMLEUm2PerS = fisherEstimate.diffusionUm2PerS;
    diffusion.diffusionFisherSigmaUm2PerS = fisherEstimate.sigmaUm2PerS;
    diffusion.diffusionFisherInformation = fisherEstimate.information;
    diffusion.meanGateAcceptanceProbability = fisherEstimate.meanAcceptance;
    diffusion.gateCensoringCorrectionApplied = opts.gateCensoringCorrection;
    diffusion.gateAcceptanceModel = ['Exact anisotropic Mahalanobis-gate ' ...
        'conditioning with Ito mean; focus selection, localisation ' ...
        'acceptance, and association mistakes are not conditioned.'];
    diffusion.nIndependentStepsForFisher = fisherEstimate.nSteps;
    diffusion.nLag1StepsForFisher = fisherEstimate.nLag1Steps;
    diffusion.msdTable = table(lagFrames, lagTimeS, nPairs, msd, noise, ...
        msdCorrected, 'VariableNames', {'lagFrames','lagTimeS','nPairs', ...
        'surfaceMsdUm2','meanFisherNoiseUm2','surfaceMsdCorrectedUm2'});
    diffusion.nTracksTotal = numel(trackIds);
    diffusion.nTracksUsed = numel(usedTrackIds);
    diffusion.truthValidation = truthDiffusionDiagnostics( ...
        trajectories, trackTable, modes, opts);
end
