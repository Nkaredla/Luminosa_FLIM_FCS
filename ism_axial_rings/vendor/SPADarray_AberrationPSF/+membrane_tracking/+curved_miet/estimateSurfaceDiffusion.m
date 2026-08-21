function diffusion = estimateSurfaceDiffusion(trackTable, trajectories, ...
        curvature, opts)
    import membrane_tracking.curved_miet.*

    if height(trackTable) == 0 || ~curvature.fitSucceeded || ...
            ~isfinite(curvature.curvaturePerUm)
        diffusion = emptyDiffusionResult(opts);
        diffusion.truthValidation = truthDiffusionDiagnostics( ...
            trajectories, trackTable, opts.curvaturePerUm, opts);
        return;
    end

    reliable = trackTable.associationProbability >= ...
        opts.minAssociationProbabilityForDiffusion;
    allTrackIds = unique(trackTable.trackId).';
    usedTrackIds = zeros(0,1);
    for trackId = allTrackIds
        inTrack = trackTable.trackId == trackId;
        if sum(reliable & inTrack) >= opts.minTrackLength
            usedTrackIds(end+1,1) = trackId; %#ok<AGROW>
        end
    end

    if isempty(usedTrackIds)
        diffusion = emptyDiffusionResult(opts);
        diffusion.nTracksTotal = numel(allTrackIds);
        diffusion.truthValidation = truthDiffusionDiagnostics( ...
            trajectories, trackTable, opts.curvaturePerUm, opts);
        return;
    end

    lagFrames = (1:opts.maxMsdLag).';
    lagTimeS = lagFrames * opts.dtS;
    nPairs = zeros(opts.maxMsdLag,1);
    msdSurface = nan(opts.maxMsdLag,1);
    localizationNoiseSurface = nan(opts.maxMsdLag,1);
    msdSurfaceCorrected = nan(opts.maxMsdLag,1);
    msdProjected = nan(opts.maxMsdLag,1);

    for lag = 1:opts.maxMsdLag
        surfaceSquared = zeros(0,1);
        noiseTrace = zeros(0,1);
        projectedSquared = zeros(0,1);
        for trackIndex = 1:numel(usedTrackIds)
            trackId = usedTrackIds(trackIndex);
            rows = find(trackTable.trackId == trackId & reliable);
            [~, order] = sort(trackTable.frame(rows));
            rows = rows(order);
            for firstIndex = 1:numel(rows)
                firstRow = rows(firstIndex);
                targetFrame = trackTable.frame(firstRow) + lag;
                secondIndex = find(trackTable.frame(rows) == ...
                    targetFrame, 1, 'first');
                if isempty(secondIndex)
                    continue;
                end
                secondRow = rows(secondIndex);
                [stepTangent, stepNoise, projectedStep, valid] = ...
                    surfaceStepFromRows(trackTable, firstRow, secondRow, ...
                    curvature.curvaturePerUm);
                if valid
                    surfaceSquared(end+1,1) = ...
                        sum(stepTangent.^2); %#ok<AGROW>
                    noiseTrace(end+1,1) = trace(stepNoise); %#ok<AGROW>
                    projectedSquared(end+1,1) = ...
                        sum(projectedStep.^2); %#ok<AGROW>
                end
            end
        end

        nPairs(lag) = numel(surfaceSquared);
        if nPairs(lag) >= opts.minMsdPairs
            msdSurface(lag) = mean(surfaceSquared);
            localizationNoiseSurface(lag) = mean(noiseTrace);
            msdSurfaceCorrected(lag) = ...
                msdSurface(lag) - localizationNoiseSurface(lag);
            msdProjected(lag) = mean(projectedSquared);
        end
    end

    fitMask = isfinite(msdSurfaceCorrected) & ...
        lagFrames <= opts.maxMsdFitLag;
    [Dmsd, msdSlope] = fitMsdThroughOrigin( ...
        lagTimeS(fitMask), msdSurfaceCorrected(fitMask), nPairs(fitMask));
    projectedMask = isfinite(msdProjected) & ...
        lagFrames <= opts.maxMsdFitLag;
    [Dprojected, projectedSlope] = fitMsdThroughOrigin( ...
        lagTimeS(projectedMask), msdProjected(projectedMask), ...
        nPairs(projectedMask));

    if isfinite(msdSurfaceCorrected(1))
        Dlag1 = max(msdSurfaceCorrected(1),0) / (4*opts.dtS);
    else
        Dlag1 = NaN;
    end

    [stepTangent, stepNoise, stepDt, driftPerDiffusion] = ...
        collectIndependentSurfaceSteps(trackTable, reliable, ...
        usedTrackIds, curvature.curvaturePerUm, opts);
    fisherEstimate = surfaceDiffusionFisherMLE(stepTangent, stepNoise, ...
        stepDt, driftPerDiffusion, Dmsd, opts.dtS);

    msdTable = table(lagFrames, lagTimeS, nPairs, msdSurface, ...
        localizationNoiseSurface, msdSurfaceCorrected, msdProjected, ...
        'VariableNames', {'lagFrames', 'lagTimeS', 'nPairs', ...
        'surfaceMsdUm2', 'meanFisherLocalizationNoiseSurfaceUm2', ...
        'surfaceMsdCorrectedUm2', 'projectedMsdUm2'});

    usedRows = ismember(trackTable.trackId, usedTrackIds) & reliable;
    localizationSigma = sqrt(0.5 * ...
        (trackTable.crbVarXUm2(usedRows) + ...
        trackTable.crbVarYUm2(usedRows)));

    diffusion = struct();
    diffusion.geometry = ...
        'local tangent metric from the estimated paraboloid curvature';
    diffusion.uncertaintyConditioning = ...
        'Diffusion Fisher information conditional on the fitted curvature';
    diffusion.trueDiffusionUm2PerS = opts.diffusionUm2PerS;
    diffusion.curvatureUsedPerUm = curvature.curvaturePerUm;
    diffusion.DLag1CorrectedUm2PerS = Dlag1;
    diffusion.DMsdCorrectedUm2PerS = Dmsd;
    diffusion.DProjectedMsdUm2PerS = Dprojected;
    diffusion.DFisherMLEUm2PerS = fisherEstimate.diffusionUm2PerS;
    diffusion.diffusionFisherSigmaUm2PerS = fisherEstimate.sigmaUm2PerS;
    diffusion.diffusionFisherInformation = fisherEstimate.information;
    diffusion.nIndependentStepsForFisher = fisherEstimate.nSteps;
    diffusion.nLag1StepsForFisher = fisherEstimate.nLag1Steps;
    diffusion.msdSlopeSurfaceCorrected = msdSlope;
    diffusion.msdSlopeProjected = projectedSlope;
    diffusion.msdTable = msdTable;
    diffusion.nTracksTotal = numel(allTrackIds);
    diffusion.nTracksUsed = numel(usedTrackIds);
    diffusion.usedTrackIds = usedTrackIds;
    diffusion.nReliableLocalizations = sum(usedRows);
    diffusion.meanLocalizationSigmaUm = meanOrNan(localizationSigma);
    diffusion.truthValidation = truthDiffusionDiagnostics( ...
        trajectories, trackTable, opts.curvaturePerUm, opts);
end
