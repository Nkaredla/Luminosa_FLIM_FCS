function validation = estimateTruthDiffusionDiagnostics( ...
        trajectories, trackTable, photophysics, opts)
    import membrane_tracking.focused_ism.*

    positions = trajectories.positionsUm;
    nFrames = size(positions, 3);

    fullR2 = zeros(0, 1);
    for t = 1:(nFrames-1)
        delta = positions(:,:,t+1) - positions(:,:,t);
        delta = applyMinimumImage(delta, opts);
        fullR2 = [fullR2; sum(delta.^2, 2)]; %#ok<AGROW>
    end

    photonCounts = photophysics.detectedPhotonCountsInFov;
    focusR2 = zeros(0, 1);
    for t = 1:(nFrames-1)
        detectable = photonCounts(:,t) >= opts.minPhotonsPerLocalization & ...
            photonCounts(:,t+1) >= opts.minPhotonsPerLocalization;
        if any(detectable)
            delta = positions(detectable,:,t+1) - ...
                positions(detectable,:,t);
            delta = applyMinimumImage(delta, opts);
            focusR2 = [focusR2; sum(delta.^2, 2)]; %#ok<AGROW>
        end
    end

    reliable = rowsUsableForDiffusion(trackTable, opts);
    acceptedTruthR2 = zeros(0, 1);
    acceptedLocalizedR2 = zeros(0, 1);
    nAcceptedAdjacentSteps = 0;
    nIdentitySwitches = 0;
    for trackId = unique(trackTable.trackId).'
        tr = trackTable(trackTable.trackId == trackId, :);
        tr = sortrows(tr, 'frame');
        trReliable = reliable(trackTable.trackId == trackId);
        [~, trOrder] = sort(trackTable.frame(trackTable.trackId == trackId));
        trReliable = trReliable(trOrder);
        for i = 1:(height(tr)-1)
            if tr.frame(i+1) - tr.frame(i) ~= 1 || ...
                    ~trReliable(i) || ~trReliable(i+1)
                continue;
            end
            nAcceptedAdjacentSteps = nAcceptedAdjacentSteps + 1;
            moleculeI = tr.nearestTrueMolecule(i);
            moleculeJ = tr.nearestTrueMolecule(i+1);
            sameTruth = isfinite(moleculeI) && isfinite(moleculeJ) && ...
                moleculeI == moleculeJ;
            if ~sameTruth
                nIdentitySwitches = nIdentitySwitches + 1;
                continue;
            end

            molecule = round(moleculeI);
            trueDelta = positions(molecule,:,tr.frame(i+1)) - ...
                positions(molecule,:,tr.frame(i));
            trueDelta = applyMinimumImage(trueDelta, opts);
            acceptedTruthR2(end+1,1) = sum(trueDelta.^2); %#ok<AGROW>

            localizedDelta = [tr.xUm(i+1)-tr.xUm(i), ...
                tr.yUm(i+1)-tr.yUm(i)];
            localizedDelta = applyMinimumImage(localizedDelta, opts);
            acceptedLocalizedR2(end+1,1) = ...
                sum(localizedDelta.^2); %#ok<AGROW>
        end
    end

    localizationErrorR2 = zeros(0, 1);
    fisherTrace = zeros(0, 1);
    reliableRows = find(reliable);
    for index = 1:numel(reliableRows)
        row = reliableRows(index);
        molecule = trackTable.nearestTrueMolecule(row);
        if ~isfinite(molecule)
            continue;
        end
        molecule = round(molecule);
        frame = trackTable.frame(row);
        errorXY = [trackTable.xUm(row), trackTable.yUm(row)] - ...
            positions(molecule,:,frame);
        errorXY = applyMinimumImage(errorXY, opts);
        localizationErrorR2(end+1,1) = sum(errorXY.^2); %#ok<AGROW>
        fisherTrace(end+1,1) = trackTable.crbVarXUm2(row) + ...
            trackTable.crbVarYUm2(row); %#ok<AGROW>
    end

    validation = struct();
    validation.configuredDiffusionUm2PerS = opts.diffusionUm2PerS;
    validation.fullTrajectoryLag1DiffusionUm2PerS = ...
        radialStepsToDiffusion(fullR2, opts.dtS);
    validation.focusDetectableTruthLag1DiffusionUm2PerS = ...
        radialStepsToDiffusion(focusR2, opts.dtS);
    validation.acceptedTrackTruthLag1DiffusionUm2PerS = ...
        radialStepsToDiffusion(acceptedTruthR2, opts.dtS);
    validation.acceptedTrackLocalizedRawLag1DiffusionUm2PerS = ...
        radialStepsToDiffusion(acceptedLocalizedR2, opts.dtS);
    validation.nFullTrajectorySteps = numel(fullR2);
    validation.nFocusDetectableTruthSteps = numel(focusR2);
    validation.nAcceptedTruthConsistentTrackSteps = ...
        numel(acceptedTruthR2);
    validation.nAcceptedAdjacentTrackSteps = nAcceptedAdjacentSteps;
    validation.truthIdentitySwitchFraction = nIdentitySwitches / ...
        max(nAcceptedAdjacentSteps, 1);
    validation.empiricalLocalizationMseRadialUm2 = ...
        meanOrNan(localizationErrorR2);
    validation.meanFisherLocalizationVarianceTraceUm2 = ...
        meanOrNan(fisherTrace);
    validation.localizationVarianceInflationOverFisher = ...
        meanOrNan(localizationErrorR2) / meanOrNan(fisherTrace);
end
