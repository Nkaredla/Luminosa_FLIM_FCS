function diagnostics = summarizeDensityDiagnostics(localizationTable, ...
        trackTable, detector, diffusion, opts)
    import membrane_tracking.focused_ism.*

    effectivePsfAreaUm2 = detector.effectiveExcitationAreaUm2;
    occupancy = opts.actualMoleculeConcentrationPerUm2 * ...
        effectivePsfAreaUm2;
    if occupancy < 0.05
        regime = 'sparse';
    elseif occupancy < 0.25
        regime = 'moderate-overlap';
    else
        regime = 'dense-information-limited';
    end

    if height(localizationTable) > 0
        [~, uniqueRegionRows] = unique(localizationTable.fitRegionId, ...
            'stable');
        regionCounts = localizationTable.jointEmitterCount(uniqueRegionRows);
        regionProbabilities = ...
            localizationTable.emitterCountPosteriorProbability( ...
            uniqueRegionRows);
        fractionMultiEmitter = mean(regionCounts > 1);
        fractionAtEmitterLimit = mean( ...
            regionCounts >= opts.maxEmittersPerFitRegion);
        meanOrderProbability = mean(regionProbabilities);
    else
        fractionMultiEmitter = NaN;
        fractionAtEmitterLimit = NaN;
        meanOrderProbability = NaN;
    end

    requestedConcentration = opts.moleculeConcentrationPerUm2;
    if isempty(requestedConcentration)
        requestedConcentration = NaN;
    end

    diagnostics = struct();
    diagnostics.requestedConcentrationPerUm2 = requestedConcentration;
    diagnostics.actualConcentrationPerUm2 = ...
        opts.actualMoleculeConcentrationPerUm2;
    diagnostics.expectedMoleculesInDetectorFov = occupancy;
    diagnostics.expectedMoleculesInEffectiveFocus = occupancy;
    diagnostics.effectivePsfAreaUm2 = effectivePsfAreaUm2;
    diagnostics.meanMoleculesPerEffectivePsfArea = occupancy;
    diagnostics.excitationWaistUm = opts.laserWaistUm;
    diagnostics.detectionPsfSigmaUm = opts.psfSigmaUm;
    diagnostics.detectorPitchUm = detector.detectorPitchUm;
    diagnostics.analysisRegime = regime;
    diagnostics.meanLocalizationsPerFrame = ...
        height(localizationTable) / opts.nFrames;
    diagnostics.fractionMultiEmitterFitRegions = fractionMultiEmitter;
    diagnostics.fractionFitRegionsAtEmitterLimit = ...
        fractionAtEmitterLimit;
    diagnostics.meanSelectedEmitterCountPosteriorProbability = ...
        meanOrderProbability;
    diagnostics.fractionAmbiguousAssociations = ...
        diffusion.fractionAmbiguousAssociations;
    diagnostics.fractionUncertainEmitterCounts = ...
        diffusion.fractionUncertainEmitterCounts;
    diagnostics.fractionRejectedForDiffusion = ...
        diffusion.fractionRejectedForDiffusion;
    diagnostics.meanAssociationProbability = ...
        diffusion.meanAssociationProbability;
    diagnostics.meanEmitterCountPosteriorProbability = ...
        diffusion.meanEmitterCountPosteriorProbability;
    diagnostics.meanHardGateAcceptanceProbability = ...
        diffusion.meanHardGateAcceptanceProbability;
    diagnostics.hardGateCensoringCorrectionApplied = ...
        diffusion.hardGateCensoringCorrectionApplied;
    diagnostics.fractionNonidentifiablePositions = ...
        diffusion.fractionNonidentifiablePositions;
    diagnostics.fractionTrackRowsAtAssociationHypothesisLimit = meanOrNan( ...
        double(trackTable.associationHypothesisCount >= ...
        opts.maxAssociationHypotheses));
    diagnostics.nTracksTotal = diffusion.nTracksTotal;
    diagnostics.nTracksUsedForDiffusion = diffusion.nTracksUsed;
    diagnostics.localizationMethod = opts.localizationMethod;
    diagnostics.trackingMethod = opts.trackingMethod;
    diagnostics.detectorLayout = detector.layout;
    diagnostics.nDetectorChannels = detector.nChannels;
    diagnostics.maxIdentifiableEmittersFromChannelCount = ...
        floor((detector.nChannels - 1) / 3);
    diagnostics.maxEmittersUsedForModelSelection = ...
        opts.maxEmittersPerFitRegion;
    diagnostics.trackingPriorDiffusionUm2PerS = ...
        opts.trackingPriorDiffusionUm2PerS;
    diagnostics.nTrackingRefinementIterations = ...
        diffusion.nTrackingRefinementIterations;
    diagnostics.nTrackRows = height(trackTable);
    diagnostics.fullTrajectoryTruthDiffusionUm2PerS = ...
        diffusion.truthValidation.fullTrajectoryLag1DiffusionUm2PerS;
    diagnostics.focusDetectableTruthDiffusionUm2PerS = ...
        diffusion.truthValidation.focusDetectableTruthLag1DiffusionUm2PerS;
    diagnostics.acceptedTrackTruthDiffusionUm2PerS = ...
        diffusion.truthValidation.acceptedTrackTruthLag1DiffusionUm2PerS;
    diagnostics.truthIdentitySwitchFraction = ...
        diffusion.truthValidation.truthIdentitySwitchFraction;
    diagnostics.localizationVarianceInflationOverFisher = ...
        diffusion.truthValidation.localizationVarianceInflationOverFisher;
end
