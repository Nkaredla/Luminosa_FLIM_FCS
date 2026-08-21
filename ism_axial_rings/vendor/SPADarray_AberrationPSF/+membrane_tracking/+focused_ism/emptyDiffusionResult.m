function diffusion = emptyDiffusionResult(opts)
    import membrane_tracking.focused_ism.*

    lagVec = (1:opts.maxMsdLag).';
    lagTimeS = lagVec * opts.dtS;
    nPairs = zeros(opts.maxMsdLag, 1);
    nanVec = nan(opts.maxMsdLag, 1);
    msdTable = table(lagVec, lagTimeS, nPairs, nanVec, nanVec, nanVec, ...
        'VariableNames', {'lagFrames', 'lagTimeS', 'nPairs', 'msdUm2', ...
        'meanFisherLocalizationNoiseUm2', 'msdCorrectedUm2'});

    diffusion = struct();
    diffusion.trueDiffusionUm2PerS = opts.diffusionUm2PerS;
    diffusion.DLag1CorrectedUm2PerS = NaN;
    diffusion.DLag1CorrectedUnconstrainedUm2PerS = NaN;
    diffusion.DMsdCorrectedUm2PerS = NaN;
    diffusion.DMsdCorrectedFreeInterceptUm2PerS = NaN;
    diffusion.DMsdRawUm2PerS = NaN;
    diffusion.DFisherMLEUm2PerS = NaN;
    diffusion.DFisherCRBUm2PerS = NaN;
    diffusion.diffusionFisherInformation = NaN;
    diffusion.meanHardGateAcceptanceProbability = NaN;
    diffusion.hardGateCensoringCorrectionApplied = false;
    diffusion.meanGateAcceptanceProbability = NaN;
    diffusion.gateCensoringCorrectionApplied = false;
    diffusion.gateAcceptanceModel = ...
        ['Exact angular quadrature for the Euclidean/Mahalanobis hard-gate ' ...
        'intersection; JPDA confidence and focus selection are not included'];
    diffusion.nIndependentStepsForFisher = 0;
    diffusion.nLag1PairsForFisher = 0;
    diffusion.msdInterceptCorrectedUm2 = NaN;
    diffusion.msdInterceptCorrectedFreeUm2 = NaN;
    diffusion.msdInterceptRawUm2 = NaN;
    diffusion.msdSlopeCorrected = [NaN; NaN];
    diffusion.msdSlopeCorrectedFree = [NaN; NaN];
    diffusion.msdSlopeRaw = [NaN; NaN];
    diffusion.msdTable = msdTable;
    diffusion.nTracksTotal = 0;
    diffusion.nTracksUsed = 0;
    diffusion.usedTrackIds = [];
    diffusion.meanLocalizationSigmaUm = NaN;
    diffusion.medianLocalizationSigmaUm = NaN;
    diffusion.nLocalizationsTotal = 0;
    diffusion.nAssociationReliableLocalizations = 0;
    diffusion.nEmitterCountReliableLocalizations = 0;
    diffusion.nReliableLocalizationsForDiffusion = 0;
    diffusion.fractionAmbiguousAssociations = NaN;
    diffusion.fractionUncertainEmitterCounts = NaN;
    diffusion.fractionRejectedForDiffusion = NaN;
    diffusion.meanAssociationProbability = NaN;
    diffusion.meanEmitterCountPosteriorProbability = NaN;
    diffusion.nPositionIdentifiableLocalizations = 0;
    diffusion.fractionNonidentifiablePositions = NaN;
    diffusion.nTrackingRefinementIterations = 0;
    diffusion.trackingPriorHistoryUm2PerS = zeros(0, 1);
    diffusion.trackingEstimateHistoryUm2PerS = zeros(0, 1);
    diffusion.trackingFixedPointRelativeChange = NaN;
end
