function diffusion = emptyDiffusionResult(opts)
    import membrane_tracking.fluctuating_miet.*

    lagFrames = (1:opts.maxMsdLag).';
    emptyValues = nan(opts.maxMsdLag,1);
    diffusion = struct();
    diffusion.geometry = [opts.diffusionGeometry ...
        ' tangent metric from the estimated membrane'];
    diffusion.trueDiffusionUm2PerS = opts.diffusionUm2PerS;
    diffusion.curvatureUsedPerUm = NaN;
    diffusion.DMsdCorrectedUm2PerS = NaN;
    diffusion.DFisherMLEUm2PerS = NaN;
    diffusion.diffusionFisherSigmaUm2PerS = NaN;
    diffusion.diffusionFisherInformation = NaN;
    diffusion.meanGateAcceptanceProbability = NaN;
    diffusion.gateCensoringCorrectionApplied = false;
    diffusion.gateAcceptanceModel = ['Exact anisotropic Mahalanobis-gate ' ...
        'conditioning with Ito mean; other selection is not conditioned.'];
    diffusion.nIndependentStepsForFisher = 0;
    diffusion.nLag1StepsForFisher = 0;
    diffusion.msdTable = table(lagFrames, lagFrames*opts.dtS, ...
        zeros(opts.maxMsdLag,1), emptyValues, emptyValues, emptyValues, ...
        'VariableNames', {'lagFrames','lagTimeS','nPairs', ...
        'surfaceMsdUm2','meanFisherNoiseUm2','surfaceMsdCorrectedUm2'});
    diffusion.nTracksTotal = 0;
    diffusion.nTracksUsed = 0;
    diffusion.nTrackingRefinementIterations = 0;
    diffusion.trackingPriorHistoryUm2PerS = zeros(0,1);
    diffusion.trackingEstimateHistoryUm2PerS = zeros(0,1);
    diffusion.trackingFixedPointRelativeChange = NaN;
end
