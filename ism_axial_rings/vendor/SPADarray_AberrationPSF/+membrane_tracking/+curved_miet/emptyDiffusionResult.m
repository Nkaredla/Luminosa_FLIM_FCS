function diffusion = emptyDiffusionResult(opts)
    import membrane_tracking.curved_miet.*

    lagFrames = (1:opts.maxMsdLag).';
    lagTimeS = lagFrames * opts.dtS;
    nPairs = zeros(opts.maxMsdLag,1);
    emptyValues = nan(opts.maxMsdLag,1);
    diffusion = struct();
    diffusion.geometry = ...
        'local tangent metric from the estimated paraboloid curvature';
    diffusion.uncertaintyConditioning = ...
        'Diffusion Fisher information conditional on the fitted curvature';
    diffusion.trueDiffusionUm2PerS = opts.diffusionUm2PerS;
    diffusion.curvatureUsedPerUm = NaN;
    diffusion.DLag1CorrectedUm2PerS = NaN;
    diffusion.DMsdCorrectedUm2PerS = NaN;
    diffusion.DProjectedMsdUm2PerS = NaN;
    diffusion.DFisherMLEUm2PerS = NaN;
    diffusion.diffusionFisherSigmaUm2PerS = NaN;
    diffusion.diffusionFisherInformation = NaN;
    diffusion.nIndependentStepsForFisher = 0;
    diffusion.nLag1StepsForFisher = 0;
    diffusion.msdSlopeSurfaceCorrected = NaN;
    diffusion.msdSlopeProjected = NaN;
    diffusion.msdTable = table(lagFrames, lagTimeS, nPairs, ...
        emptyValues, emptyValues, emptyValues, emptyValues, ...
        'VariableNames', {'lagFrames', 'lagTimeS', 'nPairs', ...
        'surfaceMsdUm2', 'meanFisherLocalizationNoiseSurfaceUm2', ...
        'surfaceMsdCorrectedUm2', 'projectedMsdUm2'});
    diffusion.nTracksTotal = 0;
    diffusion.nTracksUsed = 0;
    diffusion.usedTrackIds = zeros(0,1);
    diffusion.nReliableLocalizations = 0;
    diffusion.meanLocalizationSigmaUm = NaN;
end
