function result = emptyExperimentalBnpResult(opts)
    import membrane_tracking.focused_ism.*

    % The sampler is a finite beta-Bernoulli approximation. It is deliberately
    % labelled experimental: the single-time Brownian-bridge kernel is not the
    % block trajectory update required for reliable dense-data inference.
    result = struct();
    result.wasRun = false;
    result.estimateAvailable = false;
    result.countEstimateAvailable = false;
    result.passedApplicabilityChecks = false;
    result.isValidatedForUse = false;
    result.status = 'disabled';
    result.method = ['experimental finite beta-Bernoulli raw-count ' ...
        'trajectory sampler'];
    result.posteriorMeanD = NaN;
    result.posteriorMedianD = NaN;
    result.credibleIntervalD = [NaN NaN];
    result.effectiveSampleSizeD = NaN;
    result.splitRhatD = NaN;
    result.posteriorMedianActiveEmitters = NaN;
    result.posteriorMeanBrightEmittersPerFrame = NaN;
    result.posteriorMeanBackgroundPerChannel = NaN;
    result.posteriorMeanBrightnessScale = NaN;
    result.posteriorMeanMolecularBrightnessPhotons = NaN;
    result.posteriorCenterBrightnessPhotons = NaN;
    result.brightnessMode = opts.bnpBrightnessMode;
    result.calibratedBrightnessScale = ...
        opts.bnpCalibratedBrightnessScale;
    result.brightnessCalibrationRelativeStd = ...
        opts.bnpBrightnessCalibrationRelativeStd;
    result.meanPositionAcceptance = NaN;
    result.fractionPoolSaturated = NaN;
    result.medianInformativeSteps = NaN;
    result.elapsedSeconds = NaN;
    result.poolSize = opts.bnpPoolSize;
    result.trueContributingMolecules = NaN;
    result.photophysicsModelCompatible = false;
    result.partialDiffusionPosterior = true;
    result.periodicBridgeApproximation = true;
    result.warningMessages = {};
    result.recommendation = ['Disabled. Set runExperimentalBnpSampler=true ' ...
        'to run the experimental joint sampler.'];
    emptyColumn = zeros(0, 1);
    result.trace = table(emptyColumn, emptyColumn, emptyColumn, emptyColumn, ...
        emptyColumn, emptyColumn, emptyColumn, emptyColumn, emptyColumn, ...
        emptyColumn, emptyColumn, 'VariableNames', {'diffusionUm2PerS', ...
        'activeEmitterCount', 'meanBrightEmittersPerFrame', ...
        'backgroundCountsPerChannel', 'brightnessScale', ...
        'positionAcceptance', 'indicatorChanges', ...
        'informativeStepCount', 'backgroundAccepted', ...
        'brightnessAccepted', 'meanMolecularBrightnessPhotons'});
    result.finalTrajectoriesUm = zeros(0, 0, 2);
    result.finalActiveIndicators = false(0, 1);
    result.emitterOccupancyProbability = zeros(0, 1);
    result.frameEmitterCountProbability = zeros(0, opts.bnpPoolSize + 1);
    emptyFrame = zeros(0, 1);
    result.frameEmitterCountPosterior = table(emptyFrame, emptyFrame, ...
        emptyFrame, emptyFrame, emptyFrame, emptyFrame, emptyFrame, ...
        'VariableNames', {'frame', 'posteriorMeanCount', 'mapCount', ...
        'decisionProbability', 'entropyBits', ...
        'credibleLowerCount', 'credibleUpperCount'});
    result.posteriorMeanEmittersInFocus = NaN;
    result.medianFrameCountDecisionProbability = NaN;
    result.fractionFramesCountDecisionAboveThreshold = NaN;
    result.activeEmitterCountProbability = ...
        zeros(1, opts.bnpPoolSize + 1);
    result.posteriorMapActiveEmitterCount = NaN;
    result.activeEmitterCountDecisionProbability = NaN;
    result.activeEmitterCountEntropyBits = NaN;
    result.activeEmitterCountCredibleInterval = [NaN NaN];
    result.adjacentCountSeparationSigma = NaN;
    result.instantaneousSpatialParameterLimit = NaN;
    result.effectiveSampleSizeActiveEmitterCount = NaN;
    result.splitRhatActiveEmitterCount = NaN;
    result.effectiveSampleSizeMeanInFocusCount = NaN;
    result.splitRhatMeanInFocusCount = NaN;
    result.countDiagnosticsPassed = false;
    result.countIsExternallyValidated = false;
    result.countWarningMessages = {};
    result.diffusionDiagnosticsPassed = false;
    result.diffusionWarningMessages = {};
end
