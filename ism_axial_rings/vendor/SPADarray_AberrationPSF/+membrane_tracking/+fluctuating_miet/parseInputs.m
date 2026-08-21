function opts = parseInputs(varargin)
    import membrane_tracking.fluctuating_miet.*

    p = inputParser;
    p.FunctionName = 'simulateFluctuatingMIETMembraneTracking';

    % Mean surface and dynamics.
    addParameter(p, 'rngSeed', 23);
    addParameter(p, 'nMolecules', 1);
    addParameter(p, 'membraneRadiusUm', 0.65);
    addParameter(p, 'tipHeightUm', 0.18);
    addParameter(p, 'curvaturePerUm', -0.60);
    addParameter(p, 'diffusionUm2PerS', 0.10);
    addParameter(p, 'dtS', 1e-3);
    addParameter(p, 'nFrames', 1500);
    addParameter(p, 'substepsPerFrame', 5);

    % Fluctuating field. kT = 1 units.
    addParameter(p, 'enableFluctuations', true);
    addParameter(p, 'bendingModulusKT', 20);
    addParameter(p, 'tensionKTPerUm2', 300);
    addParameter(p, 'nModeShells', 6);
    addParameter(p, 'fluctuationRelaxationTimeS', 0.05);
    addParameter(p, 'maxFieldRmsUm', 0.05);

    % Excitation and detection.
    addParameter(p, 'laserPowerW', 1e-5);
    addParameter(p, 'excitationRatePerW', 5e11);
    addParameter(p, 'saturationRateS', 5e7);
    addParameter(p, 'quantumYield', 0.70);
    addParameter(p, 'collectionEfficiency', 0.12);
    addParameter(p, 'detectorQuantumEfficiency', 0.60);
    addParameter(p, 'laserWaistUm', 0.38);

    % MIET calibration and TCSPC.
    addParameter(p, 'lifetimeAtSubstrateNs', 1.0);
    addParameter(p, 'lifetimeSlopeNsPerUm', 6.0);
    addParameter(p, 'laserRepetitionRateHz', 80e6);
    addParameter(p, 'lifetimeBoundsNs', [0.15 8.0]);

    % ISM detector.
    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detectorPitchUm', 0.18);
    addParameter(p, 'detectorFillRatio', 1.0);
    addParameter(p, 'detectorSubsamples', 5);
    addParameter(p, 'psfSigmaUm', 0.14);
    addParameter(p, 'backgroundPhotonsPerDetectorS', 15);
    addParameter(p, 'darkCountsPerDetectorS', 2);

    % Localisation.
    addParameter(p, 'minPhotonsPerLocalization', 25);
    addParameter(p, 'maxLocalizationRadiusUm', 0.56);
    addParameter(p, 'maxReducedPoissonDeviance', 3.0);
    addParameter(p, 'minExpectedCount', 1e-10);

    % Surface estimation.
    addParameter(p, 'fluctuationKernel', 'helfrich');
    addParameter(p, 'estimateFieldParameters', true);
    addParameter(p, 'momentCorrectCurvature', true);
    addParameter(p, 'maxGpObservations', 300);
    addParameter(p, 'nGpNigpIterations', 3);
    addParameter(p, 'topographyGridSize', 31);
    addParameter(p, 'nTopographySnapshots', 3);
    addParameter(p, 'curvatureProfilePoints', 13);
    addParameter(p, 'curvatureProfileHalfWidth', 3.0);

    % Tracking and diffusion.
    addParameter(p, 'trackingPriorDiffusionUm2PerS', []);
    addParameter(p, 'trackingGateChi2', 13.8155);
    addParameter(p, 'minMahalanobisAcceptance', 0);
    addParameter(p, 'gateCensoringCorrection', true);
    addParameter(p, 'nTrackingRefinementIterations', 3);
    addParameter(p, 'trackingRefinementTolerance', 0.02);
    addParameter(p, 'gateAcceptanceQuadraturePoints', 64);
    addParameter(p, 'crbCorrectionScale', 1);
    addParameter(p, 'diffusionGeometry', 'posteriorMean');
    addParameter(p, 'maxGapFrames', 3);
    addParameter(p, 'minTrackLength', 8);
    addParameter(p, 'maxMsdLag', 5);
    addParameter(p, 'maxMsdFitLag', 3);
    addParameter(p, 'minMsdPairs', 3);

    % Output.
    addParameter(p, 'makeFigure', true);
    addParameter(p, 'saveFigures', false);
    addParameter(p, 'writeOutputs', false);
    addParameter(p, 'outputDir', '');
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;

    integerNames = {'nMolecules', 1; 'nFrames', 2; 'substepsPerFrame', 1; ...
        'detectorSubsamples', 2; 'minPhotonsPerLocalization', 1; ...
        'nModeShells', 1; 'maxGapFrames', 0; 'minTrackLength', 2; ...
        'maxMsdLag', 1; 'maxMsdFitLag', 1; 'minMsdPairs', 1; ...
        'maxGpObservations', 20; 'nGpNigpIterations', 1; ...
        'topographyGridSize', 11; 'nTopographySnapshots', 1; ...
        'curvatureProfilePoints', 5; ...
        'nTrackingRefinementIterations', 1; ...
        'gateAcceptanceQuadraturePoints', 16};
    for k = 1:size(integerNames, 1)
        opts.(integerNames{k,1}) = roundScalar(opts.(integerNames{k,1}), ...
            integerNames{k,1}, integerNames{k,2});
    end
    if ~isempty(opts.rngSeed)
        opts.rngSeed = roundScalar(opts.rngSeed, 'rngSeed', 0);
    end

    positiveNames = {'membraneRadiusUm', 'diffusionUm2PerS', 'dtS', ...
        'laserPowerW', 'excitationRatePerW', 'saturationRateS', ...
        'laserWaistUm', 'lifetimeSlopeNsPerUm', 'laserRepetitionRateHz', ...
        'detectorPitchUm', 'detectorFillRatio', 'psfSigmaUm', ...
        'maxLocalizationRadiusUm', 'maxReducedPoissonDeviance', ...
        'trackingGateChi2', 'bendingModulusKT', ...
        'fluctuationRelaxationTimeS', 'maxFieldRmsUm', ...
        'curvatureProfileHalfWidth', 'trackingRefinementTolerance', ...
        'crbCorrectionScale'};
    for k = 1:numel(positiveNames)
        opts.(positiveNames{k}) = validateScalar( ...
            opts.(positiveNames{k}), positiveNames{k}, 0, false);
    end
    nonNegativeNames = {'tipHeightUm', 'backgroundPhotonsPerDetectorS', ...
        'darkCountsPerDetectorS', 'tensionKTPerUm2', ...
        'minMahalanobisAcceptance'};
    for k = 1:numel(nonNegativeNames)
        opts.(nonNegativeNames{k}) = validateScalar( ...
            opts.(nonNegativeNames{k}), nonNegativeNames{k}, 0, true);
    end
    probabilityNames = {'quantumYield', 'collectionEfficiency', ...
        'detectorQuantumEfficiency'};
    for k = 1:numel(probabilityNames)
        opts.(probabilityNames{k}) = validateProbability( ...
            opts.(probabilityNames{k}), probabilityNames{k});
    end
    logicalNames = {'makeFigure', 'saveFigures', 'writeOutputs', ...
        'verbose', 'enableFluctuations', 'estimateFieldParameters', ...
        'momentCorrectCurvature', 'gateCensoringCorrection'};
    for k = 1:numel(logicalNames)
        opts.(logicalNames{k}) = validateLogicalScalar( ...
            opts.(logicalNames{k}), logicalNames{k});
    end

    opts.curvaturePerUm = validateFiniteScalar(opts.curvaturePerUm, ...
        'curvaturePerUm');
    opts.lifetimeAtSubstrateNs = validateFiniteScalar( ...
        opts.lifetimeAtSubstrateNs, 'lifetimeAtSubstrateNs');
    opts.lifetimeBoundsNs = validateIncreasingPair(opts.lifetimeBoundsNs, ...
        'lifetimeBoundsNs');
    opts.fluctuationKernel = normalizeChoice(opts.fluctuationKernel, ...
        {'sqexp', 'helfrich', 'none'}, 'fluctuationKernel');
    opts.diffusionGeometry = normalizeChoice(opts.diffusionGeometry, ...
        {'meanSurface', 'posteriorMean'}, 'diffusionGeometry');

    opts.repetitionPeriodNs = 1e9 / opts.laserRepetitionRateHz;
    opts.maxMsdFitLag = min(opts.maxMsdFitLag, opts.maxMsdLag);
    opts.maxLocalizationRadiusUm = min(opts.maxLocalizationRadiusUm, ...
        opts.membraneRadiusUm);
    if isempty(opts.trackingPriorDiffusionUm2PerS)
        opts.trackingPriorDiffusionUm2PerS = opts.diffusionUm2PerS;
    else
        opts.trackingPriorDiffusionUm2PerS = validateScalar( ...
            opts.trackingPriorDiffusionUm2PerS, ...
            'trackingPriorDiffusionUm2PerS', 0, false);
    end
    if opts.minMahalanobisAcceptance > 0
        warning('simulateFluctuatingMIETMembraneTracking:HardGate', ...
            ['minMahalanobisAcceptance thresholds a Gaussian kernel ' ...
             'value, not a probability. A threshold of %.2f corresponds ' ...
             'to a Mahalanobis gate of %.2f sigma, which censors long ' ...
             'steps and biases D low.'], opts.minMahalanobisAcceptance, ...
            sqrt(max(-2*log(opts.minMahalanobisAcceptance), 0)));
    end

    heights = opts.tipHeightUm + 0.5 * opts.curvaturePerUm * ...
        [0, opts.membraneRadiusUm^2];
    if min(heights) < 0
        error('simulateFluctuatingMIETMembraneTracking:SurfaceBelowSubstrate', ...
            'The mean paraboloid crosses z=0 inside the membrane disk.');
    end
    opts.meanHeightRangeUm = [min(heights) max(heights)];

    if isempty(opts.outputDir)
        projectDirectory = fileparts(fileparts(fileparts( ...
            mfilename('fullpath'))));
        opts.outputDir = fullfile(projectDirectory, ...
            'output_fluctuating_miet');
    else
        opts.outputDir = char(opts.outputDir);
    end
end
