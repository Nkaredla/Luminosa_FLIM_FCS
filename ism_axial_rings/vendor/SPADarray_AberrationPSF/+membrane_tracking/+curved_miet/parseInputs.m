function opts = parseInputs(varargin)
    import membrane_tracking.curved_miet.*

    p = inputParser;
    p.FunctionName = 'simulateMIETCurvedMembraneTracking';

    % Surface and intrinsic Brownian dynamics.
    addParameter(p, 'rngSeed', 19);
    addParameter(p, 'nMolecules', 1);
    addParameter(p, 'membraneRadiusUm', 0.65);
    addParameter(p, 'tipHeightUm', 0.18);
    addParameter(p, 'curvaturePerUm', -0.60);
    addParameter(p, 'diffusionUm2PerS', 0.10);
    addParameter(p, 'dtS', 0.001);
    addParameter(p, 'nFrames', 1500);
    addParameter(p, 'substepsPerFrame', 5);
    addParameter(p, 'initialPositionsUm', []);

    % Focused excitation and detected photon rate.
    addParameter(p, 'laserPowerW', 1e-5);
    addParameter(p, 'excitationRatePerW', 5e11);
    addParameter(p, 'saturationRateS', 5e7);
    addParameter(p, 'quantumYield', 0.70);
    addParameter(p, 'collectionEfficiency', 0.12);
    addParameter(p, 'detectorQuantumEfficiency', 0.60);
    addParameter(p, 'laserWaistUm', 0.38);

    % Linear MIET calibration and ideal pulsed TCSPC.
    addParameter(p, 'lifetimeAtSubstrateNs', 1.0);
    addParameter(p, 'lifetimeSlopeNsPerUm', 6.0);
    addParameter(p, 'laserRepetitionRateHz', 80e6);
    addParameter(p, 'lifetimeBoundsNs', [0.15 8.0]);

    % Focused finite-channel ISM detector.
    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detectorPitchUm', 0.18);
    addParameter(p, 'detectorFillRatio', 1.0);
    addParameter(p, 'detectorSubsamples', 5);
    addParameter(p, 'psfSigmaUm', 0.14);
    addParameter(p, 'backgroundPhotonsPerDetectorS', 15);
    addParameter(p, 'darkCountsPerDetectorS', 2);

    % Single-emitter localization and lifetime inference.
    addParameter(p, 'minPhotonsPerLocalization', 25);
    addParameter(p, 'maxLocalizationRadiusUm', 0.56);
    addParameter(p, 'maxReducedPoissonDeviance', 3.0);
    addParameter(p, 'minExpectedCount', 1e-10);

    % Curvature fit. The tip height is fitted jointly with signed curvature.
    addParameter(p, 'tipHeightFitBoundsUm', [0.001 0.80]);
    addParameter(p, 'curvatureFitBoundsPerUm', [-3.0 3.0]);
    addParameter(p, 'minCurvatureLocalizations', 20);

    % Low-density tracking and surface diffusion estimation.
    addParameter(p, 'trackingPriorDiffusionUm2PerS', 0.10);
    addParameter(p, 'trackingGateChi2', 13.8155);
    addParameter(p, 'maxGapFrames', 3);
    addParameter(p, 'minTrackLength', 8);
    addParameter(p, 'minAssociationProbabilityForDiffusion', 0.20);
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

    if ~isempty(opts.rngSeed)
        opts.rngSeed = roundScalar(opts.rngSeed, 'rngSeed', 0);
    end
    opts.nMolecules = roundScalar(opts.nMolecules, 'nMolecules', 1);
    opts.nFrames = roundScalar(opts.nFrames, 'nFrames', 2);
    opts.substepsPerFrame = roundScalar(opts.substepsPerFrame, ...
        'substepsPerFrame', 1);
    opts.detectorSubsamples = roundScalar(opts.detectorSubsamples, ...
        'detectorSubsamples', 2);
    opts.minPhotonsPerLocalization = roundScalar( ...
        opts.minPhotonsPerLocalization, 'minPhotonsPerLocalization', 1);
    opts.minCurvatureLocalizations = roundScalar( ...
        opts.minCurvatureLocalizations, 'minCurvatureLocalizations', 3);
    opts.maxGapFrames = roundScalar(opts.maxGapFrames, 'maxGapFrames', 0);
    opts.minTrackLength = roundScalar(opts.minTrackLength, ...
        'minTrackLength', 2);
    opts.maxMsdLag = roundScalar(opts.maxMsdLag, 'maxMsdLag', 1);
    opts.maxMsdFitLag = roundScalar(opts.maxMsdFitLag, ...
        'maxMsdFitLag', 1);
    opts.minMsdPairs = roundScalar(opts.minMsdPairs, 'minMsdPairs', 1);

    positiveNames = {'membraneRadiusUm', 'diffusionUm2PerS', 'dtS', ...
        'laserPowerW', 'excitationRatePerW', 'saturationRateS', ...
        'laserWaistUm', 'lifetimeSlopeNsPerUm', ...
        'laserRepetitionRateHz', 'detectorPitchUm', ...
        'detectorFillRatio', 'psfSigmaUm', ...
        'maxLocalizationRadiusUm', 'maxReducedPoissonDeviance', ...
        'trackingPriorDiffusionUm2PerS', 'trackingGateChi2'};
    for k = 1:numel(positiveNames)
        name = positiveNames{k};
        opts.(name) = validateScalar(opts.(name), name, 0, false);
    end

    nonnegativeNames = {'tipHeightUm', 'backgroundPhotonsPerDetectorS', ...
        'darkCountsPerDetectorS'};
    for k = 1:numel(nonnegativeNames)
        name = nonnegativeNames{k};
        opts.(name) = validateScalar(opts.(name), name, 0, true);
    end
    opts.curvaturePerUm = validateFiniteScalar(opts.curvaturePerUm, ...
        'curvaturePerUm');
    opts.lifetimeAtSubstrateNs = validateFiniteScalar( ...
        opts.lifetimeAtSubstrateNs, 'lifetimeAtSubstrateNs');

    probabilityNames = {'quantumYield', 'collectionEfficiency', ...
        'detectorQuantumEfficiency'};
    for k = 1:numel(probabilityNames)
        name = probabilityNames{k};
        opts.(name) = validateProbability(opts.(name), name);
    end

    opts.lifetimeBoundsNs = validateIncreasingPair( ...
        opts.lifetimeBoundsNs, 'lifetimeBoundsNs');
    opts.tipHeightFitBoundsUm = validateIncreasingPair( ...
        opts.tipHeightFitBoundsUm, 'tipHeightFitBoundsUm');
    opts.curvatureFitBoundsPerUm = validateIncreasingPair( ...
        opts.curvatureFitBoundsPerUm, 'curvatureFitBoundsPerUm');
    opts.minAssociationProbabilityForDiffusion = validateProbability( ...
        opts.minAssociationProbabilityForDiffusion, ...
        'minAssociationProbabilityForDiffusion');

    logicalNames = {'makeFigure', 'saveFigures', 'writeOutputs', 'verbose'};
    for k = 1:numel(logicalNames)
        name = logicalNames{k};
        opts.(name) = validateLogicalScalar(opts.(name), name);
    end

    if ~(ischar(opts.detectorLayout) || ...
            (isstring(opts.detectorLayout) && isscalar(opts.detectorLayout)))
        error('simulateMIETCurvedMembraneTracking:BadOption', ...
            'detectorLayout must be a character vector or scalar string.');
    end
    opts.detectorLayout = char(opts.detectorLayout);
    opts.repetitionPeriodNs = 1e9 / opts.laserRepetitionRateHz;
    if opts.lifetimeBoundsNs(2) >= 0.95 * opts.repetitionPeriodNs
        warning('simulateMIETCurvedMembraneTracking:LongLifetimeBound', ...
            ['The upper lifetime bound is close to the TCSPC repetition ' ...
            'period; lifetime information can become weak.']);
    end

    surfaceHeights = opts.tipHeightUm + 0.5 * opts.curvaturePerUm * ...
        [0, opts.membraneRadiusUm^2];
    minimumHeight = min(surfaceHeights);
    maximumHeight = max(surfaceHeights);
    if minimumHeight < 0
        error('simulateMIETCurvedMembraneTracking:SurfaceBelowSubstrate', ...
            ['The configured paraboloid crosses below z=0. Increase ' ...
            'tipHeightUm, reduce membraneRadiusUm, or reduce the magnitude ' ...
            'of negative curvaturePerUm.']);
    end
    trueLifetimeRange = opts.lifetimeAtSubstrateNs + ...
        opts.lifetimeSlopeNsPerUm * [minimumHeight maximumHeight];
    if trueLifetimeRange(1) <= 0 || ...
            trueLifetimeRange(1) < opts.lifetimeBoundsNs(1) || ...
            trueLifetimeRange(2) > opts.lifetimeBoundsNs(2)
        error('simulateMIETCurvedMembraneTracking:LifetimeOutOfBounds', ...
            ['The true lifetime range [%.4g %.4g] ns is outside the ' ...
            'configured positive lifetimeBoundsNs.'], ...
            trueLifetimeRange(1), trueLifetimeRange(2));
    end
    opts.membraneHeightRangeUm = [minimumHeight maximumHeight];
    opts.trueLifetimeRangeNs = trueLifetimeRange;
    if opts.tipHeightUm <= opts.tipHeightFitBoundsUm(1) || ...
            opts.tipHeightUm >= opts.tipHeightFitBoundsUm(2)
        warning('simulateMIETCurvedMembraneTracking:TipOutsideFitBounds', ...
            'The simulated tipHeightUm lies outside its fit bounds.');
    end
    if opts.curvaturePerUm <= opts.curvatureFitBoundsPerUm(1) || ...
            opts.curvaturePerUm >= opts.curvatureFitBoundsPerUm(2)
        warning('simulateMIETCurvedMembraneTracking:CurvatureOutsideFitBounds', ...
            'The simulated curvaturePerUm lies outside its fit bounds.');
    end

    if opts.maxLocalizationRadiusUm > opts.membraneRadiusUm
        opts.maxLocalizationRadiusUm = opts.membraneRadiusUm;
    end
    if opts.maxMsdFitLag > opts.maxMsdLag
        opts.maxMsdFitLag = opts.maxMsdLag;
    end

    if isempty(opts.outputDir)
        projectDirectory = fileparts(fileparts(fileparts( ...
            mfilename('fullpath'))));
        opts.outputDir = fullfile(projectDirectory, ...
            'output_miet_curved_membrane');
    elseif ~(ischar(opts.outputDir) || ...
            (isstring(opts.outputDir) && isscalar(opts.outputDir)))
        error('simulateMIETCurvedMembraneTracking:BadOption', ...
            'outputDir must be a character vector or scalar string.');
    else
        opts.outputDir = char(opts.outputDir);
    end
end
