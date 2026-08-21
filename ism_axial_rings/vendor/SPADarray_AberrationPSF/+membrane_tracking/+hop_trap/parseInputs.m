function opts = parseInputs(varargin)
%PARSEINPUTS Parse motion-model options and focused-ISM options.

    import membrane_tracking.hop_trap.*

    p = inputParser;
    p.FunctionName = 'simulateHopTrapMembraneTracking';
    p.KeepUnmatched = true;

    addParameter(p, 'diffusionMode', 'hop');
    addParameter(p, 'referencePreset', 'auto');
    addParameter(p, 'freeDiffusionUm2PerS', 0.20);
    addParameter(p, 'trappedDiffusionUm2PerS', 1e-10);
    addParameter(p, 'hopProbability', 0.01);
    addParameter(p, 'compartmentSizeUm', 0.10);
    addParameter(p, 'meanFreeTimeS', 0.030);
    addParameter(p, 'meanTrappedTimeS', 0.010);
    addParameter(p, 'startAtFocus', true);

    addParameter(p, 'minimumInferenceSteps', 30);
    addParameter(p, 'minimumMsdPairs', 12);
    addParameter(p, 'msdMaximumLagFrames', 60);
    addParameter(p, 'bicEvidenceThreshold', 10);
    addParameter(p, 'maximumTrapDiffusionRatio', 0.35);
    addParameter(p, 'minimumTrapStateOccupancyFraction', 0.05);
    addParameter(p, 'minimumHopDiffusionRatio', 1.5);
    addParameter(p, ...
        'maximumCompartmentFractionOfLocalizationRadius', 0.75);
    addParameter(p, 'trapHmmStarts', 6);
    addParameter(p, 'trapHmmMaximumIterations', 150);

    addParameter(p, 'runStickyBnp', true);
    addParameter(p, 'stickyBnpMaximumStates', 6);
    addParameter(p, 'stickyBnpBurnIn', 150);
    addParameter(p, 'stickyBnpSamples', 150);
    addParameter(p, 'stickyBnpThinning', 2);
    addParameter(p, 'stickyBnpAlpha', 2);
    addParameter(p, 'stickyBnpGamma', 1);
    addParameter(p, 'stickyBnpKappa', 20);
    addParameter(p, 'stickyBnpMinimumOccupancyFraction', 0.03);
    addParameter(p, 'stickyBnpDiffusionGridSize', 100);
    addParameter(p, 'stickyBnpLogDiffusionPriorSigma', 2.5);

    parse(p, varargin{:});
    custom = applyReferencePreset(p.Results, p.UsingDefaults);
    unmatched = p.Unmatched;

    defaults = struct();
    defaults.nMolecules = 1;
    defaults.boxSizeUm = [1.2 1.2];
    defaults.dtS = 0.001;
    defaults.nFrames = 2000;
    defaults.substepsPerFrame = 10;
    defaults.boundaryCondition = 'periodic';
    defaults.initialOnFraction = 1;
    defaults.blinkOffRateS = 0;
    defaults.blinkOnRateS = 0;
    defaults.bleachProbabilityPerAbsorbedPhoton = 0;
    defaults.trackingMethod = 'nearestNeighbor';
    defaults.trackingMaxStepUm = 0.40;
    defaults.trackingPriorDiffusionUm2PerS = custom.freeDiffusionUm2PerS;
    defaults.maxGapFrames = 3;
    defaults.minTrackLength = 12;
    defaults.maxMsdLag = custom.msdMaximumLagFrames;
    defaults.maxMsdFitLag = min(8, custom.msdMaximumLagFrames);
    defaults.minMsdPairs = custom.minimumMsdPairs;
    defaults.maxEmittersPerFitRegion = 1;
    defaults.maxLocalizationsPerFrame = 1;
    defaults.nTrackingRefinementIterations = 1;
    defaults.makeFigure = true;
    defaults.verbose = true;
    projectDirectory = fileparts(fileparts(fileparts( ...
        mfilename('fullpath'))));
    defaults.outputDir = fullfile(projectDirectory, 'output_matlab', ...
        'hop_trap_tracking');

    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~hasFieldIgnoreCase(unmatched, names{k})
            unmatched.(names{k}) = defaults.(names{k});
        end
    end
    focusedArguments = structToNameValue(unmatched);
    opts = membrane_tracking.focused_ism.parseInputs(focusedArguments{:});

    customNames = fieldnames(custom);
    for k = 1:numel(customNames)
        opts.(customNames{k}) = custom.(customNames{k});
    end
    opts.diffusionMode = normalizeChoice(opts.diffusionMode, ...
        {'free', 'hop', 'trap', 'mixed'}, 'diffusionMode');
    opts.freeDiffusionUm2PerS = positiveScalar( ...
        opts.freeDiffusionUm2PerS, 'freeDiffusionUm2PerS', true);
    opts.trappedDiffusionUm2PerS = positiveScalar( ...
        opts.trappedDiffusionUm2PerS, ...
        'trappedDiffusionUm2PerS', true);
    opts.hopProbability = probabilityScalar( ...
        opts.hopProbability, 'hopProbability');
    opts.compartmentSizeUm = positiveScalar( ...
        opts.compartmentSizeUm, 'compartmentSizeUm', false);
    opts.meanFreeTimeS = positiveScalar( ...
        opts.meanFreeTimeS, 'meanFreeTimeS', false);
    opts.meanTrappedTimeS = positiveScalar( ...
        opts.meanTrappedTimeS, 'meanTrappedTimeS', false);
    opts.startAtFocus = logicalScalar(opts.startAtFocus, 'startAtFocus');

    opts.minimumInferenceSteps = integerScalar( ...
        opts.minimumInferenceSteps, 'minimumInferenceSteps', 5);
    opts.minimumMsdPairs = integerScalar( ...
        opts.minimumMsdPairs, 'minimumMsdPairs', 3);
    opts.msdMaximumLagFrames = integerScalar( ...
        opts.msdMaximumLagFrames, 'msdMaximumLagFrames', 2);
    opts.bicEvidenceThreshold = positiveScalar( ...
        opts.bicEvidenceThreshold, 'bicEvidenceThreshold', true);
    opts.maximumTrapDiffusionRatio = positiveScalar( ...
        opts.maximumTrapDiffusionRatio, ...
        'maximumTrapDiffusionRatio', false);
    opts.minimumTrapStateOccupancyFraction = positiveScalar( ...
        opts.minimumTrapStateOccupancyFraction, ...
        'minimumTrapStateOccupancyFraction', false);
    if opts.minimumTrapStateOccupancyFraction >= 0.5
        error('simulateHopTrapMembraneTracking:BadOption', ...
            'minimumTrapStateOccupancyFraction must be less than 0.5.');
    end
    opts.minimumHopDiffusionRatio = positiveScalar( ...
        opts.minimumHopDiffusionRatio, ...
        'minimumHopDiffusionRatio', false);
    opts.maximumCompartmentFractionOfLocalizationRadius = positiveScalar( ...
        opts.maximumCompartmentFractionOfLocalizationRadius, ...
        'maximumCompartmentFractionOfLocalizationRadius', false);
    opts.trapHmmStarts = integerScalar(opts.trapHmmStarts, ...
        'trapHmmStarts', 1);
    opts.trapHmmMaximumIterations = integerScalar( ...
        opts.trapHmmMaximumIterations, 'trapHmmMaximumIterations', 10);

    opts.runStickyBnp = logicalScalar(opts.runStickyBnp, 'runStickyBnp');
    opts.stickyBnpMaximumStates = integerScalar( ...
        opts.stickyBnpMaximumStates, 'stickyBnpMaximumStates', 2);
    opts.stickyBnpBurnIn = integerScalar( ...
        opts.stickyBnpBurnIn, 'stickyBnpBurnIn', 0);
    opts.stickyBnpSamples = integerScalar( ...
        opts.stickyBnpSamples, 'stickyBnpSamples', 10);
    opts.stickyBnpThinning = integerScalar( ...
        opts.stickyBnpThinning, 'stickyBnpThinning', 1);
    opts.stickyBnpAlpha = positiveScalar( ...
        opts.stickyBnpAlpha, 'stickyBnpAlpha', false);
    opts.stickyBnpGamma = positiveScalar( ...
        opts.stickyBnpGamma, 'stickyBnpGamma', false);
    opts.stickyBnpKappa = positiveScalar( ...
        opts.stickyBnpKappa, 'stickyBnpKappa', true);
    opts.stickyBnpMinimumOccupancyFraction = positiveScalar( ...
        opts.stickyBnpMinimumOccupancyFraction, ...
        'stickyBnpMinimumOccupancyFraction', false);
    if opts.stickyBnpMinimumOccupancyFraction >= 0.5
        error('simulateHopTrapMembraneTracking:BadOption', ...
            'stickyBnpMinimumOccupancyFraction must be less than 0.5.');
    end
    opts.stickyBnpDiffusionGridSize = integerScalar( ...
        opts.stickyBnpDiffusionGridSize, ...
        'stickyBnpDiffusionGridSize', 30);
    opts.stickyBnpLogDiffusionPriorSigma = positiveScalar( ...
        opts.stickyBnpLogDiffusionPriorSigma, ...
        'stickyBnpLogDiffusionPriorSigma', false);

    % The focused-ISM model uses this field for tracking defaults and audit
    % output; the heterogeneous simulator uses the explicit free/trapped D.
    opts.diffusionUm2PerS = opts.freeDiffusionUm2PerS;
    if opts.nMolecules / prod(opts.boxSizeUm) > 2
        warning('simulateHopTrapMembraneTracking:HighDensity', ...
            ['The state inference is intended for low concentration. ' ...
            'Current density is %.3g molecules/um^2.'], ...
            opts.nMolecules / prod(opts.boxSizeUm));
    end
end
