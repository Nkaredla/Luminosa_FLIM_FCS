function out = runSimulation(varargin)
%RUNSIMULATION Orchestrate flat-membrane hop/trap simulation and inference.

    projectDirectory = fileparts(fileparts(fileparts( ...
        mfilename('fullpath'))));
    addpath(projectDirectory);
    import membrane_tracking.hop_trap.*

    opts = parseInputs(varargin{:});
    if ~isempty(opts.rngSeed)
        rng(opts.rngSeed);
    end

    detector = membrane_tracking.focused_ism.makeDetector(opts);
    maxIdentifiableEmitters = floor((detector.nChannels - 1) / 3);
    opts.requestedMaxEmittersPerFitRegion = ...
        opts.maxEmittersPerFitRegion;
    opts.maxEmittersPerFitRegion = min( ...
        opts.maxEmittersPerFitRegion, maxIdentifiableEmitters);

    mesh = makePeriodicVoronoiMesh(opts);
    trajectories = simulateHeterogeneousTrajectories(mesh, opts);
    [frames, photonEvents, photophysics] = ...
        membrane_tracking.focused_ism.simulatePhotonMovie( ...
        trajectories.positionsUm, detector, opts);
    localizationTable = membrane_tracking.focused_ism.localizeMovie( ...
        frames, detector, opts, trajectories.positionsUm, ...
        photophysics.detectedPhotonCountsInFov);
    [trackTable, baselineDiffusion, opts] = ...
        membrane_tracking.focused_ism.refineTrackingAndDiffusion( ...
        localizationTable, opts);

    stepTable = collectIndependentTrackSteps( ...
        trackTable, trajectories, opts);
    msdTable = buildMsdCurve(trackTable, opts);
    modelComparison = fitDiffusionNatureModels( ...
        stepTable, msdTable, opts);
    stickyBnp = runStickyHdpHmm(stepTable, ...
        modelComparison.freeStep, opts);
    classification = classifyDiffusionNature( ...
        modelComparison.freeStep, modelComparison.trapHmm, ...
        modelComparison.hopMsd, stickyBnp, opts);
    modelComparison.classification = classification;

    figureSummary = [];
    if opts.makeFigure
        figureSummary = plotSummary(trajectories, mesh, frames, detector, ...
            localizationTable, trackTable, stepTable, msdTable, ...
            modelComparison, stickyBnp, classification, opts);
    end

    out = struct();
    out.options = opts;
    out.detector = detector;
    out.mesh = mesh;
    out.trajectories = trajectories;
    out.frames = frames;
    out.photonEvents = photonEvents;
    out.photophysics = photophysics;
    out.localizationTable = localizationTable;
    out.trackTable = trackTable;
    out.baselineDiffusion = baselineDiffusion;
    out.stepTable = stepTable;
    out.msdTable = msdTable;
    out.modelComparison = modelComparison;
    out.stickyBnp = stickyBnp;
    out.figureSummary = figureSummary;

    if opts.writeOutputs || opts.saveFigures
        if ~exist(opts.outputDir, 'dir')
            mkdir(opts.outputDir);
        end
        if opts.writeOutputs
            writeOutputs(opts.outputDir, out);
        end
        if opts.saveFigures && ~isempty(figureSummary)
            membrane_tracking.focused_ism.safeExportFigure( ...
                figureSummary, fullfile(opts.outputDir, ...
                'hop_trap_tracking_summary.png'), 180);
        end
    end

    if opts.verbose
        printSummary(out);
    end
end
