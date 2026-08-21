function out = runSimulation(varargin)
%RUNSIMULATION Fluctuating-membrane MIET/ISM particle tracking.
%
%   This file is the readable orchestration layer. Numerical details
%   are implemented by the neighboring package functions.

    projectDirectory = fileparts(fileparts(fileparts( ...
        mfilename('fullpath'))));
    addpath(projectDirectory);
    import membrane_tracking.fluctuating_miet.*
    opts = parseInputs(varargin{:});
    if ~isempty(opts.rngSeed)
        rng(opts.rngSeed);
    end

    detector = makeDetector(opts);
    modes = makeFluctuationModes(opts);
    trajectories = simulateSurfaceTrajectories(modes, opts);
    [frames, photonEvents, photophysics] = simulatePhotonMovie( ...
        trajectories, detector, opts);
    localizationTable = localizeMovie(frames, photonEvents, ...
        photophysics, trajectories, detector, opts);
    surfaceModel = estimateFluctuatingSurfaceModel( ...
        localizationTable, trajectories, modes, opts);
    [trackTable, diffusion, opts] = refineTrackingAndDiffusion( ...
        localizationTable, trajectories, surfaceModel, modes, opts);

    figureSummary = [];
    if opts.makeFigure
        figureSummary = plotSummary(trajectories, frames, detector, ...
            localizationTable, trackTable, surfaceModel, diffusion, opts);
    end
    if opts.writeOutputs || opts.saveFigures
        ensureDir(opts.outputDir);
        if opts.writeOutputs
            writeAllOutputs(opts, detector, modes, trajectories, frames, ...
                photonEvents, photophysics, localizationTable, ...
                surfaceModel, trackTable, diffusion);
        end
        if opts.saveFigures && ~isempty(figureSummary)
            safeExportFigure(figureSummary, fullfile(opts.outputDir, ...
                'fluctuating_miet_summary.png'), 180);
        end
    end
    if opts.verbose
        reportToCommandWindow(surfaceModel, diffusion, photophysics, opts);
    end

    out = struct();
    out.options = opts;
    out.detector = detector;
    out.modes = modes;
    out.trajectories = trajectories;
    out.frames = frames;
    out.photonEvents = photonEvents;
    out.photophysics = photophysics;
    out.localizationTable = localizationTable;
    out.surfaceModel = surfaceModel;
    out.trackTable = trackTable;
    out.diffusion = diffusion;
    out.figureSummary = figureSummary;
end
