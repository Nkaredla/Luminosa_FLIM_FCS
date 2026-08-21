function out = runSimulation(varargin)
%RUNSIMULATION Static curved-membrane MIET/ISM particle tracking.
%
%   This file is the readable orchestration layer. Numerical details
%   are implemented by the neighboring package functions.

    projectDirectory = fileparts(fileparts(fileparts( ...
        mfilename('fullpath'))));
    addpath(projectDirectory);
    import membrane_tracking.curved_miet.*

    opts = parseInputs(varargin{:});
    if ~isempty(opts.rngSeed)
        rng(opts.rngSeed);
    end

    detector = makeDetector(opts);
    trajectories = simulateSurfaceTrajectories(opts);
    [frames, photonEvents, photophysics] = simulateMIETPhotonMovie( ...
        trajectories, detector, opts);
    localizationTable = localizeMIETMovie(frames, photonEvents, ...
        photophysics, trajectories, detector, opts);
    [curvature, localizationTable] = estimateCurvatureFromTCSPC( ...
        localizationTable, photonEvents, photophysics, detector, opts);
    trackTable = linkLowDensityLocalizations(localizationTable, opts);
    diffusion = estimateSurfaceDiffusion(trackTable, trajectories, ...
        curvature, opts);

    figureSummary = [];
    if opts.makeFigure
        figureSummary = plotSummary(trajectories, frames, detector, ...
            localizationTable, trackTable, curvature, diffusion, opts);
    end

    if opts.writeOutputs || opts.saveFigures
        ensureDir(opts.outputDir);
        if opts.writeOutputs
            writeOutputs(opts.outputDir, opts, detector, trajectories, ...
                frames, photonEvents, photophysics, localizationTable, ...
                trackTable, curvature, diffusion);
        end
        if opts.saveFigures && ~isempty(figureSummary)
            safeExportFigure(figureSummary, fullfile(opts.outputDir, ...
                'miet_curved_membrane_tracking_summary.png'), 180);
        end
    end

    if opts.verbose
        fprintf(['[curved MIET membrane] true apex curvature = %.4g 1/um, ' ...
            'estimated = %.4g +/- %.3g 1/um; true tip height = %.4g um, ' ...
            'estimated = %.4g +/- %.3g um.\n'], ...
            opts.curvaturePerUm, curvature.curvaturePerUm, ...
            curvature.curvatureSigmaPerUm, opts.tipHeightUm, ...
            curvature.tipHeightUm, curvature.tipHeightSigmaUm);
        fprintf(['[surface diffusion] true D = %.4g um^2/s, ' ...
            'corrected-MSD D = %.4g um^2/s, Fisher/MLE D = %.4g +/- ' ...
            '%.3g um^2/s, sampled truth = %.4g um^2/s; %d localizations, ' ...
            '%d tracks.\n'], opts.diffusionUm2PerS, ...
            diffusion.DMsdCorrectedUm2PerS, ...
            diffusion.DFisherMLEUm2PerS, ...
            diffusion.diffusionFisherSigmaUm2PerS, ...
            diffusion.truthValidation.fullSurfaceLag1DiffusionUm2PerS, ...
            height(localizationTable), diffusion.nTracksUsed);
        if photophysics.fractionFramesWithMultipleBrightMolecules > 0
            fprintf(['[low-density audit] %.2f%% of frames contained more ' ...
                'than one molecule above the localization photon threshold; ' ...
                'those frames can violate the single-emitter analysis.\n'], ...
                100 * photophysics.fractionFramesWithMultipleBrightMolecules);
        end
    end

    out = struct();
    out.options = opts;
    out.detector = detector;
    out.trajectories = trajectories;
    out.frames = frames;
    out.photonEvents = photonEvents;
    out.photophysics = photophysics;
    out.localizationTable = localizationTable;
    out.curvature = curvature;
    out.trackTable = trackTable;
    out.diffusion = diffusion;
    out.figureSummary = figureSummary;
end
