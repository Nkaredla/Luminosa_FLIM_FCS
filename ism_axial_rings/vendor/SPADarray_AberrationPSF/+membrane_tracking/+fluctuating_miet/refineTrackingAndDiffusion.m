function [trackTable, diffusion, opts] = refineTrackingAndDiffusion( ...
        localizationTable, trajectories, surfaceModel, modes, opts)
    import membrane_tracking.fluctuating_miet.*

%REFINETRACKINGANDDIFFUSION Iterate the motion prior to a fixed point.
    linkOpts = opts;
    maximumIterations = opts.nTrackingRefinementIterations;
    priorHistory = nan(maximumIterations,1);
    estimateHistory = nan(maximumIterations,1);
    relativeChange = NaN;
    trackTable = table();
    diffusion = emptyDiffusionResult(opts);

    for iteration = 1:maximumIterations
        priorHistory(iteration) = ...
            linkOpts.trackingPriorDiffusionUm2PerS;
        trackTable = linkLocalizations(localizationTable, linkOpts);
        diffusion = estimateSurfaceDiffusion(trackTable, trajectories, ...
            surfaceModel, modes, linkOpts);
        estimate = diffusion.DFisherMLEUm2PerS;
        estimateHistory(iteration) = estimate;
        if ~isfinite(estimate) || estimate <= 0
            break;
        end
        relativeChange = abs(estimate- ...
            linkOpts.trackingPriorDiffusionUm2PerS) / ...
            max(estimate,realmin);
        if relativeChange < opts.trackingRefinementTolerance || ...
                iteration == maximumIterations
            break;
        end
        linkOpts.trackingPriorDiffusionUm2PerS = estimate;
    end

    iterationsUsed = iteration;
    opts.trackingPriorDiffusionUm2PerS = priorHistory(iterationsUsed);
    diffusion.nTrackingRefinementIterations = iterationsUsed;
    diffusion.trackingPriorHistoryUm2PerS = ...
        priorHistory(1:iterationsUsed);
    diffusion.trackingEstimateHistoryUm2PerS = ...
        estimateHistory(1:iterationsUsed);
    diffusion.trackingFixedPointRelativeChange = relativeChange;
end
