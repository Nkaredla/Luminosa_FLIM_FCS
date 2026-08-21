function [trackTable, diffusion, opts] = refineTrackingAndDiffusion( ...
        localizationTable, opts)
    import membrane_tracking.focused_ism.*

    % The motion prior sets both JPDA likelihoods and hard gates. Iterate the
    % link/estimate cycle so the final result is not tied to a badly chosen
    % one-pass prior. Each reported result is kept with the exact prior and
    % gate that generated its track table.
    linkOpts = opts;
    nMaximum = opts.nTrackingRefinementIterations;
    priorHistory = nan(nMaximum, 1);
    estimateHistory = nan(nMaximum, 1);
    relativeChange = NaN;
    trackTable = table();
    diffusion = emptyDiffusionResult(opts);

    for iteration = 1:nMaximum
        priorHistory(iteration) = ...
            linkOpts.trackingPriorDiffusionUm2PerS;
        trackTable = linkLocalizations(localizationTable, linkOpts);
        diffusion = estimateDiffusionFromTracks(trackTable, linkOpts);
        estimate = diffusion.DFisherMLEUm2PerS;
        estimateHistory(iteration) = estimate;
        if ~isfinite(estimate) || estimate <= 0
            break;
        end

        relativeChange = abs(estimate - ...
            linkOpts.trackingPriorDiffusionUm2PerS) / ...
            max(estimate, realmin);
        if relativeChange < opts.trackingRefinementTolerance || ...
                iteration == nMaximum
            break;
        end

        linkOpts.trackingPriorDiffusionUm2PerS = estimate;
        if linkOpts.trackingMaxStepUmAuto
            brownianRmsStep = sqrt(4 * estimate * linkOpts.dtS);
            linkOpts.trackingMaxStepUm = max(4 * brownianRmsStep, ...
                0.5 * linkOpts.detectorPitchUm);
        end
    end

    iterationsUsed = iteration;
    opts.trackingPriorDiffusionUm2PerS = ...
        priorHistory(iterationsUsed);
    opts.trackingMaxStepUm = linkOpts.trackingMaxStepUm;
    diffusion.nTrackingRefinementIterations = iterationsUsed;
    diffusion.trackingPriorHistoryUm2PerS = ...
        priorHistory(1:iterationsUsed);
    diffusion.trackingEstimateHistoryUm2PerS = ...
        estimateHistory(1:iterationsUsed);
    diffusion.trackingFixedPointRelativeChange = relativeChange;
end
