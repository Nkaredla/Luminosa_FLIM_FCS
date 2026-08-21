function [state, accepted] = bnpUpdateBrightness(state, opts)
    import membrane_tracking.focused_ism.*

    if strcmp(opts.bnpBrightnessMode, 'fixed')
        accepted = NaN;
        return;
    end

    oldLogScale = log(max(state.brightnessScale, 1e-12));
    proposedLogScale = oldLogScale + ...
        opts.bnpLogNuisanceProposalSigma * randn;
    proposedLogScale = max(min(proposedLogScale, 20), -20);
    scaleRatio = exp(proposedLogScale - oldLogScale);
    signalMean = state.meanCounts - state.backgroundCountsPerChannel;
    proposedMean = state.backgroundCountsPerChannel + ...
        scaleRatio * signalMean;

    if strcmp(opts.bnpBrightnessMode, 'calibrated')
        priorMean = log(opts.bnpCalibratedBrightnessScale);
        priorSigma = sqrt(log(1 + ...
            opts.bnpBrightnessCalibrationRelativeStd^2));
    else
        priorMean = log(opts.bnpCalibratedBrightnessScale);
        priorSigma = opts.bnpLogBrightnessPriorSigma;
    end
    priorVariance = priorSigma^2;
    logPriorRatio = -0.5 * ...
        ((proposedLogScale - priorMean)^2 - ...
        (oldLogScale - priorMean)^2) / priorVariance;
    logRatio = bnpPoissonLogLikelihood( ...
        state.data, proposedMean, opts) - ...
        bnpPoissonLogLikelihood(state.data, state.meanCounts, opts) + ...
        logPriorRatio;
    accepted = log(rand) < min(0, logRatio);
    if accepted
        state.brightnessScale = exp(proposedLogScale);
        state.contributions = scaleRatio * state.contributions;
        state.meanCounts = proposedMean;
    end
end
