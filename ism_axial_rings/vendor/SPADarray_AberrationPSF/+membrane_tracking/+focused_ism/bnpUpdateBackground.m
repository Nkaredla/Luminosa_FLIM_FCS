function [state, accepted] = bnpUpdateBackground(state, opts)
    import membrane_tracking.focused_ism.*

    oldBackground = max(state.backgroundCountsPerChannel, 1e-12);
    proposedLogBackground = log(oldBackground) + ...
        opts.bnpLogNuisanceProposalSigma * randn;
    proposedBackground = exp(max(min(proposedLogBackground, 30), -30));
    proposedMean = state.meanCounts + ...
        (proposedBackground - oldBackground);
    logRatio = bnpPoissonLogLikelihood( ...
        state.data, proposedMean, opts) - ...
        bnpPoissonLogLikelihood(state.data, state.meanCounts, opts);
    accepted = log(rand) < min(0, logRatio);
    if accepted
        state.backgroundCountsPerChannel = proposedBackground;
        state.meanCounts = proposedMean;
    end
end
