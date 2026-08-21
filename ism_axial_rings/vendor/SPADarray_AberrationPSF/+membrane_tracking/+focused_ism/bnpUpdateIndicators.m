function [state, nChanges] = bnpUpdateIndicators(state, opts)
    import membrane_tracking.focused_ism.*

    nChanges = 0;
    for candidate = randperm(state.nCandidates)
        contribution = reshape(state.contributions(candidate,:,:), ...
            state.nFrames, state.nChannels);
        if state.active(candidate)
            meanOff = state.meanCounts - contribution;
        else
            meanOff = state.meanCounts;
        end
        meanOn = meanOff + contribution;

        nOtherActive = sum(state.active) - state.active(candidate);
        logOn = bnpPoissonLogLikelihood(state.data, meanOn, opts) + ...
            log(opts.bnpBetaPriorA + nOtherActive);
        logOff = bnpPoissonLogLikelihood(state.data, meanOff, opts) + ...
            log(opts.bnpBetaPriorB + state.nCandidates - 1 - nOtherActive);
        logOdds = logOn - logOff;
        if logOdds >= 0
            probabilityOn = 1 / (1 + exp(-logOdds));
        else
            odds = exp(logOdds);
            probabilityOn = odds / (1 + odds);
        end
        newState = rand < probabilityOn;
        if newState ~= state.active(candidate)
            nChanges = nChanges + 1;
        end
        state.active(candidate) = newState;
        if newState
            state.meanCounts = meanOn;
        else
            state.meanCounts = meanOff;
        end
    end
end
