function meanCounts = bnpRebuildMeanCounts(state)
    import membrane_tracking.focused_ism.*

    meanCounts = state.backgroundCountsPerChannel * ...
        ones(state.nFrames, state.nChannels);
    if any(state.active)
        meanCounts = meanCounts + reshape(sum( ...
            state.contributions(state.active,:,:), 1), ...
            state.nFrames, state.nChannels);
    end
end
