function count = bnpMeanBrightEmitterCount(state, opts)
    import membrane_tracking.focused_ism.*

    countsPerFrame = bnpBrightEmitterCounts(state, opts);
    count = mean(countsPerFrame);
end
