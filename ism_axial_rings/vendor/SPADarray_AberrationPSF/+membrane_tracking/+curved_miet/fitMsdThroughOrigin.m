function [diffusion, slope] = fitMsdThroughOrigin(lagTime, msd, nPairs)
    import membrane_tracking.curved_miet.*

    lagTime = lagTime(:);
    msd = msd(:);
    nPairs = nPairs(:);
    valid = isfinite(lagTime) & isfinite(msd) & nPairs > 0;
    lagTime = lagTime(valid);
    msd = msd(valid);
    nPairs = nPairs(valid);
    if isempty(lagTime)
        diffusion = NaN;
        slope = NaN;
        return;
    end
    denominator = sum(nPairs .* lagTime.^2);
    if denominator <= 0
        diffusion = NaN;
        slope = NaN;
        return;
    end
    slope = sum(nPairs .* lagTime .* msd) / denominator;
    slope = max(slope, 0);
    diffusion = slope / 4;
end
