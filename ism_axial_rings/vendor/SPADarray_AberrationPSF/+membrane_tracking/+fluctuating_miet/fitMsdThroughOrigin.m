function D = fitMsdThroughOrigin(lagTime, msd, nPairs)
    import membrane_tracking.fluctuating_miet.*

    lagTime = lagTime(:); msd = msd(:); nPairs = nPairs(:);
    valid = isfinite(lagTime) & isfinite(msd) & nPairs > 0;
    lagTime = lagTime(valid); msd = msd(valid); nPairs = nPairs(valid);
    denominator = sum(nPairs .* lagTime.^2);
    if isempty(lagTime) || denominator <= 0
        D = NaN;
        return;
    end
    D = max(sum(nPairs .* lagTime .* msd) / denominator, 0) / 4;
end
