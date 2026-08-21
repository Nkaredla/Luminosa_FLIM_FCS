function rhat = bnpSplitRhat(values)
    import membrane_tracking.focused_ism.*

    values = values(isfinite(values));
    halfLength = floor(numel(values) / 2);
    if halfLength < 2
        rhat = NaN;
        return;
    end
    chains = [values(1:halfLength), ...
        values(end-halfLength+1:end)];
    withinVariance = mean(var(chains, 0, 1));
    betweenVariance = halfLength * var(mean(chains, 1), 0);
    if withinVariance <= realmin
        if betweenVariance <= realmin
            rhat = 1;
        else
            rhat = Inf;
        end
        return;
    end
    varianceEstimate = ((halfLength - 1) / halfLength) * ...
        withinVariance + betweenVariance / halfLength;
    rhat = sqrt(max(varianceEstimate / withinVariance, 0));
end
