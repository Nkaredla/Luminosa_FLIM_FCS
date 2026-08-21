function quantiles = bnpEmpiricalQuantile(values, probabilities)
    import membrane_tracking.focused_ism.*

    values = sort(values(isfinite(values)));
    quantiles = nan(size(probabilities));
    if isempty(values)
        return;
    end
    n = numel(values);
    for k = 1:numel(probabilities)
        index = 1 + (n - 1) * probabilities(k);
        lower = floor(index);
        upper = ceil(index);
        fraction = index - lower;
        quantiles(k) = (1 - fraction) * values(lower) + ...
            fraction * values(upper);
    end
end
