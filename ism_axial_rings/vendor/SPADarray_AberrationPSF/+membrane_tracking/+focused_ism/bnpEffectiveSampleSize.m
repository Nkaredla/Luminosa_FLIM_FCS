function effectiveSampleSize = bnpEffectiveSampleSize(values)
    import membrane_tracking.focused_ism.*

    values = values(isfinite(values));
    n = numel(values);
    if n < 4
        effectiveSampleSize = NaN;
        return;
    end
    centered = values - mean(values);
    variance = sum(centered.^2) / n;
    if variance <= realmin
        effectiveSampleSize = 1;
        return;
    end

    autocorrelationSum = 0;
    for lag = 1:2:(n-2)
        rho1 = sum(centered(1:n-lag) .* centered(1+lag:n)) / ...
            ((n-lag) * variance);
        rho2 = sum(centered(1:n-lag-1) .* centered(2+lag:n)) / ...
            ((n-lag-1) * variance);
        pairSum = rho1 + rho2;
        if pairSum <= 0
            break;
        end
        autocorrelationSum = autocorrelationSum + pairSum;
    end
    effectiveSampleSize = min(n, max(1, ...
        n / (1 + 2 * autocorrelationSum)));
end
