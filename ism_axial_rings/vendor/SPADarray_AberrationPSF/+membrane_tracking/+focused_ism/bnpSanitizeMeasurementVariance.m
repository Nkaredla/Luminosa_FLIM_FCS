function variance = bnpSanitizeMeasurementVariance(variance)
    import membrane_tracking.focused_ism.*

    variance = double(variance(:));
    valid = isfinite(variance) & variance > 0;
    if any(valid)
        fallback = median(variance(valid));
    else
        fallback = 0.05^2;
    end
    variance(~valid) = fallback;
    variance = max(variance, 1e-8);
end
