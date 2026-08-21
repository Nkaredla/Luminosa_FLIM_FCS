function value = meanOrNan(values)
    import membrane_tracking.fluctuating_miet.*

    values = values(isfinite(values));
    if isempty(values)
        value = NaN;
    else
        value = mean(values);
    end
end
