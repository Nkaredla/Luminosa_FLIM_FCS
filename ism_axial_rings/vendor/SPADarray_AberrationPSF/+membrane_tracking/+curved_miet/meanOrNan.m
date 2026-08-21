function value = meanOrNan(values)
    import membrane_tracking.curved_miet.*

    values = values(isfinite(values));
    if isempty(values)
        value = NaN;
    else
        value = mean(values);
    end
end
