function value = medianOrNan(x)
    import membrane_tracking.focused_ism.*

    if isempty(x)
        value = NaN;
    else
        value = median(x);
    end
end
