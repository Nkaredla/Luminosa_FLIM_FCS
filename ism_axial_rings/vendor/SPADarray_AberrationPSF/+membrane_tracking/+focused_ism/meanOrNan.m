function value = meanOrNan(x)
    import membrane_tracking.focused_ism.*

    if isempty(x)
        value = NaN;
    else
        value = mean(x);
    end
end
