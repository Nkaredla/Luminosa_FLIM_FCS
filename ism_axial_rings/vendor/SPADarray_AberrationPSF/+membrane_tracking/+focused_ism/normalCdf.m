function y = normalCdf(x)
    import membrane_tracking.focused_ism.*

    y = 0.5 * erfc(-x / sqrt(2));
end
