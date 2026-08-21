function y = expClamped(x)
    import membrane_tracking.focused_ism.*

    y = exp(min(max(x, -40), 40));
end
