function value = expClamped(value)
    import membrane_tracking.curved_miet.*

    value = exp(min(max(value, -50), 50));
end
