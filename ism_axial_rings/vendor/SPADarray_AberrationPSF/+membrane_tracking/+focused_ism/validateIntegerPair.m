function pair = validateIntegerPair(value, name)
    import membrane_tracking.focused_ism.*

    pair = validatePair(value, name);
    rounded = round(pair);
    if any(pair ~= rounded)
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must contain integer values.', name);
    end
    pair = rounded;
end
