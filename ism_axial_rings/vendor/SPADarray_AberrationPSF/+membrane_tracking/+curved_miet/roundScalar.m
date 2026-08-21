function value = roundScalar(value, name, minimum)
    import membrane_tracking.curved_miet.*

    value = validateFiniteScalar(value, name);
    value = round(value);
    if value < minimum
        error('simulateMIETCurvedMembraneTracking:BadOption', ...
            '%s must be an integer >= %d.', name, minimum);
    end
end
