function value = roundScalar(value, name, minimum)
    import membrane_tracking.fluctuating_miet.*

    value = round(validateFiniteScalar(value, name));
    if value < minimum
        error('simulateFluctuatingMIETMembraneTracking:BadOption', ...
            '%s must be an integer >= %d.', name, minimum);
    end
end
