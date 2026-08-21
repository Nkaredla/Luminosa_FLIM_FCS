function value = validateScalar(value, name, bound, inclusive)
    import membrane_tracking.fluctuating_miet.*

    value = validateFiniteScalar(value, name);
    if (inclusive && value < bound) || (~inclusive && value <= bound)
        error('simulateFluctuatingMIETMembraneTracking:BadOption', ...
            '%s must be greater than%s %.4g.', name, ...
            repmat(' or equal to', 1, inclusive), bound);
    end
end
