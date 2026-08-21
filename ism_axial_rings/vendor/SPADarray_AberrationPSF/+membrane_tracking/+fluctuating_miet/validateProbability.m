function value = validateProbability(value, name)
    import membrane_tracking.fluctuating_miet.*

    value = validateFiniteScalar(value, name);
    if value < 0 || value > 1
        error('simulateFluctuatingMIETMembraneTracking:BadOption', ...
            '%s must lie in [0,1].', name);
    end
end
