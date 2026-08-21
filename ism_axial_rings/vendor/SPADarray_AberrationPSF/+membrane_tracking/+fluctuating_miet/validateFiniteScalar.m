function value = validateFiniteScalar(value, name)
    import membrane_tracking.fluctuating_miet.*

    value = double(value);
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value)
        error('simulateFluctuatingMIETMembraneTracking:BadOption', ...
            '%s must be a finite real scalar.', name);
    end
end
