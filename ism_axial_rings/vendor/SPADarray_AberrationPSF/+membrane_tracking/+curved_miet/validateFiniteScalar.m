function value = validateFiniteScalar(value, name)
    import membrane_tracking.curved_miet.*

    value = double(value);
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value)
        error('simulateMIETCurvedMembraneTracking:BadOption', ...
            '%s must be a finite real scalar.', name);
    end
end
