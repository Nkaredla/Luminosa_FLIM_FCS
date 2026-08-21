function value = validateLogicalScalar(value, name)
    import membrane_tracking.curved_miet.*

    if ~(isscalar(value) && (islogical(value) || isnumeric(value)))
        error('simulateMIETCurvedMembraneTracking:BadOption', ...
            '%s must be a logical scalar.', name);
    end
    value = logical(value);
end
