function value = validateLogicalScalar(value, name)
    import membrane_tracking.fluctuating_miet.*

    if ~(isscalar(value) && (islogical(value) || isnumeric(value)))
        error('simulateFluctuatingMIETMembraneTracking:BadOption', ...
            '%s must be a logical scalar.', name);
    end
    value = logical(value);
end
