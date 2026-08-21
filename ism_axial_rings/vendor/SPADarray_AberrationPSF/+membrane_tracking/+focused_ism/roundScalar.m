function value = roundScalar(value, name, minValue)
    import membrane_tracking.focused_ism.*

    if ~isnumeric(value)
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must be a numeric scalar integer >= %d.', name, minValue);
    end
    value = double(value);
    if ~isscalar(value) || ~isfinite(value) || ~isreal(value) || ...
            value ~= round(value) || value < minValue
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must be a scalar integer >= %d.', name, minValue);
    end
    value = round(value);
end
