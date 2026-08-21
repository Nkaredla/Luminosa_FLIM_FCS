function value = validateLogicalScalar(value, name)
    import membrane_tracking.focused_ism.*

    if ~(islogical(value) || isnumeric(value)) || ~isscalar(value) || ...
            ~isfinite(double(value)) || ~ismember(double(value), [0 1])
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must be a logical scalar.', name);
    end
    value = logical(value);
end
