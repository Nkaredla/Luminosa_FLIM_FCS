function value = validateProbability(value, name)
    import membrane_tracking.focused_ism.*

    if ~isnumeric(value)
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must be numeric.', name);
    end
    value = double(value);
    if ~isscalar(value) || ~isfinite(value) || ~isreal(value) || ...
            value < 0 || value > 1
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must be a finite real scalar in [0, 1].', name);
    end
end
