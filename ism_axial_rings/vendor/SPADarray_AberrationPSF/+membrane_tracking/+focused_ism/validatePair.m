function pair = validatePair(value, name)
    import membrane_tracking.focused_ism.*

    if ~isnumeric(value)
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must contain numeric values.', name);
    end
    pair = double(value(:)).';
    if numel(pair) ~= 2 || ~isreal(pair) || any(~isfinite(pair)) || ...
            any(pair <= 0)
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must contain two positive finite values.', name);
    end
end
