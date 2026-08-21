function choice = normalizeChoice(value, allowed, name)
    import membrane_tracking.focused_ism.*

    if ~(ischar(value) && isrow(value)) && ...
            ~(isstring(value) && isscalar(value))
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must be a character vector or scalar string.', name);
    end
    choice = char(value);
    matches = strcmpi(choice, allowed);
    if ~any(matches)
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must be one of: %s.', name, strjoin(allowed, ', '));
    end
    choice = allowed{find(matches, 1, 'first')};
end
