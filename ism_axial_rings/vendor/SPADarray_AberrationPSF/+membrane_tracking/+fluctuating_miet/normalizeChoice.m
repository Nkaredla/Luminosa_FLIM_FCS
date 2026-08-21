function choice = normalizeChoice(value, allowed, name)
    import membrane_tracking.fluctuating_miet.*

    choice = char(value);
    match = strcmpi(choice, allowed);
    if ~any(match)
        error('simulateFluctuatingMIETMembraneTracking:BadOption', ...
            '%s must be one of: %s.', name, strjoin(allowed, ', '));
    end
    choice = allowed{find(match, 1, 'first')};
end
