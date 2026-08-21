function value = normalizeChoice(value, choices, name)
%NORMALIZECHOICE Validate and normalize a text enumeration.

    value = char(string(value));
    match = find(strcmpi(value, choices), 1);
    if isempty(match)
        error('simulateHopTrapMembraneTracking:BadOption', ...
            '%s must be one of: %s.', name, strjoin(choices, ', '));
    end
    value = choices{match};
end
