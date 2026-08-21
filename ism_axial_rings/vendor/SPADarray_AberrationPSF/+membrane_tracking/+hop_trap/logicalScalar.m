function value = logicalScalar(value, name)
%LOGICALSCALAR Validate a logical scalar option.

    if ~(isscalar(value) && (islogical(value) || isnumeric(value)))
        error('simulateHopTrapMembraneTracking:BadOption', ...
            '%s must be a logical scalar.', name);
    end
    value = logical(value);
end
