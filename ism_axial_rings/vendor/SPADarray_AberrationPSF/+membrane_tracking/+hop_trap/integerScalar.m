function value = integerScalar(value, name, minimum)
%INTEGERSCALAR Validate and round a finite integer-valued option.

    if ~(isnumeric(value) && isscalar(value) && isfinite(value))
        error('simulateHopTrapMembraneTracking:BadOption', ...
            '%s must be a finite integer scalar.', name);
    end
    value = round(double(value));
    if value < minimum
        error('simulateHopTrapMembraneTracking:BadOption', ...
            '%s must be at least %d.', name, minimum);
    end
end
