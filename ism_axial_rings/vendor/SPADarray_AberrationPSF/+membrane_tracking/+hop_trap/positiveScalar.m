function value = positiveScalar(value, name, allowZero)
%POSITIVESCALAR Validate a finite positive or nonnegative scalar.

    if ~(isnumeric(value) && isscalar(value) && isfinite(value))
        error('simulateHopTrapMembraneTracking:BadOption', ...
            '%s must be a finite numeric scalar.', name);
    end
    if (allowZero && value < 0) || (~allowZero && value <= 0)
        error('simulateHopTrapMembraneTracking:BadOption', ...
            '%s is outside its allowed range.', name);
    end
    value = double(value);
end
