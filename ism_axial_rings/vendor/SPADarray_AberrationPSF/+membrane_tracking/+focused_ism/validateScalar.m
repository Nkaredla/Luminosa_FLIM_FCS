function value = validateScalar(value, name, lowerBound, inclusive)
    import membrane_tracking.focused_ism.*

    if ~isnumeric(value)
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must be numeric.', name);
    end
    value = double(value);
    if ~isscalar(value) || ~isfinite(value) || ~isreal(value)
        isAboveBound = false;
    else
        isAboveBound = value > lowerBound || ...
            (inclusive && value == lowerBound);
    end
    if ~isAboveBound
        if inclusive
            relation = '>=';
        else
            relation = '>';
        end
        error('simulateMembraneDiffusionParticleTracking:BadOption', ...
            '%s must be a finite real scalar %s %.6g.', ...
            name, relation, lowerBound);
    end
end
