function value = validateScalar(value, name, lowerBound, inclusive)
    import membrane_tracking.curved_miet.*

    value = validateFiniteScalar(value, name);
    if (inclusive && value < lowerBound) || ...
            (~inclusive && value <= lowerBound)
        comparator = 'greater than';
        if inclusive
            comparator = 'greater than or equal to';
        end
        error('simulateMIETCurvedMembraneTracking:BadOption', ...
            '%s must be %s %.4g.', name, comparator, lowerBound);
    end
end
