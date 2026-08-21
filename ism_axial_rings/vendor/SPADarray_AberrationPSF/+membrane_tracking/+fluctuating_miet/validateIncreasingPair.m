function pair = validateIncreasingPair(value, name)
    import membrane_tracking.fluctuating_miet.*

    pair = double(value(:).');
    if numel(pair) ~= 2 || any(~isfinite(pair)) || pair(2) <= pair(1)
        error('simulateFluctuatingMIETMembraneTracking:BadOption', ...
            '%s must be a finite increasing two-element vector.', name);
    end
end
