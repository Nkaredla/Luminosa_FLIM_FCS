function value = probabilityScalar(value, name)
%PROBABILITYSCALAR Validate a probability in the closed unit interval.

    import membrane_tracking.hop_trap.*

    value = positiveScalar(value, name, true);
    if value > 1
        error('simulateHopTrapMembraneTracking:BadOption', ...
            '%s must lie between zero and one.', name);
    end
end
