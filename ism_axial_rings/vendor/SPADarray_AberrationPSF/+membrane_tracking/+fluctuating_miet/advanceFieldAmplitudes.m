function [amplitudeA, amplitudeB] = advanceFieldAmplitudes( ...
        amplitudeA, amplitudeB, modes, dt)
    import membrane_tracking.fluctuating_miet.*

% Exact Ornstein-Uhlenbeck transition, valid for any dt relative to tau.
    decay = exp(-dt ./ max(modes.relaxationTimeS, realmin));
    innovation = sqrt(max(modes.variance .* (1 - decay.^2), 0));
    amplitudeA = decay .* amplitudeA + innovation .* randn(modes.nModes, 1);
    amplitudeB = decay .* amplitudeB + innovation .* randn(modes.nModes, 1);
end
