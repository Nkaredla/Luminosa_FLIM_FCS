function [amplitudeA, amplitudeB] = initialFieldAmplitudes(modes)
    import membrane_tracking.fluctuating_miet.*

    amplitudeA = sqrt(modes.variance) .* randn(modes.nModes, 1);
    amplitudeB = sqrt(modes.variance) .* randn(modes.nModes, 1);
end
