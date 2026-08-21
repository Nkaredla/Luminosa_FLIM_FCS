function [probability, dPdx, dPdy] = ismDetectorChannelProbability( ...
        positionUm, detector, opts)
    import membrane_tracking.focused_ism.*

    [raw, dRawDx, dRawDy] = ismDetectorRawResponse( ...
        positionUm, detector, opts);
    normalization = sum(raw);
    if ~(normalization > realmin)
        probability = ones(detector.nChannels, 1) / detector.nChannels;
        dPdx = zeros(detector.nChannels, 1);
        dPdy = zeros(detector.nChannels, 1);
        return;
    end

    % Condition on a photon reaching the array. The fitted amplitude is then
    % the total captured signal, while this vector describes its channel
    % allocation.
    probability = raw / normalization;
    if nargout > 1
        derivativeNormalizationX = sum(dRawDx);
        derivativeNormalizationY = sum(dRawDy);
        dPdx = (dRawDx * normalization - ...
            raw * derivativeNormalizationX) / normalization^2;
        dPdy = (dRawDy * normalization - ...
            raw * derivativeNormalizationY) / normalization^2;
    end
end
