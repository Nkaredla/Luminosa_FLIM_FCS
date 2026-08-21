function [probability, dPdx, dPdy] = ismDetectorChannelProbability( ...
        positionXY, detector, opts)
    import membrane_tracking.curved_miet.*

    [raw, dRawDx, dRawDy] = ismDetectorRawResponse( ...
        positionXY, detector, opts);
    normalization = sum(raw);
    if ~(normalization > realmin)
        probability = ones(detector.nChannels,1) / detector.nChannels;
        dPdx = zeros(detector.nChannels,1);
        dPdy = zeros(detector.nChannels,1);
        return;
    end

    probability = raw / normalization;
    if nargout > 1
        dNormalizationX = sum(dRawDx);
        dNormalizationY = sum(dRawDy);
        dPdx = (dRawDx * normalization - raw * dNormalizationX) / ...
            normalization^2;
        dPdy = (dRawDy * normalization - raw * dNormalizationY) / ...
            normalization^2;
    end
end
