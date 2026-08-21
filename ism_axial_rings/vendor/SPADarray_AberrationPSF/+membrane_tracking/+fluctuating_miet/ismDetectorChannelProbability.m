function [probability, dPdx, dPdy] = ismDetectorChannelProbability( ...
        position, detector, opts)
    import membrane_tracking.fluctuating_miet.*

    [raw, dRawDx, dRawDy] = ismDetectorRawResponse(position, detector, opts);
    normalization = sum(raw);
    if ~(normalization > realmin)
        probability = ones(detector.nChannels,1) / detector.nChannels;
        dPdx = zeros(detector.nChannels,1);
        dPdy = zeros(detector.nChannels,1);
        return;
    end
    probability = raw / normalization;
    if nargout > 1
        dPdx = (dRawDx*normalization - raw*sum(dRawDx)) / normalization^2;
        dPdy = (dRawDy*normalization - raw*sum(dRawDy)) / normalization^2;
    end
end
