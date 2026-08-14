function curve = curveForMask(photonBins, pixelIndices, mask, gateLength)
%CURVEFORMASK Sum a gated photon list over one binary spatial mask.

    if isempty(photonBins)
        curve = zeros(gateLength, 1);
        return;
    end
    inside = logical(mask(round(pixelIndices)));
    curve = accumarray(round(photonBins(inside)), 1, [gateLength 1], @sum, 0);
end
