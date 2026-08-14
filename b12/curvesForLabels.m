function [curves, photonCounts] = curvesForLabels(photonBins, pixelIndices, labels, gateLength)
%CURVESFORLABELS Build one TCSPC curve per connected labelled region.

    regionCount = max(labels(:));
    curves = zeros(gateLength, regionCount);
    photonCounts = zeros(regionCount, 1);
    if regionCount == 0 || isempty(photonBins)
        return;
    end

    regionId = double(labels(round(pixelIndices)));
    inside = regionId > 0;
    curves = accumarray([round(photonBins(inside)), regionId(inside)], 1, ...
        [gateLength regionCount], @sum, 0);
    photonCounts = sum(curves, 1).';
end
