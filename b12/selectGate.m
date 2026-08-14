function [relativeBins, pixelIndices] = selectGate(photonBins, photonChannels, ...
        photonColumns, photonLines, detectorId, gateStart, gateLength, ...
        periodBins, nx, ny)
%SELECTGATE Select photons from one detector and one PIE laser window.
% Modulo timing makes the function safe if a gate crosses the TCSPC period.

    timeBin = round(double(photonBins(:)));
    channel = double(photonChannels(:));
    x = round(double(photonColumns(:)));
    y = round(double(photonLines(:)));
    nPhotons = min([numel(timeBin), numel(channel), numel(x), numel(y)]);
    timeBin = timeBin(1:nPhotons);
    channel = channel(1:nPhotons);
    x = x(1:nPhotons);
    y = y(1:nPhotons);

    relative = mod(timeBin - gateStart, periodBins) + 1;
    keep = channel == detectorId & relative <= gateLength & ...
        x >= 1 & x <= nx & y >= 1 & y <= ny;
    relativeBins = relative(keep);
    pixelIndices = x(keep) + (y(keep) - 1) * nx;
end
