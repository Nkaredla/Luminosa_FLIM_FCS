function histogram = detectorHistogram(photonBins, photonChannels, detectorIds, nBins)
%DETECTORHISTOGRAM Build one global TCSPC histogram per raw detector.

    timeBin = round(double(photonBins(:)));
    channel = double(photonChannels(:));
    nPhotons = min(numel(timeBin), numel(channel));
    timeBin = timeBin(1:nPhotons);
    channel = channel(1:nPhotons);

    histogram = zeros(nBins, numel(detectorIds));
    for detector = 1:numel(detectorIds)
        use = channel == detectorIds(detector) & timeBin >= 1 & timeBin <= nBins;
        histogram(:, detector) = accumarray(timeBin(use), 1, [nBins 1], @sum, 0);
    end
end
