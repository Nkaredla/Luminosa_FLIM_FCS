function histogram = detectorHistogram(photonBins, photonChannels, ...
        detectorIds, nBins, blockSize)
%DETECTORHISTOGRAM Build one global TCSPC histogram per raw detector.

    if nargin < 5 || isempty(blockSize)
        blockSize = 1e6;
    end
    nPhotons = min(numel(photonBins), numel(photonChannels));
    histogram = zeros(nBins, numel(detectorIds));
    channelLut = zeros(256, 1, 'uint16');
    channelLut(double(detectorIds(:)) + 1) = uint16(1:numel(detectorIds));
    for first = 1:blockSize:nPhotons
        last = min(first + blockSize - 1, nPhotons);
        indices = first:last;
        timeBin = round(double(photonBins(indices)));
        localChannel = double(channelLut(double(photonChannels(indices)) + 1));
        valid = timeBin >= 1 & timeBin <= nBins & localChannel >= 1;
        if ~any(valid)
            continue;
        end
        linear = timeBin(valid) + (localChannel(valid) - 1) * nBins;
        histogram = histogram + reshape(accumarray(linear, 1, ...
            [nBins * numel(detectorIds) 1], @sum, 0), ...
            [nBins numel(detectorIds)]);
    end
end
