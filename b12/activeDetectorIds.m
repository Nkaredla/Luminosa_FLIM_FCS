function detectorIds = activeDetectorIds(photonChannels, discoveredIds, blockSize)
%ACTIVEDETECTORIDS Return sorted raw detector IDs from the full photon list.
% Using the complete list avoids missing a detector that was absent from the
% short discovery block read by PTU_FLIM_GPU.

    if nargin < 3 || isempty(blockSize)
        blockSize = 1e6;
    end
    present = false(256, 1);
    for first = 1:blockSize:numel(photonChannels)
        last = min(first + blockSize - 1, numel(photonChannels));
        ids = unique(uint8(photonChannels(first:last)));
        present(double(ids) + 1) = true;
    end
    discoveredIds = double(discoveredIds(:));
    discoveredIds = discoveredIds(isfinite(discoveredIds) & ...
        discoveredIds >= 0 & discoveredIds <= 255 & ...
        discoveredIds == round(discoveredIds));
    present(discoveredIds + 1) = true;
    detectorIds = find(present(1:32)) - 1;
end
