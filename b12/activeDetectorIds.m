function detectorIds = activeDetectorIds(photonChannels, discoveredIds)
%ACTIVEDETECTORIDS Return sorted raw detector IDs from the full photon list.
% Using the complete list avoids missing a detector that was absent from the
% short discovery block read by PTU_FLIM_GPU.

    detectorIds = unique([double(photonChannels(:)); double(discoveredIds(:))]);
    detectorIds = detectorIds(isfinite(detectorIds) & detectorIds >= 0 & detectorIds < 32);
    detectorIds = sort(detectorIds(:));
end
