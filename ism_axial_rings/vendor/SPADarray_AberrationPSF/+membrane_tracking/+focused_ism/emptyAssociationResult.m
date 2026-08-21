function association = emptyAssociationResult(nTracks, nDetections)
    import membrane_tracking.focused_ism.*

    association = struct();
    association.mapDetectionForTrack = zeros(nTracks, 1);
    association.mapAssociationProbability = ones(nTracks, 1);
    association.newDetectionProbability = ones(nDetections, 1);
    association.entropyForTrack = zeros(nTracks, 1);
    association.entropyForDetection = zeros(nDetections, 1);
    association.hypothesisCountForTrack = ones(nTracks, 1);
    association.hypothesisCountForDetection = ones(nDetections, 1);
end
