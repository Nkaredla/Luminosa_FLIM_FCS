function association = solveJpdaAssociations(logRatio, logMiss, ...
        maxHypotheses)
    import membrane_tracking.focused_ism.*

    [nTracks, nDetections] = size(logRatio);
    association = emptyAssociationResult(nTracks, nDetections);
    if nTracks == 0 || nDetections == 0
        return;
    end

    adjacency = isfinite(logRatio);
    components = bipartiteAssociationComponents(adjacency);
    for c = 1:numel(components)
        trackIndices = components(c).tracks;
        detectionIndices = components(c).detections;
        if isempty(trackIndices) || isempty(detectionIndices)
            continue;
        end

        componentResult = enumerateAssociationHypotheses( ...
            logRatio(trackIndices, detectionIndices), logMiss, ...
            maxHypotheses);
        association.mapDetectionForTrack(trackIndices) = ...
            mapLocalToGlobal(componentResult.mapDetectionForTrack, ...
            detectionIndices);
        association.mapAssociationProbability(trackIndices) = ...
            componentResult.mapAssociationProbability;
        association.entropyForTrack(trackIndices) = componentResult.entropy;
        association.hypothesisCountForTrack(trackIndices) = ...
            componentResult.nHypotheses;

        association.newDetectionProbability(detectionIndices) = ...
            componentResult.newDetectionProbability;
        association.entropyForDetection(detectionIndices) = ...
            componentResult.entropy;
        association.hypothesisCountForDetection(detectionIndices) = ...
            componentResult.nHypotheses;
    end
end
