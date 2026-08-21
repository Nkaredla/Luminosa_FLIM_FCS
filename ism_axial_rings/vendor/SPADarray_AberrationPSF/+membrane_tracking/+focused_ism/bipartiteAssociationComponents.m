function components = bipartiteAssociationComponents(adjacency)
    import membrane_tracking.focused_ism.*

    [nTracks, nDetections] = size(adjacency);
    unusedTracks = any(adjacency, 2);
    unusedDetections = any(adjacency, 1).';
    components = struct('tracks', {}, 'detections', {});

    while any(unusedTracks)
        firstTrack = find(unusedTracks, 1, 'first');
        trackQueue = firstTrack;
        detectionQueue = zeros(0, 1);
        componentTracks = zeros(0, 1);
        componentDetections = zeros(0, 1);
        unusedTracks(firstTrack) = false;

        while ~isempty(trackQueue) || ~isempty(detectionQueue)
            while ~isempty(trackQueue)
                track = trackQueue(1);
                trackQueue(1) = [];
                componentTracks(end+1, 1) = track; %#ok<AGROW>
                neighbors = find(adjacency(track,:)).';
                neighbors = neighbors(unusedDetections(neighbors));
                unusedDetections(neighbors) = false;
                detectionQueue = [detectionQueue; neighbors(:)]; %#ok<AGROW>
            end
            while ~isempty(detectionQueue)
                detection = detectionQueue(1);
                detectionQueue(1) = [];
                componentDetections(end+1, 1) = detection; %#ok<AGROW>
                neighbors = find(adjacency(:,detection));
                neighbors = neighbors(unusedTracks(neighbors));
                unusedTracks(neighbors) = false;
                trackQueue = [trackQueue; neighbors]; %#ok<AGROW>
            end
        end

        components(end+1).tracks = componentTracks; %#ok<AGROW>
        components(end).detections = componentDetections;
    end

    if nTracks == 0 || nDetections == 0
        components = struct('tracks', {}, 'detections', {});
    end
end
