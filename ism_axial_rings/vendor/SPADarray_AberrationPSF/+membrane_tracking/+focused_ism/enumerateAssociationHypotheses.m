function result = enumerateAssociationHypotheses(logRatio, logMiss, ...
        maxHypotheses)
    import membrane_tracking.focused_ism.*

    [nTracks, nDetections] = size(logRatio);
    assignments = zeros(1, 0);
    logWeights = 0;

    for track = 1:nTracks
        newAssignments = zeros(0, track);
        newLogWeights = zeros(0, 1);
        for h = 1:size(assignments, 1)
            usedDetections = assignments(h, assignments(h,:) > 0);
            options = find(isfinite(logRatio(track,:)));
            options = options(~ismember(options, usedDetections));
            options = options(:).';

            newAssignments(end+1, 1:track) = ...
                [assignments(h,:), 0]; %#ok<AGROW>
            newLogWeights(end+1, 1) = ...
                logWeights(h) + logMiss; %#ok<AGROW>
            for optionIndex = 1:numel(options)
                option = options(optionIndex);
                newAssignments(end+1, 1:track) = ...
                    [assignments(h,:), option]; %#ok<AGROW>
                newLogWeights(end+1, 1) = ...
                    logWeights(h) + logRatio(track,option); %#ok<AGROW>
            end
        end

        if numel(newLogWeights) > maxHypotheses
            [~, keep] = maxk(newLogWeights, maxHypotheses);
            newAssignments = newAssignments(keep,:);
            newLogWeights = newLogWeights(keep);
        end
        assignments = newAssignments;
        logWeights = newLogWeights;
    end

    shiftedLogWeights = logWeights - max(logWeights);
    weights = exp(shiftedLogWeights);
    weights = weights / sum(weights);
    [~, mapIndex] = max(weights);
    mapAssignment = assignments(mapIndex,:).';

    marginal = zeros(nTracks, nDetections);
    for h = 1:numel(weights)
        for track = 1:nTracks
            detection = assignments(h,track);
            if detection > 0
                marginal(track,detection) = ...
                    marginal(track,detection) + weights(h);
            end
        end
    end

    mapProbability = zeros(nTracks, 1);
    for track = 1:nTracks
        if mapAssignment(track) > 0
            mapProbability(track) = ...
                marginal(track,mapAssignment(track));
        else
            mapProbability(track) = 1 - sum(marginal(track,:));
        end
    end
    newProbability = max(1 - sum(marginal, 1).', 0);
    entropyValue = -sum(weights .* log(max(weights, realmin)));

    result = struct();
    result.mapDetectionForTrack = mapAssignment;
    result.mapAssociationProbability = mapProbability;
    result.newDetectionProbability = newProbability;
    result.entropy = entropyValue;
    result.nHypotheses = numel(weights);
end
