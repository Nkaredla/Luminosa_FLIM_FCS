function trackTable = linkLowDensityLocalizations(localizationTable, opts)
    import membrane_tracking.curved_miet.*

    trackTable = sortrows(localizationTable, 'frame');
    nLocalizations = height(trackTable);
    trackTable.trackId = zeros(nLocalizations,1);
    trackTable.isNewTrack = false(nLocalizations,1);
    trackTable.associationProbability = nan(nLocalizations,1);
    trackTable.frameGapFromPrevious = nan(nLocalizations,1);
    if nLocalizations == 0
        return;
    end

    activeId = zeros(0,1);
    activeLastRow = zeros(0,1);
    nextTrackId = 1;
    uniqueFrames = unique(trackTable.frame).';
    for frame = uniqueFrames
        detectionRows = find(trackTable.frame == frame);
        if ~isempty(activeId)
            frameGap = frame - trackTable.frame(activeLastRow);
            retain = frameGap <= opts.maxGapFrames + 1;
            activeId = activeId(retain);
            activeLastRow = activeLastRow(retain);
        end

        nActive = numel(activeId);
        nDetections = numel(detectionRows);
        cost = inf(nActive, nDetections);
        probability = zeros(nActive, nDetections);
        gapMatrix = nan(nActive, nDetections);
        for active = 1:nActive
            previousRow = activeLastRow(active);
            previousPosition = [trackTable.xUm(previousRow), ...
                trackTable.yUm(previousRow)];
            previousCovariance = localizationCovariance( ...
                trackTable, previousRow);
            for detection = 1:nDetections
                row = detectionRows(detection);
                gap = frame - trackTable.frame(previousRow);
                deltaTime = gap * opts.dtS;
                currentCovariance = localizationCovariance(trackTable, row);
                transitionCovariance = previousCovariance + ...
                    currentCovariance + ...
                    2 * opts.trackingPriorDiffusionUm2PerS * ...
                    deltaTime * eye(2);
                [inverseCovariance, covarianceValid] = ...
                    invertPositiveDefinite(transitionCovariance);
                if ~covarianceValid
                    continue;
                end
                delta = [trackTable.xUm(row), trackTable.yUm(row)] - ...
                    previousPosition;
                mahalanobisSquared = delta * inverseCovariance * delta.';
                if mahalanobisSquared <= opts.trackingGateChi2
                    cost(active,detection) = mahalanobisSquared;
                    probability(active,detection) = ...
                        exp(-0.5 * mahalanobisSquared);
                    gapMatrix(active,detection) = gap;
                end
            end
        end

        assignedActive = false(nActive,1);
        assignedDetection = false(nDetections,1);
        while any(isfinite(cost(:)))
            [~, linearIndex] = min(cost(:));
            [active, detection] = ind2sub(size(cost), linearIndex);
            row = detectionRows(detection);
            trackTable.trackId(row) = activeId(active);
            trackTable.associationProbability(row) = ...
                probability(active,detection);
            trackTable.frameGapFromPrevious(row) = ...
                gapMatrix(active,detection);
            activeLastRow(active) = row;
            assignedActive(active) = true;
            assignedDetection(detection) = true;
            cost(active,:) = Inf;
            cost(:,detection) = Inf;
        end

        for detection = find(~assignedDetection).'
            row = detectionRows(detection);
            trackTable.trackId(row) = nextTrackId;
            trackTable.isNewTrack(row) = true;
            trackTable.associationProbability(row) = 1;
            trackTable.frameGapFromPrevious(row) = NaN;
            activeId(end+1,1) = nextTrackId; %#ok<AGROW>
            activeLastRow(end+1,1) = row; %#ok<AGROW>
            nextTrackId = nextTrackId + 1;
        end

        % Assigned tracks already had their last rows updated. Unassigned
        % active tracks remain available until maxGapFrames is exceeded.
        if any(assignedActive)
            activeId = activeId(:);
            activeLastRow = activeLastRow(:);
        end
    end
end
