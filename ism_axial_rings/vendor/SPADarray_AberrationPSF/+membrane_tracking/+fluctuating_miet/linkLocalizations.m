function trackTable = linkLocalizations(localizationTable, opts)
    import membrane_tracking.fluctuating_miet.*

    trackTable = sortrows(localizationTable, 'frame');
    n = height(trackTable);
    trackTable.trackId = zeros(n,1);
    trackTable.isNewTrack = false(n,1);
    trackTable.mahalanobisSquared = nan(n,1);
    trackTable.gateRadiusUm = nan(n,1);
    trackTable.frameGapFromPrevious = nan(n,1);
    if n == 0
        return;
    end
    activeId = zeros(0,1);
    activeRow = zeros(0,1);
    nextTrackId = 1;
    for frame = unique(trackTable.frame).'
        detectionRows = find(trackTable.frame == frame);
        if ~isempty(activeId)
            keep = frame - trackTable.frame(activeRow) <= opts.maxGapFrames+1;
            activeId = activeId(keep);
            activeRow = activeRow(keep);
        end
        nActive = numel(activeId);
        nDetections = numel(detectionRows);
        cost = inf(nActive, nDetections);
        gateRadius = nan(nActive, nDetections);
        gapMatrix = nan(nActive, nDetections);
        for a = 1:nActive
            previousRow = activeRow(a);
            previousPosition = [trackTable.xUm(previousRow), ...
                trackTable.yUm(previousRow)];
            previousCovariance = localizationCovariance(trackTable, previousRow);
            for d = 1:nDetections
                row = detectionRows(d);
                gap = frame - trackTable.frame(previousRow);
                innovation = previousCovariance + ...
                    localizationCovariance(trackTable, row) + ...
                    2*opts.trackingPriorDiffusionUm2PerS*gap*opts.dtS*eye(2);
                innovation = 0.5*(innovation + innovation.');
                [inverseCovariance, ok] = invertPositiveDefinite(innovation);
                if ~ok
                    continue;
                end
                delta = [trackTable.xUm(row), trackTable.yUm(row)] - ...
                    previousPosition;
                m2 = delta * inverseCovariance * delta.';
                if m2 <= opts.trackingGateChi2
                    cost(a,d) = m2;
                    gapMatrix(a,d) = gap;
                    % Isotropic-equivalent acceptance radius for this step.
                    gateRadius(a,d) = sqrt(opts.trackingGateChi2 * ...
                        0.5*trace(innovation));
                end
            end
        end
        assignedDetection = false(nDetections,1);
        while any(isfinite(cost(:)))
            [~, linearIndex] = min(cost(:));
            [a, d] = ind2sub(size(cost), linearIndex);
            row = detectionRows(d);
            trackTable.trackId(row) = activeId(a);
            trackTable.mahalanobisSquared(row) = cost(a,d);
            trackTable.gateRadiusUm(row) = gateRadius(a,d);
            trackTable.frameGapFromPrevious(row) = gapMatrix(a,d);
            activeRow(a) = row;
            assignedDetection(d) = true;
            cost(a,:) = Inf;
            cost(:,d) = Inf;
        end
        for d = find(~assignedDetection).'
            row = detectionRows(d);
            trackTable.trackId(row) = nextTrackId;
            trackTable.isNewTrack(row) = true;
            activeId(end+1,1) = nextTrackId; %#ok<AGROW>
            activeRow(end+1,1) = row; %#ok<AGROW>
            nextTrackId = nextTrackId + 1;
        end
    end
    trackTable = sortrows(trackTable, {'trackId','frame'});
end
