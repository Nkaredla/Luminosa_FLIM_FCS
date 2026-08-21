function trackTable = linkLocalizationsNearestNeighbor(localizationTable, opts)
    import membrane_tracking.focused_ism.*

    trackTable = localizationTable;
    if height(trackTable) > 0
        [~, order] = sortrows([trackTable.frame, -trackTable.signalPhotons]);
        trackTable = trackTable(order, :);
    end
    trackTable.trackId = zeros(height(trackTable), 1);
    trackTable.associationProbability = ones(height(trackTable), 1);
    trackTable.associationEntropy = zeros(height(trackTable), 1);
    trackTable.associationHypothesisCount = ones(height(trackTable), 1);
    trackTable.isNewTrack = false(height(trackTable), 1);

    if height(trackTable) == 0
        return;
    end

    frames = unique(trackTable.frame).';
    nextTrackId = 1;
    tracks = struct('id', {}, 'lastFrame', {}, 'x', {}, 'y', {});

    for frame = frames
        locIdx = find(trackTable.frame == frame);
        nLoc = numel(locIdx);
        assignedLoc = false(nLoc, 1);

        maxFrameSeparation = opts.maxGapFrames + 1;
        activeIdx = find([tracks.lastFrame] >= frame - maxFrameSeparation & ...
            [tracks.lastFrame] < frame);

        if ~isempty(activeIdx) && nLoc > 0
            dist = inf(numel(activeIdx), nLoc);
            for a = 1:numel(activeIdx)
                tr = tracks(activeIdx(a));
                delta = [trackTable.xUm(locIdx) - tr.x, ...
                    trackTable.yUm(locIdx) - tr.y];
                delta = applyMinimumImage(delta, opts);
                dist(a, :) = hypot(delta(:,1), delta(:,2)).';
                frameSeparation = frame - tr.lastFrame;
                gate = opts.trackingMaxStepUm * sqrt(frameSeparation);
                dist(a, dist(a,:) > gate) = Inf;
            end

            usedTrack = false(numel(activeIdx), 1);
            while true
                dist(usedTrack, :) = inf;
                dist(:, assignedLoc) = inf;
                [bestDist, flatIdx] = min(dist(:));
                if isempty(bestDist) || ~isfinite(bestDist)
                    break;
                end

                [a, l] = ind2sub(size(dist), flatIdx);
                trIdx = activeIdx(a);
                row = locIdx(l);
                trackTable.trackId(row) = tracks(trIdx).id;
                tracks(trIdx).lastFrame = frame;
                tracks(trIdx).x = trackTable.xUm(row);
                tracks(trIdx).y = trackTable.yUm(row);
                usedTrack(a) = true;
                assignedLoc(l) = true;
            end
        end

        for k = 1:nLoc
            if assignedLoc(k)
                continue;
            end
            row = locIdx(k);
            trackTable.trackId(row) = nextTrackId;
            trackTable.isNewTrack(row) = true;
            tracks(end+1).id = nextTrackId; %#ok<AGROW>
            tracks(end).lastFrame = frame;
            tracks(end).x = trackTable.xUm(row);
            tracks(end).y = trackTable.yUm(row);
            nextTrackId = nextTrackId + 1;
        end
    end

    trackTable = sortrows(trackTable, {'trackId', 'frame'});
end
