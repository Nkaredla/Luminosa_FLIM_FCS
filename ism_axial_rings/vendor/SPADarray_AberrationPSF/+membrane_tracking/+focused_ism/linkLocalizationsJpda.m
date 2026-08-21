function trackTable = linkLocalizationsJpda(localizationTable, opts)
    import membrane_tracking.focused_ism.*

    trackTable = localizationTable;
    if height(trackTable) > 0
        [~, order] = sortrows([trackTable.frame, -trackTable.signalPhotons]);
        trackTable = trackTable(order, :);
    end
    nRows = height(trackTable);
    trackTable.trackId = zeros(nRows, 1);
    trackTable.associationProbability = nan(nRows, 1);
    trackTable.associationEntropy = nan(nRows, 1);
    trackTable.associationHypothesisCount = zeros(nRows, 1);
    trackTable.isNewTrack = false(nRows, 1);

    if nRows == 0
        return;
    end

    frames = unique(trackTable.frame).';
    nextTrackId = 1;
    tracks = struct('id', {}, 'lastFrame', {}, 'x', {}, 'y', {}, ...
        'covariance', {});
    % Clutter density converts each Brownian motion likelihood into a
    % dimensionless association likelihood ratio.
    fovAreaUm2 = pi * opts.maxLocalizationRadiusUm^2;
    clutterDensity = opts.associationClutterRatePerFrame / fovAreaUm2;
    logMiss = log(1 - opts.associationDetectionProbability);

    for frame = frames
        locIdx = find(trackTable.frame == frame);
        nLoc = numel(locIdx);
        assignedLoc = false(nLoc, 1);
        maxFrameSeparation = opts.maxGapFrames + 1;
        activeIdx = find([tracks.lastFrame] >= frame - maxFrameSeparation & ...
            [tracks.lastFrame] < frame);

        if ~isempty(activeIdx) && nLoc > 0
            logRatio = -inf(numel(activeIdx), nLoc);
            for a = 1:numel(activeIdx)
                tr = tracks(activeIdx(a));
                frameSeparation = frame - tr.lastFrame;
                motionVariance = 2 * opts.trackingPriorDiffusionUm2PerS * ...
                    frameSeparation * opts.dtS;
                euclideanGate = opts.trackingMaxStepUm * sqrt(frameSeparation);

                for l = 1:nLoc
                    row = locIdx(l);
                    delta = [trackTable.xUm(row) - tr.x, ...
                        trackTable.yUm(row) - tr.y];
                    delta = applyMinimumImage(delta, opts);
                    if hypot(delta(1), delta(2)) > euclideanGate
                        continue;
                    end

                    detectionCov = localizationCovarianceFromTable( ...
                        trackTable, row);
                    innovationCov = tr.covariance + detectionCov + ...
                        motionVariance * eye(2);
                    innovationCov = 0.5 * ...
                        (innovationCov + innovationCov.');
                    [inverseCov, covarianceValid] = ...
                        invertPositiveDefinite(innovationCov);
                    if ~covarianceValid
                        continue;
                    end
                    mahalanobisSquared = delta * inverseCov * delta.';
                    if mahalanobisSquared > opts.associationGateChi2
                        continue;
                    end

                    logDeterminant = log(max(det(innovationCov), realmin));
                    logMotionDensity = -log(2*pi) - ...
                        0.5 * logDeterminant - 0.5 * mahalanobisSquared;
                    % Compare an existing-track detection against the
                    % alternative that the localization is clutter/new.
                    logRatio(a,l) = log( ...
                        opts.associationDetectionProbability) + ...
                        logMotionDensity - log(clutterDensity);
                end
            end

            association = solveJpdaAssociations(logRatio, logMiss, ...
                opts.maxAssociationHypotheses);
            for a = 1:numel(activeIdx)
                l = association.mapDetectionForTrack(a);
                if l <= 0
                    continue;
                end
                row = locIdx(l);
                trIdx = activeIdx(a);
                trackTable.trackId(row) = tracks(trIdx).id;
                trackTable.associationProbability(row) = ...
                    association.mapAssociationProbability(a);
                trackTable.associationEntropy(row) = ...
                    association.entropyForTrack(a);
                trackTable.associationHypothesisCount(row) = ...
                    association.hypothesisCountForTrack(a);
                tracks(trIdx).lastFrame = frame;
                tracks(trIdx).x = trackTable.xUm(row);
                tracks(trIdx).y = trackTable.yUm(row);
                tracks(trIdx).covariance = ...
                    localizationCovarianceFromTable(trackTable, row);
                assignedLoc(l) = true;
            end
        else
            association = emptyAssociationResult(numel(activeIdx), nLoc);
        end

        for l = 1:nLoc
            if assignedLoc(l)
                continue;
            end
            row = locIdx(l);
            trackTable.trackId(row) = nextTrackId;
            if l <= numel(association.newDetectionProbability)
                newProbability = association.newDetectionProbability(l);
                entropyValue = association.entropyForDetection(l);
                hypothesisCount = association.hypothesisCountForDetection(l);
            else
                newProbability = 1;
                entropyValue = 0;
                hypothesisCount = 1;
            end
            trackTable.associationProbability(row) = newProbability;
            trackTable.associationEntropy(row) = entropyValue;
            trackTable.associationHypothesisCount(row) = hypothesisCount;
            trackTable.isNewTrack(row) = true;

            tracks(end+1).id = nextTrackId; %#ok<AGROW>
            tracks(end).lastFrame = frame;
            tracks(end).x = trackTable.xUm(row);
            tracks(end).y = trackTable.yUm(row);
            tracks(end).covariance = ...
                localizationCovarianceFromTable(trackTable, row);
            nextTrackId = nextTrackId + 1;
        end
    end

    trackTable = sortrows(trackTable, {'trackId', 'frame'});
end
