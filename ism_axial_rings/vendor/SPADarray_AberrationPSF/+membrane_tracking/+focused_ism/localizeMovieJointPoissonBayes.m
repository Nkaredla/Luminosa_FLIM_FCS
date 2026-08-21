function localizationTable = localizeMovieJointPoissonBayes( ...
        frames, detector, opts, positionsUm, photonCountsInFov)
    import membrane_tracking.focused_ism.*

    frameVec = zeros(0, 1);
    spotVec = zeros(0, 1);
    xVec = zeros(0, 1);
    yVec = zeros(0, 1);
    xPixVec = zeros(0, 1);
    yPixVec = zeros(0, 1);
    signalVec = zeros(0, 1);
    windowCountsVec = zeros(0, 1);
    bgVec = zeros(0, 1);
    crbXVec = zeros(0, 1);
    crbYVec = zeros(0, 1);
    crbXYVec = zeros(0, 1);
    crbSigmaVec = zeros(0, 1);
    nllVec = zeros(0, 1);
    trueMolVec = zeros(0, 1);
    truthErrVec = zeros(0, 1);
    regionVec = zeros(0, 1);
    jointCountVec = zeros(0, 1);
    orderProbVec = zeros(0, 1);
    orderEntropyVec = zeros(0, 1);

    nFrames = size(frames, 3);
    nextRegionId = 1;
    for t = 1:nFrames
        image = double(frames(:,:,t));
        candidates = findSpotCandidates(image, detector, opts);
        if isempty(candidates)
            continue;
        end

        groups = groupCandidateRegions(candidates(:,1:2), ...
            opts.multiEmitterGroupingRadiusPx);
        spotInFrame = 0;
        for g = 1:numel(groups)
            if spotInFrame >= opts.maxLocalizationsPerFrame
                break;
            end
            candidateRows = candidates(groups{g}, :);
            [locs, modelOrder] = fitJointEmitterRegion( ...
                image, candidateRows, detector, opts);
            if isempty(locs)
                nextRegionId = nextRegionId + 1;
                continue;
            end

            remainingSlots = opts.maxLocalizationsPerFrame - spotInFrame;
            if isfinite(remainingSlots) && numel(locs) > remainingSlots
                [~, signalOrder] = sort([locs.signalPhotons], 'descend');
                locs = locs(signalOrder(1:remainingSlots));
            end

            for e = 1:numel(locs)
                loc = locs(e);
                spotInFrame = spotInFrame + 1;
                [nearestMol, truthError] = nearestTruthPosition( ...
                    positionsUm(:,:,t), [loc.xUm loc.yUm], ...
                    photonCountsInFov(:,t), opts);

                frameVec(end+1, 1) = t; %#ok<AGROW>
                spotVec(end+1, 1) = spotInFrame; %#ok<AGROW>
                xVec(end+1, 1) = loc.xUm; %#ok<AGROW>
                yVec(end+1, 1) = loc.yUm; %#ok<AGROW>
                xPixVec(end+1, 1) = loc.xPixel; %#ok<AGROW>
                yPixVec(end+1, 1) = loc.yPixel; %#ok<AGROW>
                signalVec(end+1, 1) = loc.signalPhotons; %#ok<AGROW>
                windowCountsVec(end+1, 1) = loc.windowCounts; %#ok<AGROW>
                bgVec(end+1, 1) = loc.backgroundCountsPerPixel; %#ok<AGROW>
                crbXVec(end+1, 1) = loc.crbVarXUm2; %#ok<AGROW>
                crbYVec(end+1, 1) = loc.crbVarYUm2; %#ok<AGROW>
                crbXYVec(end+1, 1) = loc.crbCovXYUm2; %#ok<AGROW>
                crbSigmaVec(end+1, 1) = loc.crbSigmaUm; %#ok<AGROW>
                nllVec(end+1, 1) = loc.negLogLikelihood; %#ok<AGROW>
                trueMolVec(end+1, 1) = nearestMol; %#ok<AGROW>
                truthErrVec(end+1, 1) = truthError; %#ok<AGROW>
                regionVec(end+1, 1) = nextRegionId; %#ok<AGROW>
                jointCountVec(end+1, 1) = modelOrder.nEmitters; %#ok<AGROW>
                orderProbVec(end+1, 1) = modelOrder.posteriorProbability; %#ok<AGROW>
                orderEntropyVec(end+1, 1) = modelOrder.entropy; %#ok<AGROW>
            end
            nextRegionId = nextRegionId + 1;
        end
    end

    localizationTable = table(frameVec, spotVec, xVec, yVec, xPixVec, ...
        yPixVec, signalVec, windowCountsVec, bgVec, crbXVec, crbYVec, ...
        crbXYVec, crbSigmaVec, nllVec, trueMolVec, truthErrVec, ...
        regionVec, jointCountVec, orderProbVec, orderEntropyVec, ...
        'VariableNames', {'frame', 'spotIndexInFrame', 'xUm', 'yUm', ...
        'xDetectorPixel', 'yDetectorPixel', 'signalPhotons', ...
        'windowCounts', 'backgroundCountsPerPixel', 'crbVarXUm2', ...
        'crbVarYUm2', 'crbCovXYUm2', 'crbSigmaUm', ...
        'negLogLikelihood', 'nearestTrueMolecule', 'truthErrorUm', ...
        'fitRegionId', 'jointEmitterCount', ...
        'emitterCountPosteriorProbability', 'modelOrderEntropy'});
end
