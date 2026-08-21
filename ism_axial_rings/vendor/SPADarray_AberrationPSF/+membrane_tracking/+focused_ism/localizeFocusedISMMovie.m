function localizationTable = localizeFocusedISMMovie( ...
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
    identVec = false(0, 1);

    nFrames = size(frames, 2);
    for t = 1:nFrames
        data = double(frames(:,t));
        % Reject frames that cannot support even one fitted monomer before
        % running the more expensive model-order search.
        signalEstimate = sum(data) - ...
            detector.nChannels * detector.backgroundMeanCounts;
        if signalEstimate < opts.minPhotonsPerLocalization
            continue;
        end

        [locs, modelOrder] = fitFocusedISMFrame(data, detector, opts);
        if isempty(locs)
            continue;
        end
        if numel(locs) > opts.maxLocalizationsPerFrame
            [~, signalOrder] = sort([locs.signalPhotons], 'descend');
            locs = locs(signalOrder(1:opts.maxLocalizationsPerFrame));
        end

        for e = 1:numel(locs)
            loc = locs(e);
            [nearestMol, truthError] = nearestTruthPosition( ...
                positionsUm(:,:,t), [loc.xUm loc.yUm], ...
                photonCountsInFov(:,t), opts);

            frameVec(end+1, 1) = t; %#ok<AGROW>
            spotVec(end+1, 1) = e; %#ok<AGROW>
            xVec(end+1, 1) = loc.xUm; %#ok<AGROW>
            yVec(end+1, 1) = loc.yUm; %#ok<AGROW>
            xPixVec(end+1, 1) = -loc.xUm; %#ok<AGROW>
            yPixVec(end+1, 1) = -loc.yUm; %#ok<AGROW>
            signalVec(end+1, 1) = loc.signalPhotons; %#ok<AGROW>
            windowCountsVec(end+1, 1) = sum(data); %#ok<AGROW>
            bgVec(end+1, 1) = loc.backgroundCountsPerDetector; %#ok<AGROW>
            crbXVec(end+1, 1) = loc.crbVarXUm2; %#ok<AGROW>
            crbYVec(end+1, 1) = loc.crbVarYUm2; %#ok<AGROW>
            crbXYVec(end+1, 1) = loc.crbCovXYUm2; %#ok<AGROW>
            crbSigmaVec(end+1, 1) = loc.crbSigmaUm; %#ok<AGROW>
            nllVec(end+1, 1) = loc.negLogLikelihood; %#ok<AGROW>
            trueMolVec(end+1, 1) = nearestMol; %#ok<AGROW>
            truthErrVec(end+1, 1) = truthError; %#ok<AGROW>
            regionVec(end+1, 1) = t; %#ok<AGROW>
            jointCountVec(end+1, 1) = modelOrder.nEmitters; %#ok<AGROW>
            orderProbVec(end+1, 1) = ...
                modelOrder.posteriorProbability; %#ok<AGROW>
            orderEntropyVec(end+1, 1) = modelOrder.entropy; %#ok<AGROW>
            identVec(end+1, 1) = loc.positionIdentifiable; %#ok<AGROW>
        end
    end

    localizationTable = table(frameVec, spotVec, xVec, yVec, xPixVec, ...
        yPixVec, signalVec, windowCountsVec, bgVec, crbXVec, crbYVec, ...
        crbXYVec, crbSigmaVec, nllVec, trueMolVec, truthErrVec, ...
        regionVec, jointCountVec, orderProbVec, orderEntropyVec, identVec, ...
        'VariableNames', {'frame', 'spotIndexInFrame', 'xUm', 'yUm', ...
        'xDetectorImageUm', 'yDetectorImageUm', 'signalPhotons', ...
        'windowCounts', 'backgroundCountsPerDetector', 'crbVarXUm2', ...
        'crbVarYUm2', 'crbCovXYUm2', 'crbSigmaUm', ...
        'negLogLikelihood', 'nearestTrueMolecule', 'truthErrorUm', ...
        'fitRegionId', 'jointEmitterCount', ...
        'emitterCountPosteriorProbability', 'modelOrderEntropy', ...
        'positionIdentifiable'});
end
