function localizationTable = localizeMovieSingleEmitter( ...
        frames, detector, opts, positionsUm, ...
        photonCountsInFov)
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

    nFrames = size(frames, 3);
    for t = 1:nFrames
        image = double(frames(:,:,t));
        candidates = findSpotCandidates(image, detector, opts);
        if isempty(candidates)
            continue;
        end

        spotInFrame = 0;
        acceptedXY = zeros(0, 2);
        for k = 1:size(candidates, 1)
            if spotInFrame >= opts.maxLocalizationsPerFrame
                break;
            end

            row0 = candidates(k, 1);
            col0 = candidates(k, 2);
            [loc, ok] = fitGaussianSpotMLE(image, row0, col0, detector, opts);
            if ~ok
                continue;
            end

            if ~isempty(acceptedXY)
                d = hypot(acceptedXY(:,1) - loc.xUm, ...
                    acceptedXY(:,2) - loc.yUm);
                if any(d < opts.minPeakSeparationPx * detector.pixelSizeUm)
                    continue;
                end
            end

            spotInFrame = spotInFrame + 1;
            acceptedXY(end+1, :) = [loc.xUm loc.yUm]; %#ok<AGROW>

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
        end
    end

    localizationTable = table(frameVec, spotVec, xVec, yVec, xPixVec, ...
        yPixVec, signalVec, windowCountsVec, bgVec, crbXVec, crbYVec, ...
        crbXYVec, crbSigmaVec, nllVec, trueMolVec, truthErrVec, ...
        'VariableNames', {'frame', 'spotIndexInFrame', 'xUm', 'yUm', ...
        'xDetectorPixel', 'yDetectorPixel', 'signalPhotons', ...
        'windowCounts', 'backgroundCountsPerPixel', 'crbVarXUm2', ...
        'crbVarYUm2', 'crbCovXYUm2', 'crbSigmaUm', ...
        'negLogLikelihood', 'nearestTrueMolecule', 'truthErrorUm'});
    localizationTable.fitRegionId = (1:height(localizationTable)).';
    localizationTable.jointEmitterCount = ones(height(localizationTable), 1);
    localizationTable.emitterCountPosteriorProbability = ...
        ones(height(localizationTable), 1);
    localizationTable.modelOrderEntropy = zeros(height(localizationTable), 1);
end
