function localizationTable = localizeMIETMovie(frames, photonEvents, ...
        photophysics, trajectories, detector, opts)
    import membrane_tracking.curved_miet.*

    nFrames = size(frames, 2);

    frameCol = zeros(nFrames,1);
    timeCol = zeros(nFrames,1);
    xCol = zeros(nFrames,1);
    yCol = zeros(nFrames,1);
    amplitudeCol = zeros(nFrames,1);
    backgroundCol = zeros(nFrames,1);
    nllCol = zeros(nFrames,1);
    devianceCol = zeros(nFrames,1);
    crbVarXCol = zeros(nFrames,1);
    crbVarYCol = zeros(nFrames,1);
    crbCovXYCol = zeros(nFrames,1);
    lifetimeCol = zeros(nFrames,1);
    lifetimeSigmaCol = zeros(nFrames,1);
    lifetimeFisherLogCol = zeros(nFrames,1);
    heightCol = zeros(nFrames,1);
    heightSigmaCol = zeros(nFrames,1);
    lifetimePhotonCountCol = zeros(nFrames,1);
    brightMoleculeCountCol = zeros(nFrames,1);
    nearestMoleculeCol = zeros(nFrames,1);
    truthErrorCol = zeros(nFrames,1);
    trueXCol = zeros(nFrames,1);
    trueYCol = zeros(nFrames,1);
    trueZCol = zeros(nFrames,1);
    trueLifetimeCol = zeros(nFrames,1);

    nAccepted = 0;
    for frame = 1:nFrames
        data = double(frames(:,frame));
        estimatedSignal = sum(data) - ...
            detector.nChannels * detector.backgroundMeanCounts;
        if estimatedSignal < opts.minPhotonsPerLocalization
            continue;
        end

        [model, fisher, ok] = fitSingleEmitterISM(data, detector, opts);
        if ~ok || model.reducedPoissonDeviance > ...
                opts.maxReducedPoissonDeviance
            continue;
        end

        firstRow = photophysics.eventStartRow(frame);
        lastRow = photophysics.eventEndRow(frame);
        if firstRow <= 0 || lastRow < firstRow
            continue;
        end
        rows = firstRow:lastRow;
        [lifetimeEstimate, lifetimeSigma, lifetimeFisherLog, ...
            lifetimeNll, lifetimeOk] = fitFrameLifetime( ...
            photonEvents.microtimeNs(rows), ...
            photonEvents.detectorChannel(rows), model, detector, opts);
        if ~lifetimeOk
            continue;
        end

        truthXY = trajectories.positionsXYUm(:,:,frame);
        delta = truthXY - [model.xUm model.yUm];
        [truthError, nearestMolecule] = min(hypot(delta(:,1), delta(:,2)));

        nAccepted = nAccepted + 1;
        frameCol(nAccepted) = frame;
        timeCol(nAccepted) = (frame-1) * opts.dtS;
        xCol(nAccepted) = model.xUm;
        yCol(nAccepted) = model.yUm;
        amplitudeCol(nAccepted) = model.amplitude;
        backgroundCol(nAccepted) = model.background;
        nllCol(nAccepted) = model.negLogLikelihood + lifetimeNll;
        devianceCol(nAccepted) = model.reducedPoissonDeviance;
        crbVarXCol(nAccepted) = fisher.covarianceXY(1,1);
        crbVarYCol(nAccepted) = fisher.covarianceXY(2,2);
        crbCovXYCol(nAccepted) = fisher.covarianceXY(1,2);
        lifetimeCol(nAccepted) = lifetimeEstimate;
        lifetimeSigmaCol(nAccepted) = lifetimeSigma;
        lifetimeFisherLogCol(nAccepted) = lifetimeFisherLog;
        heightCol(nAccepted) = (lifetimeEstimate - ...
            opts.lifetimeAtSubstrateNs) / opts.lifetimeSlopeNsPerUm;
        heightSigmaCol(nAccepted) = ...
            lifetimeSigma / opts.lifetimeSlopeNsPerUm;
        lifetimePhotonCountCol(nAccepted) = numel(rows);
        brightMoleculeCountCol(nAccepted) = ...
            photophysics.brightMoleculeCountPerFrame(frame);
        nearestMoleculeCol(nAccepted) = nearestMolecule;
        truthErrorCol(nAccepted) = truthError;
        trueXCol(nAccepted) = truthXY(nearestMolecule,1);
        trueYCol(nAccepted) = truthXY(nearestMolecule,2);
        trueZCol(nAccepted) = ...
            trajectories.positionsZUm(nearestMolecule,frame);
        trueLifetimeCol(nAccepted) = ...
            trajectories.lifetimeNs(nearestMolecule,frame);
    end

    keep = 1:nAccepted;
    localizationTable = table(frameCol(keep), timeCol(keep), xCol(keep), ...
        yCol(keep), amplitudeCol(keep), backgroundCol(keep), nllCol(keep), ...
        devianceCol(keep), crbVarXCol(keep), crbVarYCol(keep), ...
        crbCovXYCol(keep), sqrt(crbVarXCol(keep)), ...
        sqrt(crbVarYCol(keep)), lifetimeCol(keep), ...
        lifetimeSigmaCol(keep), lifetimeFisherLogCol(keep), ...
        heightCol(keep), heightSigmaCol(keep), ...
        lifetimePhotonCountCol(keep), brightMoleculeCountCol(keep), ...
        nearestMoleculeCol(keep), truthErrorCol(keep), trueXCol(keep), ...
        trueYCol(keep), trueZCol(keep), trueLifetimeCol(keep), ...
        'VariableNames', {'frame', 'timeS', 'xUm', 'yUm', ...
        'signalPhotons', 'backgroundPerChannel', 'negativeLogLikelihood', ...
        'reducedPoissonDeviance', 'crbVarXUm2', 'crbVarYUm2', ...
        'crbCovXYUm2', 'crbSigmaXUm', 'crbSigmaYUm', 'lifetimeNs', ...
        'lifetimeSigmaNs', 'lifetimeObservedFisherLogTau', 'heightUm', ...
        'heightSigmaUm', 'nPhotonsForLifetime', ...
        'trueBrightMoleculeCount', 'nearestTrueMolecule', ...
        'truthPositionErrorUm', 'trueXUm', 'trueYUm', 'trueZUm', ...
        'trueLifetimeNs'});
end
