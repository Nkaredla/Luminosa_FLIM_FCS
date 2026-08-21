function [frames, photonEvents, photophysics] = simulatePhotonMovie( ...
        trajectories, detector, opts)
    import membrane_tracking.fluctuating_miet.*

    positionsXY = trajectories.positionsXYUm;
    positionsZ = trajectories.positionsZUm;
    lifetimes = trajectories.lifetimeNs;
    nMolecules = size(positionsXY, 1);
    nFrames = size(positionsXY, 3);

    meanSignal = zeros(nMolecules, nFrames);
    for frame = 1:nFrames
        position = positionsXY(:,:,frame);
        relativeIllumination = exp(-2 * sum(position.^2, 2) / ...
            opts.laserWaistUm^2);
        excitationRate = opts.laserPowerW * opts.excitationRatePerW * ...
            relativeIllumination;
        cycleRate = excitationRate ./ (1 + excitationRate / ...
            opts.saturationRateS);
        capture = zeros(nMolecules,1);
        for molecule = 1:nMolecules
            capture(molecule) = sum(ismDetectorRawResponse( ...
                position(molecule,:), detector, opts)) / ...
                detector.centerCaptureNormalization;
        end
        meanSignal(:,frame) = opts.quantumYield * cycleRate * ...
            opts.collectionEfficiency * opts.detectorQuantumEfficiency .* ...
            max(capture,0) * opts.dtS;
    end

    signalCounts = samplePoisson(meanSignal);
    backgroundCounts = samplePoisson(detector.backgroundMeanCounts * ...
        ones(detector.nChannels, nFrames));
    frames = backgroundCounts;

    total = sum(signalCounts(:)) + sum(backgroundCounts(:));
    frameCol = zeros(total,1);      moleculeCol = zeros(total,1);
    channelCol = zeros(total,1);    microtimeCol = zeros(total,1);
    trueZCol = nan(total,1);        trueTauCol = nan(total,1);
    isBackgroundCol = false(total,1);
    eventStartRow = zeros(nFrames,1);
    eventEndRow = zeros(nFrames,1);
    cursor = 0;

    for frame = 1:nFrames
        eventStartRow(frame) = cursor + 1;
        for molecule = 1:nMolecules
            nPhotons = signalCounts(molecule,frame);
            if nPhotons <= 0
                continue;
            end
            probability = ismDetectorChannelProbability( ...
                positionsXY(molecule,:,frame), detector, opts);
            channels = sampleCategorical(probability, nPhotons);
            idx = cursor + (1:nPhotons);
            frameCol(idx) = frame;
            moleculeCol(idx) = molecule;
            channelCol(idx) = channels;
            trueZCol(idx) = positionsZ(molecule,frame);
            trueTauCol(idx) = lifetimes(molecule,frame);
            microtimeCol(idx) = sampleTruncatedExponential( ...
                lifetimes(molecule,frame), opts.repetitionPeriodNs, nPhotons);
            frames(:,frame) = frames(:,frame) + ...
                accumarray(channels, 1, [detector.nChannels 1]);
            cursor = cursor + nPhotons;
        end
        backgroundChannels = repelem((1:detector.nChannels).', ...
            backgroundCounts(:,frame));
        nBackground = numel(backgroundChannels);
        if nBackground > 0
            idx = cursor + (1:nBackground);
            frameCol(idx) = frame;
            channelCol(idx) = backgroundChannels;
            microtimeCol(idx) = opts.repetitionPeriodNs * ...
                rand(nBackground,1);
            isBackgroundCol(idx) = true;
            cursor = cursor + nBackground;
        end
        eventEndRow(frame) = cursor;
        if eventEndRow(frame) < eventStartRow(frame)
            eventStartRow(frame) = 0;
            eventEndRow(frame) = 0;
        end
    end

    keep = 1:cursor;
    photonEvents = table(frameCol(keep), moleculeCol(keep), ...
        channelCol(keep), microtimeCol(keep), trueZCol(keep), ...
        trueTauCol(keep), isBackgroundCol(keep), 'VariableNames', ...
        {'frame','molecule','detectorChannel','microtimeNs', ...
         'trueZUm','trueLifetimeNs','isBackground'});

    brightCount = sum(signalCounts >= opts.minPhotonsPerLocalization, 1);
    photophysics = struct();
    photophysics.meanSignalPhotons = meanSignal;
    photophysics.signalPhotonCounts = signalCounts;
    photophysics.backgroundCounts = backgroundCounts;
    photophysics.eventStartRow = eventStartRow;
    photophysics.eventEndRow = eventEndRow;
    photophysics.brightMoleculeCountPerFrame = brightCount;
    photophysics.fractionFramesWithMultipleBrightMolecules = ...
        mean(brightCount > 1);
    photophysics.totalSignalPhotons = sum(signalCounts(:));
end
