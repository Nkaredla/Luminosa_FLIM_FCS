function [frames, photonEvents, photophysics] = simulateMIETPhotonMovie( ...
        trajectories, detector, opts)
    import membrane_tracking.curved_miet.*

    positionsXY = trajectories.positionsXYUm;
    lifetimes = trajectories.lifetimeNs;
    positionsZ = trajectories.positionsZUm;
    nMolecules = size(positionsXY, 1);
    nFrames = size(positionsXY, 3);

    meanSignalPhotons = zeros(nMolecules, nFrames);
    detectedRateS = zeros(nMolecules, nFrames);
    for frame = 1:nFrames
        position = positionsXY(:,:,frame);
        radialSquared = sum(position.^2, 2);
        relativeIllumination = exp(-2 * radialSquared / ...
            opts.laserWaistUm^2);
        excitationRate = opts.laserPowerW * opts.excitationRatePerW * ...
            relativeIllumination;
        excitationCycleRate = excitationRate ./ ...
            (1 + excitationRate / opts.saturationRateS);
        emissionRate = opts.quantumYield * excitationCycleRate;
        relativeCapture = ismArrayCaptureRelative(position, detector, opts);
        detectedRate = emissionRate * opts.collectionEfficiency * ...
            opts.detectorQuantumEfficiency .* relativeCapture;
        detectedRateS(:,frame) = detectedRate;
        meanSignalPhotons(:,frame) = detectedRate * opts.dtS;
    end

    signalPhotonCounts = samplePoisson(meanSignalPhotons);
    backgroundCounts = samplePoisson(detector.backgroundMeanCounts * ...
        ones(detector.nChannels, nFrames));
    frames = backgroundCounts;

    totalEvents = sum(signalPhotonCounts(:)) + sum(backgroundCounts(:));
    frameCol = zeros(totalEvents, 1);
    moleculeCol = zeros(totalEvents, 1);
    channelCol = zeros(totalEvents, 1);
    detectorXCol = zeros(totalEvents, 1);
    detectorYCol = zeros(totalEvents, 1);
    trueXCol = nan(totalEvents, 1);
    trueYCol = nan(totalEvents, 1);
    trueZCol = nan(totalEvents, 1);
    trueLifetimeCol = nan(totalEvents, 1);
    microtimeCol = zeros(totalEvents, 1);
    arrivalTimeCol = zeros(totalEvents, 1);
    isBackgroundCol = false(totalEvents, 1);
    eventStartRow = zeros(nFrames, 1);
    eventEndRow = zeros(nFrames, 1);

    cursor = 0;
    for frame = 1:nFrames
        eventStartRow(frame) = cursor + 1;
        for molecule = 1:nMolecules
            nPhotons = signalPhotonCounts(molecule,frame);
            if nPhotons <= 0
                continue;
            end

            position = positionsXY(molecule,:,frame);
            channelProbability = ismDetectorChannelProbability( ...
                position, detector, opts);
            channels = sampleCategorical(channelProbability, nPhotons);
            indices = cursor + (1:nPhotons);

            frameCol(indices) = frame;
            moleculeCol(indices) = molecule;
            channelCol(indices) = channels;
            detectorXCol(indices) = detector.detXY(channels,1);
            detectorYCol(indices) = detector.detXY(channels,2);
            trueXCol(indices) = position(1);
            trueYCol(indices) = position(2);
            trueZCol(indices) = positionsZ(molecule,frame);
            trueLifetimeCol(indices) = lifetimes(molecule,frame);
            microtimeCol(indices) = sampleTruncatedExponential( ...
                lifetimes(molecule,frame), opts.repetitionPeriodNs, nPhotons);
            arrivalTimeCol(indices) = ((frame-1) + rand(nPhotons,1)) * ...
                opts.dtS;

            frames(:,frame) = frames(:,frame) + ...
                accumarray(channels, 1, [detector.nChannels 1]);
            cursor = cursor + nPhotons;
        end

        backgroundChannels = repelem((1:detector.nChannels).', ...
            backgroundCounts(:,frame));
        nBackground = numel(backgroundChannels);
        if nBackground > 0
            indices = cursor + (1:nBackground);
            frameCol(indices) = frame;
            moleculeCol(indices) = 0;
            channelCol(indices) = backgroundChannels;
            detectorXCol(indices) = detector.detXY(backgroundChannels,1);
            detectorYCol(indices) = detector.detXY(backgroundChannels,2);
            microtimeCol(indices) = opts.repetitionPeriodNs * ...
                rand(nBackground,1);
            arrivalTimeCol(indices) = ((frame-1) + rand(nBackground,1)) * ...
                opts.dtS;
            isBackgroundCol(indices) = true;
            cursor = cursor + nBackground;
        end
        eventEndRow(frame) = cursor;
        if eventEndRow(frame) < eventStartRow(frame)
            eventStartRow(frame) = 0;
            eventEndRow(frame) = 0;
        end
    end

    photonEvents = table(frameCol(1:cursor), moleculeCol(1:cursor), ...
        channelCol(1:cursor), detectorXCol(1:cursor), ...
        detectorYCol(1:cursor), trueXCol(1:cursor), trueYCol(1:cursor), ...
        trueZCol(1:cursor), trueLifetimeCol(1:cursor), ...
        microtimeCol(1:cursor), arrivalTimeCol(1:cursor), ...
        isBackgroundCol(1:cursor), ...
        'VariableNames', {'frame', 'molecule', 'detectorChannel', ...
        'xDetectorUm', 'yDetectorUm', 'trueXUm', 'trueYUm', 'trueZUm', ...
        'trueLifetimeNs', 'microtimeNs', 'arrivalTimeS', 'isBackground'});

    brightMoleculeCount = sum(signalPhotonCounts >= ...
        opts.minPhotonsPerLocalization, 1);
    photophysics = struct();
    photophysics.meanSignalPhotons = meanSignalPhotons;
    photophysics.signalPhotonCounts = signalPhotonCounts;
    photophysics.backgroundCounts = backgroundCounts;
    photophysics.detectedRateS = detectedRateS;
    photophysics.eventStartRow = eventStartRow;
    photophysics.eventEndRow = eventEndRow;
    photophysics.brightMoleculeCountPerFrame = brightMoleculeCount;
    photophysics.fractionFramesWithMultipleBrightMolecules = ...
        mean(brightMoleculeCount > 1);
    photophysics.totalSignalPhotons = sum(signalPhotonCounts(:));
    photophysics.totalBackgroundPhotons = sum(backgroundCounts(:));
end
