function [frames, photonEvents, photonCountsInFov] = ...
        projectMoleculePhotons(frames, positionsUm, photonCounts, detector, opts)
    import membrane_tracking.focused_ism.*

    totalPhotons = sum(photonCounts(:));

    if opts.storePhotonEvents
        frameCol = zeros(totalPhotons, 1);
        moleculeCol = zeros(totalPhotons, 1);
        channelCol = zeros(totalPhotons, 1);
        detectorXCol = zeros(totalPhotons, 1);
        detectorYCol = zeros(totalPhotons, 1);
        trueXCol = zeros(totalPhotons, 1);
        trueYCol = zeros(totalPhotons, 1);
        arrivalTimeCol = zeros(totalPhotons, 1);
    else
        frameCol = zeros(0, 1);
        moleculeCol = zeros(0, 1);
        channelCol = zeros(0, 1);
        detectorXCol = zeros(0, 1);
        detectorYCol = zeros(0, 1);
        trueXCol = zeros(0, 1);
        trueYCol = zeros(0, 1);
        arrivalTimeCol = zeros(0, 1);
    end

    cursor = 0;
    nMol = size(positionsUm, 1);
    nFrames = size(positionsUm, 3);
    photonCountsInFov = zeros(nMol, nFrames);

    for t = 1:nFrames
        for m = 1:nMol
            nPhot = photonCounts(m,t);
            if nPhot <= 0
                continue;
            end

            pos = positionsUm(m,:,t);
            channelProbability = ismDetectorChannelProbability( ...
                pos, detector, opts);
            channels = sampleCategorical(channelProbability, nPhot);
            % Poisson total counts followed by categorical channel
            % allocation is equivalent to independent Poisson channel
            % counts, while retaining individual photon identities.
            frames(:,t) = frames(:,t) + accumarray(channels, 1, ...
                [detector.nChannels 1]);
            photonCountsInFov(m,t) = nPhot;

            if opts.storePhotonEvents
                idx = cursor + (1:nPhot);
                frameCol(idx) = t;
                moleculeCol(idx) = m;
                channelCol(idx) = channels;
                detectorXCol(idx) = detector.detXY(channels,1);
                detectorYCol(idx) = detector.detXY(channels,2);
                trueXCol(idx) = pos(1);
                trueYCol(idx) = pos(2);
                arrivalTimeCol(idx) = ((t - 1) + rand(nPhot,1)) * opts.dtS;
                cursor = cursor + nPhot;
            end
        end
    end

    frameCol = frameCol(1:cursor);
    moleculeCol = moleculeCol(1:cursor);
    channelCol = channelCol(1:cursor);
    detectorXCol = detectorXCol(1:cursor);
    detectorYCol = detectorYCol(1:cursor);
    trueXCol = trueXCol(1:cursor);
    trueYCol = trueYCol(1:cursor);
    arrivalTimeCol = arrivalTimeCol(1:cursor);

    photonEvents = table(frameCol, moleculeCol, channelCol, ...
        detectorXCol, detectorYCol, trueXCol, trueYCol, arrivalTimeCol, ...
        'VariableNames', {'frame', 'molecule', 'detectorChannel', ...
        'xDetectorUm', 'yDetectorUm', 'trueXUm', 'trueYUm', ...
        'arrivalTimeS'});
end
