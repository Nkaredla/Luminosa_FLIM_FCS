function photonData = buildCurvaturePhotonData(localizationTable, ...
        radialSquaredCorrected, photonEvents, photophysics, detector, opts)
    import membrane_tracking.curved_miet.*

    nLocalizations = height(localizationTable);
    estimatedCount = 0;
    for row = 1:nLocalizations
        frame = localizationTable.frame(row);
        first = photophysics.eventStartRow(frame);
        last = photophysics.eventEndRow(frame);
        if first > 0 && last >= first
            estimatedCount = estimatedCount + last-first+1;
        end
    end

    microtime = zeros(estimatedCount,1);
    localizationIndex = zeros(estimatedCount,1);
    signalScale = zeros(estimatedCount,1);
    backgroundScale = zeros(estimatedCount,1);
    cursor = 0;
    for row = 1:nLocalizations
        frame = localizationTable.frame(row);
        first = photophysics.eventStartRow(frame);
        last = photophysics.eventEndRow(frame);
        if first <= 0 || last < first
            continue;
        end
        eventRows = first:last;
        nEvents = numel(eventRows);
        indices = cursor + (1:nEvents);
        probability = ismDetectorChannelProbability( ...
            [localizationTable.xUm(row), localizationTable.yUm(row)], ...
            detector, opts);
        channels = photonEvents.detectorChannel(eventRows);
        signalFraction = localizationTable.signalPhotons(row) / ...
            (localizationTable.signalPhotons(row) + ...
            detector.nChannels * ...
            localizationTable.backgroundPerChannel(row));
        signalFraction = min(max(signalFraction, 1e-6), 1-1e-6);

        microtime(indices) = photonEvents.microtimeNs(eventRows);
        localizationIndex(indices) = row;
        signalScale(indices) = signalFraction * probability(channels);
        backgroundScale(indices) = (1-signalFraction) / ...
            (detector.nChannels * opts.repetitionPeriodNs);
        cursor = cursor + nEvents;
    end

    photonData = struct();
    photonData.microtimeNs = microtime(1:cursor);
    photonData.localizationIndex = localizationIndex(1:cursor);
    photonData.signalScale = signalScale(1:cursor);
    photonData.backgroundScale = backgroundScale(1:cursor);
    photonData.radialSquaredCorrectedUm2 = radialSquaredCorrected;
    photonData.nPhotons = cursor;
end
