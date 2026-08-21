function [frames, photonEvents, photophysics] = simulatePhotonMovie( ...
        positionsUm, detector, opts)
    import membrane_tracking.focused_ism.*

    nMol = size(positionsUm, 1);
    nFrames = size(positionsUm, 3);

    activeOn = rand(nMol, 1) < opts.initialOnFraction;
    bleached = false(nMol, 1);

    onStateAtFrameStart = false(nMol, nFrames);
    onStateAtFrameEnd = false(nMol, nFrames);
    bleachedStateAtFrameStart = false(nMol, nFrames);
    bleachedStateAtFrameEnd = false(nMol, nFrames);
    excitationRateS = zeros(nMol, nFrames);
    excitationCycleRateS = zeros(nMol, nFrames);
    detectedRateS = zeros(nMol, nFrames);
    emissionActiveTimeS = zeros(nMol, nFrames);
    photonCounts = zeros(nMol, nFrames);

    for t = 1:nFrames
        if t > 1
            activeOn = updateBlinking(activeOn, bleached, opts);
        end

        pos = positionsUm(:,:,t);
        relIllum = illuminationRelative(pos, opts);
        excRate = opts.laserPowerW * opts.excitationRatePerW * relIllum;
        excitationCycleRate = excRate ./ ...
            (1 + excRate / opts.saturationRateS);
        emissionRate = opts.quantumYield * excitationCycleRate;
        % Focused excitation controls how often the molecule emits; finite
        % detector coverage separately controls how many emitted photons
        % reach the array.
        captureRelative = ismArrayCaptureRelative(pos, detector, opts);
        detRate = emissionRate * opts.collectionEfficiency * ...
            opts.detectorQuantumEfficiency .* captureRelative;

        canEmit = activeOn & ~bleached;
        activeTime = opts.dtS * double(canEmit);
        bleachHazard = opts.bleachProbabilityPerAbsorbedPhoton * ...
            excitationCycleRate;
        timeToBleach = inf(nMol, 1);
        canBleach = canEmit & bleachHazard > 0;
        timeToBleach(canBleach) = -log(max(rand(sum(canBleach), 1), ...
            realmin)) ./ bleachHazard(canBleach);
        newlyBleached = canBleach & timeToBleach < opts.dtS;
        activeTime(newlyBleached) = timeToBleach(newlyBleached);

        meanPhotons = detRate .* activeTime;
        counts = samplePoisson(meanPhotons);

        onStateAtFrameStart(:,t) = activeOn;
        bleachedStateAtFrameStart(:,t) = bleached;
        excitationRateS(:,t) = excRate;
        excitationCycleRateS(:,t) = excitationCycleRate .* ...
            (activeTime / opts.dtS);
        detectedRateS(:,t) = detRate .* (activeTime / opts.dtS);
        emissionActiveTimeS(:,t) = activeTime;
        photonCounts(:,t) = counts;

        bleached(newlyBleached) = true;
        activeOn(newlyBleached) = false;
        onStateAtFrameEnd(:,t) = activeOn;
        bleachedStateAtFrameEnd(:,t) = bleached;
    end

    frames = samplePoisson(detector.backgroundMeanCounts * ...
        ones(detector.nChannels, nFrames));
    [frames, photonEvents, photonCountsInFov] = projectMoleculePhotons( ...
        frames, positionsUm, photonCounts, detector, opts);

    photophysics = struct();
    photophysics.onStateAtFrameStart = onStateAtFrameStart;
    photophysics.onStateAtFrameEnd = onStateAtFrameEnd;
    photophysics.bleachedStateAtFrameStart = bleachedStateAtFrameStart;
    photophysics.bleachedStateAtFrameEnd = bleachedStateAtFrameEnd;
    photophysics.onState = onStateAtFrameStart;
    photophysics.bleachedState = bleachedStateAtFrameStart;
    photophysics.excitationRateS = excitationRateS;
    photophysics.excitationCycleRateS = excitationCycleRateS;
    photophysics.detectedRateS = detectedRateS;
    photophysics.emissionActiveTimeS = emissionActiveTimeS;
    photophysics.detectedPhotonCounts = photonCounts;
    photophysics.detectedPhotonCountsInFov = photonCountsInFov;
    photophysics.totalDetectedMoleculePhotonsBeforeFov = sum(photonCounts(:));
    photophysics.totalDetectedMoleculePhotonsInFov = ...
        sum(photonCountsInFov(:));
    photophysics.totalDetectedMoleculePhotons = ...
        photophysics.totalDetectedMoleculePhotonsInFov;
    photophysics.backgroundMeanCountsPerDetector = ...
        detector.backgroundMeanCounts;
end
