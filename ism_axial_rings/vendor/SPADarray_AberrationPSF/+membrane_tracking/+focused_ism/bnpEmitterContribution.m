function contribution = bnpEmitterContribution(positionsUm, detector, ...
        opts, brightnessScale)
    import membrane_tracking.focused_ism.*

    positionsUm = reshape(positionsUm, [], 2);
    nPositions = size(positionsUm, 1);
    raw = zeros(nPositions, detector.nChannels);
    nQuadrature = numel(detector.quadratureOffsetXUm);
    sigma2 = opts.psfSigmaUm^2;
    for q = 1:nQuadrature
        x = positionsUm(:,1) + detector.detXY(:,1).' + ...
            detector.quadratureOffsetXUm(q);
        y = positionsUm(:,2) + detector.detXY(:,2).' + ...
            detector.quadratureOffsetYUm(q);
        raw = raw + exp(-0.5 * (x.^2 + y.^2) / sigma2);
    end
    raw = raw / nQuadrature;

    relativeIllumination = illuminationRelative(positionsUm, opts);
    excitationRate = opts.laserPowerW * opts.excitationRatePerW * ...
        relativeIllumination;
    cycleRate = excitationRate ./ ...
        (1 + excitationRate / opts.saturationRateS);
    countScale = brightnessScale * opts.quantumYield * ...
        opts.collectionEfficiency * opts.detectorQuantumEfficiency * ...
        opts.dtS / detector.centerCaptureNormalization;
    contribution = countScale * cycleRate .* raw;
end
