function trajectories = simulateSurfaceTrajectories(modes, opts)
    import membrane_tracking.fluctuating_miet.*

    nMolecules = opts.nMolecules;
    nFrames = opts.nFrames;
    positionsXY = zeros(nMolecules, 2, nFrames);
    positionsZ = zeros(nMolecules, nFrames);
    fieldAtMolecule = zeros(nMolecules, nFrames);
    boundaryHit = false(nMolecules, nFrames);

    if nMolecules == 1
        positionsXY(:,:,1) = [0 0];
    else
        radius = 0.45 * opts.membraneRadiusUm * sqrt(rand(nMolecules,1));
        angle = 2*pi*rand(nMolecules,1);
        positionsXY(:,:,1) = [radius.*cos(angle), radius.*sin(angle)];
    end

    [amplitudeA, amplitudeB] = initialFieldAmplitudes(modes);
    amplitudeHistoryA = zeros(max(modes.nModes,1), nFrames);
    amplitudeHistoryB = zeros(max(modes.nModes,1), nFrames);
    dtSub = opts.dtS / opts.substepsPerFrame;

    for frame = 1:nFrames
        if frame > 1
            position = positionsXY(:,:,frame-1);
            hit = false(nMolecules, 1);
            for substep = 1:opts.substepsPerFrame
                [amplitudeA, amplitudeB] = advanceFieldAmplitudes( ...
                    amplitudeA, amplitudeB, modes, dtSub);
                [position, stepHit] = surfaceBrownianStep(position, ...
                    modes, amplitudeA, amplitudeB, dtSub, opts);
                hit = hit | stepHit;
            end
            positionsXY(:,:,frame) = position;
            boundaryHit(:,frame) = hit;
        end
        if modes.nModes > 0
            amplitudeHistoryA(:,frame) = amplitudeA;
            amplitudeHistoryB(:,frame) = amplitudeB;
        end
        for molecule = 1:nMolecules
            [z, ~, ~] = totalSurfaceAt(positionsXY(molecule,:,frame), ...
                modes, amplitudeA, amplitudeB, opts);
            positionsZ(molecule,frame) = z;
            fieldAtMolecule(molecule,frame) = z - opts.tipHeightUm - ...
                0.5 * opts.curvaturePerUm * ...
                sum(positionsXY(molecule,:,frame).^2);
        end
    end

    lifetimes = opts.lifetimeAtSubstrateNs + ...
        opts.lifetimeSlopeNsPerUm * positionsZ;
    if any(positionsZ(:) < 0)
        error('simulateFluctuatingMIETMembraneTracking:SurfaceBelowSubstrate', ...
            ['A realised fluctuation placed the membrane below z=0. Raise ' ...
             'tipHeightUm, stiffen the membrane, or reduce maxFieldRmsUm.']);
    end
    if opts.nMolecules > 1
        warning('simulateFluctuatingMIETMembraneTracking:SingleEmitterFit', ...
            ['The photon likelihood fits one emitter per frame. Frames with ' ...
             'two bright molecules are model violations; inspect ' ...
             'photophysics.fractionFramesWithMultipleBrightMolecules.']);
    end
    if any(lifetimes(:) <= opts.lifetimeBoundsNs(1)) || ...
            any(lifetimes(:) >= opts.lifetimeBoundsNs(2))
        error('simulateFluctuatingMIETMembraneTracking:LifetimeOutOfBounds', ...
            ['The realised MIET lifetimes extend outside lifetimeBoundsNs. ' ...
             'Widen the bounds or change the surface/calibration settings.']);
    end

    trajectories = struct();
    trajectories.positionsXYUm = positionsXY;
    trajectories.positionsZUm = positionsZ;
    trajectories.fieldHeightUm = fieldAtMolecule;
    trajectories.lifetimeNs = lifetimes;
    trajectories.timeS = (0:nFrames-1).' * opts.dtS;
    trajectories.boundaryHit = boundaryHit;
    trajectories.amplitudeHistoryA = amplitudeHistoryA;
    trajectories.amplitudeHistoryB = amplitudeHistoryB;
    trajectories.realisedFieldRmsUm = sqrt(mean(fieldAtMolecule(:).^2));
end
