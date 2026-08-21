function trajectories = simulateTrajectories(opts)
    import membrane_tracking.focused_ism.*

    nMol = opts.nMolecules;
    nFrames = opts.nFrames;
    positions = zeros(nMol, 2, nFrames);

    if isempty(opts.initialPositionsUm)
        positions(:,:,1) = (rand(nMol, 2) - 0.5) .* opts.boxSizeUm;
    else
        init = double(opts.initialPositionsUm);
        if ~isequal(size(init), [nMol 2])
            error('simulateMembraneDiffusionParticleTracking:BadInitialPositions', ...
                'initialPositionsUm must have size [nMolecules x 2].');
        end
        if any(~isfinite(init(:))) || ~isreal(init)
            error('simulateMembraneDiffusionParticleTracking:BadInitialPositions', ...
                'initialPositionsUm must contain finite real values.');
        end
        positions(:,:,1) = applyBoundary(init, opts.boxSizeUm, ...
            opts.boundaryCondition);
    end

    dtSub = opts.dtS / opts.substepsPerFrame;
    stepSigma = sqrt(2 * opts.diffusionUm2PerS * dtSub);

    for t = 2:nFrames
        pos = positions(:,:,t-1);
        for k = 1:opts.substepsPerFrame
            pos = pos + stepSigma * randn(nMol, 2);
            pos = applyBoundary(pos, opts.boxSizeUm, opts.boundaryCondition);
        end
        positions(:,:,t) = pos;
    end

    trajectories = struct();
    trajectories.positionsUm = positions;
    trajectories.timeS = (0:nFrames-1).' * opts.dtS;
    trajectories.boxSizeUm = opts.boxSizeUm;
    trajectories.diffusionUm2PerS = opts.diffusionUm2PerS;
    trajectories.boundaryCondition = opts.boundaryCondition;
end
