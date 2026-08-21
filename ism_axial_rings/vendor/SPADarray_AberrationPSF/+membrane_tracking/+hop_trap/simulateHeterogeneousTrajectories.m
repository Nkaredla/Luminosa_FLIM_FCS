function trajectories = simulateHeterogeneousTrajectories(mesh, opts)
%SIMULATEHETEROGENEOUSTRAJECTORIES Simulate free, hop, and trap motion.

    import membrane_tracking.hop_trap.*

    nMolecule = opts.nMolecules;
    nFrame = opts.nFrames;
    nSubstep = opts.substepsPerFrame;
    dtSub = opts.dtS / nSubstep;

    positions = zeros(nMolecule, 2, nFrame);
    if ~isempty(opts.initialPositionsUm)
        initial = opts.initialPositionsUm;
    else
        initial = (rand(nMolecule, 2) - 0.5) .* opts.boxSizeUm;
        if opts.startAtFocus
            initial(1,:) = [0 0];
        end
    end
    positions(:,:,1) = initial;

    trapState = ones(nMolecule, nFrame, 'uint8');
    compartmentId = ones(nMolecule, nFrame);
    compartmentId(:,1) = periodicVoronoiCell(initial, mesh, opts.boxSizeUm);
    boundaryAttempts = zeros(nMolecule, 1);
    acceptedHops = zeros(nMolecule, 1);
    trapEntries = zeros(nMolecule, 1);
    trapExits = zeros(nMolecule, 1);

    hasTraps = ismember(opts.diffusionMode, {'trap', 'mixed'});
    hasHops = ismember(opts.diffusionMode, {'hop', 'mixed'});
    pEnter = 1 - exp(-dtSub / opts.meanFreeTimeS);
    pExit = 1 - exp(-dtSub / opts.meanTrappedTimeS);
    currentState = ones(nMolecule, 1, 'uint8');
    currentPosition = initial;
    currentCell = compartmentId(:,1);

    for frame = 2:nFrame
        for substep = 1:nSubstep
            if hasTraps
                enter = currentState == 1 & rand(nMolecule, 1) < pEnter;
                leave = currentState == 2 & rand(nMolecule, 1) < pExit;
                currentState(enter) = 2;
                currentState(leave) = 1;
                trapEntries = trapEntries + double(enter);
                trapExits = trapExits + double(leave);
            end

            diffusion = opts.freeDiffusionUm2PerS * ones(nMolecule, 1);
            diffusion(currentState == 2) = opts.trappedDiffusionUm2PerS;
            proposal = currentPosition + ...
                sqrt(2 * diffusion * dtSub) .* randn(nMolecule, 2);
            proposal = applyPeriodicBoundary(proposal, opts.boxSizeUm);

            if hasHops
                proposedCell = periodicVoronoiCell( ...
                    proposal, mesh, opts.boxSizeUm);
                crossed = proposedCell ~= currentCell;
                accepted = ~crossed | ...
                    rand(nMolecule, 1) < opts.hopProbability;
                boundaryAttempts = boundaryAttempts + double(crossed);
                acceptedHops = acceptedHops + double(crossed & accepted);
                currentPosition(accepted,:) = proposal(accepted,:);
                currentCell(accepted) = proposedCell(accepted);
            else
                currentPosition = proposal;
            end
        end
        positions(:,:,frame) = currentPosition;
        trapState(:,frame) = currentState;
        compartmentId(:,frame) = currentCell;
    end

    trajectories = struct();
    trajectories.positionsUm = positions;
    trajectories.trapState = trapState;
    trajectories.compartmentId = compartmentId;
    trajectories.boxSizeUm = opts.boxSizeUm;
    trajectories.dtS = opts.dtS;
    trajectories.boundaryAttemptsPerMolecule = boundaryAttempts;
    trajectories.acceptedHopsPerMolecule = acceptedHops;
    trajectories.trapEntriesPerMolecule = trapEntries;
    trajectories.trapExitsPerMolecule = trapExits;
    trajectories.hopAcceptanceFraction = sum(acceptedHops) / ...
        max(sum(boundaryAttempts), 1);
    trajectories.trappedFrameFraction = mean(trapState(:) == 2);
end
