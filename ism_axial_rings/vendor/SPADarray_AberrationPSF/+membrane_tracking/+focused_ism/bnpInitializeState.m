function state = bnpInitializeState(frames, detector, trackTable, ...
        initialD, opts)
    import membrane_tracking.focused_ism.*

    data = double(frames).';
    [nFrames, nChannels] = size(data);
    nCandidates = opts.bnpPoolSize;

    positions = zeros(nCandidates, nFrames, 2);
    firstPositions = (rand(nCandidates, 2) - 0.5) .* opts.boxSizeUm;
    positions(:,1,:) = reshape(firstPositions, nCandidates, 1, 2);
    stepSigma = sqrt(2 * initialD * opts.dtS);
    for frame = 2:nFrames
        previous = reshape(positions(:,frame-1,:), nCandidates, 2);
        current = previous + stepSigma * randn(nCandidates, 2);
        current = applyBoundary(current, opts.boxSizeUm, ...
            opts.boundaryCondition);
        positions(:,frame,:) = reshape(current, nCandidates, 1, 2);
    end

    active = false(nCandidates, 1);
    if height(trackTable) > 0
        usableRows = isfinite(trackTable.xUm) & isfinite(trackTable.yUm) & ...
            trackTable.trackId > 0;
        if ismember('positionIdentifiable', ...
                trackTable.Properties.VariableNames)
            usableRows = usableRows & logical(trackTable.positionIdentifiable);
        end
        seedTable = trackTable(usableRows,:);
        trackIds = unique(seedTable.trackId);
        trackLengths = zeros(numel(trackIds), 1);
        for k = 1:numel(trackIds)
            trackLengths(k) = sum(seedTable.trackId == trackIds(k));
        end
        [~, order] = sort(trackLengths, 'descend');
        nSeeds = min(nCandidates, numel(order));
        for seed = 1:nSeeds
            tr = seedTable(seedTable.trackId == trackIds(order(seed)), :);
            positions(seed,:,:) = reshape( ...
                bnpTrackSeedPath(tr, initialD, nFrames, opts), ...
                1, nFrames, 2);
            active(seed) = true;
        end
    end

    contributions = zeros(nCandidates, nFrames, nChannels);
    initialBrightnessScale = opts.bnpCalibratedBrightnessScale;
    for candidate = 1:nCandidates
        path = reshape(positions(candidate,:,:), nFrames, 2);
        contributions(candidate,:,:) = reshape( ...
            bnpEmitterContribution(path, detector, opts, ...
            initialBrightnessScale), ...
            1, nFrames, nChannels);
    end

    state = struct();
    state.data = data;
    state.nFrames = nFrames;
    state.nChannels = nChannels;
    state.nCandidates = nCandidates;
    state.positionsUm = positions;
    state.active = active;
    state.contributions = contributions;
    state.diffusionUm2PerS = initialD;
    state.backgroundCountsPerChannel = ...
        max(detector.backgroundMeanCounts, 1e-9);
    state.brightnessScale = initialBrightnessScale;
    state.meanCounts = bnpRebuildMeanCounts(state);
end
