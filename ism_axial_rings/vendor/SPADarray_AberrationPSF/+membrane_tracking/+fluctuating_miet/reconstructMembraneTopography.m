function topography = reconstructMembraneTopography( ...
        surfaceModel, trajectories, modes, opts)
    import membrane_tracking.fluctuating_miet.*

%RECONSTRUCTMEMBRANETOPOGRAPHY Posterior membrane maps at selected times.
    gridVector = linspace(-opts.membraneRadiusUm, ...
        opts.membraneRadiusUm, opts.topographyGridSize);
    [gridX, gridY] = meshgrid(gridVector, gridVector);
    inside = gridX.^2+gridY.^2 <= opts.membraneRadiusUm^2;
    queryPosition = [gridX(inside), gridY(inside)];
    nFrames = numel(trajectories.timeS);
    snapshotFrame = unique(round(linspace(1, nFrames, ...
        opts.nTopographySnapshots))).';
    nSnapshots = numel(snapshotFrame);
    outputSize = [size(gridX), nSnapshots];
    meanHeight = nan(outputSize);
    posteriorField = nan(outputSize);
    posteriorHeight = nan(outputSize);
    posteriorSigma = nan(outputSize);
    truthHeight = nan(outputSize);

    meanAtQuery = surfaceModel.tipHeightUm + ...
        0.5*surfaceModel.curvaturePerUm*sum(queryPosition.^2,2);
    hasGp = isfield(surfaceModel.gp,'isValid') && surfaceModel.gp.isValid && ...
        isfield(surfaceModel.gp,'trainingAlpha');
    for snapshot = 1:nSnapshots
        frame = snapshotFrame(snapshot);
        queryTime = trajectories.timeS(frame)*ones(size(queryPosition,1),1);
        fieldMean = zeros(size(meanAtQuery));
        fieldSigma = nan(size(meanAtQuery));
        if hasGp
            crossCovariance = fluctuationKernelCross(queryPosition, ...
                queryTime, surfaceModel.gp.trainingPositionUm, ...
                surfaceModel.gp.trainingTimeS, ...
                surfaceModel.gp.logParameters, opts, modes);
            fieldMean = crossCovariance*surfaceModel.gp.trainingAlpha;
            solved = surfaceModel.gp.trainingCholesky \ ...
                crossCovariance.';
            priorVariance = exp(2*surfaceModel.gp.logParameters(1));
            fieldSigma = sqrt(max(priorVariance-sum(solved.^2,1).',0));
        end

        truthAtQuery = zeros(size(meanAtQuery));
        amplitudeA = trajectories.amplitudeHistoryA(1:modes.nModes,frame);
        amplitudeB = trajectories.amplitudeHistoryB(1:modes.nModes,frame);
        for point = 1:size(queryPosition,1)
            truthAtQuery(point) = totalSurfaceAt(queryPosition(point,:), ...
                modes, amplitudeA, amplitudeB, opts);
        end

        layer = nan(size(gridX));
        layer(inside) = meanAtQuery;
        meanHeight(:,:,snapshot) = layer;
        layer(inside) = fieldMean;
        posteriorField(:,:,snapshot) = layer;
        layer(inside) = meanAtQuery+fieldMean;
        posteriorHeight(:,:,snapshot) = layer;
        layer(inside) = fieldSigma;
        posteriorSigma(:,:,snapshot) = layer;
        layer(inside) = truthAtQuery;
        truthHeight(:,:,snapshot) = layer;
    end

    difference = posteriorHeight-truthHeight;
    topography = struct();
    topography.isValid = surfaceModel.fitSucceeded;
    topography.xGridUm = gridX;
    topography.yGridUm = gridY;
    topography.insideMembrane = inside;
    topography.snapshotFrame = snapshotFrame;
    topography.snapshotTimeS = trajectories.timeS(snapshotFrame);
    topography.meanHeightUm = meanHeight;
    topography.posteriorFieldHeightUm = posteriorField;
    topography.posteriorHeightUm = posteriorHeight;
    topography.posteriorSigmaUm = posteriorSigma;
    topography.trueHeightUm = truthHeight;
    topography.rmseUm = sqrt(mean(difference(isfinite(difference)).^2));
    topography.kernel = opts.fluctuationKernel;
    topography.note = ['Posterior maps are conditional reconstructions from ' ...
        'heights sampled only along the particle trajectory. Away from that ' ...
        'trajectory the posterior reverts to the fitted mean surface.'];
end
