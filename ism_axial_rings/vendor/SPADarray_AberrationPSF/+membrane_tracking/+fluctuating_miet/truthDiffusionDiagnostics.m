function validation = truthDiffusionDiagnostics( ...
        trajectories, trackTable, modes, opts)
    import membrane_tracking.fluctuating_miet.*

    positions = trajectories.positionsXYUm;
    boundaryHit = trajectories.boundaryHit;
    [nMolecules, ~, nFrames] = size(positions);
    fullSquared = zeros(0,1);
    meanSurfaceSquared = zeros(0,1);
    for frame = 1:(nFrames-1)
        amplitudeA = 0.5*(trajectories.amplitudeHistoryA( ...
            1:modes.nModes,frame)+trajectories.amplitudeHistoryA( ...
            1:modes.nModes,frame+1));
        amplitudeB = 0.5*(trajectories.amplitudeHistoryB( ...
            1:modes.nModes,frame)+trajectories.amplitudeHistoryB( ...
            1:modes.nModes,frame+1));
        for molecule = 1:nMolecules
            if boundaryHit(molecule,frame+1)
                continue;
            end
            first = positions(molecule,:,frame);
            second = positions(molecule,:,frame+1);
            midpoint = 0.5*(first+second);
            [~, gradient, hessian] = totalSurfaceAt(midpoint, modes, ...
                amplitudeA, amplitudeB, opts);
            [metricScale, ~] = surfaceGeometry(gradient,hessian);
            transform = tangentTransformFromGradient(gradient,metricScale);
            fullSquared(end+1,1) = ...
                sum((transform*(second-first).').^2); %#ok<AGROW>
            meanGradient = opts.curvaturePerUm*midpoint;
            meanMetricScale = 1+sum(meanGradient.^2);
            meanTransform = tangentTransformFromGradient( ...
                meanGradient,meanMetricScale);
            meanSurfaceSquared(end+1,1) = sum( ...
                (meanTransform*(second-first).').^2); %#ok<AGROW>
        end
    end
    validation = struct();
    validation.configuredDiffusionUm2PerS = opts.diffusionUm2PerS;
    validation.fullSurfaceLag1DiffusionUm2PerS = ...
        meanOrNan(fullSquared) / (4*opts.dtS);
    validation.meanSurfaceLag1DiffusionUm2PerS = ...
        meanOrNan(meanSurfaceSquared) / (4*opts.dtS);
    validation.nFullSurfaceSteps = numel(fullSquared);
    validation.nBoundaryStepsExcluded = sum(boundaryHit(:));
    validation.realisedFieldRmsUm = trajectories.realisedFieldRmsUm;
    if height(trackTable) > 0
        validation.lateralLocalizationRmseUm = ...
            sqrt(mean(trackTable.truthPositionErrorUm.^2));
        validation.heightRmseUm = ...
            sqrt(mean((trackTable.heightUm - trackTable.trueZUm).^2));
        validation.heightRmseAgainstMeanSurfaceUm = sqrt(mean( ...
            (trackTable.heightUm - (trackTable.trueZUm - ...
             trackTable.trueFieldHeightUm)).^2));
    else
        validation.lateralLocalizationRmseUm = NaN;
        validation.heightRmseUm = NaN;
        validation.heightRmseAgainstMeanSurfaceUm = NaN;
    end
end
