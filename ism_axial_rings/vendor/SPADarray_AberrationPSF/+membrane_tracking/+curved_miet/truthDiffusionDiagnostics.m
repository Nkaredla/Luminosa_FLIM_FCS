function validation = truthDiffusionDiagnostics(trajectories, trackTable, ...
        trueCurvature, opts)
    import membrane_tracking.curved_miet.*

    positions = trajectories.positionsXYUm;
    boundaryHit = trajectories.boundaryHit;
    nMolecules = size(positions,1);
    nFrames = size(positions,3);

    surfaceSquared = zeros(0,1);
    projectedSquared = zeros(0,1);
    for frame = 1:(nFrames-1)
        for molecule = 1:nMolecules
            if boundaryHit(molecule,frame+1)
                continue;
            end
            first = positions(molecule,:,frame);
            second = positions(molecule,:,frame+1);
            delta = second-first;
            [transform, ~] = surfaceGeometryAt( ...
                0.5*(first+second), trueCurvature);
            tangent = (transform*delta.').';
            surfaceSquared(end+1,1) = sum(tangent.^2); %#ok<AGROW>
            projectedSquared(end+1,1) = sum(delta.^2); %#ok<AGROW>
        end
    end

    acceptedTruthSquared = zeros(0,1);
    identitySwitches = 0;
    adjacentTrackSteps = 0;
    if height(trackTable) > 0
        reliable = trackTable.associationProbability >= ...
            opts.minAssociationProbabilityForDiffusion;
        trackIds = unique(trackTable.trackId).';
        for trackId = trackIds
            rows = find(trackTable.trackId == trackId);
            [~, order] = sort(trackTable.frame(rows));
            rows = rows(order);
            for index = 1:(numel(rows)-1)
                firstRow = rows(index);
                secondRow = rows(index+1);
                if ~reliable(firstRow) || ~reliable(secondRow) || ...
                        trackTable.frame(secondRow) - ...
                        trackTable.frame(firstRow) ~= 1
                    continue;
                end
                adjacentTrackSteps = adjacentTrackSteps + 1;
                moleculeFirst = trackTable.nearestTrueMolecule(firstRow);
                moleculeSecond = trackTable.nearestTrueMolecule(secondRow);
                if moleculeFirst ~= moleculeSecond
                    identitySwitches = identitySwitches + 1;
                    continue;
                end
                molecule = round(moleculeFirst);
                frameFirst = trackTable.frame(firstRow);
                frameSecond = trackTable.frame(secondRow);
                if boundaryHit(molecule,frameSecond)
                    continue;
                end
                first = positions(molecule,:,frameFirst);
                second = positions(molecule,:,frameSecond);
                [transform, ~] = surfaceGeometryAt( ...
                    0.5*(first+second), trueCurvature);
                tangent = transform * (second-first).';
                acceptedTruthSquared(end+1,1) = ...
                    sum(tangent.^2); %#ok<AGROW>
            end
        end
    end

    validation = struct();
    validation.configuredDiffusionUm2PerS = opts.diffusionUm2PerS;
    validation.fullSurfaceLag1DiffusionUm2PerS = ...
        meanOrNan(surfaceSquared) / (4*opts.dtS);
    validation.fullProjectedLag1DiffusionUm2PerS = ...
        meanOrNan(projectedSquared) / (4*opts.dtS);
    validation.acceptedTrackTruthSurfaceLag1DiffusionUm2PerS = ...
        meanOrNan(acceptedTruthSquared) / (4*opts.dtS);
    validation.nFullSurfaceSteps = numel(surfaceSquared);
    validation.nAcceptedTruthSteps = numel(acceptedTruthSquared);
    validation.nBoundaryAffectedStepsExcluded = sum(boundaryHit(:));
    validation.identitySwitchFraction = identitySwitches / ...
        max(adjacentTrackSteps,1);
    if height(trackTable) > 0
        validation.lateralLocalizationRmseUm = sqrt(mean( ...
            trackTable.truthPositionErrorUm.^2));
        validation.lifetimeRmseNs = sqrt(mean( ...
            (trackTable.lifetimeNs-trackTable.trueLifetimeNs).^2));
        validation.heightRmseUm = sqrt(mean( ...
            (trackTable.heightUm-trackTable.trueZUm).^2));
    else
        validation.lateralLocalizationRmseUm = NaN;
        validation.lifetimeRmseNs = NaN;
        validation.heightRmseUm = NaN;
    end
end
