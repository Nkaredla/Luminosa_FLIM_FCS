function [stepTangent, noiseCovariance, stepDt, driftPerDiffusion] = ...
        collectIndependentSurfaceSteps(trackTable, reliable, usedTrackIds, ...
        curvature, opts)
    import membrane_tracking.curved_miet.*

    stepTangent = zeros(0,2);
    noiseCovariance = zeros(2,2,0);
    stepDt = zeros(0,1);
    driftPerDiffusion = zeros(0,2);

    for trackIndex = 1:numel(usedTrackIds)
        rows = find(trackTable.trackId == usedTrackIds(trackIndex));
        [~, order] = sort(trackTable.frame(rows));
        rows = rows(order);
        index = 1;
        while index < numel(rows)
            firstRow = rows(index);
            secondRow = rows(index+1);
            if ~reliable(firstRow) || ~reliable(secondRow)
                index = index + 1;
                continue;
            end

            frameGap = trackTable.frame(secondRow) - ...
                trackTable.frame(firstRow);
            if frameGap <= 0
                index = index + 1;
                continue;
            end
            [step, noise, ~, valid] = surfaceStepFromRows( ...
                trackTable, firstRow, secondRow, curvature);
            if valid
                firstPosition = [trackTable.xUm(firstRow), ...
                    trackTable.yUm(firstRow)];
                secondPosition = [trackTable.xUm(secondRow), ...
                    trackTable.yUm(secondRow)];
                midpoint = 0.5 * (firstPosition + secondPosition);
                [transform, projectedDriftPerDiffusion] = ...
                    surfaceGeometryAt(midpoint, curvature);
                dt = frameGap * opts.dtS;
                stepTangent(end+1,:) = step; %#ok<AGROW>
                noiseCovariance(:,:,end+1) = noise; %#ok<AGROW>
                stepDt(end+1,1) = dt; %#ok<AGROW>
                driftPerDiffusion(end+1,:) = ...
                    (transform * projectedDriftPerDiffusion.').'*dt; %#ok<AGROW>
            end

            % Adjacent displacement pairs share one localization error.
            index = index + 2;
        end
    end
end
