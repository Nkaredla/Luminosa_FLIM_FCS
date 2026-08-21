function [stepXY, noiseCovariance, stepDt, driftTerm, ...
        gateRadialLimit] = ...
        collectIndependentSteps(trackTable, usable, usedTrackIds, ...
        surfaceModel, modes, opts)
    import membrane_tracking.fluctuating_miet.*

    stepXY = zeros(0,2);
    noiseCovariance = zeros(2,2,0);
    stepDt = zeros(0,1);
    driftTerm = zeros(0,2);
    gateRadialLimit = zeros(0,opts.gateAcceptanceQuadraturePoints);
    for k = 1:numel(usedTrackIds)
        rows = find(trackTable.trackId == usedTrackIds(k));
        [~, order] = sort(trackTable.frame(rows));
        rows = rows(order);
        i = 1;
        while i < numel(rows)
            j = i + 1;
            gap = trackTable.frame(rows(j)) - trackTable.frame(rows(i));
            if ~usable(rows(i)) || ~usable(rows(j)) || gap <= 0
                i = i + 1;
                continue;
            end
            [step, stepNoise, transform, drift, ok] = tangentStep( ...
                trackTable, rows(i), rows(j), surfaceModel, modes, opts);
            if ok
                dt = gap * opts.dtS;
                stepXY(end+1,:) = step; %#ok<AGROW>
                noiseCovariance(:,:,end+1) = stepNoise; %#ok<AGROW>
                stepDt(end+1,1) = dt; %#ok<AGROW>
                driftTerm(end+1,:) = (transform*drift.').' * dt; %#ok<AGROW>
                rawPositionCovariance = ...
                    localizationCovariance(trackTable,rows(i)) + ...
                    localizationCovariance(trackTable,rows(j));
                innovation = rawPositionCovariance + ...
                    2*opts.trackingPriorDiffusionUm2PerS*dt*eye(2);
                effectiveGateChi2 = opts.trackingGateChi2;
                if opts.minMahalanobisAcceptance > 0
                    effectiveGateChi2 = min(effectiveGateChi2, ...
                        -2*log(opts.minMahalanobisAcceptance));
                end
                radialLimit = transformedGateRadialLimit( ...
                    transform, innovation, effectiveGateChi2, opts);
                gateRadialLimit(end+1,:) = radialLimit; %#ok<AGROW>
            end
            % Adjacent displacements share one localisation error; skipping
            % both endpoints keeps the retained steps independent.
            i = i + 2;
        end
    end
end
