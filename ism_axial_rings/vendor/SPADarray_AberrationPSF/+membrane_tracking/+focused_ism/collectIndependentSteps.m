function [stepXY, noiseCov, stepDtS, gateRadialLimitUm] = ...
        collectIndependentSteps(trackTable, usedTrackIds, opts)
    import membrane_tracking.focused_ism.*

    stepXY = zeros(0, 2);
    noiseCov = zeros(2, 2, 0);
    stepDtS = zeros(0, 1);
    gateRadialLimitUm = zeros(0, opts.gateAcceptanceQuadraturePoints);

    for k = 1:numel(usedTrackIds)
        tr = trackTable(trackTable.trackId == usedTrackIds(k), :);
        tr = sortrows(tr, 'frame');
        usable = rowsUsableForDiffusion(tr, opts);
        i = 1;
        while i < height(tr)
            j = i + 1;
            if ~(usable(i) && usable(j))
                i = i + 1;
                continue;
            end
            delta = [tr.xUm(j) - tr.xUm(i), ...
                tr.yUm(j) - tr.yUm(i)];
            delta = applyMinimumImage(delta, opts);

            covI = [tr.crbVarXUm2(i), tr.crbCovXYUm2(i); ...
                tr.crbCovXYUm2(i), tr.crbVarYUm2(i)];
            covJ = [tr.crbVarXUm2(j), tr.crbCovXYUm2(j); ...
                tr.crbCovXYUm2(j), tr.crbVarYUm2(j)];
            covStep = 0.5 * (covI + covJ + (covI + covJ).');
            [~, covarianceValid] = invertPositiveDefinite(covStep);

            frameSeparation = tr.frame(j) - tr.frame(i);
            if covarianceValid && frameSeparation > 0
                deltaTime = frameSeparation * opts.dtS;
                innovationCovariance = covStep + ...
                    2 * opts.trackingPriorDiffusionUm2PerS * ...
                    deltaTime * eye(2);
                radialLimit = hardGateRadialLimit( ...
                    opts.trackingMaxStepUm * sqrt(frameSeparation), ...
                    innovationCovariance, opts);
                stepXY(end+1, :) = delta; %#ok<AGROW>
                noiseCov(:,:,end+1) = ...
                    opts.crbCorrectionScale * covStep; %#ok<AGROW>
                stepDtS(end+1, 1) = deltaTime; %#ok<AGROW>
                gateRadialLimitUm(end+1,:) = radialLimit; %#ok<AGROW>
            end

            % Consecutive displacement pairs share a localization error.
            % Skipping both endpoints makes the retained steps independent.
            i = i + 2;
        end
    end
end
