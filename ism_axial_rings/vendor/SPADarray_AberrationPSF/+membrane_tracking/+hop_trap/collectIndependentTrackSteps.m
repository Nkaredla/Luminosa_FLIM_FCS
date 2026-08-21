function stepTable = collectIndependentTrackSteps( ...
        trackTable, trajectories, opts)
%COLLECTINDEPENDENTTRACKSTEPS Build ordered steps without shared endpoints.

    variableNames = {'sequenceId', 'trackId', 'frameStart', 'frameEnd', ...
        'timeCenterS', 'stepDtS', 'dxUm', 'dyUm', 'noiseVarXUm2', ...
        'noiseVarYUm2', 'noiseCovXYUm2', 'trueTrapState'};
    variableTypes = {'double', 'double', 'double', 'double', 'double', ...
        'double', 'double', 'double', 'double', 'double', 'double', ...
        'double'};
    stepTable = table('Size', [0 numel(variableNames)], ...
        'VariableTypes', variableTypes, 'VariableNames', variableNames);
    if isempty(trackTable) || height(trackTable) < 2
        return;
    end

    nextSequence = 1;
    trackIds = unique(trackTable.trackId).';
    for trackId = trackIds
        tr = sortrows(trackTable(trackTable.trackId == trackId, :), 'frame');
        if height(tr) < opts.minTrackLength
            continue;
        end
        usable = membrane_tracking.focused_ism.rowsUsableForDiffusion( ...
            tr, opts);
        rows = find(usable);
        pairIndex = 1;
        previousEnd = NaN;
        currentSequence = nextSequence;
        contributed = false;
        while pairIndex + 1 <= numel(rows)
            i = rows(pairIndex);
            j = rows(pairIndex + 1);
            frameSeparation = tr.frame(j) - tr.frame(i);
            pairIndex = pairIndex + 2;
            if frameSeparation < 1 || ...
                    frameSeparation > opts.maxGapFrames + 1
                previousEnd = NaN;
                continue;
            end

            if isfinite(previousEnd) && ...
                    tr.frame(i) - previousEnd > opts.maxGapFrames + 1
                currentSequence = nextSequence + 1;
                nextSequence = currentSequence;
            end

            delta = [tr.xUm(j) - tr.xUm(i), tr.yUm(j) - tr.yUm(i)];
            delta = membrane_tracking.focused_ism.applyMinimumImage( ...
                delta, opts);
            covarianceI = [tr.crbVarXUm2(i), tr.crbCovXYUm2(i); ...
                tr.crbCovXYUm2(i), tr.crbVarYUm2(i)];
            covarianceJ = [tr.crbVarXUm2(j), tr.crbCovXYUm2(j); ...
                tr.crbCovXYUm2(j), tr.crbVarYUm2(j)];
            covariance = opts.crbCorrectionScale * ...
                0.5 * (covarianceI + covarianceJ + ...
                (covarianceI + covarianceJ).');
            if any(~isfinite(covariance), 'all') || ...
                    min(eig(covariance)) < 0
                continue;
            end

            truthState = NaN;
            molecule = round(tr.nearestTrueMolecule(j));
            frameTruth = round(tr.frame(j));
            if isfinite(molecule) && molecule >= 1 && ...
                    molecule <= size(trajectories.trapState, 1) && ...
                    frameTruth >= 1 && ...
                    frameTruth <= size(trajectories.trapState, 2)
                truthState = double( ...
                    trajectories.trapState(molecule, frameTruth));
            end

            newRow = {currentSequence, trackId, tr.frame(i), tr.frame(j), ...
                0.5 * (tr.frame(i) + tr.frame(j) - 2) * opts.dtS, ...
                frameSeparation * opts.dtS, delta(1), delta(2), ...
                covariance(1,1), covariance(2,2), covariance(1,2), ...
                truthState};
            stepTable(end+1,:) = newRow; %#ok<AGROW>
            previousEnd = tr.frame(j);
            contributed = true;
        end
        if contributed
            nextSequence = nextSequence + 1;
        end
    end

    if height(stepTable) > 0
        [~, ~, compactSequence] = unique(stepTable.sequenceId, 'stable');
        stepTable.sequenceId = compactSequence;
    end
end
