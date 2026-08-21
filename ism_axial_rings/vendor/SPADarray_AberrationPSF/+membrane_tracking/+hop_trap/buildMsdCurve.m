function msdTable = buildMsdCurve(trackTable, opts)
%BUILDMSDCURVE Compute Fisher-corrected ensemble MSD over exact frame lags.

    lagFrame = (1:opts.msdMaximumLagFrames).';
    lagTimeS = lagFrame * opts.dtS;
    rawMsdUm2 = nan(size(lagFrame));
    correctedMsdUm2 = nan(size(lagFrame));
    standardErrorUm2 = nan(size(lagFrame));
    nPairs = zeros(size(lagFrame));

    if isempty(trackTable)
        msdTable = table(lagFrame, lagTimeS, rawMsdUm2, ...
            correctedMsdUm2, standardErrorUm2, nPairs);
        return;
    end

    trackIds = unique(trackTable.trackId).';
    for lag = 1:opts.msdMaximumLagFrames
        rawValues = zeros(0, 1);
        correctedValues = zeros(0, 1);
        for trackId = trackIds
            tr = sortrows(trackTable(trackTable.trackId == trackId, :), ...
                'frame');
            if height(tr) < opts.minTrackLength
                continue;
            end
            usable = membrane_tracking.focused_ism.rowsUsableForDiffusion( ...
                tr, opts);
            [present, targetRows] = ismember(tr.frame + lag, tr.frame);
            startRows = find(present & usable);
            targetRows = targetRows(startRows);
            valid = usable(targetRows);
            startRows = startRows(valid);
            targetRows = targetRows(valid);
            if isempty(startRows)
                continue;
            end
            delta = [tr.xUm(targetRows) - tr.xUm(startRows), ...
                tr.yUm(targetRows) - tr.yUm(startRows)];
            delta = membrane_tracking.focused_ism.applyMinimumImage( ...
                delta, opts);
            squared = sum(delta.^2, 2);
            noiseTrace = opts.crbCorrectionScale * ...
                (tr.crbVarXUm2(startRows) + tr.crbVarYUm2(startRows) + ...
                tr.crbVarXUm2(targetRows) + tr.crbVarYUm2(targetRows));
            rawValues = [rawValues; squared]; %#ok<AGROW>
            correctedValues = [correctedValues; squared - noiseTrace]; ...
                %#ok<AGROW>
        end

        nPairs(lag) = numel(correctedValues);
        if nPairs(lag) >= opts.minimumMsdPairs
            rawMsdUm2(lag) = mean(rawValues);
            correctedMsdUm2(lag) = mean(correctedValues);
            if nPairs(lag) > 1
                % Overlapping lag-L displacements are correlated.
                effectivePairs = max(nPairs(lag) / lag, 1);
                standardErrorUm2(lag) = ...
                    std(correctedValues) / sqrt(effectivePairs);
            end
        end
    end

    msdTable = table(lagFrame, lagTimeS, rawMsdUm2, correctedMsdUm2, ...
        standardErrorUm2, nPairs);
end
