function diffusion = estimateDiffusionFromTracks(trackTable, opts)
    import membrane_tracking.focused_ism.*

    if height(trackTable) == 0
        diffusion = emptyDiffusionResult(opts);
        return;
    end

    trackIds = unique(trackTable.trackId);
    trackLengths = zeros(numel(trackIds), 1);
    reliableAssociation = trackTable.associationProbability >= ...
        opts.minAssociationProbabilityForDiffusion;
    reliableEmitterCount = ...
        trackTable.emitterCountPosteriorProbability >= ...
        opts.minEmitterCountPosteriorProbabilityForDiffusion;
    % A displacement enters diffusion inference only when both the emitter
    % count and the temporal identity are sufficiently supported.
    reliableLocalization = rowsUsableForDiffusion(trackTable, opts);
    for k = 1:numel(trackIds)
        trackLengths(k) = sum(trackTable.trackId == trackIds(k) & ...
            reliableLocalization);
    end
    usedTrackIds = trackIds(trackLengths >= opts.minTrackLength);

    lagVec = (1:opts.maxMsdLag).';
    nPairs = zeros(opts.maxMsdLag, 1);
    msd = nan(opts.maxMsdLag, 1);
    locNoise = nan(opts.maxMsdLag, 1);
    msdCorrected = nan(opts.maxMsdLag, 1);

    for lag = 1:opts.maxMsdLag
        r2All = [];
        noiseAll = [];

        for k = 1:numel(usedTrackIds)
            tr = trackTable(trackTable.trackId == usedTrackIds(k), :);
            tr = sortrows(tr, 'frame');
            n = height(tr);
            reliableInTrack = rowsUsableForDiffusion(tr, opts);
            for i = 1:(n-1)
                for j = (i+1):n
                    frameLag = tr.frame(j) - tr.frame(i);
                    if frameLag > lag
                        break;
                    end
                    if frameLag == lag
                        if ~all(reliableInTrack(i:j))
                            continue;
                        end
                        delta = [tr.xUm(j) - tr.xUm(i), ...
                            tr.yUm(j) - tr.yUm(i)];
                        delta = applyMinimumImage(delta, opts);
                        dx = delta(1);
                        dy = delta(2);
                        noiseX = opts.crbCorrectionScale * ...
                            (tr.crbVarXUm2(i) + tr.crbVarXUm2(j));
                        noiseY = opts.crbCorrectionScale * ...
                            (tr.crbVarYUm2(i) + tr.crbVarYUm2(j));
                        r2All(end+1, 1) = dx^2 + dy^2; %#ok<AGROW>
                        % E[|dr|^2] = 4*D*lagTime plus the trace of the two
                        % independent localization covariance matrices.
                        noiseAll(end+1, 1) = noiseX + noiseY; %#ok<AGROW>
                    end
                end
            end
        end

        nPairs(lag) = numel(r2All);
        if nPairs(lag) > 0
            msd(lag) = mean(r2All);
            locNoise(lag) = mean(noiseAll);
            msdCorrected(lag) = mean(r2All - noiseAll);
        end
    end

    lagTimeS = lagVec * opts.dtS;
    msdTable = table(lagVec, lagTimeS, nPairs, msd, locNoise, ...
        msdCorrected, 'VariableNames', {'lagFrames', 'lagTimeS', ...
        'nPairs', 'msdUm2', 'meanFisherLocalizationNoiseUm2', ...
        'msdCorrectedUm2'});

    valid = isfinite(msdCorrected) & nPairs >= opts.minMsdPairs;
    validFit = valid & lagVec <= opts.maxMsdFitLag;
    [DcorrFree, interceptCorrFree, betaCorrFree] = fitMsdSlope( ...
        lagTimeS(validFit), msdCorrected(validFit), nPairs(validFit));
    [Dcorr, betaCorr] = fitMsdThroughOrigin(lagTimeS(validFit), ...
        msdCorrected(validFit), nPairs(validFit));
    [Draw, interceptRaw, betaRaw] = fitMsdSlope(lagTimeS(validFit), ...
        msd(validFit), nPairs(validFit));

    lagOne = find(lagVec == 1 & isfinite(msdCorrected), 1, 'first');
    if isempty(lagOne)
        DLag1 = NaN;
        DLag1Unconstrained = NaN;
    else
        DLag1Unconstrained = msdCorrected(lagOne) / (4 * opts.dtS);
        DLag1 = max(DLag1Unconstrained, 0);
    end
    [stepXY, stepNoiseCov, stepDtS, gateRadialLimitUm] = ...
        collectIndependentSteps( ...
        trackTable, usedTrackIds, opts);
    fisherD = fisherDiffusionMLE(stepXY, stepNoiseCov, stepDtS, DLag1, ...
        opts.dtS, gateRadialLimitUm, opts);

    locSigma = sqrt(max(0.5 * (trackTable.crbVarXUm2 + ...
        trackTable.crbVarYUm2), 0));
    usedMask = ismember(trackTable.trackId, usedTrackIds);
    usedSigma = locSigma(usedMask & reliableLocalization & ...
        isfinite(locSigma));

    diffusion = struct();
    diffusion.trueDiffusionUm2PerS = opts.diffusionUm2PerS;
    diffusion.DLag1CorrectedUm2PerS = DLag1;
    diffusion.DLag1CorrectedUnconstrainedUm2PerS = ...
        DLag1Unconstrained;
    diffusion.DMsdCorrectedUm2PerS = Dcorr;
    diffusion.DMsdCorrectedFreeInterceptUm2PerS = DcorrFree;
    diffusion.DMsdRawUm2PerS = Draw;
    diffusion.DFisherMLEUm2PerS = fisherD.diffusionUm2PerS;
    diffusion.DFisherCRBUm2PerS = fisherD.crbUm2PerS;
    diffusion.diffusionFisherInformation = fisherD.information;
    diffusion.meanHardGateAcceptanceProbability = ...
        fisherD.meanHardGateAcceptanceProbability;
    diffusion.hardGateCensoringCorrectionApplied = ...
        opts.gateCensoringCorrection;
    % Backward-compatible aliases used by the bias-decomposition utility.
    diffusion.meanGateAcceptanceProbability = ...
        diffusion.meanHardGateAcceptanceProbability;
    diffusion.gateCensoringCorrectionApplied = ...
        diffusion.hardGateCensoringCorrectionApplied;
    diffusion.gateAcceptanceModel = ...
        ['Exact angular quadrature for the Euclidean/Mahalanobis hard-gate ' ...
        'intersection; JPDA confidence and focus selection are not included'];
    diffusion.nIndependentStepsForFisher = fisherD.nSteps;
    diffusion.nLag1PairsForFisher = fisherD.nLag1Steps;
    diffusion.msdInterceptCorrectedUm2 = 0;
    diffusion.msdInterceptCorrectedFreeUm2 = interceptCorrFree;
    diffusion.msdInterceptRawUm2 = interceptRaw;
    diffusion.msdSlopeCorrected = betaCorr;
    diffusion.msdSlopeCorrectedFree = betaCorrFree;
    diffusion.msdSlopeRaw = betaRaw;
    diffusion.msdTable = msdTable;
    diffusion.nTracksTotal = numel(trackIds);
    diffusion.nTracksUsed = numel(usedTrackIds);
    diffusion.usedTrackIds = usedTrackIds;
    diffusion.meanLocalizationSigmaUm = meanOrNan(usedSigma);
    diffusion.medianLocalizationSigmaUm = medianOrNan(usedSigma);
    diffusion.nLocalizationsTotal = height(trackTable);
    diffusion.nAssociationReliableLocalizations = ...
        sum(reliableAssociation);
    diffusion.nEmitterCountReliableLocalizations = ...
        sum(reliableEmitterCount);
    diffusion.nReliableLocalizationsForDiffusion = ...
        sum(reliableLocalization);
    diffusion.fractionAmbiguousAssociations = ...
        1 - sum(reliableAssociation) / max(height(trackTable), 1);
    diffusion.fractionUncertainEmitterCounts = ...
        1 - sum(reliableEmitterCount) / max(height(trackTable), 1);
    diffusion.fractionRejectedForDiffusion = ...
        1 - sum(reliableLocalization) / max(height(trackTable), 1);
    diffusion.meanAssociationProbability = meanOrNan( ...
        trackTable.associationProbability( ...
        isfinite(trackTable.associationProbability)));
    diffusion.meanEmitterCountPosteriorProbability = meanOrNan( ...
        trackTable.emitterCountPosteriorProbability( ...
        isfinite(trackTable.emitterCountPosteriorProbability)));
    if ismember('positionIdentifiable', ...
            trackTable.Properties.VariableNames)
        diffusion.nPositionIdentifiableLocalizations = ...
            sum(trackTable.positionIdentifiable);
        diffusion.fractionNonidentifiablePositions = ...
            1 - mean(trackTable.positionIdentifiable);
    else
        diffusion.nPositionIdentifiableLocalizations = height(trackTable);
        diffusion.fractionNonidentifiablePositions = 0;
    end
end
