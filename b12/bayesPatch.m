function [floatPatch, modelPatch, posteriorPatch, xRange, yRange, ok, ...
        message, photonCount, bayesTau, modelSummary] = bayesPatch( ...
        cube, regionMask, seedTau, seedTauMean, seedFraction, irf, ...
        pulsePeriodNs, dtNs, useGpu, batchSize, includeBackground, ...
        minimumPixelPhotons, minimumRegionPhotons, signalGrid, fractionGrid, ...
        singleTauGrid, shiftBounds, minimumSeparationFraction, ...
        minimumAmplitudeFraction)
%BAYESPATCH Calculate one regional Bayesian FLIM patch without structures.
% floatPatch layers:
%   1 tau mean, 2 tau std, 3 tau MAP, 4 signal fraction,
%   5-6 component lifetimes, 7-8 amplitudes, 9-10 fractions,
%   11 background amplitude.
% Pixels below the posterior photon threshold retain the selected regional
% fit, so the filled cell domain has no lifetime-map holes. posteriorPatch
% records which pixels actually received a pixel-level Bayesian posterior.

    floatPatch = [];
    modelPatch = [];
    posteriorPatch = [];
    xRange = [];
    yRange = [];
    ok = false;
    message = '';
    photonCount = 0;
    bayesTau = nan(1, 2);
    modelSummary = [0 0 0]; % requested states, used states, collapsed flag

    [x, y] = find(regionMask);
    if isempty(x)
        message = 'Empty region.';
        return;
    end
    xRange = min(x):max(x);
    yRange = min(y):max(y);
    mask = regionMask(xRange, yRange);
    croppedCube = cube(xRange, yRange, :);
    for bin = 1:size(croppedCube, 3)
        slice = croppedCube(:, :, bin);
        slice(~mask) = 0;
        croppedCube(:, :, bin) = slice;
    end
    intensity = sum(single(croppedCube), 3);
    photonCount = double(sum(intensity(:)));

    floatPatch = nan([size(mask), 11], 'single');
    modelPatch = zeros(size(mask), 'uint8');
    posteriorPatch = false(size(mask));
    seedTau = double(seedTau(isfinite(seedTau) & seedTau > 0));
    seedStateCount = min(numel(seedTau), 2);
    modelSummary = [seedStateCount, seedStateCount, 0];

    % First fill the complete boundary with regional estimates. Bayesian
    % values overwrite these estimates only where enough photons are present.
    if seedStateCount > 0
        if ~isfinite(seedTauMean)
            validSeedFraction = seedFraction(1:seedStateCount);
            validSeedFraction = validSeedFraction ./ max(sum(validSeedFraction), eps);
            seedTauMean = sum(validSeedFraction .* seedTau(:).');
        end
        floatPatch(:, :, 1) = maskedConstant(mask, seedTauMean);
        floatPatch(:, :, 3) = maskedConstant(mask, seedTauMean);
        floatPatch(:, :, 4) = maskedConstant(mask, 1);
        modelPatch(mask) = uint8(seedStateCount);
        for state = 1:seedStateCount
            fraction = seedFraction(state);
            if ~isfinite(fraction)
                fraction = 1 / seedStateCount;
            end
            floatPatch(:, :, 4 + state) = maskedConstant(mask, seedTau(state));
            amplitude = intensity * single(fraction);
            amplitude(~mask) = NaN;
            floatPatch(:, :, 6 + state) = amplitude;
            floatPatch(:, :, 8 + state) = maskedConstant(mask, fraction);
        end
        background = zeros(size(mask), 'single');
        background(~mask) = NaN;
        floatPatch(:, :, 11) = background;
    end

    if seedStateCount == 0 || photonCount < minimumRegionPhotons
        message = sprintf('Posterior skipped: %d photons or no valid seed.', round(photonCount));
        return;
    end

    % The existing Bayesian routine has a structure-based public API. It is
    % unpacked here and never passed to another B12 helper function.
    options = struct('useGPU', useGpu, 'batchSize', batchSize, ...
        'includeBackground', includeBackground, 'optimizeTau', false, ...
        'signalGrid', signalGrid, 'fractionGrid', fractionGrid, ...
        'singleExpTauGrid', singleTauGrid, 'shiftBounds', shiftBounds, ...
        'minStateSeparationFraction', minimumSeparationFraction, ...
        'minStateAmplitudeFraction', minimumAmplitudeFraction);
    try
        bayes = flim_bayes_lowphoton(croppedCube, irf, pulsePeriodNs, ...
            dtNs, seedTau, options);
    catch exception
        message = sprintf('Bayesian posterior failed: %s', exception.message);
        warning('b12:BayesPatchFailed', '%s', message);
        return;
    end

    eligible = mask & intensity >= minimumPixelPhotons;
    posteriorPatch(eligible) = true;
    floatPatch(:, :, 1) = overwriteMasked(floatPatch(:, :, 1), bayes.tauMeanArithmetic, eligible);
    floatPatch(:, :, 2) = overwriteMasked(floatPatch(:, :, 2), bayes.tauPosteriorStd, eligible);
    floatPatch(:, :, 3) = overwriteMasked(floatPatch(:, :, 3), bayes.tauMap, eligible);
    floatPatch(:, :, 4) = overwriteMasked(floatPatch(:, :, 4), bayes.signalFractionMean, eligible);

    [fittedTau, order] = sort(double(bayes.globalFit.tauFit(:)), 'ascend');
    usedStateCount = min(numel(fittedTau), 2);
    bayesTau(1:usedStateCount) = fittedTau(1:usedStateCount);
    modelPatch(eligible) = uint8(usedStateCount);
    signalCounts = intensity .* single(bayes.signalFractionMean);
    if usedStateCount == 1
        floatPatch(:, :, 5) = overwriteConstant(floatPatch(:, :, 5), eligible, fittedTau(1));
        floatPatch(:, :, 7) = overwriteMasked(floatPatch(:, :, 7), signalCounts, eligible);
        floatPatch(:, :, 9) = overwriteConstant(floatPatch(:, :, 9), eligible, 1);
        floatPatch(:, :, 6) = clearMasked(floatPatch(:, :, 6), eligible);
        floatPatch(:, :, 8) = clearMasked(floatPatch(:, :, 8), eligible);
        floatPatch(:, :, 10) = clearMasked(floatPatch(:, :, 10), eligible);
    else
        fraction1 = single(bayes.stateFractionMean);
        if order(1) ~= 1
            fraction1 = 1 - fraction1;
        end
        floatPatch(:, :, 5) = overwriteConstant(floatPatch(:, :, 5), eligible, fittedTau(1));
        floatPatch(:, :, 6) = overwriteConstant(floatPatch(:, :, 6), eligible, fittedTau(2));
        floatPatch(:, :, 7) = overwriteMasked(floatPatch(:, :, 7), signalCounts .* fraction1, eligible);
        floatPatch(:, :, 8) = overwriteMasked(floatPatch(:, :, 8), signalCounts .* (1 - fraction1), eligible);
        floatPatch(:, :, 9) = overwriteMasked(floatPatch(:, :, 9), fraction1, eligible);
        floatPatch(:, :, 10) = overwriteMasked(floatPatch(:, :, 10), 1 - fraction1, eligible);
    end
    floatPatch(:, :, 11) = overwriteMasked(floatPatch(:, :, 11), intensity - signalCounts, eligible);

    collapsed = 0;
    if isfield(bayes, 'modelSelection') && isfield(bayes.modelSelection, 'collapsed')
        collapsed = double(bayes.modelSelection.collapsed);
    end
    modelSummary = [seedStateCount, usedStateCount, collapsed];
    ok = true;
end
