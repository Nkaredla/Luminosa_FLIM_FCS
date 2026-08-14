function [tauNs, tauMeanNs, coefficients, speciesFraction, fitCounts, ...
        deviance, bic, reducedChiSquare] = fitFixedOrder(counts, irf, dtNs, ...
        exponentialCount, oneStateSeed, twoStateSeed, lifetimeBounds, ...
        includeBackground, maximumIterations)
%FITFIXEDORDER Fit one specified number of IRF-convolved exponentials.
% Outputs have fixed sizes: two lifetime/fraction slots and three amplitude
% slots [background, component 1, component 2]. Unused slots are NaN or zero.

    counts = max(double(counts(:)), 0);
    lowerBound = max(lifetimeBounds(1), 0.01);
    upperBound = max(lifetimeBounds(2), lowerBound + 0.01);
    if exponentialCount == 1
        seedRows = [oneStateSeed; 0.6 * oneStateSeed; 1.6 * oneStateSeed];
    else
        seedRows = [twoStateSeed(:).'; 0.7 * twoStateSeed(:).'; ...
            1.4 * twoStateSeed(:).'; 0.3 3.0; 0.8 4.0];
    end

    bestValue = inf;
    bestLifetime = seedRows(1, :);
    options = optimset('Display', 'off', 'MaxIter', maximumIterations, ...
        'MaxFunEvals', 5 * maximumIterations, 'TolX', 1e-5, 'TolFun', 1e-5);
    for seed = 1:size(seedRows, 1)
        initialLifetime = min(max(seedRows(seed, :), lowerBound), upperBound);
        objective = @(logTau) fitObjective(logTau, counts, irf, dtNs, ...
            includeBackground, lowerBound, upperBound);
        try
            fittedLogLifetime = fminsearch(objective, log(initialLifetime), options);
            lifetime = min(max(exp(fittedLogLifetime), lowerBound), upperBound);
        catch
            lifetime = initialLifetime;
        end
        value = objective(log(lifetime));
        if value < bestValue
            bestValue = value;
            bestLifetime = lifetime;
        end
    end

    fittedLifetime = sort(bestLifetime(:), 'ascend');
    model = decayModel(irf, dtNs, fittedLifetime, includeBackground);
    rawCoefficients = lsqnonneg(model, counts);
    fitCounts = max(model * rawCoefficients, eps);
    rawSpecies = rawCoefficients(1 + double(includeBackground):end);
    rawFraction = rawSpecies ./ max(sum(rawSpecies), eps);

    tauNs = nan(1, 2);
    speciesFraction = nan(1, 2);
    coefficients = zeros(1, 3);
    tauNs(1:exponentialCount) = fittedLifetime;
    speciesFraction(1:exponentialCount) = rawFraction;
    if includeBackground
        coefficients(1) = rawCoefficients(1);
    end
    coefficients(2:1 + exponentialCount) = rawSpecies;
    tauMeanNs = sum(rawFraction .* fittedLifetime);

    deviance = poissonDeviance(counts, fitCounts);
    parameterCount = 2 * exponentialCount + double(includeBackground);
    bic = deviance + parameterCount * log(max(numel(counts), 2));
    chiSquare = sum((counts - fitCounts).^2 ./ max(fitCounts, 1));
    reducedChiSquare = chiSquare / max(numel(counts) - parameterCount, 1);
end
