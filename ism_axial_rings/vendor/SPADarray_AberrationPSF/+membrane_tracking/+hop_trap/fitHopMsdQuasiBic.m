function comparison = fitHopMsdQuasiBic(msdTable, opts)
%FITHOPMSDQUASIBIC Compare free and compartmentalized ensemble MSD models.

    import membrane_tracking.hop_trap.*

    comparison = emptyHopMsdResult();
    valid = isfinite(msdTable.correctedMsdUm2) & ...
        isfinite(msdTable.standardErrorUm2) & ...
        msdTable.standardErrorUm2 > 0 & ...
        msdTable.nPairs >= opts.minimumMsdPairs;
    tau = msdTable.lagTimeS(valid);
    observed = msdTable.correctedMsdUm2(valid);
    sigma = msdTable.standardErrorUm2(valid);
    comparison.nLagPoints = numel(tau);
    if numel(tau) < 6
        comparison.status = 'insufficient MSD lags';
        return;
    end

    sigmaFloor = max(median(sigma) * 0.1, 1e-10);
    sigma = max(sigma, sigmaFloor);
    weightedSlope = sum((4*tau) .* observed ./ sigma.^2) / ...
        max(sum((4*tau).^2 ./ sigma.^2), realmin);
    freeD = max(weightedSlope, 1e-8);
    freePrediction = 4 * freeD * tau;
    freeLogLikelihood = msdGaussianQuasiLogLikelihood( ...
        observed, freePrediction, sigma);

    shortD = max(observed(1) / max(4 * tau(1), realmin), freeD);
    amplitudeGuess = max(max(observed - 4*freeD*tau), ...
        opts.compartmentSizeUm^2 / 12);
    timeGuess = max(opts.dtS, ...
        opts.compartmentSizeUm^2 / max(4 * shortD, realmin));
    initial = [log(freeD), log(amplitudeGuess), log(timeGuess)];
    settings = optimset('Display', 'off', 'MaxIter', 1000, ...
        'MaxFunEvals', 3000, 'TolX', 1e-8, 'TolFun', 1e-8);
    objective = @(theta) -msdGaussianQuasiLogLikelihood(observed, ...
        hopMsdModel(exp(min(max(theta, -30), 30)), tau), sigma);

    bestObjective = Inf;
    bestTheta = initial;
    scales = [0.25 0.7 1.5 4];
    for scale = scales
        trial = initial + [log(scale), -log(scale), 0.5*log(scale)];
        [theta, value] = fminsearch(objective, trial, settings);
        if isfinite(value) && value < bestObjective
            bestObjective = value;
            bestTheta = theta;
        end
    end

    parameters = exp(min(max(bestTheta, -30), 30));
    hopPrediction = hopMsdModel(parameters, tau);
    hopLogLikelihood = msdGaussianQuasiLogLikelihood( ...
        observed, hopPrediction, sigma);
    n = numel(tau);
    freeBic = log(n) - 2 * freeLogLikelihood;
    hopBic = 3 * log(n) - 2 * hopLogLikelihood;

    comparison.status = 'ok';
    comparison.estimateAvailable = true;
    comparison.validRows = find(valid);
    comparison.lagTimeS = tau;
    comparison.observedMsdUm2 = observed;
    comparison.standardErrorUm2 = sigma;
    comparison.freeDiffusionUm2PerS = freeD;
    comparison.freePredictionUm2 = freePrediction;
    comparison.freeLogLikelihood = freeLogLikelihood;
    comparison.freeQuasiBic = freeBic;
    comparison.hopParameters = parameters;
    comparison.hopPredictionUm2 = hopPrediction;
    comparison.hopLogLikelihood = hopLogLikelihood;
    comparison.hopQuasiBic = hopBic;
    comparison.deltaHopQuasiBic = freeBic - hopBic;
    comparison.longRangeDiffusionUm2PerS = parameters(1);
    comparison.shortRangeDiffusionUm2PerS = ...
        parameters(1) + parameters(2) / (4 * parameters(3));
    comparison.shortToLongDiffusionRatio = ...
        comparison.shortRangeDiffusionUm2PerS / ...
        max(comparison.longRangeDiffusionUm2PerS, realmin);
    % For a uniformly explored square, the 2-D MSD plateau is L^2/3.
    comparison.effectiveCompartmentSizeUm = sqrt(3 * parameters(2));
    comparison.confinementTimeS = parameters(3);
end
