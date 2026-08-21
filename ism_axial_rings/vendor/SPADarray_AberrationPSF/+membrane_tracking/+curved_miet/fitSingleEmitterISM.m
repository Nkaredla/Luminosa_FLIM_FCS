function [model, fisher, ok] = fitSingleEmitterISM(data, detector, opts)
    import membrane_tracking.curved_miet.*

    background0 = max(detector.backgroundMeanCounts, 1e-4);
    residual = max(data - background0, 0);
    residualSum = sum(residual);
    if residualSum > 0
        x0 = -sum(residual .* detector.detXY(:,1)) / residualSum;
        y0 = -sum(residual .* detector.detXY(:,2)) / residualSum;
    else
        x0 = 0;
        y0 = 0;
    end
    radius0 = hypot(x0, y0);
    if radius0 > 0.85 * opts.maxLocalizationRadiusUm
        scale = 0.85 * opts.maxLocalizationRadiusUm / radius0;
        x0 = x0 * scale;
        y0 = y0 * scale;
    end

    amplitude0 = max(sum(data) - ...
        detector.nChannels * background0, 1);
    theta0 = [x0, y0, log(amplitude0), log(background0)];
    objective = @(theta) singleEmitterPoissonNll( ...
        theta, data, detector, opts);
    fitOptions = optimset('Display', 'off', 'MaxIter', 450, ...
        'MaxFunEvals', 1800, 'TolX', 1e-7, 'TolFun', 1e-7);
    [thetaHat, nll, exitFlag] = fminsearch(objective, theta0, fitOptions);

    model = struct();
    fisher = struct('information', nan(4), ...
        'covariance', nan(4), 'covarianceXY', nan(2), 'isValid', false);
    ok = exitFlag > 0 && all(isfinite(thetaHat)) && isfinite(nll);
    if ~ok
        return;
    end

    model.xUm = thetaHat(1);
    model.yUm = thetaHat(2);
    model.amplitude = expClamped(thetaHat(3));
    model.background = expClamped(thetaHat(4));
    model.probability = ismDetectorChannelProbability( ...
        [model.xUm model.yUm], detector, opts);
    model.mu = model.amplitude * model.probability + model.background;
    model.negLogLikelihood = nll;
    model.poissonDeviance = poissonDeviance(data, model.mu);
    model.reducedPoissonDeviance = model.poissonDeviance / ...
        max(detector.nChannels - 4, 1);

    inside = hypot(model.xUm, model.yUm) <= opts.maxLocalizationRadiusUm;
    ok = inside && model.amplitude >= opts.minPhotonsPerLocalization;
    if ~ok
        return;
    end

    fisher = singleEmitterFisher(model, detector, opts);
    ok = fisher.isValid && all(diag(fisher.covarianceXY) > 0) && ...
        all(isfinite(fisher.covarianceXY(:)));
end
