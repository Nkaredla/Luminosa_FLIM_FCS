function [model, fisher, ok] = fitSingleEmitterISM(data, detector, opts)
    import membrane_tracking.fluctuating_miet.*

    background0 = max(detector.backgroundMeanCounts, 1e-4);
    residual = max(data - background0, 0);
    residualSum = sum(residual);
    if residualSum > 0
        x0 = -sum(residual .* detector.detXY(:,1)) / residualSum;
        y0 = -sum(residual .* detector.detXY(:,2)) / residualSum;
    else
        x0 = 0; y0 = 0;
    end
    radius0 = hypot(x0, y0);
    if radius0 > 0.85 * opts.maxLocalizationRadiusUm
        scale = 0.85 * opts.maxLocalizationRadiusUm / radius0;
        x0 = x0*scale; y0 = y0*scale;
    end
    amplitude0 = max(sum(data) - detector.nChannels*background0, 1);
    theta0 = [x0, y0, log(amplitude0), log(background0)];

    objective = @(theta) singleEmitterPoissonNll(theta, data, detector, opts);
    fitOptions = optimset('Display','off','MaxIter',450, ...
        'MaxFunEvals',1800,'TolX',1e-7,'TolFun',1e-7);
    [thetaHat, nll, exitFlag] = fminsearch(objective, theta0, fitOptions);

    model = struct();
    fisher = struct('covarianceXY', nan(2), 'isValid', false);
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
    model.reducedPoissonDeviance = poissonDeviance(data, model.mu) / ...
        max(detector.nChannels - 4, 1);

    ok = hypot(model.xUm, model.yUm) <= opts.maxLocalizationRadiusUm && ...
        model.amplitude >= opts.minPhotonsPerLocalization;
    if ~ok
        return;
    end
    [probability, dPdx, dPdy] = ismDetectorChannelProbability( ...
        [model.xUm model.yUm], detector, opts);
    mu = max(model.amplitude*probability + model.background, ...
        opts.minExpectedCount);
    jacobian = [model.amplitude*dPdx, model.amplitude*dPdy, ...
        probability, ones(detector.nChannels,1)];
    information = jacobian.' * bsxfun(@rdivide, jacobian, mu);
    information = 0.5*(information + information.');
    [covariance, valid] = invertPositiveDefinite(information);
    if valid
        fisher.covarianceXY = covariance(1:2,1:2);
        fisher.isValid = true;
    end
    ok = fisher.isValid && all(diag(fisher.covarianceXY) > 0);
end
