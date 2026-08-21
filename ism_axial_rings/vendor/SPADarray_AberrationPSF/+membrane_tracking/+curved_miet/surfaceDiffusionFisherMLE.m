function estimate = surfaceDiffusionFisherMLE(stepTangent, noiseCovariance, ...
        stepDt, driftPerDiffusion, initialD, frameDt)
    import membrane_tracking.curved_miet.*

    if isempty(stepTangent)
        estimate = struct('diffusionUm2PerS', NaN, ...
            'sigmaUm2PerS', NaN, 'information', NaN, ...
            'nSteps', 0, 'nLag1Steps', 0);
        return;
    end

    if ~isfinite(initialD) || initialD <= 0
        noiseTrace = reshape(noiseCovariance(1,1,:) + ...
            noiseCovariance(2,2,:), [], 1);
        corrected = sum(stepTangent.^2,2) - noiseTrace;
        initialD = max(mean(corrected ./ (4*stepDt)), 1e-6);
    end
    theta0 = log(max(initialD, 1e-9));
    objective = @(theta) surfaceDiffusionStepNll(theta, stepTangent, ...
        noiseCovariance, stepDt, driftPerDiffusion);
    fitOptions = optimset('Display', 'off', 'MaxIter', 180, ...
        'MaxFunEvals', 500, 'TolX', 1e-8, 'TolFun', 1e-8);
    [thetaHat, ~, exitFlag] = fminsearch(objective, theta0, fitOptions);

    if exitFlag <= 0 || ~isfinite(thetaHat)
        estimate = struct('diffusionUm2PerS', NaN, ...
            'sigmaUm2PerS', NaN, 'information', NaN, ...
            'nSteps', size(stepTangent,1), ...
            'nLag1Steps', sum(abs(stepDt-frameDt) <= ...
            10*eps(max(frameDt,1))));
        return;
    end

    diffusion = exp(min(max(thetaHat, -30), 30));
    information = expectedSurfaceDiffusionFisher(diffusion, ...
        noiseCovariance, stepDt, driftPerDiffusion);
    if information > 0 && isfinite(information)
        sigma = sqrt(1/information);
    else
        sigma = NaN;
    end

    estimate = struct();
    estimate.diffusionUm2PerS = diffusion;
    estimate.sigmaUm2PerS = sigma;
    estimate.information = information;
    estimate.nSteps = size(stepTangent,1);
    estimate.nLag1Steps = sum(abs(stepDt-frameDt) <= ...
        10*eps(max(frameDt,1)));
end
