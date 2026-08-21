function estimate = fisherDiffusionMLE(stepXY, noiseCov, stepDtS, initialD, ...
        frameDtS, gateRadialLimitUm, opts)
    import membrane_tracking.focused_ism.*

    if isempty(stepXY)
        estimate = struct('diffusionUm2PerS', NaN, 'crbUm2PerS', NaN, ...
            'information', NaN, 'nSteps', 0, 'nLag1Steps', 0, ...
            'meanHardGateAcceptanceProbability', NaN);
        return;
    end

    if ~isfinite(initialD) || initialD <= 0
        radialNoise = squeeze(noiseCov(1,1,:) + noiseCov(2,2,:));
        radialStep = sum(stepXY.^2, 2) - radialNoise(:);
        initialD = max(mean(radialStep ./ (4 * stepDtS)), 1e-6);
    end

    theta0 = log(max(initialD, 1e-9));
    if ~opts.gateCensoringCorrection
        likelihoodGateLimit = zeros(0, ...
            opts.gateAcceptanceQuadraturePoints);
    else
        likelihoodGateLimit = gateRadialLimitUm;
    end
    objective = @(theta) diffusionStepNll(theta, stepXY, noiseCov, ...
        stepDtS, likelihoodGateLimit);
    fminOpts = optimset('Display', 'off', 'MaxIter', 200, ...
        'MaxFunEvals', 500, 'TolX', 1e-8, 'TolFun', 1e-8);
    [thetaHat, ~, exitFlag] = fminsearch(objective, theta0, fminOpts);
    if exitFlag <= 0 || ~isfinite(thetaHat)
        estimate = struct('diffusionUm2PerS', NaN, 'crbUm2PerS', NaN, ...
            'information', NaN, 'nSteps', size(stepXY, 1), ...
            'nLag1Steps', sum(abs(stepDtS - frameDtS) <= ...
            10 * eps(max(frameDtS, 1))), ...
            'meanHardGateAcceptanceProbability', NaN);
        return;
    end
    Dhat = exp(min(max(thetaHat, -30), 30));

    % Evaluate curvature with respect to D, not log(D), using the same
    % conditioned objective that was optimized.
    h = min(max(0.01 * Dhat, eps(max(Dhat,1))), 0.25 * Dhat);
    fMinus = diffusionStepNll(log(Dhat-h), stepXY, noiseCov, ...
        stepDtS, likelihoodGateLimit);
    fZero = diffusionStepNll(log(Dhat), stepXY, noiseCov, ...
        stepDtS, likelihoodGateLimit);
    fPlus = diffusionStepNll(log(Dhat+h), stepXY, noiseCov, ...
        stepDtS, likelihoodGateLimit);
    information = (fPlus - 2*fZero + fMinus) / h^2;
    if information > 0 && isfinite(information)
        crb = sqrt(1 / information);
    else
        crb = NaN;
    end

    estimate = struct();
    estimate.diffusionUm2PerS = Dhat;
    estimate.crbUm2PerS = crb;
    estimate.information = information;
    estimate.nSteps = size(stepXY, 1);
    estimate.nLag1Steps = sum(abs(stepDtS - frameDtS) <= ...
        10 * eps(max(frameDtS, 1)));
    gateAcceptance = hardGateAcceptanceProbability( ...
        Dhat, noiseCov, stepDtS, gateRadialLimitUm);
    estimate.meanHardGateAcceptanceProbability = ...
        meanOrNan(gateAcceptance);
end
