function estimate = surfaceDiffusionMLE(stepXY, noiseCovariance, stepDt, ...
        driftTerm, gateRadialLimit, initialD, frameDt)
    import membrane_tracking.fluctuating_miet.*

    estimate = struct('diffusionUm2PerS', NaN, 'sigmaUm2PerS', NaN, ...
        'information', NaN, 'nSteps', size(stepXY,1), ...
        'nLag1Steps', sum(abs(stepDt-frameDt) <= ...
        10*eps(max(frameDt,1))), 'meanAcceptance', NaN);
    if isempty(stepXY)
        estimate.nSteps = 0;
        return;
    end
    if ~isfinite(initialD) || initialD <= 0
        noiseTrace = reshape(noiseCovariance(1,1,:)+noiseCovariance(2,2,:), [],1);
        initialD = max(mean((sum(stepXY.^2,2)-noiseTrace)./(4*stepDt)), 1e-6);
    end
    objective = @(theta) surfaceDiffusionNll(theta, stepXY, ...
        noiseCovariance, stepDt, driftTerm, gateRadialLimit);
    fitOptions = optimset('Display','off','MaxIter',200, ...
        'MaxFunEvals',400,'TolX',1e-8,'TolFun',1e-8);
    [thetaHat, ~, exitFlag] = fminsearch(objective, ...
        log(max(initialD,1e-9)), fitOptions);
    if exitFlag <= 0 || ~isfinite(thetaHat)
        return;
    end
    Dhat = exp(min(max(thetaHat,-30),30));
    h = min(max(0.01*Dhat,eps(max(Dhat,1))),0.25*Dhat);
    information = (objective(log(Dhat+h)) - ...
        2*objective(log(Dhat)) + objective(log(Dhat-h))) / h^2;
    estimate.diffusionUm2PerS = Dhat;
    estimate.information = information;
    if information > 0 && isfinite(information)
        estimate.sigmaUm2PerS = sqrt(1/information);
    end
    estimate.nSteps = size(stepXY,1);
    acceptance = surfaceGateAcceptanceProbability(Dhat, ...
        noiseCovariance, stepDt, driftTerm, gateRadialLimit);
    estimate.meanAcceptance = meanOrNan(acceptance);
end
