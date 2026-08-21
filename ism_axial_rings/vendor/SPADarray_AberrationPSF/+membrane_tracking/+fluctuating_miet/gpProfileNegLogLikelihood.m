function [nll, beta, betaCovariance, ok] = gpProfileNegLogLikelihood( ...
        logParameters, gpData, opts, fixedCurvature, noiseCurvature)
    import membrane_tracking.fluctuating_miet.*

%GPPROFILENEGLOGLIKELIHOOD Marginal likelihood with the mean profiled out.
    if nargin < 4
        fixedCurvature = [];
    end
    if nargin < 5 || isempty(noiseCurvature)
        noiseCurvature = 0;
    end
    nll = realmax('double')/10;
    beta = [NaN; NaN];
    betaCovariance = nan(2);
    ok = false;
    if any(~isfinite(logParameters))
        return;
    end
    K = fluctuationKernelMatrix(gpData, logParameters, opts);
    noiseDiagonal = effectiveNoiseDiagonal(gpData, noiseCurvature);
    S = K + diag(noiseDiagonal);
    jitter = max(1e-12,1e-10*mean(diag(S)));
    S = S+jitter*eye(size(S));
    [L, flag] = chol(S, 'lower');
    if flag ~= 0
        return;
    end
    logDeterminant = 2*sum(log(diag(L)));
    if isempty(fixedCurvature)
        X = gpData.design;
        y = gpData.z;
    else
        X = gpData.design(:,1);
        measurementSlope = gpData.attenuation * fixedCurvature;
        y = gpData.z - measurementSlope * gpData.design(:,2);
    end
    solveX = L \ X;
    solveY = L \ y;
    information = solveX.' * solveX;
    [betaCovariance, invertible] = invertPositiveDefinite(information);
    if ~invertible
        return;
    end
    betaHat = betaCovariance * (solveX.' * solveY);
    residual = solveY - solveX*betaHat;
    nll = 0.5*(logDeterminant + residual.'*residual);
    if isempty(fixedCurvature)
        beta = betaHat;
    else
        beta = [betaHat(1); fixedCurvature];
        betaCovariance = blkdiag(betaCovariance, 0);
        betaCovariance = betaCovariance(1:2,1:2);
    end
    ok = isfinite(nll);
end
