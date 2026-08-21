function value = surfaceDiffusionNll(logD, stepXY, noiseCovariance, ...
        stepDt, driftTerm, gateRadialLimit)
    import membrane_tracking.fluctuating_miet.*

    D = exp(min(max(logD,-30),30));
    value = 0;
    for s = 1:size(stepXY,1)
        covariance = noiseCovariance(:,:,s) + 2*D*stepDt(s)*eye(2);
        determinant = covariance(1,1)*covariance(2,2) - covariance(1,2)^2;
        if determinant <= 0 || ~isfinite(determinant)
            value = realmax('double')/10;
            return;
        end
        inverse = [covariance(2,2), -covariance(1,2); ...
            -covariance(1,2), covariance(1,1)] / determinant;
        residual = stepXY(s,:) - D*driftTerm(s,:);
        value = value + 0.5*(log(determinant) + residual*inverse*residual.');
    end
    if ~isempty(gateRadialLimit)
        % Every retained step passed a displacement gate, so the likelihood
        % is conditioned on acceptance. The probability uses the actual
        % anisotropic gate and retains the small non-zero Ito drift.
        acceptance = surfaceGateAcceptanceProbability(D, ...
            noiseCovariance, stepDt, driftTerm, gateRadialLimit);
        value = value+sum(log(max(acceptance,realmin)));
    end
    if ~isfinite(value)
        value = realmax('double')/10;
    end
end
