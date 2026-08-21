function value = surfaceDiffusionStepNll(logDiffusion, stepTangent, ...
        noiseCovariance, stepDt, driftPerDiffusion)
    import membrane_tracking.curved_miet.*

    diffusion = exp(min(max(logDiffusion, -30), 30));
    value = 0;
    for step = 1:size(stepTangent,1)
        covariance = noiseCovariance(:,:,step) + ...
            2 * diffusion * stepDt(step) * eye(2);
        covariance = 0.5 * (covariance + covariance.');
        determinant = covariance(1,1)*covariance(2,2) - ...
            covariance(1,2)^2;
        if determinant <= 0 || ~isfinite(determinant)
            value = realmax('double') / 10;
            return;
        end
        inverse = [covariance(2,2), -covariance(1,2); ...
            -covariance(1,2), covariance(1,1)] / determinant;
        residual = stepTangent(step,:) - ...
            diffusion * driftPerDiffusion(step,:);
        value = value + 0.5 * ...
            (log(determinant) + residual*inverse*residual.');
    end
    if ~isfinite(value)
        value = realmax('double') / 10;
    end
end
