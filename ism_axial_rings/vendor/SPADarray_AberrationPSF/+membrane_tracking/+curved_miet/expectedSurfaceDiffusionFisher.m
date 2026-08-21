function information = expectedSurfaceDiffusionFisher(diffusion, ...
        noiseCovariance, stepDt, driftPerDiffusion)
    import membrane_tracking.curved_miet.*

    information = 0;
    for step = 1:numel(stepDt)
        covariance = noiseCovariance(:,:,step) + ...
            2 * diffusion * stepDt(step) * eye(2);
        determinant = covariance(1,1)*covariance(2,2) - ...
            covariance(1,2)^2;
        if determinant <= 0 || ~isfinite(determinant)
            information = NaN;
            return;
        end
        inverse = [covariance(2,2), -covariance(1,2); ...
            -covariance(1,2), covariance(1,1)] / determinant;
        derivativeMean = driftPerDiffusion(step,:).';
        meanTerm = derivativeMean.' * inverse * derivativeMean;
        covarianceTerm = 0.5 * (2*stepDt(step))^2 * ...
            trace(inverse * inverse);
        information = information + meanTerm + covarianceTerm;
    end
end
