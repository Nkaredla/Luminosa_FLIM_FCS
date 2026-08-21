function logDensity = logGaussianStepDensity(stepTable, diffusionUm2PerS)
%LOGGAUSSIANSTEPDENSITY Heteroscedastic 2-D Brownian step log density.

    diffusionUm2PerS = max(double(diffusionUm2PerS), 1e-12);
    processVariance = 2 * diffusionUm2PerS .* stepTable.stepDtS;
    vxx = processVariance + stepTable.noiseVarXUm2;
    vyy = processVariance + stepTable.noiseVarYUm2;
    vxy = stepTable.noiseCovXYUm2;
    determinant = max(vxx .* vyy - vxy.^2, realmin);
    quadratic = (vyy .* stepTable.dxUm.^2 - ...
        2 * vxy .* stepTable.dxUm .* stepTable.dyUm + ...
        vxx .* stepTable.dyUm.^2) ./ determinant;
    logDensity = -log(2*pi) - 0.5 * (log(determinant) + quadratic);
    logDensity(~isfinite(logDensity)) = -realmax('double') / 100;
end
