function probability = hardGateAcceptanceProbability(D, noiseCov, ...
        stepDtS, gateRadialLimitUm)
    import membrane_tracking.focused_ism.*

    if isempty(gateRadialLimitUm)
        probability = ones(size(stepDtS));
        return;
    end

    nAngles = size(gateRadialLimitUm, 2);
    angles = ((0:nAngles-1) + 0.5) * (2*pi/nAngles);
    ux = cos(angles);
    uy = sin(angles);
    [vxx, vyy, vxy, determinant] = diffusionStepCovariance( ...
        D, noiseCov, stepDtS);
    probability = zeros(size(stepDtS));
    for step = 1:numel(stepDtS)
        inverseXX = vyy(step) / determinant(step);
        inverseYY = vxx(step) / determinant(step);
        inverseXY = -vxy(step) / determinant(step);
        directionalPrecision = inverseXX*ux.^2 + ...
            2*inverseXY*ux.*uy + inverseYY*uy.^2;
        exponent = -0.5 * directionalPrecision .* ...
            gateRadialLimitUm(step,:).^2;
        radialIntegral = -expm1(exponent) ./ ...
            max(directionalPrecision, realmin);
        probability(step) = mean(radialIntegral) / ...
            sqrt(determinant(step));
    end
    probability = min(max(probability, realmin), 1);
end
