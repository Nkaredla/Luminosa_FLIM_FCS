function value = diffusionStepNll(theta, stepXY, noiseCov, stepDtS, ...
        gateRadialLimitUm)
    import membrane_tracking.focused_ism.*

    D = exp(min(max(theta, -30), 30));
    [vxx, vyy, vxy, determinant] = diffusionStepCovariance( ...
        D, noiseCov, stepDtS);
    dx = stepXY(:,1);
    dy = stepXY(:,2);
    quadratic = (vyy .* dx.^2 - 2 * vxy .* dx .* dy + ...
        vxx .* dy.^2) ./ determinant;
    value = 0.5 * sum(log(determinant) + quadratic);
    if nargin >= 5 && ~isempty(gateRadialLimitUm)
        acceptance = hardGateAcceptanceProbability( ...
            D, noiseCov, stepDtS, gateRadialLimitUm);
        % -log[p(step)/P(accepted)] = -log p(step) + log P(accepted).
        value = value + sum(log(max(acceptance, realmin)));
    end
    if ~isfinite(value)
        value = realmax('double') / 10;
    end
end
