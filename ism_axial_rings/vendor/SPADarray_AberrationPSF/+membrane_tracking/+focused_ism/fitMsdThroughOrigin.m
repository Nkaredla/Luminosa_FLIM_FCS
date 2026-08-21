function [D, beta] = fitMsdThroughOrigin(tau, y, nPairs)
    import membrane_tracking.focused_ism.*

    tau = tau(:);
    y = y(:);
    nPairs = nPairs(:);
    if isempty(tau)
        D = NaN;
        beta = [0; NaN];
        return;
    end

    % Fisher localization noise has already been subtracted, so the
    % Brownian model is MSD_corrected(tau)=4*D*tau. Enforcing its known zero
    % intercept avoids fitting noise in a small set of short-lag points.
    weights = max(nPairs, 1);
    denominator = sum(weights .* tau.^2);
    if denominator <= 0
        D = NaN;
        beta = [0; NaN];
        return;
    end
    slope = sum(weights .* tau .* y) / denominator;
    slope = max(slope, 0);
    D = slope / 4;
    beta = [0; slope];
end
