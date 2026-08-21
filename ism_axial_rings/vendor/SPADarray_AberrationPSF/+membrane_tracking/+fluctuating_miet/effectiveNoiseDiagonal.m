function noiseDiagonal = effectiveNoiseDiagonal(gpData, curvature)
    import membrane_tracking.fluctuating_miet.*

%EFFECTIVENOISEDIAGONAL Measurement noise plus propagated input noise.
%
%   NIGP: a Gaussian process with noisy inputs is, to first order, a
%   process with exact inputs and an inflated noise variance
%   grad(m)' * C * grad(m). For m = h0 + 0.5*kappa*|q|^2 the gradient is
%   kappa*q. The full lateral Fisher covariance is retained.
%   Without this term the GP treats crowded, poorly localised frames as if
%   they carried full information about the surface.
    inputTerm = zeros(gpData.nObservations,1);
    for row = 1:gpData.nObservations
        gradient = curvature * gpData.position(row,:).';
        inputTerm(row) = gradient.' * ...
            gpData.positionCovariance(:,:,row) * gradient;
    end
    noiseDiagonal = gpData.measurementVariance + inputTerm;
end
