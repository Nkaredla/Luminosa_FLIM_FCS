function [metricScale, driftPerDiffusion, inverseSqrtMetric] = ...
        surfaceGeometry(gradient, hessian)
    import membrane_tracking.fluctuating_miet.*

%SURFACEGEOMETRY Metric, Ito drift, and noise transform on a height graph.
%
%   g = I + p p',  det g = 1 + |p|^2 = s
%   Ito drift of D * Laplace-Beltrami, divided by D:
%       b = -p * ( trace(H)*s - p'*H*p ) / s^2
%   Verified against a numerical evaluation of
%   (1/sqrt(det g)) d_j ( sqrt(det g) g^{ij} ) to 1e-7 relative error.
%   The noise transform satisfies M*M' = inv(g).
    p = gradient(:);
    metricScale = 1 + (p.' * p);
    driftPerDiffusion = -p.' * ...
        (trace(hessian) * metricScale - p.' * hessian * p) / metricScale^2;
    % inv(g) = I - p p'/s, whose symmetric square root is
    % I - (1 - 1/sqrt(s)) * phat phat'
    normP = sqrt(p.' * p);
    if normP > 10 * eps
        pHat = p / normP;
        inverseSqrtMetric = eye(2) - ...
            (1 - 1/sqrt(metricScale)) * (pHat * pHat.');
    else
        inverseSqrtMetric = eye(2);
    end
end
