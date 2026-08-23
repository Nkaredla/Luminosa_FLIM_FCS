function [value, amplitude] = ring_constrained_deviance(tau, height, ...
        counts, irf, timeNs, periodNs, weights)
%RING_CONSTRAINED_DEVIANCE Deviance of the height-constrained model.
%
% [value, amplitude] = ring_constrained_deviance(tau, height, counts, irf, ...
%                                               timeNs, periodNs, weights)
%
% Takes lifetimes and heights in PHYSICAL units, not in any search
% parameterisation, so it can be used to profile the objective directly - fix
% every parameter but one, scan it, and see whether the objective actually has
% curvature there. That is the only honest way to tell an identified parameter
% from one the optimiser simply left near its starting value.

    tau = tau(:)';
    height = height(:)';
    if any(~isfinite(tau)) || any(tau <= 0) || any(~isfinite(height))
        value = 1e12; amplitude = zeros(1, numel(tau)); return;
    end
    M = ring_build_constrained_design(tau, height, irf, timeNs, periodNs, ...
        weights);
    amplitude = lsqnonneg(M, counts(:));
    value = ring_poisson_deviance(counts(:), M * amplitude);
    if ~isfinite(value); value = 1e12; end
    amplitude = amplitude(:)';
end
