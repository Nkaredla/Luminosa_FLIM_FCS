function M = ring_build_constrained_design(tau, height, irf, timeNs, ...
        periodNs, weights)
%RING_BUILD_CONSTRAINED_DESIGN Design matrix for the height-constrained model.
%
% M = ring_build_constrained_design(tau, height, irf, timeNs, periodNs, weights)
%
% Column j is the outer product w(:, height(j)) * pattern(tau(j))', flattened
% in the same [ring, bin] order as counts(:). So the model for ring r, bin t is
%
%   sum_j a_j * w_r(z_j) * pattern_j(t)
%
% One amplitude per component, not one per ring: the ring dependence is fixed
% by the optics. That is the whole point of the constrained model - a component
% must claim a single height, and a spurious component has no height that
% explains its ring pattern.
%
% M is [nRing*nBins x numel(tau)].

    nRing = size(weights.table, 1);
    nBins = numel(timeNs);
    if numel(height) ~= numel(tau)
        error('ring_build_constrained_design:Size', ...
            'tau and height must have the same number of elements.');
    end
    M = zeros(nRing * nBins, numel(tau));
    for j = 1:numel(tau)
        w = ring_interpolate_weights(weights, height(j));
        pattern = ring_periodic_decay(irf, timeNs, periodNs, tau(j));
        M(:, j) = reshape(w(:) * pattern(:)', [], 1);
    end
end
