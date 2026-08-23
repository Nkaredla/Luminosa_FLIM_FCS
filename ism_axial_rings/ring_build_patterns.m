function P = ring_build_patterns(tau, irf, timeNs, periodNs)
%RING_BUILD_PATTERNS Decay pattern matrix, one unit-sum column per lifetime.
%
% P = ring_build_patterns(tau, irf, timeNs, periodNs)   % [nBins x numel(tau)]

    P = zeros(numel(timeNs), numel(tau));
    for j = 1:numel(tau)
        P(:, j) = ring_periodic_decay(irf, timeNs, periodNs, tau(j));
    end
end
