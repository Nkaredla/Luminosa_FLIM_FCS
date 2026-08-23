function tau = ring_seed_tau(opts, nComponent, start)
%RING_SEED_TAU Starting lifetimes for restart number `start`.
%
% tau = ring_seed_tau(opts, nComponent, start)
%
% Start 1 uses opts.seedTauNs unchanged, so a single-start fit is exactly
% reproducible. Later starts draw log-uniformly across the allowed range, which
% is the right scale for lifetimes because the deviance surface is far closer to
% symmetric in log(tau) than in tau.
%
% Restarts exist to answer a specific question: is a reported difference
% between two models real, or is one of them simply landing in a worse local
% minimum from a fixed seed?

    if start <= 1
        tau = opts.seedTauNs(:)';
        return;
    end
    lo = log(max(opts.tauBoundsNs(1), 0.1));
    hi = log(min(opts.tauBoundsNs(2), 8));
    tau = sort(exp(lo + (hi - lo) * rand(1, nComponent)));
end
