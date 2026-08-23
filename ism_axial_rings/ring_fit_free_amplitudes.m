function fit = ring_fit_free_amplitudes(counts, irf, timeNs, periodNs, ...
        nComponent, opts)
%RING_FIT_FREE_AMPLITUDES Shared lifetimes, free non-negative amplitude per ring.
%
% fit = ring_fit_free_amplitudes(counts, irf, timeNs, periodNs, nComponent)
% fit = ring_fit_free_amplitudes(..., opts)
%
% counts  [nRing x nBins] photon counts. With a SINGLE row this is the ordinary
%         summed-decay multiexponential fit, so the same function serves both
%         the "summed" and the "ring-free" model and no code path differs
%         between them.
%
% Variable projection: the amplitudes enter linearly, so only the
% log-lifetimes are searched and the amplitudes come from an inner
% non-negative least squares. That removes the amplitude-lifetime degeneracy
% which makes naive multiexponential fitting unstable.
%
% opts fields (all optional)
%   restarts        number of starts; > 1 randomises the seed and keeps the
%                   best deviance. Default 1.
%   seedTauNs       explicit starting lifetimes for the first start.
%                   Default logspace(0.3, 4, nComponent).
%   evalsPerParam   fminsearch MaxFunEvals per searched parameter, default 400.
%                   Set evalBudget instead to give every model the same total.
%   evalBudget      absolute MaxFunEvals, overrides evalsPerParam.
%   tolFun, tolX    default 1e-4 each.
%   tauBoundsNs     [lo hi] hard reject outside, default [0.05 20].
%
% OUTPUT fit
%   .tau            sorted lifetimes, ns
%   .amplitude      [nComponent x nRing]
%   .deviance       Poisson deviance at the optimum
%   .funcCount      total objective evaluations used, summed over restarts
%   .converged      true if the winning start met fminsearch's tolerances
%   .nRing, .nComponent

    if nargin < 6 || isempty(opts); opts = struct(); end
    opts = ring_fit_defaults(opts, nComponent);

    counts = double(counts);
    if isvector(counts); counts = reshape(counts, 1, []); end
    nRing = size(counts, 1);

    objective = @(p) freeDeviance(p, counts, irf, timeNs, periodNs, ...
        opts.tauBoundsNs);

    best = [];
    bestValue = inf;
    bestConverged = false;
    totalEvals = 0;
    for start = 1:opts.restarts
        seed = log(ring_seed_tau(opts, nComponent, start));
        [candidate, value, exitflag, output] = fminsearch(objective, seed, ...
            ring_search_options(opts, numel(seed)));
        totalEvals = totalEvals + output.funcCount;
        if value < bestValue
            bestValue = value;
            best = candidate;
            bestConverged = exitflag == 1;
        end
    end

    tau = sort(exp(best));
    P = ring_build_patterns(tau, irf, timeNs, periodNs);
    amplitude = zeros(nComponent, nRing);
    for r = 1:nRing
        amplitude(:, r) = lsqnonneg(P, counts(r, :)');
    end
    fit = struct('tau', tau(:)', 'amplitude', amplitude, ...
        'deviance', bestValue, 'funcCount', totalEvals, ...
        'converged', bestConverged, 'nRing', nRing, ...
        'nComponent', nComponent);
end

function value = freeDeviance(logTau, counts, irf, timeNs, periodNs, bounds)
    tau = exp(logTau);
    if any(~isfinite(tau)) || any(tau <= bounds(1)) || any(tau > bounds(2))
        value = 1e12; return;
    end
    P = ring_build_patterns(tau, irf, timeNs, periodNs);
    value = 0;
    for r = 1:size(counts, 1)
        y = counts(r, :)';
        a = lsqnonneg(P, y);
        value = value + ring_poisson_deviance(y, P * a);
    end
    if ~isfinite(value); value = 1e12; end
end
