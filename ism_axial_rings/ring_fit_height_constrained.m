function fit = ring_fit_height_constrained(counts, irf, timeNs, periodNs, ...
        nComponent, weights, opts)
%RING_FIT_HEIGHT_CONSTRAINED Shared lifetimes, ring pattern forced to a*w(:,z).
%
% fit = ring_fit_height_constrained(counts, irf, timeNs, periodNs, ...
%                                   nComponent, weights)
% fit = ring_fit_height_constrained(..., opts)
%
% Each component costs ONE height instead of nRing free amplitudes, so the ring
% dimension constrains the fit rather than merely enlarging it. A component
% that does not correspond to a real population at some height cannot be
% accommodated, which is the property the free-amplitude model lacks.
%
% Heights are searched on a logit scale spanning the tabulated height range, so
% the optimiser cannot walk outside the weight table where w(:, z) is
% undefined. Note the consequence: a fitted height at either endpoint of the
% table means "at or past the edge of the model", never a measurement. Callers
% should check for it - see .atBound.
%
% opts: as ring_fit_defaults, plus
%   seedHeightUm   starting heights, one per component
%   heightGroups   [1 x nComponent] group label per component; components
%                  sharing a label share ONE searched height. Default
%                  1:nComponent (every component free). Set to
%                  ones(1, nComponent) for the "everything sits at one
%                  height" model, which is the natural null for asking
%                  whether any lifetime component is axially displaced from
%                  the others. Labels apply to components in SEED ORDER,
%                  i.e. increasing seed lifetime.
%
% OUTPUT fit
%   .tau, .height   sorted by tau, with heights permuted to match
%   .amplitude      one per component
%   .deviance, .funcCount, .converged
%   .atBound        true where a fitted height is within 1% of the table edge

    if nargin < 7 || isempty(opts); opts = struct(); end
    opts = ring_fit_defaults(opts, nComponent);

    counts = double(counts);
    zLow = weights.heights(1);
    zHigh = weights.heights(end);
    if isempty(opts.seedHeightUm)
        opts.seedHeightUm = linspace(zLow + 0.01, min(zHigh, 0.5), nComponent);
    end
    if ~isfield(opts, 'heightGroups') || isempty(opts.heightGroups)
        opts.heightGroups = 1:nComponent;
    end
    groups = round(opts.heightGroups(:)');
    if numel(groups) ~= nComponent || any(groups < 1) || ...
            ~isequal(unique(groups), 1:max(groups))
        error('ring_fit_height_constrained:HeightGroups', ...
            'heightGroups must label components 1..k with no gaps.');
    end
    nHeightParam = max(groups);

    objective = @(p) transformedDeviance(p, counts, irf, timeNs, periodNs, ...
        weights, nComponent, zLow, zHigh, opts.tauBoundsNs, groups);

    best = [];
    bestValue = inf;
    bestConverged = false;
    totalEvals = 0;
    for start = 1:opts.restarts
        seedTau = ring_seed_tau(opts, nComponent, start);
        if start <= 1
            seedZ = opts.seedHeightUm;
        else
            seedZ = zLow + (zHigh - zLow) * rand(1, nComponent);
        end
        % One seed per searched height, taking the mean of the components
        % sharing that height.
        groupSeed = accumarray(groups(:), seedZ(:), [nHeightParam 1], @mean)';
        seed = [log(seedTau), toSearch(groupSeed, zLow, zHigh)];
        [candidate, value, exitflag, output] = fminsearch(objective, seed, ...
            ring_search_options(opts, numel(seed)));
        totalEvals = totalEvals + output.funcCount;
        if value < bestValue
            bestValue = value;
            best = candidate;
            bestConverged = exitflag == 1;
        end
    end

    tau = exp(best(1:nComponent));
    groupHeight = fromSearch(best(nComponent+1:end), zLow, zHigh);
    height = groupHeight(groups);
    [tau, order] = sort(tau);
    height = height(order);
    [~, amplitude] = ring_constrained_deviance(tau, height, counts, irf, ...
        timeNs, periodNs, weights);
    span = zHigh - zLow;
    fit = struct('tau', tau(:)', 'height', height(:)', ...
        'amplitude', amplitude(:)', 'deviance', bestValue, ...
        'funcCount', totalEvals, 'converged', bestConverged, ...
        'atBound', (height(:)' - zLow < 0.01 * span) | ...
                   (zHigh - height(:)' < 0.01 * span), ...
        'heightGroups', groups, 'nHeightParam', nHeightParam, ...
        'nRing', size(counts, 1), 'nComponent', nComponent);
end

function value = transformedDeviance(p, counts, irf, timeNs, periodNs, ...
        weights, nComponent, zLow, zHigh, tauBounds, groups)
    tau = exp(p(1:nComponent));
    if any(~isfinite(tau)) || any(tau <= tauBounds(1)) || any(tau > tauBounds(2))
        value = 1e12; return;
    end
    groupHeight = fromSearch(p(nComponent+1:end), zLow, zHigh);
    height = groupHeight(groups);
    value = ring_constrained_deviance(tau, height, counts, irf, timeNs, ...
        periodNs, weights);
end

function u = toSearch(z, low, high)
    frac = min(max((z - low) ./ (high - low), 1e-4), 1 - 1e-4);
    u = log(frac ./ (1 - frac));
end

function z = fromSearch(u, low, high)
    z = low + (high - low) ./ (1 + exp(-u));
end
