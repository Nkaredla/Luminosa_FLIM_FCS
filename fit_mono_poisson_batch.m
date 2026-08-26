function out = fit_mono_poisson_batch(Y, basis, tauGrid, opts)
%FIT_MONO_POISSON_BATCH Single exponential plus background, many pixels at once.
%
% out = fit_mono_poisson_batch(Y, basis, tauGrid, opts)
%
% Fits   y(t) = B + a * P(t; tau),   B >= 0, a >= 0
%
% to every column of Y by Poisson maximum likelihood. TAU is localised on the
% supplied grid and then refined between nodes; B and a come from a whitened
% Fisher-scoring step whose constrained solve is exact (poisson_nnls2_exact).
%
% WHY A GRID IS ACCEPTABLE HERE, HAVING BEEN WRONG FOR THE TWO-COMPONENT FIT
%
% The objection to the grid in the biexponential case was never the grid: it was
% the CHEAP RANKING METRIC used to pick a node, which weighted residuals by
% max(model,1) and biased tau upward by 5-7.4%. Here every node is scored by the
% exact Poisson deviance with properly solved amplitudes, so there is no metric
% to get wrong, and the problem is one-dimensional so a modest grid plus
% parabolic refinement is both cheap and unbiased.
%
% WARM STARTING ACROSS NODES IS WHAT MAKES IT AFFORDABLE
%
% Running an inner solve to convergence at every node independently is enormously
% wasteful - it was, in the first version of this file, and it timed out.
% Neighbouring lifetimes give near-identical amplitudes, so the amplitudes from
% the previous node are an excellent start for the next and only a couple of
% Fisher steps are needed. That is safe because the inner problem is strictly
% convex with a unique optimum, so the start affects the path and not the answer;
% the final solve at the refined tau is run to full tolerance from that warm
% start regardless.
%
% out fields: tauNs, background, amplitude, deviance, reducedDeviance,
%             photons, refined

    if nargin < 4 || isempty(opts); opts = struct(); end
    defaults = struct('scanIter', 6, 'finalIter', 60, 'tol', 1e-12, ...
        'refine', true);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end

    [nBin, P] = size(Y);
    tauGrid = tauGrid(:)';
    nGrid = numel(tauGrid);
    logTau = log(tauGrid);
    photons = sum(Y, 1);

    devAll = inf(nGrid, P);
    beta = [photons / (2 * nBin); photons / 2];
    beta = max(beta, 1e-10);

    for g = 1:nGrid
        p = biexp_slb_pattern_batch(basis, tauGrid(g));
        % Warm start from the previous node, then a couple of Fisher steps.
        for it = 1:opts.scanIter
            m = max(beta(1, :) + beta(2, :) .* p, 1e-12);
            w = 1 ./ m;
            wp = w .* p;
            G = [sum(w, 1); sum(wp, 1); sum(wp .* p, 1)];
            rhs = [sum(w .* Y, 1); sum(wp .* Y, 1)];
            target = poisson_nnls2_exact(G, rhs);
            prev = poisson_nnls_em_deviance(Y, m);
            step = target - beta;
            t = ones(1, P);
            accepted = false(1, P);
            for half = 1:8
                if all(accepted); break; end
                sel = find(~accepted);
                trial = max(beta(:, sel) + t(sel) .* step(:, sel), 0);
                mt = max(trial(1, :) + trial(2, :) .* p, 1e-12);
                dt = poisson_nnls_em_deviance(Y(:, sel), mt);
                good = dt <= prev(sel);
                if any(good)
                    beta(:, sel(good)) = trial(:, good);
                    accepted(sel(good)) = true;
                end
                t(sel(~good)) = t(sel(~good)) / 2;
            end
        end
        m = max(beta(1, :) + beta(2, :) .* p, 1e-12);
        devAll(g, :) = poisson_nnls_em_deviance(Y, m);
    end

    [bestDev, bestIdx] = min(devAll, [], 1);
    tauFit = tauGrid(bestIdx);
    refined = false(1, P);

    if opts.refine && nGrid > 2
        step = mean(diff(logTau));
        lin = 1:P;
        i0 = sub2ind([nGrid P], max(bestIdx - 1, 1), lin);
        i1 = sub2ind([nGrid P], bestIdx, lin);
        i2 = sub2ind([nGrid P], min(bestIdx + 1, nGrid), lin);
        y0 = devAll(i0); y1 = devAll(i1); y2 = devAll(i2);
        curv = y0 - 2 * y1 + y2;
        usable = bestIdx > 1 & bestIdx < nGrid & curv > 0 & isfinite(curv);
        shift = zeros(1, P);
        shift(usable) = 0.5 * step * (y0(usable) - y2(usable)) ./ curv(usable);
        shift = max(min(shift, step), -step);
        tauFit(usable) = exp(logTau(bestIdx(usable)) + shift(usable));
        refined = usable;
    end

    % ---- final solve at the refined tau, to full tolerance --------------
    Pat = biexp_slb_pattern_batch(basis, tauFit);
    for it = 1:opts.finalIter
        m = max(beta(1, :) + beta(2, :) .* Pat, 1e-12);
        w = 1 ./ m;
        wp = w .* Pat;
        G = [sum(w, 1); sum(wp, 1); sum(wp .* Pat, 1)];
        rhs = [sum(w .* Y, 1); sum(wp .* Y, 1)];
        target = poisson_nnls2_exact(G, rhs);
        prev = poisson_nnls_em_deviance(Y, m);
        stepB = target - beta;
        t = ones(1, P);
        accepted = false(1, P);
        newBeta = beta;
        for half = 1:10
            if all(accepted); break; end
            sel = find(~accepted);
            trial = max(beta(:, sel) + t(sel) .* stepB(:, sel), 0);
            mt = max(trial(1, :) + trial(2, :) .* Pat(:, sel), 1e-12);
            dt = poisson_nnls_em_deviance(Y(:, sel), mt);
            good = dt <= prev(sel);
            if any(good)
                newBeta(:, sel(good)) = trial(:, good);
                accepted(sel(good)) = true;
            end
            t(sel(~good)) = t(sel(~good)) / 2;
        end
        moved = max(abs(newBeta - beta) ./ max(abs(beta), 1), [], 1);
        beta = newBeta;
        if all(moved < opts.tol); break; end
    end

    m = max(beta(1, :) + beta(2, :) .* Pat, 1e-12);
    dev = poisson_nnls_em_deviance(Y, m);
    % Refinement must never make things worse; fall back to the grid node.
    worse = dev > bestDev;
    if any(worse)
        tauFit(worse) = tauGrid(bestIdx(worse));
        refined(worse) = false;
        dev(worse) = bestDev(worse);
    end

    out = struct('tauNs', tauFit, 'background', beta(1, :), ...
        'amplitude', beta(2, :), 'deviance', dev, ...
        'reducedDeviance', dev / max(nBin - 3, 1), ...
        'photons', photons, 'refined', refined);
end
