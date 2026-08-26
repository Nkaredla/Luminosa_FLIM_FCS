function [beta, deviance, model] = poisson_nnls_em_batch(Y, P1, P2, opts, seed)
%POISSON_NNLS_EM_BATCH Poisson EM amplitudes when EVERY pixel has its own design.
%
% [beta, deviance, model] = poisson_nnls_em_batch(Y, P1, P2, opts, seed)
%
% Y, P1, P2 are all nBin-by-P: column k is pixel k's decay and its two unit-sum
% patterns, each at that pixel's OWN lifetimes. BETA is 3-by-P as [B; a1; a2],
% DEVIANCE is 1-by-P.
%
% This exists because variable projection gives every pixel a different
% (tau1, tau2), so poisson_nnls_em - which takes one shared design - does not
% apply. The grid pipeline could share a design across a whole block; a per-pixel
% optimizer cannot, and that was the reason to think a per-pixel architecture
% could never vectorise as well as the grid.
%
% It vectorises anyway, and better, because of the unit-sum normalisation. The EM
% normaliser is sum over bins of each design column, which here is
%
%     [ sum(ones) ; sum(P1) ; sum(P2) ] = [ nBin ; 1 ; 1 ]
%
% for EVERY pixel, since the patterns are normalised to unit sum by
% construction. So the update needs no matrix product at all - just elementwise
% products and column reductions:
%
%     B  <- B  .* sum(ratio) / nBin
%     a1 <- a1 .* sum(P1 .* ratio)
%     a2 <- a2 .* sum(P2 .* ratio)
%
% with ratio = Y ./ model. All of it is nBin-by-P elementwise work.
%
% Convergence is tested PER COLUMN, and converged columns stop being updated.
% A shared all(...) flag would run every column until the slowest one finished
% and then truncate that slowest column at the cap - measured as a 3.2x waste
% that also damaged exactly the pixels needing the most iterations.
%
% SEED, when given, is a 3-by-P warm start from the previous outer iterate.
% That is safe: EM is globally convergent on a non-negative design, so the seed
% changes the path and not the fixed point. It roughly halves the iteration
% count.

    if nargin < 4 || isempty(opts); opts = struct(); end
    defaults = struct('maxIter', 400, 'tol', 1e-12, 'checkEvery', 10);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end

    [nBin, P] = size(Y);
    if nargin < 5 || isempty(seed)
        total = max(sum(Y, 1), eps);
        beta = [total / (3 * nBin); total / 3; total / 3];
    else
        beta = max(seed, 1e-12);
    end
    beta = max(beta, 1e-12);

    previous = inf(1, P);
    active = true(1, P);
    for iter = 1:opts.maxIter
        if ~any(active); break; end
        idx = active;
        b = beta(:, idx);
        model = max(b(1, :) + b(2, :) .* P1(:, idx) + ...
            b(3, :) .* P2(:, idx), 1e-12);
        ratio = Y(:, idx) ./ model;
        beta(:, idx) = max([ ...
            b(1, :) .* (sum(ratio, 1) / nBin); ...
            b(2, :) .* sum(P1(:, idx) .* ratio, 1); ...
            b(3, :) .* sum(P2(:, idx) .* ratio, 1)], 0);

        if mod(iter, opts.checkEvery) == 0 || iter == opts.maxIter
            cols = find(active);
            m = max(beta(1, cols) + beta(2, cols) .* P1(:, cols) + ...
                beta(3, cols) .* P2(:, cols), 1e-12);
            current = poisson_nnls_em_deviance(Y(:, cols), m);
            moved = (previous(cols) - current) ./ max(abs(previous(cols)), 1);
            previous(cols) = current;
            active(cols(moved < opts.tol)) = false;
        end
    end

    model = max(beta(1, :) + beta(2, :) .* P1 + beta(3, :) .* P2, 1e-12);
    deviance = poisson_nnls_em_deviance(Y, model);
end
