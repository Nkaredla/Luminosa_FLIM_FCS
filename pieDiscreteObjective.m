function [dev, beta, model, design] = pieDiscreteObjective(x, decay, basis, ...
        nPulse, K)
%PIEDISCRETEOBJECTIVE Penalty-free deviance of a discrete two-pulse PIE fit.
%
% Separate file rather than a nested function because this project keeps one
% function per file. X holds log lifetimes, K per pulse, pulse-major.

    taus = reshape(exp(x), K, nPulse)';
    if any(~isfinite(taus(:))) || any(taus(:) <= 0) || ...
            any(taus(:) > basis.periodNs)
        dev = 1e12; beta = []; model = []; design = [];
        return;
    end
    cols = cell(1, nPulse);
    for p = 1:nPulse
        cols{p} = pie_pattern_columns(basis, taus(p, :), p);
    end
    design = [ones(basis.nBin, 1), cols{:}];
    [beta, dev] = poisson_nnls_whitened(design, decay, ...
        struct('maxIter', 60, 'tol', 1e-11));
    model = max(design * beta, 1e-12);
end
