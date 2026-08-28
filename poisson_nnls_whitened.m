function [beta, deviance, iterations] = poisson_nnls_whitened(design, y, opts)
%POISSON_NNLS_WHITENED Non-negative Poisson ML fit for ANY number of components.
%
% [beta, deviance, iterations] = poisson_nnls_whitened(design, y, opts)
%
% Minimises the Poisson deviance of y against design*beta subject to beta >= 0,
% for one decay and an arbitrary number of columns. Fisher scoring with weights
% w = 1/m, each step solved as a non-negative least squares problem on the
% WHITENED system sqrt(w).*design against sqrt(w).*y.
%
% THE WHITENING IS THE POINT
%
% Unconstrained, minimising ||sqrt(w).*(A*beta - y)|| and solving the normal
% equations A'WA*beta = A'Wy are the same thing. CONSTRAINED they are not, and
% that difference is a real defect elsewhere in this project: PIRLSnonneg calls
% lsqnonneg(x'*w*x, x'*w*y), handing the GRAM MATRIX to a solver that expects a
% design matrix, which minimises the norm of the normal-equation residual
% instead of the weighted residual. It agrees with the truth while every bound
% is inactive and is wrong by a median of 28 deviance units when one is active -
% on about half the pixels in this data. Passing the whitened design to
% lsqnonneg is the correct formulation.
%
% This is for ONE decay with many columns, so lsqnonneg is used directly; the
% per-pixel fits use the exact enumerated solvers instead, because lsqnonneg
% cannot be batched across pixels.
%
% Steps are damped until the deviance actually falls, which restores the
% monotonicity that Fisher scoring lacks by construction.
%
% opts: maxIter (default 200), tol (default 1e-12), maxHalvings (default 20)

    if nargin < 3 || isempty(opts); opts = struct(); end
    defaults = struct('maxIter', 200, 'tol', 1e-12, 'maxHalvings', 20);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end

    y = double(y(:));
    design = double(design);
    if any(design(:) < 0)
        error('poisson_nnls_whitened:NegativeDesign', ...
            'The design must be non-negative.');
    end
    nCol = size(design, 2);

    % Strictly positive start so the first weight is finite; the solve is free
    % to put any component straight back on the bound.
    beta = repmat(max(sum(y), eps) / (nCol * max(numel(y), 1)), nCol, 1);
    model = max(design * beta, 1e-12);
    deviance = poisson_nnls_em_deviance(y, model);
    iterations = opts.maxIter;

    for iter = 1:opts.maxIter
        sw = 1 ./ sqrt(model);
        target = lsqnonneg(design .* sw, y .* sw);

        step = target - beta;
        accepted = false;
        t = 1;
        for h = 1:opts.maxHalvings
            trial = max(beta + t * step, 0);
            mt = max(design * trial, 1e-12);
            dt = poisson_nnls_em_deviance(y, mt);
            if dt <= deviance
                moved = (deviance - dt) / max(abs(deviance), 1);
                beta = trial;
                model = mt;
                deviance = dt;
                accepted = true;
                break;
            end
            t = t / 2;
        end
        if ~accepted || moved < opts.tol
            iterations = iter;
            break;
        end
    end
end
