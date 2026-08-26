function value = poisson_nnls_em_deviance(y, model)
%POISSON_NNLS_EM_DEVIANCE Column-wise Poisson deviance for a matrix of decays.
%
% value = poisson_nnls_em_deviance(y, model)
%
% Y and MODEL are nBin-by-nPixel. VALUE is 1-by-nPixel:
%
%     2 * sum( model - y + y .* log(y ./ model) )
%
% The y.*log(y) term is skipped where y == 0, since its limit there is zero.
% Written vectorised over columns rather than calling biexp_slb_deviance in a
% loop, because the whole point of the EM solver is that it runs on thousands of
% pixels at once; a per-column loop here would dominate its cost.

    model = max(model, 1e-12);
    good = y > 0;
    terms = model - y;
    % log only where it is needed, so empty bins cost nothing and cannot warn.
    terms(good) = terms(good) + y(good) .* log(y(good) ./ model(good));
    value = 2 * sum(terms, 1);
end
