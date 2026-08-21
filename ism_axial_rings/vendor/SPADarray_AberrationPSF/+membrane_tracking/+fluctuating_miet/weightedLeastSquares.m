function [beta, covariance] = weightedLeastSquares(design, response, weights)
    import membrane_tracking.fluctuating_miet.*

    information = design.' * bsxfun(@times, design, weights(:));
    [covariance, ok] = invertPositiveDefinite(information);
    if ~ok
        beta = [NaN; NaN];
        covariance = nan(2);
        return;
    end
    beta = covariance * (design.' * (weights(:) .* response));
end
