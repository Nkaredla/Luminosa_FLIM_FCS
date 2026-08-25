function value = biexp_slb_deviance(y, model)
%BIEXP_SLB_DEVIANCE Poisson deviance, the correct loss for photon counts.
%
% 2 * sum(model - y + y*log(y/model)), skipping the y*log(y) term where y == 0
% since its limit there is zero. Differences in this between nested models are
% the likelihood-ratio statistic, and it stays finite in empty bins - which
% matters because a single pixel's decay has many.

    model = max(model(:), 1e-12);
    y = y(:);
    good = y > 0;
    value = 2 * (sum(model - y) + sum(y(good) .* log(y(good) ./ model(good))));
end
