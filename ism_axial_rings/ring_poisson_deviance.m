function value = ring_poisson_deviance(y, model)
%RING_POISSON_DEVIANCE Poisson deviance between counts and a model.
%
% value = ring_poisson_deviance(y, model)
%
% 2 * sum(model - y + y*log(y/model)), with the y*log(y) term skipped where
% y == 0 (its limit is zero). This is twice the log-likelihood ratio against
% the saturated model, so differences in it between nested models are the
% likelihood-ratio statistic. It is the correct loss for photon counts and
% stays finite in near-empty bins, which matters here because splitting the
% same photons across five rings empties the bins fivefold.

    model = max(model(:), 1e-12);
    y = y(:);
    good = y > 0;
    value = 2 * (sum(model - y) + sum(y(good) .* log(y(good) ./ model(good))));
end
