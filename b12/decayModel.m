function model = decayModel(irf, dtNs, lifetimes, includeBackground)
%DECAYMODEL Build constant-background and IRF-convolved exponential columns.
% Each exponential pattern is normalised to unit area, so its fitted
% coefficient is interpretable as a component photon amplitude.

    irf = max(double(irf(:)), 0);
    irf = irf ./ max(sum(irf), eps);
    lifetimes = max(double(lifetimes(:)), eps);
    binCount = numel(irf);
    t0 = (0:binCount-1)' * dtNs;
    t1 = (1:binCount)' * dtNs;
    model = zeros(binCount, numel(lifetimes) + double(includeBackground));
    if includeBackground
        model(:, 1) = 1;
    end

    for component = 1:numel(lifetimes)
        tau = lifetimes(component);
        integratedDecay = tau * (exp(-t0 / tau) - exp(-t1 / tau));
        convolved = conv(irf, integratedDecay, 'full');
        pattern = max(convolved(1:binCount), 0);
        model(:, double(includeBackground) + component) = ...
            pattern ./ max(sum(pattern), eps);
    end
end
