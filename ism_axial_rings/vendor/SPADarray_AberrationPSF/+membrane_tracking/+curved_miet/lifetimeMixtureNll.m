function value = lifetimeMixtureNll(logTau, microtimes, ...
        signalChannelProbability, signalFraction, nChannels, opts)
    import membrane_tracking.curved_miet.*

    tau = exp(logTau);
    lower = opts.lifetimeBoundsNs(1);
    upper = opts.lifetimeBoundsNs(2);
    tauForDensity = min(max(tau, lower), upper);
    density = truncatedExponentialDensity( ...
        microtimes, tauForDensity, opts.repetitionPeriodNs);
    mixtureDensity = signalFraction * ...
        signalChannelProbability .* density + ...
        (1-signalFraction) / ...
        (nChannels * opts.repetitionPeriodNs);
    value = -sum(log(max(mixtureDensity, realmin)));
    if tau < lower
        value = value + 1e8 * (log(tau/lower))^2;
    elseif tau > upper
        value = value + 1e8 * (log(tau/upper))^2;
    end
    if ~isfinite(value)
        value = realmax('double') / 10;
    end
end
