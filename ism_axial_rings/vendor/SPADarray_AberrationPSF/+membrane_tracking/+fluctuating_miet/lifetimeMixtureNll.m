function value = lifetimeMixtureNll(logTau, microtimes, channelWeight, ...
        signalFraction, nChannels, opts)
    import membrane_tracking.fluctuating_miet.*

    tau = exp(logTau);
    bounded = min(max(tau, opts.lifetimeBoundsNs(1)), opts.lifetimeBoundsNs(2));
    density = truncatedExponentialDensity(microtimes, bounded, ...
        opts.repetitionPeriodNs);
    mixture = signalFraction * channelWeight .* density + ...
        (1-signalFraction) / (nChannels * opts.repetitionPeriodNs);
    value = -sum(log(max(mixture, realmin)));
    if tau < opts.lifetimeBoundsNs(1)
        value = value + 1e8*(log(tau/opts.lifetimeBoundsNs(1)))^2;
    elseif tau > opts.lifetimeBoundsNs(2)
        value = value + 1e8*(log(tau/opts.lifetimeBoundsNs(2)))^2;
    end
    if ~isfinite(value)
        value = realmax('double')/10;
    end
end
