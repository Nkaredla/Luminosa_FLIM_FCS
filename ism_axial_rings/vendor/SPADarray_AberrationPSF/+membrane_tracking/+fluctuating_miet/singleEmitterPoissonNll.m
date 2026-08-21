function value = singleEmitterPoissonNll(theta, data, detector, opts)
    import membrane_tracking.fluctuating_miet.*

    amplitude = expClamped(theta(3));
    background = expClamped(theta(4));
    probability = ismDetectorChannelProbability(theta(1:2), detector, opts);
    mu = max(amplitude*probability + background, opts.minExpectedCount);
    value = sum(mu - data.*log(mu) + gammaln(data+1));
    radius = hypot(theta(1), theta(2));
    if radius > opts.maxLocalizationRadiusUm
        value = value + 1e10*(radius - opts.maxLocalizationRadiusUm)^2;
    end
    if ~isfinite(value)
        value = realmax('double')/10;
    end
end
