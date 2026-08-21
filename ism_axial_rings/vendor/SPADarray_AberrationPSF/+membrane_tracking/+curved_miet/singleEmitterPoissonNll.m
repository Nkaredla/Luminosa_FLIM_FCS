function value = singleEmitterPoissonNll(theta, data, detector, opts)
    import membrane_tracking.curved_miet.*

    x = theta(1);
    y = theta(2);
    amplitude = expClamped(theta(3));
    background = expClamped(theta(4));
    probability = ismDetectorChannelProbability([x y], detector, opts);
    mu = max(amplitude * probability + background, ...
        opts.minExpectedCount);
    value = sum(mu - data .* log(mu) + gammaln(data + 1));

    radius = hypot(x, y);
    if radius > opts.maxLocalizationRadiusUm
        value = value + 1e10 * ...
            (radius - opts.maxLocalizationRadiusUm)^2;
    end
    if ~isfinite(value)
        value = realmax('double') / 10;
    end
end
