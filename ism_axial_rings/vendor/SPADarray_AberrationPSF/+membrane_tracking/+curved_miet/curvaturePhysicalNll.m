function value = curvaturePhysicalNll(beta, photonData, opts)
    import membrane_tracking.curved_miet.*

    heightByLocalization = beta(1) + 0.5 * beta(2) * ...
        photonData.radialSquaredCorrectedUm2;
    lifetimeByLocalization = opts.lifetimeAtSubstrateNs + ...
        opts.lifetimeSlopeNsPerUm * heightByLocalization;

    lower = opts.lifetimeBoundsNs(1);
    upper = opts.lifetimeBoundsNs(2);
    invalidLow = max(lower - lifetimeByLocalization, 0);
    invalidHigh = max(lifetimeByLocalization - upper, 0);
    negativeHeight = max(-heightByLocalization, 0);
    penalty = 1e9 * (sum(invalidLow.^2) + sum(invalidHigh.^2) + ...
        sum(negativeHeight.^2));
    lifetimeByLocalization = min(max(lifetimeByLocalization, lower), upper);
    eventLifetime = lifetimeByLocalization( ...
        photonData.localizationIndex);
    density = truncatedExponentialDensity(photonData.microtimeNs, ...
        eventLifetime, opts.repetitionPeriodNs);
    mixture = photonData.signalScale .* density + ...
        photonData.backgroundScale;
    value = -sum(log(max(mixture, realmin))) + penalty;
    if ~isfinite(value)
        value = realmax('double') / 10;
    end
end
