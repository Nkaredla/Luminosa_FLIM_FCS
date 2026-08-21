function latent = boundedToLatent(value, bounds)
    import membrane_tracking.curved_miet.*

    fraction = (value - bounds(1)) / diff(bounds);
    fraction = min(max(fraction, 1e-8), 1-1e-8);
    latent = log(fraction / (1-fraction));
end
