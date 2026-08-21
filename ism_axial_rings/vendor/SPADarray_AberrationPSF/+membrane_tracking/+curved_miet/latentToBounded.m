function bounded = latentToBounded(latent, bounds)
    import membrane_tracking.curved_miet.*

    if latent >= 0
        sigmoid = 1 / (1 + exp(-min(latent, 50)));
    else
        exponential = exp(max(latent, -50));
        sigmoid = exponential / (1 + exponential);
    end
    bounded = bounds(1) + diff(bounds) * sigmoid;
end
