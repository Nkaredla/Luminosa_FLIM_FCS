function counts = samplePoisson(mu)
    import membrane_tracking.fluctuating_miet.*

    mu = max(double(mu), 0);
    if exist('poissrnd','file') == 2
        counts = poissrnd(mu);
        return;
    end
    outputSize = size(mu);
    lambda = mu(:);
    flat = zeros(size(lambda));
    small = lambda < 30;
    if any(small)
        flat(small) = samplePoissonByInversion(lambda(small));
    end
    for index = find(~small).'
        flat(index) = samplePoissonPtrs(lambda(index));
    end
    counts = reshape(flat, outputSize);
end
