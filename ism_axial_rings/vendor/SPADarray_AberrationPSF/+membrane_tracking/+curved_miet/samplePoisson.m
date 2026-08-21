function counts = samplePoisson(mu)
    import membrane_tracking.curved_miet.*

    mu = max(double(mu), 0);
    if exist('poissrnd', 'file') == 2
        counts = poissrnd(mu);
        return;
    end

    outputSize = size(mu);
    lambda = mu(:);
    flatCounts = zeros(size(lambda));
    small = lambda < 30;
    if any(small)
        flatCounts(small) = samplePoissonByInversion(lambda(small));
    end
    largeIndices = find(~small);
    for index = largeIndices.'
        flatCounts(index) = samplePoissonPtrs(lambda(index));
    end
    counts = reshape(flatCounts, outputSize);
end
