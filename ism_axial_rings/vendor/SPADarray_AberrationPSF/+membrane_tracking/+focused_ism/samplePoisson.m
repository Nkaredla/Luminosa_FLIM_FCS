function counts = samplePoisson(mu)
    import membrane_tracking.focused_ism.*

    mu = max(double(mu), 0);
    if exist('poissrnd', 'file') == 2
        counts = poissrnd(mu);
        return;
    end

    % Exact fallback for MATLAB installations without Statistics Toolbox.
    % Inversion is efficient for low-count detector backgrounds; PTRS is
    % efficient for the larger per-molecule photon means.
    outputSize = size(mu);
    lambda = mu(:);
    flatCounts = zeros(size(lambda));
    smallMask = lambda < 30;
    if any(smallMask)
        flatCounts(smallMask) = samplePoissonByInversion(lambda(smallMask));
    end

    largeIndices = find(~smallMask);
    for k = 1:numel(largeIndices)
        index = largeIndices(k);
        flatCounts(index) = samplePoissonPtrs(lambda(index));
    end
    counts = reshape(flatCounts, outputSize);
end
