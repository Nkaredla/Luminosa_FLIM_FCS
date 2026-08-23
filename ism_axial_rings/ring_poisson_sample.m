function counts = ring_poisson_sample(lambda)
%RING_POISSON_SAMPLE Poisson counts, exact where it matters.
%
% counts = ring_poisson_sample(lambda)
%
% Uses the Statistics Toolbox poissrnd when it is licensed and available,
% otherwise falls back to inversion (Knuth product) below a mean of 100 and a
% rounded normal approximation above it.
%
% Why the fallback threshold is 100 and not something smaller: every statistic
% in this folder is a Poisson DEVIANCE, computed under an exact-Poisson
% assumption, and thresholds for those statistics are calibrated by simulation.
% If the sampler were the coarser of the two, the calibration would absorb the
% sampler's error and the reported false-positive rates would describe the
% simulator rather than the instrument. At a mean of 100 the normal
% approximation's skewness error is about 0.1 of a count, which is far below
% the Poisson standard deviation of 10.

    persistent hasPoissrnd
    if isempty(hasPoissrnd)
        hasPoissrnd = exist('poissrnd', 'file') > 0 && ...
            license('test', 'Statistics_Toolbox');
    end
    if hasPoissrnd
        counts = poissrnd(lambda);
        return;
    end

    counts = zeros(size(lambda));
    large = lambda > 100;
    if any(large(:))
        draw = lambda(large) + sqrt(lambda(large)) .* randn(nnz(large), 1);
        counts(large) = max(0, round(draw));
    end
    small = find(~large & lambda > 0);
    for k = 1:numel(small)
        limit = exp(-lambda(small(k)));
        product = rand();
        n = 0;
        while product > limit
            product = product * rand();
            n = n + 1;
        end
        counts(small(k)) = n;
    end
end
