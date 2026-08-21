function value = logSumExp(values, dimension)
%LOGSUMEXP Stable logarithm of a sum of exponentials.

    if nargin < 2
        dimension = 1;
    end
    maximum = max(values, [], dimension);
    shifted = values - maximum;
    value = maximum + log(sum(exp(shifted), dimension));
    allNegativeInfinity = ~isfinite(maximum);
    value(allNegativeInfinity) = maximum(allNegativeInfinity);
end
