function value = bnpPoissonLogLikelihood(data, meanCounts, opts)
    import membrane_tracking.focused_ism.*

    meanCounts = max(meanCounts, opts.minExpectedCount);
    value = sum(data(:) .* log(meanCounts(:)) - meanCounts(:));
end
