function smoothed = bnpKalmanSmoothSeries(measurements, measurementVariance, ...
        observedFrames, diffusionUm2PerS, dtS)
    import membrane_tracking.focused_ism.*

    n = numel(measurements);
    filteredMean = zeros(n, 1);
    filteredVariance = zeros(n, 1);
    predictedMean = zeros(n, 1);
    predictedVariance = zeros(n, 1);

    filteredMean(1) = measurements(1);
    filteredVariance(1) = measurementVariance(1);
    for k = 2:n
        frameGap = observedFrames(k) - observedFrames(k-1);
        processVariance = 2 * diffusionUm2PerS * dtS * frameGap;
        predictedMean(k) = filteredMean(k-1);
        predictedVariance(k) = filteredVariance(k-1) + processVariance;
        gain = predictedVariance(k) / ...
            (predictedVariance(k) + measurementVariance(k));
        filteredMean(k) = predictedMean(k) + ...
            gain * (measurements(k) - predictedMean(k));
        filteredVariance(k) = max((1 - gain) * ...
            predictedVariance(k), 1e-12);
    end

    smoothed = filteredMean;
    smoothedVariance = filteredVariance;
    for k = n-1:-1:1
        frameGap = observedFrames(k+1) - observedFrames(k);
        processVariance = 2 * diffusionUm2PerS * dtS * frameGap;
        nextPredictionVariance = filteredVariance(k) + processVariance;
        gain = filteredVariance(k) / max(nextPredictionVariance, 1e-12);
        smoothed(k) = filteredMean(k) + ...
            gain * (smoothed(k+1) - filteredMean(k));
        smoothedVariance(k) = filteredVariance(k) + gain^2 * ...
            (smoothedVariance(k+1) - nextPredictionVariance);
    end
end
