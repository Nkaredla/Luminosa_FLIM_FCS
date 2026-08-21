function categories = sampleCategorical(probability, nSamples)
    import membrane_tracking.curved_miet.*

    probability = max(probability(:), 0);
    probability = probability / sum(probability);
    cumulative = cumsum(probability);
    cumulative(end) = 1;
    draws = rand(nSamples,1);
    categories = ones(nSamples,1);
    for channel = 1:(numel(cumulative)-1)
        categories = categories + (draws > cumulative(channel));
    end
end
