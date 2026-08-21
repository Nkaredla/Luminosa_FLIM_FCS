function index = sampleDiscreteLogWeights(logWeights)
%SAMPLEDISCRETELOGWEIGHTS Draw one index from unnormalized log probabilities.

    logWeights = double(logWeights(:).');
    maximum = max(logWeights);
    weights = exp(logWeights - maximum);
    cumulative = cumsum(weights);
    if isempty(cumulative) || ~(cumulative(end) > 0) || ...
            ~isfinite(cumulative(end))
        index = randi(max(numel(logWeights), 1));
        return;
    end
    threshold = rand * cumulative(end);
    index = find(cumulative >= threshold, 1, 'first');
end
