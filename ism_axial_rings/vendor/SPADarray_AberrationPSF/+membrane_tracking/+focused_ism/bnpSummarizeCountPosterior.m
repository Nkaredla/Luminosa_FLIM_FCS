function summary = bnpSummarizeCountPosterior(countHistogram)
%BNPSUMMARIZECOUNTPOSTERIOR Summarize discrete emitter-count probabilities.

    rowTotals = sum(countHistogram, 2);
    probability = countHistogram ./ max(rowTotals, 1);
    countValues = 0:size(countHistogram, 2)-1;
    posteriorMean = probability * countValues.';
    [decisionProbability, mapIndex] = max(probability, [], 2);
    mapCount = mapIndex - 1;
    entropyBits = -sum(probability .* log2(max(probability, realmin)), 2);
    cumulative = cumsum(probability, 2);
    lowerCount = zeros(size(rowTotals));
    upperCount = zeros(size(rowTotals));
    for row = 1:size(probability, 1)
        lowerIndex = find(cumulative(row,:) >= 0.025, 1, 'first');
        upperIndex = find(cumulative(row,:) >= 0.975, 1, 'first');
        lowerCount(row) = lowerIndex - 1;
        upperCount(row) = upperIndex - 1;
    end

    summary = struct();
    summary.probability = probability;
    summary.posteriorMean = posteriorMean;
    summary.mapCount = mapCount;
    summary.decisionProbability = decisionProbability;
    summary.entropyBits = entropyBits;
    summary.credibleIntervalCount = [lowerCount upperCount];
end
