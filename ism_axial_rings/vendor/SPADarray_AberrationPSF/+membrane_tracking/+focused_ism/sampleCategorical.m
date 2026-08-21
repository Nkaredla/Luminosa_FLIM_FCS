function categories = sampleCategorical(probability, nSamples)
    import membrane_tracking.focused_ism.*

    probability = max(probability(:), 0);
    probability = probability / sum(probability);
    cumulative = cumsum(probability);
    cumulative(end) = 1;
    uniform = rand(nSamples, 1);
    categories = ones(nSamples, 1);
    % This toolbox-free inverse-CDF sampler returns one physical detector
    % channel for every detected molecular photon.
    for k = 1:(numel(cumulative)-1)
        categories = categories + (uniform > cumulative(k));
    end
end
