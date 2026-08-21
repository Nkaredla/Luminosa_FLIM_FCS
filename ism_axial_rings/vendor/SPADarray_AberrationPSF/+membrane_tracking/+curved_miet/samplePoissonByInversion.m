function counts = samplePoissonByInversion(lambda)
    import membrane_tracking.curved_miet.*

    lambda = lambda(:);
    uniformDraw = rand(size(lambda));
    probability = exp(-lambda);
    cumulative = probability;
    counts = zeros(size(lambda));
    active = uniformDraw > cumulative;
    while any(active)
        counts(active) = counts(active) + 1;
        probability(active) = probability(active) .* lambda(active) ./ ...
            counts(active);
        cumulative(active) = cumulative(active) + probability(active);
        active = uniformDraw > cumulative & probability > 0;
    end
end
