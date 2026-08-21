function counts = samplePoissonByInversion(lambda)
    import membrane_tracking.fluctuating_miet.*

    lambda = lambda(:);
    draw = rand(size(lambda));
    probability = exp(-lambda);
    cumulative = probability;
    counts = zeros(size(lambda));
    active = draw > cumulative;
    while any(active)
        counts(active) = counts(active) + 1;
        probability(active) = probability(active).*lambda(active)./counts(active);
        cumulative(active) = cumulative(active) + probability(active);
        active = draw > cumulative & probability > 0;
    end
end
