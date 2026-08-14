function deviance = poissonDeviance(observed, expected)
%POISSONDEVIANCE Calculate twice the Poisson log-likelihood ratio.

    observed = max(double(observed(:)), 0);
    expected = max(double(expected(:)), eps);
    term = expected - observed;
    positive = observed > 0;
    term(positive) = term(positive) + observed(positive) .* ...
        log(observed(positive) ./ expected(positive));
    deviance = 2 * sum(term);
end
