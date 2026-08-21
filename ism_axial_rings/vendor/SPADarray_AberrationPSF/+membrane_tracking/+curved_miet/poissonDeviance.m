function deviance = poissonDeviance(data, mu)
    import membrane_tracking.curved_miet.*

    data = double(data(:));
    mu = max(double(mu(:)), realmin);
    positive = data > 0;
    terms = mu - data;
    terms(positive) = terms(positive) + ...
        data(positive) .* log(data(positive) ./ mu(positive));
    deviance = 2 * sum(terms);
end
