function deviance = poissonDeviance(data, mu)
    import membrane_tracking.fluctuating_miet.*

    data = double(data(:));
    mu = max(double(mu(:)), realmin);
    terms = mu - data;
    positive = data > 0;
    terms(positive) = terms(positive) + ...
        data(positive).*log(data(positive)./mu(positive));
    deviance = 2*sum(terms);
end
