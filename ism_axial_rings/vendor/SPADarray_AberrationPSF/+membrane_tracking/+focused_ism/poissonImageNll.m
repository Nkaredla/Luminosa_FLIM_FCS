function value = poissonImageNll(data, mu, opts)
    import membrane_tracking.focused_ism.*

    mu = max(mu, opts.minExpectedCount);
    value = sum(mu(:) - data(:) .* log(mu(:)));
    if ~isfinite(value)
        value = realmax('double') / 10;
    end
end
