function value = spotNegativeLogLikelihood(theta, data, detector, rows, cols, ...
        opts, xCenter0, yCenter0, fitLimitUm)
    import membrane_tracking.focused_ism.*

    if abs(theta(1) - xCenter0) > fitLimitUm || ...
            abs(theta(2) - yCenter0) > fitLimitUm
        value = 1e20 + 1e12 * ((theta(1) - xCenter0)^2 + ...
            (theta(2) - yCenter0)^2);
        return;
    end

    A = expClamped(theta(3));
    b = expClamped(theta(4));
    prob = gaussianPixelModel(theta(1), theta(2), opts.psfSigmaUm, ...
        detector.xEdges(cols(1):(cols(end)+1)), ...
        detector.yEdges(rows(1):(rows(end)+1)));
    mu = max(b + A * prob, opts.minExpectedCount);
    value = sum(mu(:) - data(:) .* log(mu(:)));
    if ~isfinite(value)
        value = realmax('double') / 10;
    end
end
