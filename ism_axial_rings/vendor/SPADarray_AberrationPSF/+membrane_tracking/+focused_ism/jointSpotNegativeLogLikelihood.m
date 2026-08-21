function value = jointSpotNegativeLogLikelihood(theta, data, rows, cols, ...
        detector, opts, xBounds, yBounds)
    import membrane_tracking.focused_ism.*

    nEmitters = (numel(theta) - 1) / 3;
    background = expClamped(theta(end));
    mu = background * ones(size(data));
    xy = zeros(nEmitters, 2);

    for e = 1:nEmitters
        base = 3 * (e - 1);
        x = theta(base + 1);
        y = theta(base + 2);
        if x < xBounds(1) || x > xBounds(2) || ...
                y < yBounds(1) || y > yBounds(2)
            value = 1e20 + 1e12 * ( ...
                max(xBounds(1)-x, 0)^2 + max(x-xBounds(2), 0)^2 + ...
                max(yBounds(1)-y, 0)^2 + max(y-yBounds(2), 0)^2);
            return;
        end
        xy(e,:) = [x y];
        amplitude = expClamped(theta(base + 3));
        probability = gaussianPixelModel(x, y, opts.psfSigmaUm, ...
            detector.xEdges(cols(1):(cols(end)+1)), ...
            detector.yEdges(rows(1):(rows(end)+1)));
        mu = mu + amplitude * probability;
    end

    value = poissonImageNll(data, mu, opts);
    minSeparationUm = opts.minResolvedEmitterSeparationPx * ...
        detector.pixelSizeUm;
    for a = 1:nEmitters
        for b = (a+1):nEmitters
            distance = hypot(xy(a,1)-xy(b,1), xy(a,2)-xy(b,2));
            if distance < minSeparationUm
                value = value + 1e6 * ...
                    ((minSeparationUm-distance) / minSeparationUm)^2;
            end
        end
    end
end
