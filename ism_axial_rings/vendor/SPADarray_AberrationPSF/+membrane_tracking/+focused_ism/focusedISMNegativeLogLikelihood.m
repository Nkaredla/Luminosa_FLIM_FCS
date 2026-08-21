function value = focusedISMNegativeLogLikelihood(theta, data, detector, opts)
    import membrane_tracking.focused_ism.*

    nEmitters = (numel(theta) - 1) / 3;
    background = expClamped(theta(end));
    mu = background * ones(detector.nChannels, 1);
    xy = zeros(nEmitters, 2);
    value = 0;
    for e = 1:nEmitters
        base = 3 * (e - 1);
        xy(e,:) = theta(base+(1:2));
        radius = hypot(xy(e,1), xy(e,2));
        if radius > opts.maxLocalizationRadiusUm
            value = value + 1e12 * ...
                (radius - opts.maxLocalizationRadiusUm)^2;
        end
        amplitude = expClamped(theta(base+3));
        probability = ismDetectorChannelProbability( ...
            xy(e,:), detector, opts);
        mu = mu + amplitude * probability;
    end
    value = value + poissonImageNll(data, mu, opts);

    if ~strcmp(opts.emitterSeparationMode, 'none')
        minSeparationUm = opts.minResolvedEmitterSeparationPx * ...
            detector.detectorPitchUm;
        for a = 1:nEmitters
            for b = (a+1):nEmitters
                distance = hypot(xy(a,1)-xy(b,1), ...
                    xy(a,2)-xy(b,2));
                if distance < minSeparationUm
                    % This mode is retained for sensitivity analysis. The
                    % default 'none' lets the joint Fisher matrix, rather
                    % than a geometric wall, diagnose unresolved positions.
                    value = value + 1e6 * ...
                        ((minSeparationUm-distance) / minSeparationUm)^2;
                end
            end
        end
    end
end
