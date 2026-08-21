function [model, ok] = fitFocusedISMModel(data, seedXY, detector, ...
        opts, previousModel)
    import membrane_tracking.focused_ism.*

    nEmitters = size(seedXY, 1);
    theta0 = zeros(1, 3*nEmitters + 1);
    background0 = max(detector.backgroundMeanCounts, 1e-4);
    residualSignal = max(sum(data) - ...
        detector.nChannels * background0, nEmitters);

    if ~isempty(previousModel) && ...
            previousModel.nEmitters == nEmitters - 1
        for e = 1:(nEmitters-1)
            base = 3 * (e - 1);
            theta0(base+(1:3)) = [previousModel.xUm(e), ...
                previousModel.yUm(e), log(previousModel.amplitude(e))];
        end
        background0 = previousModel.background;
    end
    newBase = 3 * (nEmitters - 1);
    theta0(newBase+(1:3)) = [seedXY(end,:), ...
        log(max(residualSignal / nEmitters, 1))];
    if isempty(previousModel)
        for e = 1:(nEmitters-1)
            base = 3 * (e - 1);
            theta0(base+(1:3)) = [seedXY(e,:), ...
                log(max(residualSignal / nEmitters, 1))];
        end
    end
    theta0(end) = log(background0);

    % Log amplitudes and log background enforce positivity without requiring
    % Optimization Toolbox constraints.
    objective = @(theta) focusedISMNegativeLogLikelihood( ...
        theta, data, detector, opts);
    fminOpts = optimset('Display', 'off', ...
        'MaxIter', 300 + 120*nEmitters, ...
        'MaxFunEvals', 1000 + 600*nEmitters, ...
        'TolX', 1e-6, 'TolFun', 1e-6);
    [thetaHat, fval, exitFlag] = fminsearch(objective, theta0, fminOpts);

    ok = exitFlag > 0 && all(isfinite(thetaHat)) && isfinite(fval);
    model = struct();
    if ~ok
        return;
    end

    xUm = thetaHat(1:3:(3*nEmitters)).';
    yUm = thetaHat(2:3:(3*nEmitters)).';
    amplitude = expClamped(thetaHat(3:3:(3*nEmitters))).';
    background = expClamped(thetaHat(end));
    probability = zeros(detector.nChannels, nEmitters);
    mu = background * ones(detector.nChannels, 1);
    for e = 1:nEmitters
        probability(:,e) = ismDetectorChannelProbability( ...
            [xUm(e) yUm(e)], detector, opts);
        mu = mu + amplitude(e) * probability(:,e);
    end

    minSeparationUm = opts.minResolvedEmitterSeparationPx * ...
        detector.detectorPitchUm;
    minPairSeparationUm = Inf;
    for a = 1:nEmitters
        for b = (a+1):nEmitters
            minPairSeparationUm = min(minPairSeparationUm, hypot( ...
                xUm(a)-xUm(b), yUm(a)-yUm(b)));
        end
    end
    insideFocus = hypot(xUm, yUm) <= opts.maxLocalizationRadiusUm;
    ok = all(amplitude >= opts.minPhotonsPerLocalization) && ...
        all(insideFocus);
    if strcmp(opts.emitterSeparationMode, 'hard')
        ok = ok && minPairSeparationUm >= minSeparationUm;
    end
    if ~ok
        return;
    end

    model.nEmitters = nEmitters;
    model.xUm = xUm;
    model.yUm = yUm;
    model.amplitude = amplitude;
    model.background = background;
    model.probability = probability;
    model.mu = mu;
    model.negLogLikelihood = fval;
    model.minPairSeparationUm = minPairSeparationUm;
end
