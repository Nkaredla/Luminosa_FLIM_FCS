function [model, ok] = fitJointGaussianModel(data, rows, cols, seedXY, ...
        detector, opts, previousModel)
    import membrane_tracking.focused_ism.*

    nEmitters = size(seedXY, 1);
    background0 = max(detector.backgroundMeanCounts, 1e-4);
    signalTotal = max(sum(max(data - background0, 0), 'all'), nEmitters);

    theta0 = zeros(1, 3 * nEmitters + 1);
    if ~isempty(previousModel) && ...
            previousModel.nEmitters == nEmitters - 1
        for e = 1:(nEmitters-1)
            base = 3 * (e - 1);
            theta0(base + (1:3)) = [previousModel.xUm(e), ...
                previousModel.yUm(e), log(previousModel.amplitude(e))];
        end
        background0 = previousModel.background;
    end
    newEmitter = nEmitters;
    newBase = 3 * (newEmitter - 1);
    theta0(newBase + (1:3)) = [seedXY(newEmitter,:), ...
        log(max(signalTotal / nEmitters, 1))];
    if isempty(previousModel)
        for e = 1:(nEmitters-1)
            base = 3 * (e - 1);
            theta0(base + (1:3)) = [seedXY(e,:), ...
                log(max(signalTotal / nEmitters, 1))];
        end
    end
    theta0(end) = log(background0);

    xBounds = detector.xEdges([cols(1), cols(end)+1]);
    yBounds = detector.yEdges([rows(1), rows(end)+1]);
    objective = @(theta) jointSpotNegativeLogLikelihood(theta, data, ...
        rows, cols, detector, opts, xBounds, yBounds);
    fminOpts = optimset('Display', 'off', ...
        'MaxIter', 250 + 100 * nEmitters, ...
        'MaxFunEvals', 800 + 500 * nEmitters, ...
        'TolX', 1e-5, 'TolFun', 1e-5);
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
    probability = zeros(numel(rows), numel(cols), nEmitters);
    mu = background * ones(size(data));
    signalPhotons = zeros(nEmitters, 1);
    for e = 1:nEmitters
        probability(:,:,e) = gaussianPixelModel(xUm(e), yUm(e), ...
            opts.psfSigmaUm, ...
            detector.xEdges(cols(1):(cols(end)+1)), ...
            detector.yEdges(rows(1):(rows(end)+1)));
        mu = mu + amplitude(e) * probability(:,:,e);
        signalPhotons(e) = amplitude(e) * sum(probability(:,:,e), 'all');
    end

    minSeparationUm = opts.minResolvedEmitterSeparationPx * ...
        detector.pixelSizeUm;
    separationValid = true;
    for a = 1:nEmitters
        for b = (a+1):nEmitters
            separationValid = separationValid && ...
                hypot(xUm(a)-xUm(b), yUm(a)-yUm(b)) >= minSeparationUm;
        end
    end
    ok = all(signalPhotons >= opts.minPhotonsPerLocalization) && ...
        separationValid;
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
end
