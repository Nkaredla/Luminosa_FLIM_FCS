function [loc, ok] = fitGaussianSpotMLE(image, row0, col0, detector, opts)
    import membrane_tracking.focused_ism.*

    ok = false;
    loc = struct();

    radius = opts.fitWindowRadiusPx;
    rows = (row0-radius):(row0+radius);
    cols = (col0-radius):(col0+radius);
    rows = rows(rows >= 1 & rows <= detector.ny);
    cols = cols(cols >= 1 & cols <= detector.nx);

    data = double(image(rows, cols));
    bg0 = max(detector.backgroundMeanCounts, 1e-4);
    signal = max(data - bg0, 0);
    signalInWindow = sum(signal(:));
    if signalInWindow < opts.minPhotonsPerLocalization
        return;
    end

    [xGrid, yGrid] = meshgrid(detector.xCenters(cols), ...
        detector.yCenters(rows));
    if signalInWindow > 0
        xInit = sum(signal(:) .* xGrid(:)) / signalInWindow;
        yInit = sum(signal(:) .* yGrid(:)) / signalInWindow;
    else
        xInit = detector.xCenters(col0);
        yInit = detector.yCenters(row0);
    end

    xInit = min(max(xInit, detector.xEdges(cols(1))), ...
        detector.xEdges(cols(end)+1));
    yInit = min(max(yInit, detector.yEdges(rows(1))), ...
        detector.yEdges(rows(end)+1));

    [probInit, ~, ~] = gaussianPixelModel(xInit, yInit, opts.psfSigmaUm, ...
        detector.xEdges(cols(1):(cols(end)+1)), ...
        detector.yEdges(rows(1):(rows(end)+1)));
    A0 = max(signalInWindow / max(sum(probInit(:)), eps), 1);
    theta0 = [xInit, yInit, log(A0), log(bg0)];

    xCenter0 = detector.xCenters(col0);
    yCenter0 = detector.yCenters(row0);
    fitLimitUm = opts.maxFitOffsetPx * detector.pixelSizeUm;
    nll = @(theta) spotNegativeLogLikelihood(theta, data, detector, ...
        rows, cols, opts, xCenter0, yCenter0, fitLimitUm);

    fminOpts = optimset('Display', 'off', 'MaxIter', 300, ...
        'MaxFunEvals', 1000, 'TolX', 1e-5, 'TolFun', 1e-5);
    [thetaHat, fval, exitFlag] = fminsearch(nll, theta0, fminOpts);

    if exitFlag <= 0 || ~all(isfinite(thetaHat)) || ~isfinite(fval)
        return;
    end

    xHat = thetaHat(1);
    yHat = thetaHat(2);
    if abs(xHat - xCenter0) > fitLimitUm || abs(yHat - yCenter0) > fitLimitUm
        return;
    end

    Ahat = expClamped(thetaHat(3));
    bhat = expClamped(thetaHat(4));
    [prob, dPdx, dPdy] = gaussianPixelModel(xHat, yHat, opts.psfSigmaUm, ...
        detector.xEdges(cols(1):(cols(end)+1)), ...
        detector.yEdges(rows(1):(rows(end)+1)));
    signalPhotons = Ahat * sum(prob(:));
    if signalPhotons < opts.minPhotonsPerLocalization
        return;
    end

    fisher = localizationFisher(Ahat, bhat, prob, dPdx, dPdy, opts);
    if ~fisher.isValid
        return;
    end

    loc.xUm = xHat;
    loc.yUm = yHat;
    loc.xPixel = (xHat - detector.xEdges(1)) / detector.pixelSizeUm + 0.5;
    loc.yPixel = (yHat - detector.yEdges(1)) / detector.pixelSizeUm + 0.5;
    loc.signalPhotons = signalPhotons;
    loc.windowCounts = sum(data(:));
    loc.backgroundCountsPerPixel = bhat;
    loc.crbVarXUm2 = fisher.covarianceXY(1,1);
    loc.crbVarYUm2 = fisher.covarianceXY(2,2);
    loc.crbCovXYUm2 = fisher.covarianceXY(1,2);
    loc.crbSigmaUm = sqrt(max(0.5 * trace(fisher.covarianceXY), 0));
    loc.negLogLikelihood = fval;
    ok = true;
end
