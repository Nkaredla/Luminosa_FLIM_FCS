function [tauHat, tauSigma, fisherLogTau, nll, ok] = ...
        fitFrameLifetime(microtimes, channels, model, detector, opts)
    import membrane_tracking.curved_miet.*

    microtimes = double(microtimes(:));
    channels = double(channels(:));
    signalFraction = model.amplitude / ...
        (model.amplitude + detector.nChannels * model.background);
    signalFraction = min(max(signalFraction, 1e-6), 1-1e-6);
    signalChannelProbability = model.probability(channels);

    initialTau = min(max(mean(microtimes), opts.lifetimeBoundsNs(1)*1.1), ...
        opts.lifetimeBoundsNs(2)*0.9);
    theta0 = log(initialTau);
    objective = @(theta) lifetimeMixtureNll(theta, microtimes, ...
        signalChannelProbability, signalFraction, detector.nChannels, opts);
    fitOptions = optimset('Display', 'off', 'MaxIter', 180, ...
        'MaxFunEvals', 500, 'TolX', 1e-8, 'TolFun', 1e-8);
    [thetaHat, nll, exitFlag] = fminsearch(objective, theta0, fitOptions);
    tauHat = exp(thetaHat);

    ok = exitFlag > 0 && isfinite(tauHat) && isfinite(nll) && ...
        tauHat > opts.lifetimeBoundsNs(1) && ...
        tauHat < opts.lifetimeBoundsNs(2);
    tauSigma = NaN;
    fisherLogTau = NaN;
    if ~ok
        return;
    end

    h = 1e-3;
    fisherLogTau = (objective(thetaHat+h) - 2*nll + ...
        objective(thetaHat-h)) / h^2;
    if fisherLogTau > 0 && isfinite(fisherLogTau)
        tauSigma = tauHat / sqrt(fisherLogTau);
    else
        ok = false;
    end
end
