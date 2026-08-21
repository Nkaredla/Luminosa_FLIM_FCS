function [tauHat, tauSigma, ok] = fitFrameLifetime(microtimes, channels, ...
        model, detector, opts)
    import membrane_tracking.fluctuating_miet.*

    microtimes = double(microtimes(:));
    signalFraction = model.amplitude / ...
        (model.amplitude + detector.nChannels*model.background);
    signalFraction = min(max(signalFraction, 1e-6), 1-1e-6);
    channelWeight = model.probability(double(channels(:)));
    initialTau = min(max(mean(microtimes), 1.1*opts.lifetimeBoundsNs(1)), ...
        0.9*opts.lifetimeBoundsNs(2));

    objective = @(theta) lifetimeMixtureNll(theta, microtimes, ...
        channelWeight, signalFraction, detector.nChannels, opts);
    fitOptions = optimset('Display','off','MaxIter',180, ...
        'MaxFunEvals',500,'TolX',1e-8,'TolFun',1e-8);
    [thetaHat, nll, exitFlag] = fminsearch(objective, log(initialTau), ...
        fitOptions);
    tauHat = exp(thetaHat);
    tauSigma = NaN;
    ok = exitFlag > 0 && isfinite(tauHat) && isfinite(nll) && ...
        tauHat > opts.lifetimeBoundsNs(1) && tauHat < opts.lifetimeBoundsNs(2);
    if ~ok
        return;
    end
    h = 1e-3;
    observedInformation = (objective(thetaHat+h) - 2*nll + ...
        objective(thetaHat-h)) / h^2;
    if observedInformation > 0 && isfinite(observedInformation)
        tauSigma = tauHat / sqrt(observedInformation);
    else
        ok = false;
    end
end
