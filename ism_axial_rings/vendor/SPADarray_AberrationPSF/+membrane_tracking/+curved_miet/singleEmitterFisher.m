function fisher = singleEmitterFisher(model, detector, opts)
    import membrane_tracking.curved_miet.*

    [probability, dPdx, dPdy] = ismDetectorChannelProbability( ...
        [model.xUm model.yUm], detector, opts);
    mu = max(model.amplitude * probability + model.background, ...
        opts.minExpectedCount);
    jacobian = [model.amplitude*dPdx, model.amplitude*dPdy, ...
        probability, ones(detector.nChannels,1)];
    information = jacobian.' * bsxfun(@rdivide, jacobian, mu);
    information = 0.5 * (information + information.');
    [covariance, valid] = invertPositiveDefinite(information);

    covarianceXY = nan(2);
    if valid
        covarianceXY = covariance(1:2,1:2);
    end
    fisher = struct('information', information, ...
        'covariance', covariance, 'covarianceXY', covarianceXY, ...
        'isValid', valid);
end
