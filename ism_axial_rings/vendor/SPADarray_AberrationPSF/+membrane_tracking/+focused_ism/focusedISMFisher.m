function fisher = focusedISMFisher(model, detector, opts)
    import membrane_tracking.focused_ism.*

    nEmitters = model.nEmitters;
    nParameters = 3*nEmitters + 1;
    % Column order is all emitter positions, all emitter amplitudes, then
    % shared background. Inverting the full matrix marginalizes every
    % nuisance amplitude and neighboring emitter when extracting each 2-D
    % position covariance.
    J = zeros(detector.nChannels, nParameters);
    mu = model.background * ones(detector.nChannels, 1);
    for e = 1:nEmitters
        [probability, dPdx, dPdy] = ismDetectorChannelProbability( ...
            [model.xUm(e) model.yUm(e)], detector, opts);
        mu = mu + model.amplitude(e) * probability;
        positionBase = 2 * (e - 1);
        J(:,positionBase+1) = model.amplitude(e) * dPdx;
        J(:,positionBase+2) = model.amplitude(e) * dPdy;
        J(:,2*nEmitters+e) = probability;
    end
    J(:,end) = 1;
    mu = max(mu, opts.minExpectedCount);
    information = J.' * bsxfun(@rdivide, J, mu);
    information = 0.5 * (information + information.');
    [covariance, ok] = invertPositiveDefinite(information);

    covarianceXY = nan(2, 2, nEmitters);
    if ok
        for e = 1:nEmitters
            indices = 2*(e-1) + (1:2);
            covarianceXY(:,:,e) = covariance(indices,indices);
        end
    end
    fisher = struct('informationFull', information, ...
        'covarianceFull', covariance, 'covarianceXY', covarianceXY, ...
        'isValid', ok);
end
