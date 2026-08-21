function gpFit = fitGaussianProcessSurface(gpData, betaInitial, opts)
    import membrane_tracking.fluctuating_miet.*

    gpFit = struct('isValid', false);
    observationSpan = max(hypot(gpData.position(:,1), gpData.position(:,2)));
    timeSpan = max(gpData.timeS) - min(gpData.timeS);
    residualInitial = gpData.z - gpData.design*betaInitial;
    amplitude0 = max(std(residualInitial), 1e-4);
    if strcmp(opts.fluctuationKernel, 'helfrich')
        logParameters0 = [log(max(amplitude0, ...
            gpData.modes.fieldRmsUm/3)), 0, 0];
    else
        logParameters0 = [log(amplitude0), ...
            log(max(0.5*observationSpan, 1e-3)), ...
            log(max(0.05*timeSpan, 10*eps))];
    end

    noiseCurvature = betaInitial(2);
    if opts.estimateFieldParameters
        fitOptions = optimset('Display','off','MaxIter',300, ...
            'MaxFunEvals',900,'TolX',1e-6,'TolFun',1e-6);
        if strcmp(opts.fluctuationKernel, 'helfrich')
            objective = @(logAmplitude) firstOutput(@() ...
                gpProfileNegLogLikelihood( ...
                [logAmplitude, logParameters0(2:3)], gpData, opts, [], ...
                noiseCurvature));
            [logAmplitude, ~, exitFlag] = fminsearch(objective, ...
                logParameters0(1), fitOptions);
            logParameters = [logAmplitude, logParameters0(2:3)];
        else
            objective = @(theta) firstOutput(@() ...
                gpProfileNegLogLikelihood(theta, gpData, opts, [], ...
                noiseCurvature));
            [logParameters, ~, exitFlag] = fminsearch(objective, ...
                logParameters0, fitOptions);
        end
        if exitFlag <= 0 || any(~isfinite(logParameters))
            logParameters = logParameters0;
        end
    else
        logParameters = logParameters0;
    end

    nIterations = 0;
    beta = [NaN;NaN];
    betaCovariance = nan(2);
    nll = NaN;
    correctionTransform = [1, ...
        (1-1/gpData.attenuation)*gpData.regressorMean; ...
        0, 1/gpData.attenuation];
    for iteration = 1:opts.nGpNigpIterations
        [nll, beta, betaCovariance, ok] = gpProfileNegLogLikelihood( ...
            logParameters, gpData, opts, [], noiseCurvature);
        if ~ok
            return;
        end
        betaCorrected = correctionTransform*beta;
        nIterations = iteration;
        newNoiseCurvature = betaCorrected(2);
        relativeChange = abs(newNoiseCurvature-noiseCurvature) / ...
            max(abs(newNoiseCurvature), 1e-6);
        noiseCurvature = newNoiseCurvature;
        if relativeChange < 1e-3
            break;
        end
    end

    betaCorrected = correctionTransform*beta;
    correctedCovariance = correctionTransform*betaCovariance* ...
        correctionTransform.';
    K = fluctuationKernelMatrix(gpData, logParameters, opts);
    noiseDiagonal = effectiveNoiseDiagonal(gpData, noiseCurvature);
    S = K+diag(noiseDiagonal);
    jitter = max(1e-12,1e-10*mean(diag(S)));
    [L, flag] = chol(S+jitter*eye(size(S)), 'lower');
    if flag ~= 0
        return;
    end
    residual = gpData.z-gpData.design*beta;
    alpha = L.' \ (L \ residual);
    [correlationLength, correlationTime] = effectiveGpScales( ...
        logParameters, gpData.modes, opts);

    gpFit.isValid = true;
    gpFit.logParameters = logParameters;
    gpFit.tipHeightUm = betaCorrected(1);
    gpFit.curvaturePerUm = betaCorrected(2);
    gpFit.tipHeightSigmaUm = sqrt(max(correctedCovariance(1,1),0));
    gpFit.curvatureSigmaPerUm = sqrt(max(correctedCovariance(2,2),0));
    gpFit.fieldRmsUm = exp(logParameters(1));
    gpFit.correlationLengthUm = correlationLength;
    gpFit.correlationTimeS = correlationTime;
    gpFit.negativeLogLikelihood = nll;
    gpFit.nObservations = gpData.nObservations;
    gpFit.attenuationFactor = gpData.attenuation;
    gpFit.nNigpIterations = nIterations;
    gpFit.betaMeasurementSpace = beta;
    gpFit.betaCorrected = betaCorrected;
    gpFit.correctedCovariance = correctedCovariance;
    gpFit.trainingPositionUm = gpData.position;
    gpFit.trainingTimeS = gpData.timeS;
    gpFit.trainingAlpha = alpha;
    gpFit.trainingCholesky = L;
    gpFit.trainingNoiseVarianceUm2 = noiseDiagonal;
end
