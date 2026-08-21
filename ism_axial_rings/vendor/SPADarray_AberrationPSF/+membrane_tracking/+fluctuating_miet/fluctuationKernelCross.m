function K = fluctuationKernelCross(queryPosition, queryTime, ...
        trainingPosition, trainingTime, logParameters, opts, modes)
    import membrane_tracking.fluctuating_miet.*

    amplitudeSquared = exp(2*logParameters(1));
    queryTime = queryTime(:);
    trainingTime = trainingTime(:);
    switch opts.fluctuationKernel
        case 'sqexp'
            lengthScale = max(exp(logParameters(2)), realmin);
            timeScale = max(exp(logParameters(3)), realmin);
            dx = queryPosition(:,1) - trainingPosition(:,1).';
            dy = queryPosition(:,2) - trainingPosition(:,2).';
            dt = abs(queryTime - trainingTime.');
            K = amplitudeSquared * exp( ...
                -(dx.^2+dy.^2)/(2*lengthScale^2) - dt/timeScale);
        case 'helfrich'
            K = zeros(size(queryPosition,1), size(trainingPosition,1));
            varianceSum = sum(modes.variance);
            if modes.nModes == 0 || varianceSum <= 0
                return;
            end
            normalizedVariance = modes.variance / varianceSum;
            dt = abs(queryTime - trainingTime.');
            for mode = 1:modes.nModes
                queryPhase = queryPosition * modes.qVectors(mode,:).';
                trainingPhase = trainingPosition * ...
                    modes.qVectors(mode,:).';
                spatial = cos(queryPhase-trainingPhase.');
                temporal = exp(-dt/max( ...
                    modes.relaxationTimeS(mode), realmin));
                K = K + normalizedVariance(mode)*(spatial.*temporal);
            end
            K = amplitudeSquared*K;
        otherwise
            K = zeros(size(queryPosition,1), size(trainingPosition,1));
    end
end
