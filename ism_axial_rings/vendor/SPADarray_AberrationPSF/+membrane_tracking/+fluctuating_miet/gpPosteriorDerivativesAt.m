function [gradient, hessian] = gpPosteriorDerivativesAt( ...
        position, timeS, gpFit, modes, opts)
    import membrane_tracking.fluctuating_miet.*

%GPPOSTERIORDERIVATIVESAT Derivatives of the posterior mean fluctuation.
    trainingPosition = gpFit.trainingPositionUm;
    trainingTime = gpFit.trainingTimeS;
    alpha = gpFit.trainingAlpha;
    amplitudeSquared = exp(2*gpFit.logParameters(1));
    gradient = zeros(1,2);
    hessian = zeros(2);
    switch opts.fluctuationKernel
        case 'sqexp'
            lengthScale = max(exp(gpFit.logParameters(2)),realmin);
            delta = position(:).'-trainingPosition;
            dt = abs(timeS-trainingTime);
            kernel = amplitudeSquared*exp( ...
                -sum(delta.^2,2)/(2*lengthScale^2) - ...
                dt/max(exp(gpFit.logParameters(3)),realmin));
            derivative = -delta.*kernel/lengthScale^2;
            gradient = (derivative.'*alpha).';
            hessian(1,1) = sum(alpha.*kernel.* ...
                (delta(:,1).^2/lengthScale^4-1/lengthScale^2));
            hessian(2,2) = sum(alpha.*kernel.* ...
                (delta(:,2).^2/lengthScale^4-1/lengthScale^2));
            hessian(1,2) = sum(alpha.*kernel.* ...
                (delta(:,1).*delta(:,2)/lengthScale^4));
            hessian(2,1) = hessian(1,2);
        case 'helfrich'
            varianceSum = sum(modes.variance);
            if modes.nModes == 0 || varianceSum <= 0
                return;
            end
            normalizedVariance = modes.variance/varianceSum;
            for mode = 1:modes.nModes
                q = modes.qVectors(mode,:);
                phase = (position(:).'-trainingPosition)*q.';
                temporal = exp(-abs(timeS-trainingTime) / ...
                    max(modes.relaxationTimeS(mode),realmin));
                scale = amplitudeSquared*normalizedVariance(mode);
                weightedSine = sum(alpha.*temporal.*sin(phase));
                weightedCosine = sum(alpha.*temporal.*cos(phase));
                gradient = gradient-scale*weightedSine*q;
                hessian = hessian-scale*weightedCosine*(q.'*q);
            end
    end
end
