function [lengthScale, timeScale] = effectiveGpScales( ...
        logParameters, modes, opts)
    import membrane_tracking.fluctuating_miet.*

    if strcmp(opts.fluctuationKernel, 'sqexp')
        lengthScale = exp(logParameters(2));
        timeScale = exp(logParameters(3));
        return;
    end
    varianceSum = sum(modes.variance);
    if modes.nModes == 0 || varianceSum <= 0
        lengthScale = NaN;
        timeScale = NaN;
        return;
    end
    meanQSquared = sum(modes.variance.*modes.qMagnitudePerUm.^2) / ...
        varianceSum;
    lengthScale = sqrt(2/max(meanQSquared, realmin));
    timeScale = sum(modes.variance.*modes.relaxationTimeS) / varianceSum;
end
