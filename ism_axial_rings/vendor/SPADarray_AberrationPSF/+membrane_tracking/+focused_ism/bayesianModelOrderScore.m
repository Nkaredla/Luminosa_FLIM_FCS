function score = bayesianModelOrderScore(nll, parameterCount, nPixels, ...
        nEmitters, priorMean)
    import membrane_tracking.focused_ism.*

    % BIC is the large-sample Laplace evidence approximation; the Poisson
    % count prior supplies concentration information without fixing K.
    logPrior = -priorMean + nEmitters * log(priorMean) - ...
        gammaln(nEmitters + 1);
    score = 2 * nll + parameterCount * log(max(nPixels, 2)) - ...
        2 * logPrior;
end
