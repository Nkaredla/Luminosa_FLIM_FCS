function count = samplePoissonPtrs(lambda)
    import membrane_tracking.curved_miet.*

    sqrtLambda = sqrt(lambda);
    b = 0.931 + 2.53 * sqrtLambda;
    a = -0.059 + 0.02483 * b;
    inverseAlpha = 1.1239 + 1.1328 / (b - 3.4);
    quickAccept = 0.9277 - 3.6224 / (b - 2);
    while true
        u = rand() - 0.5;
        v = rand();
        uSymmetric = 0.5 - abs(u);
        if uSymmetric <= 0
            continue;
        end
        count = floor((2*a/uSymmetric + b)*u + lambda + 0.43);
        if uSymmetric >= 0.07 && v <= quickAccept
            return;
        end
        if count < 0 || (uSymmetric < 0.013 && v > uSymmetric)
            continue;
        end
        logAcceptance = log(v * inverseAlpha / ...
            (a/uSymmetric^2 + b));
        logProbability = -lambda + count*log(lambda) - ...
            gammaln(count+1);
        if logAcceptance <= logProbability
            return;
        end
    end
end
