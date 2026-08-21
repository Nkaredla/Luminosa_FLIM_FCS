function samples = sampleTruncatedExponential(tauNs, periodNs, nSamples)
    import membrane_tracking.curved_miet.*

    normalization = -expm1(-periodNs / tauNs);
    samples = -tauNs * log(1 - rand(nSamples,1) * normalization);
end
