function density = truncatedExponentialDensity(microtimeNs, tauNs, periodNs)
    import membrane_tracking.curved_miet.*

    normalization = -expm1(-periodNs ./ tauNs);
    density = exp(-microtimeNs ./ tauNs) ./ (tauNs .* normalization);
end
