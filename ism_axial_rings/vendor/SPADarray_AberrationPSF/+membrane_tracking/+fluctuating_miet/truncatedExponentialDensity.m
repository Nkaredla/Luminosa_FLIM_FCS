function density = truncatedExponentialDensity(microtimeNs, tauNs, periodNs)
    import membrane_tracking.fluctuating_miet.*

    normalization = -expm1(-periodNs ./ tauNs);
    density = exp(-microtimeNs ./ tauNs) ./ (tauNs .* normalization);
end
