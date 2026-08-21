function fit = packSurfaceFit(beta, covariance, n)
    import membrane_tracking.fluctuating_miet.*

    fit = struct();
    fit.tipHeightUm = beta(1);
    fit.curvaturePerUm = beta(2);
    fit.tipHeightSigmaUm = sqrt(max(covariance(1,1), 0));
    fit.curvatureSigmaPerUm = sqrt(max(covariance(2,2), 0));
    fit.covariance = covariance;
    fit.nObservations = n;
end
