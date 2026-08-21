function value = curvatureLatentNll(theta, photonData, opts)
    import membrane_tracking.curved_miet.*

    beta = [latentToBounded(theta(1), opts.tipHeightFitBoundsUm), ...
        latentToBounded(theta(2), opts.curvatureFitBoundsPerUm)];
    value = curvaturePhysicalNll(beta, photonData, opts);
end
