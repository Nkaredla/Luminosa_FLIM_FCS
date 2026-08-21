function curvature = emptyCurvatureResult(opts)
    import membrane_tracking.curved_miet.*

    curvature = struct();
    curvature.model = ...
        'z=tipHeightUm+0.5*curvaturePerUm*(x^2+y^2)';
    curvature.method = 'not estimable';
    curvature.tipHeightUm = NaN;
    curvature.curvaturePerUm = NaN;
    curvature.apexRadiusUm = NaN;
    curvature.covariance = nan(2);
    curvature.observedFisherInformation = nan(2);
    curvature.tipHeightSigmaUm = NaN;
    curvature.curvatureSigmaPerUm = NaN;
    curvature.negativeLogLikelihood = NaN;
    curvature.nLocalizations = 0;
    curvature.nPhotons = 0;
    curvature.fitSucceeded = false;
    curvature.weightedLifetimeRegression = struct();
    curvature.trueTipHeightUm = opts.tipHeightUm;
    curvature.trueCurvaturePerUm = opts.curvaturePerUm;
    curvature.positionErrorCorrection = ...
        'r^2 = xhat^2+yhat^2-tr(Cxy)';
    curvature.uncertaintyConditioning = ...
        ['Observed photon Fisher information conditional on fitted lateral ' ...
        'positions and the linear MIET calibration'];
end
