function surfaceModel = fillEmptySurfaceModel(surfaceModel)
    import membrane_tracking.fluctuating_miet.*

    surfaceModel.method = 'not estimable';
    surfaceModel.tipHeightUm = NaN;
    surfaceModel.curvaturePerUm = NaN;
    surfaceModel.curvatureSigmaPerUm = NaN;
    surfaceModel.tipHeightSigmaUm = NaN;
    surfaceModel.fieldRmsUm = NaN;
    surfaceModel.fieldCorrelationLengthUm = NaN;
    surfaceModel.fieldCorrelationTimeS = NaN;
    surfaceModel.apexRadiusUm = NaN;
    surfaceModel.fitSucceeded = false;
    surfaceModel.naive = packSurfaceFit([NaN;NaN], nan(2), 0);
    surfaceModel.moment = packSurfaceFit([NaN;NaN], nan(2), 0);
    surfaceModel.gp = struct('isValid', false);
    surfaceModel.identifiability = struct();
    surfaceModel.topography = struct('isValid', false);
end
