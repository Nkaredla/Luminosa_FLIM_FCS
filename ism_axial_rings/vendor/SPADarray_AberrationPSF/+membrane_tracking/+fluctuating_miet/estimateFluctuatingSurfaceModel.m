function surfaceModel = estimateFluctuatingSurfaceModel( ...
        localizationTable, trajectories, modes, opts)
    import membrane_tracking.fluctuating_miet.*

%ESTIMATEFLUCTUATINGSURFACEMODEL Mean paraboloid under a fluctuating field.
%
%   Three estimators are reported, in increasing order of honesty:
%
%   1. naive        weighted least squares of zhat on 0.5*rhat^2, using only
%                   the per-frame MIET height variance. Point estimate is
%                   attenuated by localisation error and error bars ignore
%                   the field entirely.
%   2. moment       regression-dilution correction. Subtracting trace(C)
%                   from rhat^2 makes the REGRESSOR unbiased but leaves the
%                   SLOPE attenuated by lambda = 1 - Vbar/Suu, where Vbar is
%                   the mean regressor noise variance. Dividing by lambda
%                   removes it. Monte Carlo at 35 nm precision: naive -28%,
%                   moment-corrected -0.3%.
%   3. gp           the fluctuation field is marginalised as a Gaussian
%                   process, so kappa is a generalised-least-squares
%                   estimate under a correlated covariance and the interval
%                   comes from a profile likelihood rather than a Fisher
%                   approximation.
%
%   Reported kappa is the GP estimate when a field kernel is fitted, and
%   the moment-corrected estimate otherwise.

    surfaceModel = struct();
    surfaceModel.trueTipHeightUm = opts.tipHeightUm;
    surfaceModel.trueCurvaturePerUm = opts.curvaturePerUm;
    n = height(localizationTable);
    surfaceModel.nLocalizations = n;
    if n < 8
        surfaceModel = fillEmptySurfaceModel(surfaceModel);
        return;
    end

    % ---- observations -------------------------------------------------
    z = localizationTable.heightUm;
    zVariance = localizationTable.heightSigmaUm.^2;
    qHat = [localizationTable.xUm, localizationTable.yUm];
    timeS = localizationTable.timeS;
    positionVarianceTrace = localizationTable.crbVarXUm2 + ...
        localizationTable.crbVarYUm2;
    positionCovariance = zeros(2, 2, n);
    positionCovariance(1,1,:) = localizationTable.crbVarXUm2;
    positionCovariance(2,2,:) = localizationTable.crbVarYUm2;
    positionCovariance(1,2,:) = localizationTable.crbCovXYUm2;
    positionCovariance(2,1,:) = localizationTable.crbCovXYUm2;

    valid = isfinite(z) & isfinite(zVariance) & zVariance > 0 & ...
        all(isfinite(qHat), 2) & isfinite(positionVarianceTrace);
    z = z(valid); zVariance = zVariance(valid); qHat = qHat(valid,:);
    timeS = timeS(valid); positionVarianceTrace = positionVarianceTrace(valid);
    positionCovariance = positionCovariance(:,:,valid);
    n = numel(z);
    surfaceModel.nUsed = n;
    if n < 8
        surfaceModel = fillEmptySurfaceModel(surfaceModel);
        return;
    end

    % Unbiased estimator of r^2. E[|qhat|^2] = r^2 + trace(C).
    radiusSquaredHat = sum(qHat.^2, 2);
    u = radiusSquaredHat - positionVarianceTrace;
    design = [ones(n,1), 0.5*u];

    % ---- 1. naive weighted least squares -------------------------------
    weights = 1 ./ zVariance;
    [betaNaive, covarianceNaive] = weightedLeastSquares(design, z, weights);
    surfaceModel.naive = packSurfaceFit(betaNaive, covarianceNaive, n);

    % ---- 2. moment (regression-dilution) correction ---------------------
    % Var(|qhat|^2 | q) = 4*sigma^2*r^2 + 4*sigma^4 for isotropic error
    % with per-axis variance sigma^2 = 0.5*trace(C).
    perAxisVariance = 0.5 * positionVarianceTrace;
    regressorNoiseVariance = zeros(n,1);
    for row = 1:n
        covariance = positionCovariance(:,:,row);
        covarianceSquaredTrace = trace(covariance*covariance);
        qCq = qHat(row,:) * covariance * qHat(row,:).' - ...
            covarianceSquaredTrace;
        regressorNoiseVariance(row) = 4*max(qCq,0) + ...
            2*covarianceSquaredTrace;
    end
    weightSum = sum(weights);
    weightedMeanU = sum(weights.*u) / max(weightSum, realmin);
    weightedMeanZ = sum(weights.*z) / max(weightSum, realmin);
    signalVariance = sum(weights.*(u-weightedMeanU).^2) / ...
        max(weightSum, realmin);
    meanRegressorNoiseVariance = sum(weights.*regressorNoiseVariance) / ...
        max(weightSum, realmin);
    attenuation = 1 - meanRegressorNoiseVariance / ...
        max(signalVariance, realmin);
    attenuation = min(max(attenuation, 1e-3), 1);
    meanRegressor = 0.5 * weightedMeanU;
    correctionTransform = [1, (1-1/attenuation)*meanRegressor; ...
        0, 1/attenuation];
    betaMoment = correctionTransform * betaNaive;
    betaMoment(1) = weightedMeanZ - betaMoment(2)*meanRegressor;
    covarianceMoment = correctionTransform * covarianceNaive * ...
        correctionTransform.';
    surfaceModel.moment = packSurfaceFit(betaMoment, covarianceMoment, n);
    surfaceModel.moment.attenuationFactor = attenuation;
    betaForEstimator = betaMoment;
    covarianceForEstimator = covarianceMoment;
    attenuationForEstimator = attenuation;
    if ~opts.momentCorrectCurvature
        betaForEstimator = betaNaive;
        covarianceForEstimator = covarianceNaive;
        attenuationForEstimator = 1;
    end

    % ---- 3. Gaussian-process marginalisation ---------------------------
    surfaceModel.gp = struct('isValid', false);
    surfaceModel.identifiability = struct();
    useGp = ~strcmp(opts.fluctuationKernel, 'none') && n >= 20 && ...
        (~strcmp(opts.fluctuationKernel, 'helfrich') || modes.nModes > 0);
    if useGp
        subset = subsampleForGp(n, opts.maxGpObservations);
        gpData = struct('z', z(subset), 'design', design(subset,:), ...
            'position', qHat(subset,:), 'timeS', timeS(subset), ...
            'measurementVariance', zVariance(subset), ...
            'perAxisVariance', perAxisVariance(subset), ...
            'positionCovariance', positionCovariance(:,:,subset), ...
            'attenuation', attenuationForEstimator, ...
            'regressorMean', meanRegressor, 'modes', modes);
        gpData.nObservations = numel(subset);
        gpFit = fitGaussianProcessSurface(gpData, betaForEstimator, opts);
        surfaceModel.gp = gpFit;
        if gpFit.isValid
            surfaceModel.identifiability = curvatureProfileLikelihood( ...
                gpData, gpFit, opts);
        end
    end

    % ---- reported values ------------------------------------------------
    if isfield(surfaceModel.gp,'isValid') && surfaceModel.gp.isValid
        surfaceModel.method = 'GP-marginalised generalised least squares';
        surfaceModel.tipHeightUm = surfaceModel.gp.tipHeightUm;
        surfaceModel.curvaturePerUm = surfaceModel.gp.curvaturePerUm;
        surfaceModel.curvatureSigmaPerUm = surfaceModel.gp.curvatureSigmaPerUm;
        surfaceModel.tipHeightSigmaUm = surfaceModel.gp.tipHeightSigmaUm;
        surfaceModel.fieldRmsUm = surfaceModel.gp.fieldRmsUm;
        surfaceModel.fieldCorrelationLengthUm = ...
            surfaceModel.gp.correlationLengthUm;
        surfaceModel.fieldCorrelationTimeS = surfaceModel.gp.correlationTimeS;
    else
        if opts.momentCorrectCurvature
            surfaceModel.method = 'moment-corrected weighted least squares';
        else
            surfaceModel.method = 'mean-corrected weighted least squares';
        end
        surfaceModel.tipHeightUm = betaForEstimator(1);
        surfaceModel.curvaturePerUm = betaForEstimator(2);
        surfaceModel.curvatureSigmaPerUm = ...
            sqrt(max(covarianceForEstimator(2,2),0));
        surfaceModel.tipHeightSigmaUm = ...
            sqrt(max(covarianceForEstimator(1,1),0));
        surfaceModel.fieldRmsUm = NaN;
        surfaceModel.fieldCorrelationLengthUm = NaN;
        surfaceModel.fieldCorrelationTimeS = NaN;
    end
    surfaceModel.apexRadiusUm = 1 / max(abs(surfaceModel.curvaturePerUm), realmin);
    surfaceModel.fitSucceeded = isfinite(surfaceModel.curvaturePerUm);
    surfaceModel.observationRadiusUm = sqrt(max(radiusSquaredHat));
    surfaceModel.confoundedWavelengthUm = ...
        2 * surfaceModel.observationRadiusUm;
    surfaceModel.note = ['Fluctuation modes with wavelength near ' ...
        'confoundedWavelengthUm are quadratic across the observed patch ' ...
        'and are not separable from kappa by any estimator.'];
    surfaceModel.topography = reconstructMembraneTopography( ...
        surfaceModel, trajectories, modes, opts);
end
