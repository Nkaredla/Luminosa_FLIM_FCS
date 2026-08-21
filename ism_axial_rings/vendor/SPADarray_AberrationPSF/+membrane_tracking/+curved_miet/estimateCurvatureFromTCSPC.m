function [curvature, localizationTable] = estimateCurvatureFromTCSPC( ...
        localizationTable, photonEvents, photophysics, detector, opts)
    import membrane_tracking.curved_miet.*

    nLocalizations = height(localizationTable);
    if nLocalizations == 0
        curvature = emptyCurvatureResult(opts);
        localizationTable.radialPositionUm = zeros(0,1);
        localizationTable.radialSquaredCorrectedUm2 = zeros(0,1);
        localizationTable.fittedSurfaceHeightUm = zeros(0,1);
        localizationTable.fittedLifetimeNs = zeros(0,1);
        localizationTable.fittedMeanCurvaturePerUm = zeros(0,1);
        return;
    end

    radialSquared = localizationTable.xUm.^2 + localizationTable.yUm.^2;
    localizationVarianceTrace = localizationTable.crbVarXUm2 + ...
        localizationTable.crbVarYUm2;
    % E[xhat^2+yhat^2-tr(C)] = x^2+y^2 for an unbiased localization.
    radialSquaredCorrected = radialSquared - localizationVarianceTrace;

    regression = weightedLifetimeCurvatureRegression( ...
        localizationTable, radialSquaredCorrected, opts);
    photonData = buildCurvaturePhotonData(localizationTable, ...
        radialSquaredCorrected, photonEvents, photophysics, detector, opts);

    usePhotonMLE = nLocalizations >= opts.minCurvatureLocalizations && ...
        photonData.nPhotons > 0;
    if usePhotonMLE
        initial = [regression.tipHeightUm, regression.curvaturePerUm];
        if any(~isfinite(initial))
            finiteHeight = localizationTable.heightUm( ...
                isfinite(localizationTable.heightUm));
            if isempty(finiteHeight)
                initialHeight = mean(opts.tipHeightFitBoundsUm);
            else
                initialHeight = median(finiteHeight);
            end
            initial = [initialHeight, 0];
        end
        initial(1) = clampInside(initial(1), opts.tipHeightFitBoundsUm);
        initial(2) = clampInside(initial(2), ...
            opts.curvatureFitBoundsPerUm);
        theta0 = [boundedToLatent(initial(1), ...
            opts.tipHeightFitBoundsUm), boundedToLatent(initial(2), ...
            opts.curvatureFitBoundsPerUm)];
        objective = @(theta) curvatureLatentNll(theta, photonData, opts);
        fitOptions = optimset('Display', 'off', 'MaxIter', 350, ...
            'MaxFunEvals', 1200, 'TolX', 1e-7, 'TolFun', 1e-7);
        [thetaHat, nll, exitFlag] = fminsearch( ...
            objective, theta0, fitOptions);
        estimate = [latentToBounded(thetaHat(1), ...
            opts.tipHeightFitBoundsUm), latentToBounded(thetaHat(2), ...
            opts.curvatureFitBoundsPerUm)];
        fitSucceeded = exitFlag > 0 && all(isfinite(estimate)) && ...
            isfinite(nll);
        if ~fitSucceeded && regression.isValid
            estimate = [regression.tipHeightUm, ...
                regression.curvaturePerUm];
            nll = NaN;
            usePhotonMLE = false;
            fitSucceeded = true;
        end
    else
        estimate = [regression.tipHeightUm, regression.curvaturePerUm];
        nll = NaN;
        fitSucceeded = regression.isValid;
    end

    covariance = nan(2);
    information = nan(2);
    covarianceValid = false;
    if fitSucceeded && usePhotonMLE
        physicalObjective = @(beta) curvaturePhysicalNll( ...
            beta, photonData, opts);
        [information, informationValid] = numericalHessian2( ...
            physicalObjective, estimate, ...
            [opts.tipHeightFitBoundsUm; opts.curvatureFitBoundsPerUm]);
        if informationValid
            [covariance, covarianceValid] = ...
                invertPositiveDefinite(information);
        end
    elseif fitSucceeded && regression.isValid
        covariance = regression.covariance;
        information = regression.information;
        covarianceValid = all(isfinite(covariance(:)));
    end

    curvature = struct();
    curvature.model = ...
        'z=tipHeightUm+0.5*curvaturePerUm*(x^2+y^2)';
    curvature.method = 'global channel/microtime mixture maximum likelihood';
    if ~usePhotonMLE
        curvature.method = 'weighted per-frame lifetime regression';
    end
    curvature.tipHeightUm = estimate(1);
    curvature.curvaturePerUm = estimate(2);
    curvature.apexRadiusUm = 1 / abs(estimate(2));
    if estimate(2) == 0
        curvature.apexRadiusUm = Inf;
    end
    curvature.covariance = covariance;
    curvature.observedFisherInformation = information;
    curvature.tipHeightSigmaUm = NaN;
    curvature.curvatureSigmaPerUm = NaN;
    if covarianceValid
        curvature.tipHeightSigmaUm = sqrt(max(covariance(1,1),0));
        curvature.curvatureSigmaPerUm = sqrt(max(covariance(2,2),0));
    end
    curvature.negativeLogLikelihood = nll;
    curvature.nLocalizations = nLocalizations;
    curvature.nPhotons = photonData.nPhotons;
    curvature.fitSucceeded = fitSucceeded;
    curvature.weightedLifetimeRegression = regression;
    curvature.trueTipHeightUm = opts.tipHeightUm;
    curvature.trueCurvaturePerUm = opts.curvaturePerUm;
    curvature.positionErrorCorrection = ...
        'r^2 = xhat^2+yhat^2-tr(Cxy)';
    curvature.uncertaintyConditioning = ...
        ['Observed photon Fisher information conditional on fitted lateral ' ...
        'positions and the linear MIET calibration'];

    localizationTable.radialPositionUm = sqrt(radialSquared);
    localizationTable.radialSquaredCorrectedUm2 = ...
        radialSquaredCorrected;
    localizationTable.fittedSurfaceHeightUm = estimate(1) + ...
        0.5 * estimate(2) * radialSquaredCorrected;
    localizationTable.fittedLifetimeNs = opts.lifetimeAtSubstrateNs + ...
        opts.lifetimeSlopeNsPerUm * ...
        localizationTable.fittedSurfaceHeightUm;
    s = 1 + estimate(2)^2 * max(radialSquaredCorrected, 0);
    radialCurvature = estimate(2) ./ s.^(3/2);
    azimuthalCurvature = estimate(2) ./ sqrt(s);
    localizationTable.fittedMeanCurvaturePerUm = ...
        0.5 * (radialCurvature + azimuthalCurvature);
end
