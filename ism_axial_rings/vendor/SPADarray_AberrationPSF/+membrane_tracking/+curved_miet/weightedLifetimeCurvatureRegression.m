function regression = weightedLifetimeCurvatureRegression( ...
        localizationTable, radialSquaredCorrected, opts)
    import membrane_tracking.curved_miet.*

    varianceZ = localizationTable.heightSigmaUm.^2;
    valid = isfinite(localizationTable.heightUm) & ...
        isfinite(radialSquaredCorrected) & ...
        isfinite(varianceZ) & varianceZ > 0;

    regression = struct('tipHeightUm', NaN, 'curvaturePerUm', NaN, ...
        'tipHeightSigmaUm', NaN, 'curvatureSigmaPerUm', NaN, ...
        'covariance', nan(2), 'information', nan(2), ...
        'nLocalizations', sum(valid), 'isValid', false);
    if sum(valid) < 3
        return;
    end

    response = localizationTable.heightUm(valid);
    design = [ones(sum(valid),1), ...
        0.5 * radialSquaredCorrected(valid)];
    weights = 1 ./ varianceZ(valid);
    medianWeight = median(weights);
    if medianWeight > 0
        weights = min(weights, 100 * medianWeight);
    end
    information = design.' * bsxfun(@times, design, weights);
    [covariance, ok] = invertPositiveDefinite(information);
    if ~ok
        return;
    end
    beta = covariance * (design.' * (weights .* response));
    beta(1) = clampInside(beta(1), opts.tipHeightFitBoundsUm);
    beta(2) = clampInside(beta(2), opts.curvatureFitBoundsPerUm);

    regression.tipHeightUm = beta(1);
    regression.curvaturePerUm = beta(2);
    regression.tipHeightSigmaUm = sqrt(max(covariance(1,1),0));
    regression.curvatureSigmaPerUm = sqrt(max(covariance(2,2),0));
    regression.covariance = covariance;
    regression.information = information;
    regression.isValid = true;
end
