function [step, stepNoise, transform, driftPerDiffusion, ok] = ...
        tangentStep(trackTable, firstRow, secondRow, ...
        surfaceModel, modes, opts)
    import membrane_tracking.fluctuating_miet.*

    first = [trackTable.xUm(firstRow), trackTable.yUm(firstRow)];
    second = [trackTable.xUm(secondRow), trackTable.yUm(secondRow)];
    midpointTime = 0.5*(trackTable.timeS(firstRow)+ ...
        trackTable.timeS(secondRow));
    [transform, driftPerDiffusion] = estimatedSurfaceGeometryAt( ...
        0.5*(first+second), midpointTime, surfaceModel, modes, opts);
    step = (transform * (second-first).').';
    covariance = localizationCovariance(trackTable, firstRow) + ...
        localizationCovariance(trackTable, secondRow);
    covariance = opts.crbCorrectionScale*0.5*(covariance+covariance.');
    stepNoise = transform*covariance*transform.';
    [~, ok] = invertPositiveDefinite(stepNoise);
    ok = ok && all(isfinite(step)) && all(isfinite(transform(:))) && ...
        all(isfinite(driftPerDiffusion));
end
