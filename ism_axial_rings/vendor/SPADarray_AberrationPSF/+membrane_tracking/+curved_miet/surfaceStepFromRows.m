function [stepTangent, stepNoise, projectedStep, valid] = ...
        surfaceStepFromRows(trackTable, firstRow, secondRow, curvature)
    import membrane_tracking.curved_miet.*

    firstPosition = [trackTable.xUm(firstRow), trackTable.yUm(firstRow)];
    secondPosition = [trackTable.xUm(secondRow), trackTable.yUm(secondRow)];
    projectedStep = secondPosition - firstPosition;
    midpoint = 0.5 * (firstPosition + secondPosition);
    [transform, ~] = surfaceGeometryAt(midpoint, curvature);
    stepTangent = (transform * projectedStep.').';

    covariance = localizationCovariance(trackTable, firstRow) + ...
        localizationCovariance(trackTable, secondRow);
    covariance = 0.5 * (covariance + covariance.');
    stepNoise = transform * covariance * transform.';
    [~, valid] = invertPositiveDefinite(stepNoise);
    valid = valid && all(isfinite(stepTangent));
end
