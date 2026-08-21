function [transform, driftPerDiffusion] = estimatedSurfaceGeometryAt( ...
        position, timeS, surfaceModel, modes, opts)
    import membrane_tracking.fluctuating_miet.*

%ESTIMATEDSURFACEGEOMETRYAT Metric from the fitted mean or GP posterior.
    gradient = surfaceModel.curvaturePerUm * position(:).';
    hessian = surfaceModel.curvaturePerUm * eye(2);
    usePosterior = strcmp(opts.diffusionGeometry,'posteriorMean') && ...
        isfield(surfaceModel.gp,'isValid') && surfaceModel.gp.isValid && ...
        isfield(surfaceModel.gp,'trainingAlpha');
    if usePosterior
        [fieldGradient, fieldHessian] = gpPosteriorDerivativesAt( ...
            position, timeS, surfaceModel.gp, modes, opts);
        gradient = gradient+fieldGradient;
        hessian = hessian+fieldHessian;
    end
    [metricScale, driftPerDiffusion] = surfaceGeometry(gradient, hessian);
    transform = tangentTransformFromGradient(gradient,metricScale);
end
