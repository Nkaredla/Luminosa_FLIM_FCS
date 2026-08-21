function [height, gradient, hessian] = totalSurfaceAt(position, modes, ...
        amplitudeA, amplitudeB, opts)
    import membrane_tracking.fluctuating_miet.*

%TOTALSURFACEAT Mean paraboloid plus fluctuation field.
    [fieldHeight, fieldGradient, fieldHessian] = fieldAt(position, modes, ...
        amplitudeA, amplitudeB);
    radiusSquared = position(1)^2 + position(2)^2;
    height = opts.tipHeightUm + 0.5 * opts.curvaturePerUm * radiusSquared ...
        + fieldHeight;
    gradient = opts.curvaturePerUm * position(:).' + fieldGradient;
    hessian = opts.curvaturePerUm * eye(2) + fieldHessian;
end
