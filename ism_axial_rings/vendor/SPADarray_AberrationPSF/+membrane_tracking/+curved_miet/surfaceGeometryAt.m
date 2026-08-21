function [transform, driftPerDiffusion] = surfaceGeometryAt(position, curvature)
    import membrane_tracking.curved_miet.*

    radius = hypot(position(1), position(2));
    if radius > 10*eps
        radial = position / radius;
    else
        radial = [1 0];
    end
    tangential = [-radial(2), radial(1)];
    metricScale = 1 + curvature^2 * radius^2;

    % Multiplication by this matrix maps projected coordinate increments to
    % orthonormal local surface coordinates. Thus T'*T is the surface metric.
    transform = [sqrt(metricScale)*radial; tangential];
    driftPerDiffusion = -curvature^2 * ...
        (2 + curvature^2 * radius^2) / metricScale^2 * position;
end
