function [nextPosition, boundaryHit] = surfaceBrownianStep( ...
        position, curvature, diffusion, dt, membraneRadius)
    import membrane_tracking.curved_miet.*

    n = size(position, 1);
    radius = hypot(position(:,1), position(:,2));
    nonzero = radius > 10*eps(max(membraneRadius, 1));

    radial = zeros(n, 2);
    radial(nonzero,:) = position(nonzero,:) ./ radius(nonzero);
    radial(~nonzero,:) = repmat([1 0], sum(~nonzero), 1);
    tangential = [-radial(:,2), radial(:,1)];

    metricScale = 1 + curvature^2 * radius.^2;
    eta = randn(n, 2);
    projectedNoise = radial .* (eta(:,1) ./ sqrt(metricScale)) + ...
        tangential .* eta(:,2);

    % For g=I+kappa^2*q*q', this is
    % (1/sqrt(det(g)))*div(sqrt(det(g))*inv(g)), without the factor D.
    driftPerDiffusion = -curvature^2 * ...
        (2 + curvature^2 * radius.^2) ./ metricScale.^2 .* position;
    nextPosition = position + diffusion * driftPerDiffusion * dt + ...
        sqrt(2 * diffusion * dt) * projectedNoise;

    proposedRadius = hypot(nextPosition(:,1), nextPosition(:,2));
    boundaryHit = proposedRadius > membraneRadius;
    if any(boundaryHit)
        reflectedRadius = mod(proposedRadius(boundaryHit), ...
            2 * membraneRadius);
        outside = reflectedRadius > membraneRadius;
        reflectedRadius(outside) = ...
            2 * membraneRadius - reflectedRadius(outside);
        scale = reflectedRadius ./ max(proposedRadius(boundaryHit), realmin);
        nextPosition(boundaryHit,:) = ...
            nextPosition(boundaryHit,:) .* scale;
    end
end
