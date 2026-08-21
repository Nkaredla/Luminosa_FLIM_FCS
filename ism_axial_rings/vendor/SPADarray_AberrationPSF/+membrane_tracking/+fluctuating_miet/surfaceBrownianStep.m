function [nextPosition, boundaryHit] = surfaceBrownianStep(position, ...
        modes, amplitudeA, amplitudeB, dt, opts)
    import membrane_tracking.fluctuating_miet.*

    nMolecules = size(position, 1);
    nextPosition = position;
    D = opts.diffusionUm2PerS;
    for molecule = 1:nMolecules
        [~, gradient, hessian] = totalSurfaceAt(position(molecule,:), ...
            modes, amplitudeA, amplitudeB, opts);
        [~, driftPerDiffusion, inverseSqrtMetric] = ...
            surfaceGeometry(gradient, hessian);
        noise = (inverseSqrtMetric * randn(2,1)).';
        nextPosition(molecule,:) = position(molecule,:) + ...
            D * driftPerDiffusion * dt + sqrt(2*D*dt) * noise;
    end
    proposedRadius = hypot(nextPosition(:,1), nextPosition(:,2));
    boundaryHit = proposedRadius > opts.membraneRadiusUm;
    if any(boundaryHit)
        folded = mod(proposedRadius(boundaryHit), 2*opts.membraneRadiusUm);
        outside = folded > opts.membraneRadiusUm;
        folded(outside) = 2*opts.membraneRadiusUm - folded(outside);
        scale = folded ./ max(proposedRadius(boundaryHit), realmin);
        nextPosition(boundaryHit,:) = nextPosition(boundaryHit,:) .* scale;
    end
end
