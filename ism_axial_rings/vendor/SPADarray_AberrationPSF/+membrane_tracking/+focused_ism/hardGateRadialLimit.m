function radialLimit = hardGateRadialLimit(euclideanRadius, ...
        innovationCovariance, opts)
    import membrane_tracking.focused_ism.*

    nAngles = opts.gateAcceptanceQuadraturePoints;
    angles = ((0:nAngles-1) + 0.5) * (2*pi/nAngles);
    radialLimit = euclideanRadius * ones(1, nAngles);
    if ~strcmp(opts.trackingMethod, 'jpda')
        return;
    end

    [inverseInnovation, valid] = ...
        invertPositiveDefinite(innovationCovariance);
    if ~valid
        radialLimit(:) = 0;
        return;
    end
    ux = cos(angles);
    uy = sin(angles);
    directionalPrecision = inverseInnovation(1,1)*ux.^2 + ...
        2*inverseInnovation(1,2)*ux.*uy + ...
        inverseInnovation(2,2)*uy.^2;
    mahalanobisRadius = sqrt(opts.associationGateChi2 ./ ...
        max(directionalPrecision, realmin));
    radialLimit = min(radialLimit, mahalanobisRadius);
end
