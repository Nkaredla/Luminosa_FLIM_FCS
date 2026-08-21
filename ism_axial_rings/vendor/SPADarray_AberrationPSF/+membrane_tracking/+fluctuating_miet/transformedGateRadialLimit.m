function radialLimit = transformedGateRadialLimit( ...
        transform, innovationCovariance, gateChi2, opts)
    import membrane_tracking.fluctuating_miet.*

%TRANSFORMEDGATERADIALLIMIT Gate ellipse expressed in tangent coordinates.
    nAngles = opts.gateAcceptanceQuadraturePoints;
    angles = ((0:nAngles-1)+0.5)*(2*pi/nAngles);
    radialLimit = zeros(1,nAngles);
    [inverseInnovation, valid] = ...
        invertPositiveDefinite(innovationCovariance);
    if ~valid || rcond(transform) < 1e-12
        return;
    end
    inverseTransform = transform \ eye(2);
    tangentGatePrecision = inverseTransform.'*inverseInnovation* ...
        inverseTransform;
    ux = cos(angles);
    uy = sin(angles);
    directionalPrecision = tangentGatePrecision(1,1)*ux.^2 + ...
        2*tangentGatePrecision(1,2)*ux.*uy + ...
        tangentGatePrecision(2,2)*uy.^2;
    radialLimit = sqrt(gateChi2./max(directionalPrecision,realmin));
end
