function probability = surfaceGateAcceptanceProbability(D, ...
        noiseCovariance, stepDt, driftTerm, gateRadialLimit)
    import membrane_tracking.fluctuating_miet.*

%SURFACEGATEACCEPTANCEPROBABILITY Gaussian mass inside each gate ellipse.
%
% The gate is centred on zero displacement, whereas curvature gives the
% diffusion increment a small Ito mean D*driftTerm. Angular midpoint
% quadrature and an analytic radial integral retain that non-zero mean.
    if isempty(gateRadialLimit)
        probability = ones(size(stepDt));
        return;
    end
    stepDt = stepDt(:);
    if ismatrix(noiseCovariance)
        noiseCovariance = reshape(noiseCovariance,2,2,1);
    end
    if isvector(driftTerm)
        driftTerm = reshape(driftTerm,1,2);
    end
    if isvector(gateRadialLimit)
        gateRadialLimit = reshape(gateRadialLimit,1,[]);
    end
    nAngles = size(gateRadialLimit,2);
    angles = ((0:nAngles-1)+0.5)*(2*pi/nAngles);
    direction = [cos(angles); sin(angles)];
    probability = zeros(numel(stepDt),1);
    for stepIndex = 1:numel(stepDt)
        covariance = noiseCovariance(:,:,stepIndex) + ...
            2*D*stepDt(stepIndex)*eye(2);
        [inverseCovariance, valid] = invertPositiveDefinite(covariance);
        determinant = det(covariance);
        if ~valid || determinant <= 0
            probability(stepIndex) = realmin;
            continue;
        end
        meanStep = D*driftTerm(stepIndex,:).';
        a = sum(direction.*(inverseCovariance*direction),1);
        b = direction.'*(inverseCovariance*meanStep);
        b = b.';
        c = meanStep.'*inverseCovariance*meanStep;
        radius = gateRadialLimit(stepIndex,:);
        rootTwoA = sqrt(2*a);
        erfDifference = erf((a.*radius-b)./rootTwoA) + ...
            erf(b./rootTwoA);
        noncentralRemainder = max(c-b.^2./a,0);
        gaussianIntegral = exp(-0.5*noncentralRemainder) .* ...
            sqrt(pi./(2*a)) .* erfDifference;
        endpointExponent = max(a.*radius.^2-2*b.*radius+c,0);
        radialIntegral = (b./a).*gaussianIntegral + ...
            (exp(-0.5*c)-exp(-0.5*endpointExponent))./a;
        probability(stepIndex) = mean(max(radialIntegral,0)) / ...
            sqrt(determinant);
    end
    probability = min(max(probability,realmin),1);
end
