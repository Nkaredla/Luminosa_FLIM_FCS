function [height, gradient, hessian] = fieldAt(position, modes, ...
        amplitudeA, amplitudeB)
    import membrane_tracking.fluctuating_miet.*

%FIELDAT Fluctuation height, gradient, and Hessian at one point.
    height = 0;
    gradient = zeros(1, 2);
    hessian = zeros(2, 2);
    if modes.nModes == 0
        return;
    end
    phase = modes.qVectors * position(:);
    cosPhase = cos(phase);
    sinPhase = sin(phase);
    height = sum(amplitudeA .* cosPhase + amplitudeB .* sinPhase);
    % d/dr [ a cos(q.r) + b sin(q.r) ] = q ( -a sin + b cos )
    radialWeight = -amplitudeA .* sinPhase + amplitudeB .* cosPhase;
    gradient = (radialWeight.' * modes.qVectors);
    % second derivative brings down -q q' ( a cos + b sin )
    curvatureWeight = -(amplitudeA .* cosPhase + amplitudeB .* sinPhase);
    hessian(1,1) = sum(curvatureWeight .* modes.qVectors(:,1).^2);
    hessian(2,2) = sum(curvatureWeight .* modes.qVectors(:,2).^2);
    hessian(1,2) = sum(curvatureWeight .* ...
        modes.qVectors(:,1) .* modes.qVectors(:,2));
    hessian(2,1) = hessian(1,2);
end
