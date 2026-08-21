function z = surfaceHeight(positionsXY, tipHeight, curvature)
    import membrane_tracking.curved_miet.*

    if ismatrix(positionsXY)
        radialSquared = sum(positionsXY.^2, 2);
    else
        radialSquared = positionsXY(:,1,:).^2 + positionsXY(:,2,:).^2;
        radialSquared = reshape(radialSquared, ...
            size(positionsXY,1), size(positionsXY,3));
    end
    z = tipHeight + 0.5 * curvature * radialSquared;
end
