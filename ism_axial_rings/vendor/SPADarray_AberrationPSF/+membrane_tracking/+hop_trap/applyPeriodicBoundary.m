function positionUm = applyPeriodicBoundary(positionUm, boxSizeUm)
%APPLYPERIODICBOUNDARY Wrap positions into a box centered at the origin.

    positionUm = mod(positionUm + 0.5 * boxSizeUm, boxSizeUm) - ...
        0.5 * boxSizeUm;
end
