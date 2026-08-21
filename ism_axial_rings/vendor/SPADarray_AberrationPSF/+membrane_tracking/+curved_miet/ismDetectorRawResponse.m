function [raw, dRawDx, dRawDy] = ismDetectorRawResponse( ...
        positionXY, detector, opts)
    import membrane_tracking.curved_miet.*

    % The ISM detector image is inverted relative to the sample coordinate.
    x = positionXY(1) + detector.detXY(:,1) + ...
        detector.quadratureOffsetXUm.';
    y = positionXY(2) + detector.detXY(:,2) + ...
        detector.quadratureOffsetYUm.';
    sigmaSquared = opts.psfSigmaUm^2;
    samples = exp(-0.5 * (x.^2 + y.^2) / sigmaSquared);
    raw = mean(samples, 2);
    if nargout > 1
        dRawDx = mean(-x .* samples / sigmaSquared, 2);
        dRawDy = mean(-y .* samples / sigmaSquared, 2);
    end
end
