function [raw, dRawDx, dRawDy] = ismDetectorRawResponse(position, ...
        detector, opts)
    import membrane_tracking.fluctuating_miet.*

% The ISM microimage is inverted relative to the sample coordinate.
    x = position(1) + detector.detXY(:,1) + detector.quadratureOffsetXUm.';
    y = position(2) + detector.detXY(:,2) + detector.quadratureOffsetYUm.';
    sigmaSquared = opts.psfSigmaUm^2;
    samples = exp(-0.5 * (x.^2 + y.^2) / sigmaSquared);
    raw = mean(samples, 2);
    if nargout > 1
        dRawDx = mean(-x .* samples / sigmaSquared, 2);
        dRawDy = mean(-y .* samples / sigmaSquared, 2);
    end
end
