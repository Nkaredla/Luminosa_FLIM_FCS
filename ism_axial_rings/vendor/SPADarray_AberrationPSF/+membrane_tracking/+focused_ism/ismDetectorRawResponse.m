function [raw, dRawDx, dRawDy] = ismDetectorRawResponse( ...
        positionUm, detector, opts)
    import membrane_tracking.focused_ism.*

    % The detector image is inverted: a molecule at +x shifts the emission
    % pattern toward detector channels at -x. Averaging over quadrature
    % offsets integrates the Gaussian detection PSF over each finite cell.
    x = positionUm(1) + detector.detXY(:,1) + ...
        detector.quadratureOffsetXUm.';
    y = positionUm(2) + detector.detXY(:,2) + ...
        detector.quadratureOffsetYUm.';
    sigma2 = opts.psfSigmaUm^2;
    sampleResponse = exp(-0.5 * (x.^2 + y.^2) / sigma2);
    raw = mean(sampleResponse, 2);
    if nargout > 1
        dRawDx = mean(-x .* sampleResponse / sigma2, 2);
        dRawDy = mean(-y .* sampleResponse / sigma2, 2);
    end
end
