function detector = makeDetector(opts)
    import membrane_tracking.fluctuating_miet.*

    [detXY, indexGrid, layoutInfo] = detectorLayout( ...
        opts.detectorLayout, opts.detectorPitchUm);
    model = struct();
    model.detSize = opts.detectorFillRatio * opts.detectorPitchUm;
    model.detectorPixelShape = 'hex';
    if strcmpi(opts.detectorLayout, 'square25')
        model.detectorPixelShape = 'square';
    end
    model.detectorHexRadius = model.detSize / sqrt(3);
    [du, dv] = detectorPixelOffsets(model, opts.detectorSubsamples);

    detector = struct();
    detector.layout = opts.detectorLayout;
    detector.layoutInfo = layoutInfo;
    detector.detXY = detXY;
    detector.indexGrid = indexGrid;
    detector.nChannels = size(detXY, 1);
    detector.detectorPitchUm = opts.detectorPitchUm;
    detector.detectorHexRadiusUm = model.detectorHexRadius;
    detector.quadratureOffsetXUm = du(:);
    detector.quadratureOffsetYUm = dv(:);
    detector.backgroundMeanCounts = ...
        (opts.backgroundPhotonsPerDetectorS + ...
         opts.darkCountsPerDetectorS) * opts.dtS;
    detector.centerCaptureNormalization = ...
        sum(ismDetectorRawResponse([0 0], detector, opts));
    if ~(detector.centerCaptureNormalization > 0)
        error('simulateFluctuatingMIETMembraneTracking:BadDetectorModel', ...
            'The modelled ISM detector captures no light at the focus.');
    end
end
