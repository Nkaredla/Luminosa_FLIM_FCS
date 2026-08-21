function detector = makeDetector(opts)
    import membrane_tracking.curved_miet.*

    [detXY, indexGrid, layoutInfo] = detectorLayout( ...
        opts.detectorLayout, opts.detectorPitchUm);

    detectorModel = struct();
    detectorModel.detSize = opts.detectorFillRatio * opts.detectorPitchUm;
    detectorModel.detectorPixelShape = 'hex';
    if strcmpi(opts.detectorLayout, 'square25')
        detectorModel.detectorPixelShape = 'square';
    end
    detectorModel.detectorHexRadius = detectorModel.detSize / sqrt(3);
    [du, dv] = detectorPixelOffsets(detectorModel, ...
        opts.detectorSubsamples);

    detector = struct();
    detector.layout = opts.detectorLayout;
    detector.layoutInfo = layoutInfo;
    detector.detXY = detXY;
    detector.indexGrid = indexGrid;
    detector.nChannels = size(detXY, 1);
    detector.detectorPitchUm = opts.detectorPitchUm;
    detector.detectorSizeUm = detectorModel.detSize;
    detector.detectorHexRadiusUm = detectorModel.detectorHexRadius;
    detector.quadratureOffsetXUm = du(:);
    detector.quadratureOffsetYUm = dv(:);
    detector.backgroundMeanCounts = ...
        (opts.backgroundPhotonsPerDetectorS + ...
        opts.darkCountsPerDetectorS) * opts.dtS;
    detector.focusCenterUm = [0 0];

    centerRaw = ismDetectorRawResponse([0 0], detector, opts);
    detector.centerCaptureNormalization = sum(centerRaw);
    if ~(detector.centerCaptureNormalization > 0)
        error('simulateMIETCurvedMembraneTracking:BadDetectorModel', ...
            'The modeled ISM detector captures no light at the focus.');
    end
end
