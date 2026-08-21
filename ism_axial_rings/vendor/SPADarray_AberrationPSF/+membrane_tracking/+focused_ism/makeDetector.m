function detector = makeDetector(opts)
    import membrane_tracking.focused_ism.*

    % Reuse the detector geometry used by the rest of the Luminosa code.
    [detXY, indexGrid, layoutInfo] = detectorLayout( ...
        opts.detectorLayout, opts.detectorPitchUm);
    detectorModel = struct();
    detectorModel.detSize = opts.detectorFillRatio * opts.detectorPitchUm;
    detectorModel.detectorPixelShape = 'hex';
    if strcmp(opts.detectorLayout, 'square25')
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
    detector.backgroundMeanCounts = (opts.backgroundPhotonsPerPixelS + ...
        opts.darkCountsPerPixelS) * opts.dtS;
    detector.focusCenterUm = [0 0];
    detector.effectiveExcitationAreaUm2 = pi * opts.laserWaistUm^2 / 2;

    % collectionEfficiency is defined at the center of the excitation
    % focus. Off-axis finite-array losses are expressed relative to this
    % reference and applied to the molecular photon rate.
    [centerRaw, ~, ~] = ismDetectorRawResponse([0 0], detector, opts);
    detector.centerCaptureNormalization = sum(centerRaw);
    if ~(detector.centerCaptureNormalization > 0)
        error('simulateMembraneDiffusionParticleTracking:BadDetectorModel', ...
            'The ISM detector has zero modeled capture at the focus center.');
    end
end
