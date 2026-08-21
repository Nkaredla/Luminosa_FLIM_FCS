function relativeCapture = ismArrayCaptureRelative(positionsXY, detector, opts)
    import membrane_tracking.curved_miet.*

    relativeCapture = zeros(size(positionsXY,1), 1);
    for molecule = 1:size(positionsXY,1)
        raw = ismDetectorRawResponse(positionsXY(molecule,:), detector, opts);
        relativeCapture(molecule) = sum(raw) / ...
            detector.centerCaptureNormalization;
    end
    relativeCapture = max(relativeCapture, 0);
end
