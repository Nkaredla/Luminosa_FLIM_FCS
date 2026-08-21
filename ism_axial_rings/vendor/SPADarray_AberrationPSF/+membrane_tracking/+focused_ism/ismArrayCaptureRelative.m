function captureRelative = ismArrayCaptureRelative(positionsUm, detector, opts)
    import membrane_tracking.focused_ism.*

    captureRelative = zeros(size(positionsUm, 1), 1);
    for m = 1:size(positionsUm, 1)
        raw = ismDetectorRawResponse(positionsUm(m,:), detector, opts);
        captureRelative(m) = sum(raw) / ...
            detector.centerCaptureNormalization;
    end
    captureRelative = max(captureRelative, 0);
end
