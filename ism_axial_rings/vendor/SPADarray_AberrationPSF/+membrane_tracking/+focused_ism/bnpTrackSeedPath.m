function path = bnpTrackSeedPath(trackRows, initialD, nFrames, opts)
    import membrane_tracking.focused_ism.*

    [~, order] = sort(trackRows.frame);
    trackRows = trackRows(order,:);
    [observedFrames, uniqueRows] = unique(trackRows.frame, 'stable');
    trackRows = trackRows(uniqueRows,:);
    observedFrames = double(observedFrames(:));
    observed = [trackRows.xUm, trackRows.yUm];
    observed = bnpUnwrapTrack(observed, opts);

    varianceX = bnpSanitizeMeasurementVariance(trackRows.crbVarXUm2);
    varianceY = bnpSanitizeMeasurementVariance(trackRows.crbVarYUm2);
    smoothed = [ ...
        bnpKalmanSmoothSeries(observed(:,1), varianceX, observedFrames, ...
            initialD, opts.dtS), ...
        bnpKalmanSmoothSeries(observed(:,2), varianceY, observedFrames, ...
            initialD, opts.dtS)];

    path = nan(nFrames, 2);
    firstFrame = observedFrames(1);
    lastFrame = observedFrames(end);
    span = (firstFrame:lastFrame).';
    if numel(observedFrames) == 1
        path(span,:) = repmat(smoothed(1,:), numel(span), 1);
    else
        path(span,1) = interp1(observedFrames, smoothed(:,1), ...
            span, 'linear');
        path(span,2) = interp1(observedFrames, smoothed(:,2), ...
            span, 'linear');
    end

    stepSigma = sqrt(2 * initialD * opts.dtS);
    for frame = firstFrame-1:-1:1
        path(frame,:) = path(frame+1,:) + stepSigma * randn(1, 2);
    end
    for frame = lastFrame+1:nFrames
        path(frame,:) = path(frame-1,:) + stepSigma * randn(1, 2);
    end
    path = applyBoundary(path, opts.boxSizeUm, opts.boundaryCondition);
end
