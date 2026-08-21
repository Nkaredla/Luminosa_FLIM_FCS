function candidates = findSpotCandidates(image, detector, opts)
    import membrane_tracking.focused_ism.*

    if isempty(opts.peakThresholdCounts)
        bg = detector.backgroundMeanCounts;
        threshold = bg + opts.detectionThresholdSigma * sqrt(max(bg, 1));
        threshold = max(threshold, bg + 4);
    else
        threshold = opts.peakThresholdCounts;
    end

    isMax = image > threshold;
    [ny, nx] = size(image);

    for dy = -1:1
        for dx = -1:1
            if dy == 0 && dx == 0
                continue;
            end
            shifted = -inf(ny, nx);
            srcRows = max(1, 1 - dy):min(ny, ny - dy);
            srcCols = max(1, 1 - dx):min(nx, nx - dx);
            dstRows = srcRows + dy;
            dstCols = srcCols + dx;
            shifted(dstRows, dstCols) = image(srcRows, srcCols);
            isMax = isMax & image >= shifted;
        end
    end

    edge = opts.fitWindowRadiusPx;
    isMax(1:min(edge, ny), :) = false;
    isMax(max(1, ny-edge+1):ny, :) = false;
    isMax(:, 1:min(edge, nx)) = false;
    isMax(:, max(1, nx-edge+1):nx) = false;

    [rows, cols] = find(isMax);
    if isempty(rows)
        candidates = zeros(0, 3);
        return;
    end

    values = image(isMax);
    [values, order] = sort(values, 'descend');
    rows = rows(order);
    cols = cols(order);

    accepted = false(numel(rows), 1);
    keptRows = zeros(0, 1);
    keptCols = zeros(0, 1);
    minSep = opts.minPeakSeparationPx;

    for k = 1:numel(rows)
        if isempty(keptRows)
            keep = true;
        else
            keep = all(hypot(rows(k) - keptRows, cols(k) - keptCols) >= minSep);
        end
        if keep
            accepted(k) = true;
            keptRows(end+1, 1) = rows(k); %#ok<AGROW>
            keptCols(end+1, 1) = cols(k); %#ok<AGROW>
        end
    end

    candidates = [rows(accepted), cols(accepted), values(accepted)];
end
