function seedXY = focusedISMResidualSeed(data, mu, existingXY, detector, opts)
    import membrane_tracking.focused_ism.*

    residual = max(data - mu, 0);
    if sum(residual) < opts.minPhotonsPerLocalization
        seedXY = [];
        return;
    end

    gridStep = min(detector.detectorPitchUm / 2, opts.psfSigmaUm / 2);
    nGridHalf = ceil(opts.maxLocalizationRadiusUm / gridStep);
    axisValues = (-nGridHalf:nGridHalf) * ...
        (opts.maxLocalizationRadiusUm / nGridHalf);
    [gridX, gridY] = meshgrid(axisValues, axisValues);
    candidates = [gridX(:), gridY(:)];
    inside = sum(candidates.^2, 2) <= opts.maxLocalizationRadiusUm^2;
    candidates = candidates(inside,:);

    minSeparationUm = opts.seedExclusionRadiusPx * ...
        detector.detectorPitchUm;
    if ~isempty(existingXY)
        keep = true(size(candidates,1), 1);
        for e = 1:size(existingXY,1)
            keep = keep & hypot(candidates(:,1) - existingXY(e,1), ...
                candidates(:,2) - existingXY(e,2)) >= minSeparationUm;
        end
        candidates = candidates(keep,:);
    end
    if isempty(candidates)
        seedXY = [];
        return;
    end

    % Matched-filter the positive channel residual against every candidate
    % detector pattern. The subsequent Poisson fit refines this grid seed.
    scores = -inf(size(candidates,1), 1);
    for c = 1:size(candidates,1)
        p = ismDetectorChannelProbability(candidates(c,:), detector, opts);
        scores(c) = (residual.' * p) / sqrt(sum(p.^2));
    end
    [~, best] = max(scores);
    seedXY = candidates(best,:);
end
