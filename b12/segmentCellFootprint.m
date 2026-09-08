function [mask, boundary, detection, coreThreshold, core, growthThreshold, ...
        backgroundLevel, backgroundNoise] = segmentCellFootprint(image, ...
        smoothSigma, thresholdScale, closeRadius, dilateRadius, minimumArea, ...
        growthFraction, backgroundSigma, minimumCoreArea)
%SEGMENTCELLFOOTPRINT B12 intensity cores with connected dim-edge growth.
% Adapts the lower-population background and hysteresis reconstruction in
% immune_cell_MIET_segment. B12 intensity, rather than MIET mean arrival,
% supplies the classification signal. No convex hull or border exclusion.
% Lower growthFraction includes dimmer connected edges; larger values give
% tighter boundaries. Unconnected dim background cannot seed a cell.

    validateattributes(image, {'numeric'}, {'2d', 'real', 'nonempty'});
    validateattributes(smoothSigma, {'numeric'}, {'scalar','finite','nonnegative'});
    validateattributes(thresholdScale, {'numeric'}, {'scalar','finite','positive'});
    validateattributes(growthFraction, {'numeric'}, {'scalar','finite','>',0,'<=',1});
    validateattributes(backgroundSigma, {'numeric'}, {'scalar','finite','nonnegative'});
    validateattributes(minimumCoreArea, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(minimumArea, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(closeRadius, {'numeric'}, {'scalar','finite','nonnegative'});
    validateattributes(dilateRadius, {'numeric'}, {'scalar','finite','nonnegative'});
    valid = isfinite(image) & image >= 0;
    clean = double(image);
    clean(~valid) = 0;
    detection = robustNormalise(log1p(clean));
    if smoothSigma > 0
        weight = imgaussfilt(double(valid), smoothSigma, 'Padding', 'replicate');
        detection = imgaussfilt(detection, smoothSigma, 'Padding', 'replicate') ...
            ./ max(weight, eps);
    end
    values = sort(detection(valid));
    mask = false(size(image));
    boundary = mask;
    core = mask;
    coreThreshold = inf;
    growthThreshold = inf;
    backgroundLevel = NaN;
    backgroundNoise = NaN;
    if isempty(values) || values(end) - values(1) <= eps
        return;
    end
    otsu = graythresh(values);
    % Use the global lower population: cells may touch any image edge.
    lowCeiling = min(otsu, values(max(1, round(0.45 * numel(values)))));
    low = values(values <= lowCeiling);
    if isempty(low)
        low = values(1:max(1, round(0.45 * numel(values))));
    end
    backgroundLevel = median(low);
    backgroundNoise = 1.4826 * median(abs(low - backgroundLevel));
    noiseFloor = max(backgroundNoise, 0.005 * (values(end) - values(1)));
    coreThreshold = max(otsu * thresholdScale, ...
        backgroundLevel + backgroundSigma * noiseFloor);
    growthThreshold = min(coreThreshold, backgroundLevel + max( ...
        growthFraction * (coreThreshold - backgroundLevel), ...
        0.5 * backgroundSigma * noiseFloor));
    core = valid & detection > coreThreshold;
    core = bwareaopen(core, minimumCoreArea, 8);
    permissive = valid & detection > growthThreshold;
    mask = imreconstruct(core, permissive | core, 8);
    if closeRadius > 0
        mask = imclose(mask, strel('disk', round(closeRadius), 0));
    end
    mask = imfill(mask, 'holes');
    % Reject small islands before dilation can inflate them into cells.
    mask = bwareaopen(mask, minimumArea, 8);
    if dilateRadius > 0
        mask = imdilate(mask, strel('disk', round(dilateRadius), 0));
    end
    mask = imfill(mask, 'holes');
    core = core & mask;
    boundary = bwperim(mask, 8);
end
