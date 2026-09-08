function [cellMask, cellBoundary, mitoMask, lysoMask, cytosolMask, ...
        compartmentMap, mitoLabels, lysoLabels, cellDetectionImage, ...
        cellThreshold, mitoEnhanced, mitoThreshold, lysosomeResponse, ...
        lysosomeThreshold, lysosomeBrightThreshold, cellCore, ...
        cellGrowthThreshold, cellBackgroundLevel, cellBackgroundNoise] = ...
        segmentCompartments(mitoImage, b12Image, cellSmoothSigma, ...
        cellThresholdScale, cellCloseRadius, cellDilateRadius, ...
        minimumCellArea, mitoSmoothSigma, mitoBackgroundRadius, ...
        mitoThresholdScale, mitoCloseRadius, mitoDilateRadius, minimumMitoArea, ...
        lysoSigmaSmall, lysoSigmaLarge, lysoThresholdMad, lysoMinPeakFraction, ...
        lysoBrightQuantile, lysoExcludeMitoRadius, lysoCloseRadius, ...
        lysoDilateRadius, minimumLysoArea, cellGrowthFraction, ...
        cellBackgroundSigma, minimumCellCoreArea)
%SEGMENTCOMPARTMENTS Segment cell, mitochondria, lysosomes, and cytosol.
% Mitochondria use local-background-subtracted 640/Det2 intensity.
% Lysosomes use coarse, merged bright-spot regions in B12.

    if nargin < 23 || isempty(cellGrowthFraction), cellGrowthFraction = 0.5; end
    if nargin < 24 || isempty(cellBackgroundSigma), cellBackgroundSigma = 3; end
    if nargin < 25 || isempty(minimumCellCoreArea), minimumCellCoreArea = 10; end
    [cellMask, cellBoundary, cellDetectionImage, cellThreshold, cellCore, ...
        cellGrowthThreshold, cellBackgroundLevel, cellBackgroundNoise] = ...
        segmentCellFootprint(b12Image, cellSmoothSigma, cellThresholdScale, ...
        cellCloseRadius, cellDilateRadius, minimumCellArea, ...
        cellGrowthFraction, cellBackgroundSigma, minimumCellCoreArea);
    if ~any(cellMask(:))
        error('Cell segmentation is empty. Reduce cellThresholdScale or minimumCellArea.');
    end

    % Remove slowly varying 640 nm background before intensity thresholding.
    mitoSmooth = imgaussfilt(double(mitoImage), mitoSmoothSigma);
    mitoBackground = imopen(mitoSmooth, strel('disk', round(mitoBackgroundRadius), 0));
    mitoEnhanced = max(mitoSmooth - mitoBackground, 0);
    mitoNormalised = robustNormalise(mitoEnhanced);
    mitoValues = mitoNormalised(cellMask);
    if isempty(mitoValues) || max(mitoValues) <= min(mitoValues)
        mitoThreshold = inf;
    else
        mitoThreshold = graythresh(mitoValues(:)) * mitoThresholdScale;
    end
    mitoMask = cellMask & mitoNormalised > max(mitoThreshold, eps);
    mitoMask = bwareaopen(mitoMask, minimumMitoArea, 8);
    if mitoCloseRadius > 0
        mitoMask = imclose(mitoMask, strel('disk', round(mitoCloseRadius), 0));
    end
    if mitoDilateRadius > 0
        mitoMask = imdilate(mitoMask, strel('disk', round(mitoDilateRadius), 0));
    end
    mitoMask = mitoMask & cellMask;

    % Strong smoothing suppresses isolated photon pixels. The local response
    % detects spot contrast, while the absolute bright threshold guarantees
    % that the brightest B12 regions cannot be rejected by a shape/area test.
    b12Coarse = imgaussfilt(double(b12Image), lysoSigmaSmall);
    b12Background = imgaussfilt(double(b12Image), lysoSigmaLarge);
    lysosomeResponse = max(b12Coarse - b12Background, 0);
    candidateArea = cellMask;

    responseValues = lysosomeResponse(candidateArea);
    intensityValues = sort(b12Coarse(candidateArea & b12Coarse > 0));
    if isempty(responseValues) || isempty(intensityValues)
        lysosomeThreshold = inf;
        lysosomeBrightThreshold = inf;
        lysoCore = false(size(cellMask));
    else
        medianResponse = median(responseValues);
        robustSigma = 1.4826 * median(abs(responseValues - medianResponse));
        lysosomeThreshold = max(medianResponse + lysoThresholdMad * robustSigma, ...
            lysoMinPeakFraction * max(responseValues));
        quantileIndex = max(1, min(numel(intensityValues), ...
            round(lysoBrightQuantile * numel(intensityValues))));
        lysosomeBrightThreshold = intensityValues(quantileIndex);

        contrastCore = lysosomeResponse >= lysosomeThreshold;
        brightCore = b12Coarse >= lysosomeBrightThreshold;
        lysoCore = candidateArea & (contrastCore | brightCore);
    end

    % Closing, dilation, and filling deliberately create coarse connected
    % regions instead of one label per tiny high-frequency fragment.
    lysoMask = lysoCore;
    if lysoCloseRadius > 0
        lysoMask = imclose(lysoMask, strel('disk', round(lysoCloseRadius), 0));
    end
    if lysoDilateRadius > 0
        lysoMask = imdilate(lysoMask, strel('disk', round(lysoDilateRadius), 0));
    end
    lysoMask = imfill(lysoMask, 'holes');
    lysoMask = bwareaopen(lysoMask, minimumLysoArea, 8);

    % Mitochondria keep priority, but exclusion happens after bright-spot
    % detection so neighbouring bright spots are not lost prematurely.
    mitoExclusion = mitoMask;
    if lysoExcludeMitoRadius > 0
        mitoExclusion = imdilate(mitoExclusion, strel('disk', round(lysoExcludeMitoRadius), 0));
    end
    lysoMask = lysoMask & cellMask & ~mitoExclusion;

    mitoLabels = bwlabel(mitoMask, 8);
    lysoLabels = bwlabel(lysoMask, 8);
    cytosolMask = cellMask & ~mitoMask & ~lysoMask;

    compartmentMap = zeros(size(cellMask), 'uint8');
    compartmentMap(mitoMask) = 1;
    compartmentMap(lysoMask) = 2;
    compartmentMap(cytosolMask) = 3;
end
