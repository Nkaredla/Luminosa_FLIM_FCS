function [masks, diagnostics] = immune_cell_MIET_segment(intensity, meanArrivalNs, cfg)
%IMMUNE_CELL_MIET_SEGMENT Segment a basal cell-membrane footprint from SLB.
%
% [masks, diagnostics] = immune_cell_MIET_segment(intensity, ...
%     meanArrivalNs, cfg)
%
% Segments a cell-membrane footprint from the surrounding supported lipid
% bilayer (SLB) in one reassigned FLIM channel. The classification signal is
% strictly the longer mean photon-arrival time of the labelled cell
% membrane. Photon counts are used only to reject unreliable pixels and to
% form a confidence-weighted local mean-arrival estimate; intensity is never
% thresholded or used to find the cell boundary. No figures, dialogs, or
% other GUI objects are created.
%
% Inputs
%   intensity       2-D reassigned photon-count image.
%   meanArrivalNs   2-D mean photon-arrival-time image in ns. NaNs are OK.
%   cfg             Optional scalar structure. Important fields are:
%       minPhotonsPerPixel       reliable-pixel threshold (default 5)
%       smoothSigma              weighted-FLIM smoothing, pixels (1.25)
%       slbSeedMask              known outside-SLB pixels
%                                ([] = global lower FLIM population)
%       cellSeedMask             known cell pixels ([] = automatic)
%       minLifetimeContrastNs    cell-core excess over SLB (0.15 ns)
%       cellSeparationSigma      cell-core excess in SLB MAD units (2.5)
%       minCellArea              [] chooses an image-size default
%       cellCloseRadius          closing radius, pixels (2)
%       cellDilateRadius         final dilation radius, pixels (1)
%       slbGuardRadius           exclusion around cell, pixels (3)
%       maxCellComponents        largest components retained (1)
%       slbLowPopulationQuantile upper quantile for automatic low-FLIM
%                                SLB candidates (0.45)
%       minPlausibleCellFraction small-mask recovery trigger (0.005)
%       maxPlausibleCellFraction reject near-whole-field masks (0.90)
%
% Masks
%   valid                 finite lifetime with enough detected photons.
%   cellFootprint         filled spatial footprint, including low-count gaps.
%   cellMembrane          cellFootprint restricted to valid fitting pixels.
%   cellCore              high-confidence long-arrival cell pixels.
%   cellBoundary          geometric perimeter of cellFootprint.
%   outsideSlbGeometry    all pixels outside the guarded cell footprint,
%                        including low-count and zero-count pixels.
%   outsideSlb            valid, border-connected SLB outside the guard band.
%   slbReference          stable subset recommended for pooled SLB TCSPC.
%   transitionBand        valid pixels deliberately excluded between classes.
%
% The reported SLB mean-arrival value is a segmentation statistic.  It is
% not an exponential lifetime.  The fixed SLB lifetime for Bayesian unmixing
% must be estimated by fitting the pooled slbReference TCSPC decay with the
% measured IRF.

    if nargin < 3 || isempty(cfg)
        cfg = struct();
    end
    validateInputs(intensity, meanArrivalNs, cfg);
    cfg = fillDefaults(cfg, size(intensity));
    validateConfig(cfg, size(intensity));
    requireImageProcessingToolbox();

    intensity = double(intensity);
    meanArrivalNs = double(meanArrivalNs);
    imageSize = size(intensity);

    finiteIntensity = isfinite(intensity) & intensity >= 0;
    valid = finiteIntensity & isfinite(meanArrivalNs) & ...
        meanArrivalNs >= 0 & intensity >= cfg.minPhotonsPerPixel;
    if nnz(valid) < cfg.minimumValidPixels
        error('immune_cell_MIET_segment:InsufficientData', ...
            ['Only %d pixels have finite mean arrival and at least %g photons; ' ...
             'at least %d are required.'], nnz(valid), ...
            cfg.minPhotonsPerPixel, cfg.minimumValidPixels);
    end

    [smoothTau, smoothWeight, smoothSupport, ~] = ...
        weightedSmooth(meanArrivalNs, intensity, valid, cfg);
    analysisSupport = isfinite(smoothTau) & ...
        smoothWeight >= cfg.minSmoothedPhotonWeight & ...
        smoothSupport >= cfg.minValidNeighbourFraction;
    if nnz(analysisSupport) < cfg.minimumValidPixels
        error('immune_cell_MIET_segment:InsufficientSupport', ...
            'Weighted smoothing left only %d supported pixels.', ...
            nnz(analysisSupport));
    end

    borderMask = makeBorderMask(imageSize, cfg.borderWidthPixels);
    supportedValues = smoothTau(analysisSupport);
    [otsuThreshold, otsuRange] = lifetimeOtsuThreshold( ...
        smoothTau, analysisSupport);
    lowQuantileCeiling = finiteQuantile(supportedValues, ...
        cfg.slbLowPopulationQuantile);
    lowPopulationCeiling = minFinite(otsuThreshold, lowQuantileCeiling);

    % The image border is not assumed to be SLB: a valid cell may be
    % cropped by one or more acquisition edges. Instead, initialise the SLB
    % from the equal-pixel lower FLIM population over the entire supported
    % field. Photon counts have already served their confidence role in
    % smoothTau/analysisSupport and do not weight this classification.
    if isempty(cfg.slbSeedMask)
        initialSlbSeed = analysisSupport & ...
            smoothTau <= lowPopulationCeiling;
        seedSource = 'automatic global lower FLIM population';
        if nnz(initialSlbSeed) < cfg.minimumSlbPixels
            lowPopulationCeiling = lowQuantileCeiling;
            initialSlbSeed = analysisSupport & ...
                smoothTau <= lowPopulationCeiling;
            seedSource = ...
                'automatic global lower FLIM quantile (Otsu class too small)';
        end
    else
        initialSlbSeed = logical(cfg.slbSeedMask) & analysisSupport;
        seedSource = 'cfg.slbSeedMask';
    end
    if nnz(initialSlbSeed) < cfg.minimumSlbPixels
        error('immune_cell_MIET_segment:InsufficientSLBSeed', ...
            'The SLB seed contains %d supported pixels; at least %d are required.', ...
            nnz(initialSlbSeed), cfg.minimumSlbPixels);
    end

    [slbLocation, slbScale, slbEstimateCount] = estimateLowPopulation( ...
        smoothTau, initialSlbSeed, Inf, cfg);

    cellFootprint = false(imageSize);
    cellCore = false(imageSize);
    coreThreshold = Inf;
    growthThreshold = Inf;
    iterationsCompleted = 0;
    fallbackInfo = struct('used', false, 'method', '', ...
        'coreThresholdNs', NaN, 'growthThresholdNs', NaN, ...
        'intensityFallbackUsed', false, 'attempted', false, ...
        'replacedPrimary', false, 'primaryAreaPixels', 0, ...
        'recoveredAreaPixels', 0, 'selectionReason', 'primary hysteresis');
    for iteration = 1:cfg.refinementIterations
        effectiveSlbScale = max(slbScale, cfg.minSlbSigmaNs);
        robustCoreThreshold = slbLocation + max( ...
            cfg.minLifetimeContrastNs, ...
            cfg.cellSeparationSigma * effectiveSlbScale);
        if isfinite(otsuThreshold)
            coreThreshold = max(robustCoreThreshold, otsuThreshold);
        else
            coreThreshold = robustCoreThreshold;
        end
        growthThreshold = slbLocation + max( ...
            cfg.growthMinLifetimeContrastNs, ...
            cfg.growthSeparationSigma * effectiveSlbScale);
        growthThreshold = max(growthThreshold, slbLocation + ...
            cfg.growthFractionOfCore * (coreThreshold - slbLocation));
        % Otsu supplies the conservative seed/core boundary. Growth stays
        % below that valley so hysteresis can recover the connected lower-
        % contrast part of the same cell footprint.
        growthThreshold = min(growthThreshold, coreThreshold);

        [cellFootprint, cellCore] = segmentFromThresholds( ...
            smoothTau, analysisSupport, coreThreshold, growthThreshold, ...
            cfg.cellSeedMask, cfg);
        exteriorSpatial = outsideOfCell(cellFootprint, cfg.slbGuardRadius, ...
            cfg.slbSeedMask);
        refinementUpper = slbLocation + ...
            cfg.refinementUpperSigma * effectiveSlbScale;
        refinementSeed = exteriorSpatial & analysisSupport & ...
            smoothTau <= refinementUpper;
        iterationsCompleted = iteration;

        if iteration >= cfg.refinementIterations || ...
                ~any(cellFootprint(:)) || ...
                nnz(refinementSeed) < cfg.minimumSlbPixels
            break;
        end
        [newLocation, newScale, newCount] = estimateLowPopulation( ...
            smoothTau, refinementSeed, refinementUpper, cfg);
        if isfinite(newLocation) && isfinite(newScale)
            maximumShift = cfg.refinementMaxLocationShiftSigma * ...
                effectiveSlbScale;
            slbLocation = min(max(newLocation, ...
                slbLocation - maximumShift), slbLocation + maximumShift);
            maximumScale = cfg.refinementMaxScaleFactor * ...
                effectiveSlbScale;
            slbScale = min(newScale, maximumScale);
            slbEstimateCount = newCount;
        end
    end

    % A conservative threshold can leave an empty or deceptively small
    % high-FLIM island. That must not suppress recovery of a much larger
    % coherent cell. The recovery uses only smoothTau and spatial
    % connectivity; it never falls back to intensity.
    primaryFootprint = cellFootprint;
    primaryCore = cellCore;
    primaryArea = nnz(primaryFootprint);
    primaryFraction = primaryArea / max(nnz(analysisSupport), 1);
    primaryImplausiblySmall = primaryFraction < ...
        cfg.minPlausibleCellFraction;
    if (~any(primaryFootprint(:)) || primaryImplausiblySmall) && ...
            cfg.enableCentralFootprintFallback
        [recoveredFootprint, recoveredCore, recoveryInfo] = ...
            centralFootprintFallback(smoothTau, analysisSupport, ...
            slbLocation, effectiveScaleForFallback(slbScale, cfg), ...
            otsuThreshold, cfg);
        recoveredArea = nnz(recoveredFootprint);
        recoveryInfo.attempted = true;
        recoveryInfo.primaryAreaPixels = primaryArea;
        recoveryInfo.recoveredAreaPixels = recoveredArea;
        replacePrimary = recoveredArea > 0 && (...
            primaryArea == 0 || recoveredArea >= ...
            cfg.fallbackReplacementAreaRatio * primaryArea);
        if replacePrimary
            cellFootprint = recoveredFootprint;
            cellCore = recoveredCore;
            recoveryInfo.used = true;
            recoveryInfo.replacedPrimary = true;
            recoveryInfo.selectionReason = ...
                'FLIM-only recovery replaced an empty/implausibly small primary mask';
        else
            cellFootprint = primaryFootprint;
            cellCore = primaryCore;
            recoveryInfo.used = false;
            recoveryInfo.replacedPrimary = false;
            if recoveredArea == 0
                recoveryInfo.selectionReason = ...
                    'FLIM-only recovery found no coherent candidate';
            else
                recoveryInfo.selectionReason = ...
                    'primary retained because recovery was not sufficiently larger';
            end
        end
        fallbackInfo = recoveryInfo;
    end

    % Recalculate the exterior after the last refinement pass. A guard band
    % prevents mixed edge pixels entering the pooled fixed-SLB decay.
    exteriorSpatial = outsideOfCell(cellFootprint, cfg.slbGuardRadius, ...
        cfg.slbSeedMask);
    outsideSlb = exteriorSpatial & valid;
    effectiveSlbScale = max(slbScale, cfg.minSlbSigmaNs);
    referenceRange = slbLocation + effectiveSlbScale * ...
        [-cfg.slbReferenceLowerSigma, cfg.slbReferenceUpperSigma];
    slbReference = outsideSlb & isfinite(smoothTau) & ...
        smoothTau >= referenceRange(1) & smoothTau <= referenceRange(2);

    messages = {};
    if fallbackInfo.attempted
        messages{end + 1} = fallbackInfo.selectionReason;
    end
    referenceFallbackUsed = false;
    if nnz(slbReference) < cfg.minimumSlbPixels && ...
            nnz(outsideSlb) >= cfg.minimumSlbPixels
        slbReference = outsideSlb;
        referenceFallbackUsed = true;
        messages{end + 1} = ...
            'Stable-SLB filtering was too restrictive; outsideSlb is used as slbReference.';
    end

    cellMembrane = cellFootprint & valid;
    cellBoundary = bwperim(cellFootprint, 8);
    transitionBand = valid & ~cellMembrane & ~outsideSlb;
    labelMap = zeros(imageSize, 'uint8');
    labelMap(outsideSlb) = 1;
    labelMap(transitionBand) = 2;
    labelMap(cellMembrane) = 3;

    masks = struct();
    masks.valid = valid;
    masks.cellFootprint = cellFootprint;
    masks.cellMembrane = cellMembrane;
    masks.cellCore = cellCore & cellFootprint;
    masks.cellBoundary = cellBoundary;
    masks.outsideSlbGeometry = exteriorSpatial;
    masks.outsideSlb = outsideSlb;
    masks.slbReference = slbReference;
    masks.transitionBand = transitionBand;
    masks.labelMap = labelMap;

    cellMedian = medianForMask(smoothTau, cellMembrane);
    outsideMedian = medianForMask(smoothTau, outsideSlb);
    lifetimeContrast = cellMedian - slbLocation;
    componentSummary = summariseComponents( ...
        cellFootprint, intensity, meanArrivalNs);

    touchedEdges = [any(cellFootprint(1, :)), ...
        any(cellFootprint(end, :)), any(cellFootprint(:, 1)), ...
        any(cellFootprint(:, end))];
    touchesBorder = any(touchedEdges);
    finalCellFraction = nnz(cellFootprint) / max(nnz(analysisSupport), 1);
    plausibleCellArea = finalCellFraction >= cfg.minPlausibleCellFraction && ...
        finalCellFraction <= cfg.maxPlausibleCellFraction;
    if ~any(cellFootprint(:))
        status = 'no_cell_detected';
        messages{end + 1} = ...
            'No spatially coherent long-arrival cell component passed the thresholds.';
    elseif ~plausibleCellArea
        if finalCellFraction < cfg.minPlausibleCellFraction
            status = 'implausibly_small_cell';
            messages{end + 1} = sprintf( ...
                ['The FLIM-only footprint occupies %.3g of supported pixels, ' ...
                 'below the configured %.3g minimum.'], ...
                finalCellFraction, cfg.minPlausibleCellFraction);
        else
            status = 'implausibly_large_cell';
            messages{end + 1} = sprintf( ...
                ['The FLIM-only footprint occupies %.3g of supported pixels, ' ...
                 'above the configured %.3g maximum.'], ...
                finalCellFraction, cfg.maxPlausibleCellFraction);
        end
    elseif nnz(slbReference) < cfg.minimumSlbPixels
        status = 'insufficient_slb_reference';
        messages{end + 1} = ...
            'Too few outside-SLB reference pixels remain for a stable pooled decay.';
    elseif ~isfinite(lifetimeContrast) || ...
            lifetimeContrast < cfg.minLifetimeContrastNs
        status = 'weak_lifetime_separation';
        messages{end + 1} = ...
            'The filled cell footprint has weak mean-arrival contrast relative to the SLB.';
    else
        status = 'ok';
    end
    if touchesBorder && isempty(cfg.slbSeedMask)
        messages{end + 1} = ...
            ['The FLIM footprint touches an acquisition edge and is recorded ' ...
             'as a cropped-cell candidate; automatic SLB initialization did ' ...
             'not use the image border.'];
    end

    diagnostics = struct();
    diagnostics.algorithmVersion = 2;
    diagnostics.method = ['FLIM-only global-lower-population hysteresis; ' ...
        'photon counts used only for confidence smoothing'];
    if fallbackInfo.used
        diagnostics.method = sprintf('%s; %s', diagnostics.method, fallbackInfo.method);
    end
    diagnostics.status = status;
    diagnostics.messages = messages;
    diagnostics.config = cfg;
    diagnostics.smoothedMeanArrivalNs = smoothTau;
    diagnostics.smoothedPhotonWeight = smoothWeight;
    diagnostics.smoothedValidSupport = smoothSupport;
    diagnostics.analysisSupport = analysisSupport;
    diagnostics.normalisedLogIntensity = robustNormaliseLocal(log1p(intensity));
    diagnostics.lifetimeExcessZ = (smoothTau - slbLocation) ./ effectiveSlbScale;
    diagnostics.borderMask = borderMask;
    diagnostics.initialSlbSeedMask = initialSlbSeed;
    diagnostics.exteriorSpatialMask = exteriorSpatial;
    diagnostics.thresholds = struct( ...
        'slbMeanArrivalNs', slbLocation, ...
        'slbRobustSigmaNs', slbScale, ...
        'effectiveSlbSigmaNs', effectiveSlbScale, ...
        'otsuMeanArrivalNs', otsuThreshold, ...
        'otsuRobustRangeNs', otsuRange, ...
        'slbLowQuantileCeilingNs', lowQuantileCeiling, ...
        'slbLowPopulationCeilingNs', lowPopulationCeiling, ...
        'cellCoreMeanArrivalNs', coreThreshold, ...
        'cellGrowthMeanArrivalNs', growthThreshold, ...
        'slbReferenceRangeNs', referenceRange);
    diagnostics.slb = struct( ...
        'seedSource', seedSource, ...
        'meanArrivalNs', slbLocation, ...
        'robustSigmaNs', slbScale, ...
        'outsideMedianMeanArrivalNs', outsideMedian, ...
        'initialCandidatePixelCount', nnz(initialSlbSeed), ...
        'estimatePixelCount', slbEstimateCount, ...
        'referencePixelCount', nnz(slbReference), ...
        'referencePhotonCount', sumFinite(intensity(slbReference)), ...
        'referenceFallbackUsed', referenceFallbackUsed);
    diagnostics.cell = struct( ...
        'medianMeanArrivalNs', cellMedian, ...
        'contrastOverSlbNs', lifetimeContrast, ...
        'separationInSlbSigma', lifetimeContrast / effectiveSlbScale, ...
        'touchesImageBorder', touchesBorder, ...
        'touchedEdgesTopBottomLeftRight', touchedEdges, ...
        'croppedCellCandidate', touchesBorder && nnz(outsideSlb) >= ...
            cfg.minimumSlbPixels, ...
        'footprintFractionOfSupport', finalCellFraction, ...
        'minimumPlausibleFraction', cfg.minPlausibleCellFraction, ...
        'maximumPlausibleFraction', cfg.maxPlausibleCellFraction, ...
        'plausibleArea', plausibleCellArea, ...
        'primaryAreaPixels', primaryArea, ...
        'primaryFractionOfSupport', primaryFraction, ...
        'components', componentSummary);
    diagnostics.pixelCounts = struct( ...
        'image', numel(intensity), ...
        'valid', nnz(valid), ...
        'analysisSupport', nnz(analysisSupport), ...
        'cellFootprint', nnz(cellFootprint), ...
        'cellMembrane', nnz(cellMembrane), ...
        'cellCore', nnz(masks.cellCore), ...
        'outsideSlbGeometry', nnz(exteriorSpatial), ...
        'outsideSlb', nnz(outsideSlb), ...
        'slbReference', nnz(slbReference), ...
        'transitionBand', nnz(transitionBand));
    diagnostics.refinementIterationsCompleted = iterationsCompleted;
    diagnostics.fallback = fallbackInfo;
    diagnostics.qc = struct( ...
        'strictlyFlimClassified', true, ...
        'intensityThresholdOrBoundaryUsed', false, ...
        'cellAreaPlausible', plausibleCellArea, ...
        'cellFootprintFractionOfSupport', finalCellFraction, ...
        'croppedCellCandidate', touchesBorder && ...
            nnz(outsideSlb) >= cfg.minimumSlbPixels, ...
        'touchedEdgeCount', nnz(touchedEdges), ...
        'outsideSlbFractionOfValidPixels', nnz(outsideSlb) / ...
            max(nnz(valid), 1), ...
        'recoveryAttempted', fallbackInfo.attempted, ...
        'recoverySelected', fallbackInfo.replacedPrimary);
    diagnostics.labelLegend = struct( ...
        'unclassified', uint8(0), 'outsideSlb', uint8(1), ...
        'transitionBand', uint8(2), 'cellMembrane', uint8(3));
    diagnostics.assumptions = { ...
        'Cell-membrane pixels have later mean arrival than the surrounding SLB.', ...
        'Automatic SLB initialization uses the global lower FLIM population, not an image border.', ...
        'Photon intensity affects confidence smoothing/support only, never class membership or a boundary.', ...
        'Outside SLB is connected to an image edge after the FLIM footprint is derived.', ...
        'cellFootprint denotes the filled basal-membrane analysis domain, not only its perimeter.', ...
        'SLB mean arrival is for segmentation only; fit pooled SLB TCSPC with the IRF to fix tau_SLB.'};
end

function cfg = fillDefaults(cfg, imageSize)
    pixelCount = prod(imageSize);
    defaults = struct();
    defaults.minPhotonsPerPixel = 5;
    defaults.smoothSigma = 1.25;
    defaults.minSmoothedPhotonWeight = 2.5;
    defaults.minValidNeighbourFraction = 0.20;
    defaults.intensityWeightCapQuantile = 0.99;
    defaults.borderFraction = 0.12;
    defaults.borderWidthPixels = [];
    defaults.slbSeedMask = [];
    defaults.cellSeedMask = [];
    defaults.minimumValidPixels = max(20, ceil(0.001 * pixelCount));
    defaults.minimumSlbPixels = max(20, ceil(0.001 * pixelCount));
    defaults.slbLowPopulationQuantile = 0.45;
    defaults.minSlbSigmaNs = 0.03;
    defaults.slbClipSigma = 3.0;
    defaults.cellSeparationSigma = 2.5;
    defaults.minLifetimeContrastNs = 0.15;
    defaults.growthSeparationSigma = 1.25;
    defaults.growthMinLifetimeContrastNs = 0.075;
    defaults.growthFractionOfCore = 0.50;
    defaults.minCellArea = max(20, ceil(0.002 * pixelCount));
    defaults.minCellCoreArea = [];
    defaults.cellCloseRadius = 2;
    defaults.cellDilateRadius = 1;
    defaults.slbGuardRadius = 3;
    defaults.maxCellComponents = 1;
    defaults.refinementIterations = 2;
    defaults.refinementUpperSigma = 3.0;
    defaults.refinementMaxLocationShiftSigma = 0.75;
    defaults.refinementMaxScaleFactor = 1.5;
    defaults.slbReferenceLowerSigma = 3.0;
    defaults.slbReferenceUpperSigma = 2.5;
    defaults.minPlausibleCellFraction = 0.005;
    defaults.maxPlausibleCellFraction = 0.90;
    defaults.enableCentralFootprintFallback = true;
    defaults.fallbackLifetimeQuantile = 0.75;
    % Seed recovery at the configured high-FLIM quantile, then grow down to
    % the Otsu valley. A value of one uses the full quantile-to-SLB contrast.
    defaults.fallbackContrastFraction = 1.00;
    defaults.fallbackGrowthFraction = 0.12;
    defaults.fallbackMinContrastNs = 0.03;
    defaults.fallbackCloseRadius = max(4, ceil(0.008 * min(imageSize)));
    defaults.fallbackDilateRadius = 2;
    defaults.fallbackMinArea = max(20, ceil(0.0005 * pixelCount));
    defaults.fallbackCentralityWeight = 0.35;
    defaults.fallbackReplacementAreaRatio = 2.0;
    % Retained for configuration compatibility; intensity recovery is
    % deliberately disabled regardless of this legacy field's value.
    defaults.fallbackUseIntensity = false;
    defaults.fallbackIntensityCoreLevel = 0.60;
    defaults.fallbackIntensityGrowthLevel = 0.30;

    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(cfg, names{k}) || isempty(cfg.(names{k}))
            cfg.(names{k}) = defaults.(names{k});
        end
    end
    if isempty(cfg.borderWidthPixels)
        cfg.borderWidthPixels = max(1, ceil(cfg.borderFraction * min(imageSize)));
    end
    if isempty(cfg.minCellCoreArea)
        cfg.minCellCoreArea = max(3, ceil(0.10 * cfg.minCellArea));
    end
end

function validateInputs(intensity, meanArrivalNs, cfg)
    if ~(isnumeric(intensity) || islogical(intensity)) || ...
            ~isreal(intensity) || isempty(intensity) || ~ismatrix(intensity)
        error('immune_cell_MIET_segment:IntensityInput', ...
            'intensity must be a nonempty real 2-D numeric image.');
    end
    if ~isnumeric(meanArrivalNs) || ~isreal(meanArrivalNs) || ...
            isempty(meanArrivalNs) || ~ismatrix(meanArrivalNs)
        error('immune_cell_MIET_segment:LifetimeInput', ...
            'meanArrivalNs must be a nonempty real 2-D numeric image.');
    end
    if ~isequal(size(intensity), size(meanArrivalNs))
        error('immune_cell_MIET_segment:ImageSize', ...
            'intensity and meanArrivalNs must have identical sizes.');
    end
    if any(isfinite(intensity(:)) & intensity(:) < 0)
        error('immune_cell_MIET_segment:NegativeIntensity', ...
            'intensity cannot contain finite negative photon counts.');
    end
    if ~isstruct(cfg) || ~isscalar(cfg)
        error('immune_cell_MIET_segment:ConfigInput', ...
            'cfg must be a scalar structure.');
    end
end

function validateConfig(cfg, imageSize)
    nonnegativeScalar(cfg.minPhotonsPerPixel, 'minPhotonsPerPixel');
    nonnegativeScalar(cfg.smoothSigma, 'smoothSigma');
    nonnegativeScalar(cfg.minSmoothedPhotonWeight, 'minSmoothedPhotonWeight');
    boundedScalar(cfg.minValidNeighbourFraction, 0, 1, ...
        'minValidNeighbourFraction');
    boundedScalar(cfg.intensityWeightCapQuantile, eps, 1, ...
        'intensityWeightCapQuantile');
    boundedScalar(cfg.borderFraction, eps, 0.5, 'borderFraction');
    positiveInteger(cfg.borderWidthPixels, 'borderWidthPixels');
    positiveInteger(cfg.minimumValidPixels, 'minimumValidPixels');
    positiveInteger(cfg.minimumSlbPixels, 'minimumSlbPixels');
    boundedScalar(cfg.slbLowPopulationQuantile, 0.05, 0.50, ...
        'slbLowPopulationQuantile');
    nonnegativeScalar(cfg.minSlbSigmaNs, 'minSlbSigmaNs');
    nonnegativeScalar(cfg.slbClipSigma, 'slbClipSigma');
    nonnegativeScalar(cfg.cellSeparationSigma, 'cellSeparationSigma');
    nonnegativeScalar(cfg.minLifetimeContrastNs, 'minLifetimeContrastNs');
    nonnegativeScalar(cfg.growthSeparationSigma, 'growthSeparationSigma');
    nonnegativeScalar(cfg.growthMinLifetimeContrastNs, ...
        'growthMinLifetimeContrastNs');
    boundedScalar(cfg.growthFractionOfCore, 0, 1, 'growthFractionOfCore');
    positiveInteger(cfg.minCellArea, 'minCellArea');
    positiveInteger(cfg.minCellCoreArea, 'minCellCoreArea');
    nonnegativeInteger(cfg.cellCloseRadius, 'cellCloseRadius');
    nonnegativeInteger(cfg.cellDilateRadius, 'cellDilateRadius');
    nonnegativeInteger(cfg.slbGuardRadius, 'slbGuardRadius');
    if ~(isnumeric(cfg.maxCellComponents) && isscalar(cfg.maxCellComponents) && ...
            isreal(cfg.maxCellComponents) && ...
            (isinf(cfg.maxCellComponents) || ...
             (isfinite(cfg.maxCellComponents) && cfg.maxCellComponents >= 1 && ...
              cfg.maxCellComponents == round(cfg.maxCellComponents))))
        error('immune_cell_MIET_segment:ConfigValue', ...
            'cfg.maxCellComponents must be Inf or a positive integer.');
    end
    positiveInteger(cfg.refinementIterations, 'refinementIterations');
    nonnegativeScalar(cfg.refinementUpperSigma, 'refinementUpperSigma');
    nonnegativeScalar(cfg.refinementMaxLocationShiftSigma, ...
        'refinementMaxLocationShiftSigma');
    if ~(isnumeric(cfg.refinementMaxScaleFactor) && ...
            isscalar(cfg.refinementMaxScaleFactor) && ...
            isreal(cfg.refinementMaxScaleFactor) && ...
            isfinite(cfg.refinementMaxScaleFactor) && ...
            cfg.refinementMaxScaleFactor >= 1)
        error('immune_cell_MIET_segment:ConfigValue', ...
            'cfg.refinementMaxScaleFactor must be a finite scalar >= 1.');
    end
    nonnegativeScalar(cfg.slbReferenceLowerSigma, 'slbReferenceLowerSigma');
    nonnegativeScalar(cfg.slbReferenceUpperSigma, 'slbReferenceUpperSigma');
    boundedScalar(cfg.minPlausibleCellFraction, 0, 0.5, ...
        'minPlausibleCellFraction');
    boundedScalar(cfg.maxPlausibleCellFraction, 0.5, 1, ...
        'maxPlausibleCellFraction');
    if cfg.maxPlausibleCellFraction <= cfg.minPlausibleCellFraction
        error('immune_cell_MIET_segment:ConfigValue', ...
            'cfg.maxPlausibleCellFraction must exceed cfg.minPlausibleCellFraction.');
    end
    if ~(islogical(cfg.enableCentralFootprintFallback) || ...
            (isnumeric(cfg.enableCentralFootprintFallback) && ...
             isscalar(cfg.enableCentralFootprintFallback)))
        error('immune_cell_MIET_segment:ConfigValue', ...
            'cfg.enableCentralFootprintFallback must be scalar logical/numeric.');
    end
    boundedScalar(cfg.fallbackLifetimeQuantile, 0.5, 0.99, ...
        'fallbackLifetimeQuantile');
    boundedScalar(cfg.fallbackContrastFraction, 0, 1, ...
        'fallbackContrastFraction');
    boundedScalar(cfg.fallbackGrowthFraction, 0, 1, ...
        'fallbackGrowthFraction');
    nonnegativeScalar(cfg.fallbackMinContrastNs, 'fallbackMinContrastNs');
    nonnegativeInteger(cfg.fallbackCloseRadius, 'fallbackCloseRadius');
    nonnegativeInteger(cfg.fallbackDilateRadius, 'fallbackDilateRadius');
    positiveInteger(cfg.fallbackMinArea, 'fallbackMinArea');
    boundedScalar(cfg.fallbackCentralityWeight, 0, 1, ...
        'fallbackCentralityWeight');
    if ~(isnumeric(cfg.fallbackReplacementAreaRatio) && ...
            isscalar(cfg.fallbackReplacementAreaRatio) && ...
            isreal(cfg.fallbackReplacementAreaRatio) && ...
            isfinite(cfg.fallbackReplacementAreaRatio) && ...
            cfg.fallbackReplacementAreaRatio >= 1)
        error('immune_cell_MIET_segment:ConfigValue', ...
            'cfg.fallbackReplacementAreaRatio must be a finite scalar >= 1.');
    end
    if ~(islogical(cfg.fallbackUseIntensity) || ...
            (isnumeric(cfg.fallbackUseIntensity) && ...
             isscalar(cfg.fallbackUseIntensity) && ...
             isfinite(cfg.fallbackUseIntensity)))
        error('immune_cell_MIET_segment:ConfigValue', ...
            'cfg.fallbackUseIntensity must be scalar logical/numeric.');
    end
    boundedScalar(cfg.fallbackIntensityCoreLevel, 0, 1, ...
        'fallbackIntensityCoreLevel');
    boundedScalar(cfg.fallbackIntensityGrowthLevel, 0, 1, ...
        'fallbackIntensityGrowthLevel');
    validateOptionalMask(cfg.slbSeedMask, imageSize, 'slbSeedMask');
    validateOptionalMask(cfg.cellSeedMask, imageSize, 'cellSeedMask');
end

function scale = effectiveScaleForFallback(slbScale, cfg)
    scale = max(double(slbScale), cfg.minSlbSigmaNs);
    if ~isfinite(scale)
        scale = cfg.minSlbSigmaNs;
    end
end

function [footprint, core, info] = centralFootprintFallback( ...
        smoothTau, support, slbLocation, slbScale, otsuThreshold, cfg)
    footprint = false(size(support));
    core = false(size(support));
    info = struct('used', false, 'method', '', ...
        'coreThresholdNs', NaN, 'growthThresholdNs', NaN, ...
        'intensityFallbackUsed', false, 'attempted', false, ...
        'replacedPrimary', false, 'primaryAreaPixels', 0, ...
        'recoveredAreaPixels', 0, 'selectionReason', '');
    values = smoothTau(support & isfinite(smoothTau));
    if isempty(values) || ~isfinite(slbLocation)
        return;
    end

    highLifetime = finiteQuantile(values, cfg.fallbackLifetimeQuantile);
    availableContrast = max(highLifetime - slbLocation, 0);
    coreContrast = max([cfg.fallbackMinContrastNs, ...
        cfg.fallbackContrastFraction * availableContrast, 0.5 * slbScale]);
    growthContrast = max(cfg.fallbackGrowthFraction * coreContrast, ...
        min(0.5 * cfg.fallbackMinContrastNs, coreContrast));
    coreThreshold = slbLocation + coreContrast;
    growthThreshold = slbLocation + growthContrast;
    if isfinite(otsuThreshold)
        growthThreshold = max(growthThreshold, ...
            min(otsuThreshold, coreThreshold));
    end
    core = support & smoothTau >= coreThreshold;
    growthMask = support & smoothTau >= growthThreshold;

    core = bwareaopen(core, max(3, ceil(0.05 * cfg.fallbackMinArea)), 8);
    if ~any(core(:))
        return;
    end
    growthMask = growthMask | core;
    candidate = imreconstruct(core, growthMask, 8);
    if cfg.fallbackCloseRadius > 0
        candidate = imclose(candidate, strel('disk', cfg.fallbackCloseRadius, 0));
    end
    candidate = imfill(candidate, 'holes');
    if cfg.fallbackDilateRadius > 0
        candidate = imdilate(candidate, strel('disk', cfg.fallbackDilateRadius, 0));
    end
    candidate = bwareaopen(candidate, cfg.fallbackMinArea, 8);
    footprint = selectCentralLifetimeComponents(candidate, smoothTau, ...
        slbLocation, cfg.maxCellComponents, cfg.fallbackCentralityWeight);
    if any(footprint(:))
        footprint = imfill(footprint, 'holes');
        core = core & footprint;
        info.used = true;
        info.method = 'central high-lifetime footprint fallback';
        info.coreThresholdNs = coreThreshold;
        info.growthThresholdNs = growthThreshold;
    else
        core(:) = false;
    end
end

function selected = selectCentralLifetimeComponents(candidate, tau, ...
        slbLocation, maximumCount, centralityWeight)
    selected = false(size(candidate));
    components = bwconncomp(candidate, 8);
    if components.NumObjects == 0
        return;
    end
    geometry = regionprops(components, 'Area', 'Centroid');
    centre = ([size(candidate, 2), size(candidate, 1)] + 1) / 2;
    distanceScale = max(min(size(candidate)) / 2, 1);
    score = zeros(components.NumObjects, 1);
    for k = 1:components.NumObjects
        distance = norm(double(geometry(k).Centroid) - centre) / distanceScale;
        centrality = exp(-0.5 * distance^2);
        tauValues = tau(components.PixelIdxList{k});
        tauValues = tauValues(isfinite(tauValues));
        if isempty(tauValues)
            lifetimeFactor = 1;
        else
            lifetimeFactor = max(median(tauValues) - slbLocation, eps);
        end
        centreFactor = (1 - centralityWeight) + centralityWeight * centrality;
        score(k) = double(geometry(k).Area) * centreFactor * lifetimeFactor;
    end
    [~, order] = sort(score, 'descend');
    if isfinite(maximumCount)
        keepCount = min(round(maximumCount), numel(order));
    else
        % Retain the dominant cell. Additional tiny components are more
        % likely debris or detached objects than the basal membrane.
        keepCount = 1;
    end
    for k = 1:keepCount
        selected(components.PixelIdxList{order(k)}) = true;
    end
end

function requireImageProcessingToolbox()
    required = {'imgaussfilt', 'graythresh', 'imreconstruct', 'imclose', ...
        'imfill', 'imdilate', 'strel', 'bwareaopen', 'bwconncomp', ...
        'bwperim', 'regionprops'};
    for k = 1:numel(required)
        if exist(required{k}, 'file') ~= 2
            error('immune_cell_MIET_segment:MissingImageProcessingToolbox', ...
                'Required Image Processing Toolbox function %s is unavailable.', ...
                required{k});
        end
    end
end

function [smoothTau, smoothWeight, smoothSupport, cappedWeight] = ...
        weightedSmooth(tau, intensity, valid, cfg)
    positiveIntensities = intensity(valid & intensity > 0);
    cap = finiteQuantile(positiveIntensities, cfg.intensityWeightCapQuantile);
    if ~isfinite(cap) || cap <= 0
        cap = max(intensity(valid));
    end
    cappedWeight = min(max(intensity, 0), cap);
    cappedWeight(~valid) = 0;
    tauForFilter = tau;
    tauForFilter(~valid) = 0;

    if cfg.smoothSigma > 0
        smoothWeight = imgaussfilt(cappedWeight, cfg.smoothSigma, ...
            'Padding', 'replicate');
        numerator = imgaussfilt(tauForFilter .* cappedWeight, ...
            cfg.smoothSigma, 'Padding', 'replicate');
        smoothSupport = imgaussfilt(double(valid), cfg.smoothSigma, ...
            'Padding', 'replicate');
    else
        smoothWeight = cappedWeight;
        numerator = tauForFilter .* cappedWeight;
        smoothSupport = double(valid);
    end
    smoothTau = numerator ./ max(smoothWeight, eps);
    smoothTau(smoothWeight <= eps) = NaN;
end

function mask = makeBorderMask(imageSize, width)
    width = min([round(width), floor(min(imageSize) / 2)]);
    width = max(width, 1);
    mask = false(imageSize);
    mask(1:width, :) = true;
    mask(end - width + 1:end, :) = true;
    mask(:, 1:width) = true;
    mask(:, end - width + 1:end) = true;
end

function [location, scale, count] = estimateLowPopulation( ...
        tau, mask, locationUpperBound, cfg)
    keep = mask & isfinite(tau);
    values = tau(keep);
    if isempty(values)
        location = NaN;
        scale = NaN;
        count = 0;
        return;
    end

    location = median(values);
    if isfinite(locationUpperBound)
        location = min(location, locationUpperBound);
    end
    lower = values <= location;
    if any(lower)
        scale = 1.4826 * median(location - values(lower));
    else
        scale = 0;
    end
    effectiveScale = max(scale, cfg.minSlbSigmaNs);
    clipped = values >= location - cfg.slbClipSigma * effectiveScale & ...
        values <= location + cfg.slbClipSigma * effectiveScale;
    if nnz(clipped) >= min(cfg.minimumSlbPixels, numel(values))
        values = values(clipped);
        location = median(values);
        scale = 1.4826 * median(abs(values - location));
    end
    if ~isfinite(scale)
        scale = 0;
    end
    count = numel(values);
end

function [threshold, robustRange] = lifetimeOtsuThreshold(tau, mask)
    values = tau(mask & isfinite(tau));
    if numel(values) < 2
        threshold = Inf;
        robustRange = [NaN NaN];
        return;
    end
    low = finiteQuantile(values, 0.02);
    high = finiteQuantile(values, 0.98);
    robustRange = [low high];
    if ~isfinite(low) || ~isfinite(high) || high <= low
        threshold = Inf;
        return;
    end
    normalised = min(max((values - low) ./ (high - low), 0), 1);
    threshold = low + graythresh(normalised(:)) * (high - low);
end

function [cellMask, core] = segmentFromThresholds( ...
        tau, support, coreThreshold, growthThreshold, cellSeedMask, cfg)
    core = support & tau >= coreThreshold;
    if ~isempty(cellSeedMask)
        core = core | (logical(cellSeedMask) & support);
    end
    core = bwareaopen(core, cfg.minCellCoreArea, 8);
    if ~any(core(:))
        cellMask = false(size(core));
        return;
    end

    permissive = support & tau >= growthThreshold;
    permissive = permissive | core;
    cellMask = imreconstruct(core, permissive, 8);
    if cfg.cellCloseRadius > 0
        cellMask = imclose(cellMask, ...
            strel('disk', cfg.cellCloseRadius, 0));
    end
    cellMask = imfill(cellMask, 'holes');
    if cfg.cellDilateRadius > 0
        cellMask = imdilate(cellMask, ...
            strel('disk', cfg.cellDilateRadius, 0));
    end
    cellMask = bwareaopen(cellMask, cfg.minCellArea, 8);
    cellMask = keepLargestComponents(cellMask, cfg.maxCellComponents);
    cellMask = imfill(cellMask, 'holes');
    core = core & cellMask;
end

function exterior = outsideOfCell(cellMask, guardRadius, slbSeedMask)
    guardedCell = cellMask;
    if guardRadius > 0 && any(cellMask(:))
        guardedCell = imdilate(cellMask, strel('disk', guardRadius, 0));
    end
    possibleOutside = ~guardedCell;
    if isempty(slbSeedMask)
        marker = false(size(cellMask));
        marker(1, :) = possibleOutside(1, :);
        marker(end, :) = possibleOutside(end, :);
        marker(:, 1) = possibleOutside(:, 1);
        marker(:, end) = possibleOutside(:, end);
    else
        marker = logical(slbSeedMask) & possibleOutside;
    end
    exterior = imreconstruct(marker, possibleOutside, 8);
end

function mask = keepLargestComponents(mask, maximumCount)
    if ~isfinite(maximumCount)
        return;
    end
    components = bwconncomp(mask, 8);
    if components.NumObjects <= maximumCount
        return;
    end
    areas = cellfun(@numel, components.PixelIdxList);
    [~, order] = sort(areas, 'descend');
    keep = false(size(mask));
    for k = 1:maximumCount
        keep(components.PixelIdxList{order(k)}) = true;
    end
    mask = keep;
end

function summary = summariseComponents(cellMask, intensity, tau)
    components = bwconncomp(cellMask, 8);
    if components.NumObjects == 0
        summary = struct('areaPixels', {}, 'centroidXY', {}, ...
            'boundingBox', {}, 'photonCount', {}, ...
            'medianMeanArrivalNs', {});
        return;
    end
    geometry = regionprops(components, 'Area', 'Centroid', 'BoundingBox');
    summary = repmat(struct('areaPixels', 0, 'centroidXY', [NaN NaN], ...
        'boundingBox', [NaN NaN NaN NaN], 'photonCount', 0, ...
        'medianMeanArrivalNs', NaN), components.NumObjects, 1);
    for k = 1:components.NumObjects
        indices = components.PixelIdxList{k};
        tauValues = tau(indices);
        tauValues = tauValues(isfinite(tauValues));
        summary(k).areaPixels = geometry(k).Area;
        summary(k).centroidXY = geometry(k).Centroid;
        summary(k).boundingBox = geometry(k).BoundingBox;
        summary(k).photonCount = sumFinite(intensity(indices));
        if ~isempty(tauValues)
            summary(k).medianMeanArrivalNs = median(tauValues);
        end
    end
end

function value = medianForMask(values, mask)
    keep = mask & isfinite(values);
    if ~any(keep(:))
        value = NaN;
    else
        value = median(values(keep));
    end
end

function value = finiteQuantile(values, probability)
    values = sort(double(values(isfinite(values))));
    if isempty(values)
        value = NaN;
        return;
    end
    position = 1 + (numel(values) - 1) * min(max(probability, 0), 1);
    lowerIndex = floor(position);
    upperIndex = ceil(position);
    fraction = position - lowerIndex;
    value = values(lowerIndex) * (1 - fraction) + ...
        values(upperIndex) * fraction;
end

function value = minFinite(a, b)
    candidates = double([a b]);
    candidates = candidates(isfinite(candidates));
    if isempty(candidates)
        value = NaN;
    else
        value = min(candidates);
    end
end

function image = robustNormaliseLocal(image)
    finiteValues = image(isfinite(image));
    if isempty(finiteValues)
        image(:) = 0;
        return;
    end
    low = finiteQuantile(finiteValues, 0.01);
    high = finiteQuantile(finiteValues, 0.99);
    if high <= low
        image(:) = 0;
        return;
    end
    image = min(max((image - low) ./ (high - low), 0), 1);
    image(~isfinite(image)) = 0;
end

function total = sumFinite(values)
    values = double(values(:));
    total = sum(values(isfinite(values)));
end

function validateOptionalMask(mask, imageSize, name)
    if isempty(mask)
        return;
    end
    if ~(isnumeric(mask) || islogical(mask)) || ~isreal(mask) || ...
            ~isequal(size(mask), imageSize) || any(~isfinite(double(mask(:))))
        error('immune_cell_MIET_segment:ConfigMask', ...
            'cfg.%s must be a finite numeric/logical mask matching the images.', name);
    end
end

function nonnegativeScalar(value, name)
    if ~(isnumeric(value) && isscalar(value) && isreal(value) && ...
            isfinite(value) && value >= 0)
        error('immune_cell_MIET_segment:ConfigValue', ...
            'cfg.%s must be a finite nonnegative scalar.', name);
    end
end

function boundedScalar(value, lower, upper, name)
    if ~(isnumeric(value) && isscalar(value) && isreal(value) && ...
            isfinite(value) && value >= lower && value <= upper)
        error('immune_cell_MIET_segment:ConfigValue', ...
            'cfg.%s must be in [%g, %g].', name, lower, upper);
    end
end

function positiveInteger(value, name)
    if ~(isnumeric(value) && isscalar(value) && isreal(value) && ...
            isfinite(value) && value >= 1 && value == round(value))
        error('immune_cell_MIET_segment:ConfigValue', ...
            'cfg.%s must be a positive integer.', name);
    end
end

function nonnegativeInteger(value, name)
    if ~(isnumeric(value) && isscalar(value) && isreal(value) && ...
            isfinite(value) && value >= 0 && value == round(value))
        error('immune_cell_MIET_segment:ConfigValue', ...
            'cfg.%s must be a nonnegative integer.', name);
    end
end
