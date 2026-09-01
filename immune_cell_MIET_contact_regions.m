function out = immune_cell_MIET_contact_regions(heightMapsMat, cfg)
%IMMUNE_CELL_MIET_CONTACT_REGIONS Low-height contact regions and CD58 enrichment.
%
% out = immune_cell_MIET_contact_regions(heightMapsMat)
% out = immune_cell_MIET_contact_regions(heightMapsMat, cfg)
%
% Finds extended regions where the MIET membrane height falls below a fixed
% threshold - candidate close contacts - and tests whether CD58-Atto488 is
% enriched inside them, over the whole footprint rather than along a single cut.
%
% WHY THE THRESHOLD IS ABSOLUTE AND NOT PER CELL
%
% Contact is a physical claim: the membrane is or is not within some distance of
% the substrate. A per-cell threshold - the lower mode of the height histogram,
% or a fixed quantile - always returns regions, so it would manufacture contacts
% in cells that have none, and the cells in this data set with no low pixels at
% all are a real observation that such a threshold would erase. The cut is
% therefore a number in nanometres, identical for every cell, and it is printed
% on every figure so no reader has to guess what "low" meant.
%
% WHY A PIXELWISE TEST WOULD BE WRONG
%
% The obvious test - rank-sum of CD58 over pixels inside contacts against pixels
% outside - is invalid here, and not marginally so. Sliding 4x4 TCSPC windows
% overlap by 75%, so neighbouring height pixels are not independent, and the
% CD58 puncta span hundreds of pixels each. A test over ~1e5 pixels would return
% p < 1e-10 for evidence that really rests on a handful of independent blobs.
% The effective sample size is the number of REGIONS, not the number of pixels.
%
% THE NULL: ROTATION ABOUT THE FOOTPRINT CENTROID
%
% Every contact region is kept exactly as measured - same shape, same area - and
% rotated about the footprint centroid through a random angle, with a small
% translational jitter. The CD58 image is never touched. Only the RELATIVE
% position of the contacts and the ligand pattern is randomised, so the spatial
% structure of each is preserved, and the observed statistic is read against the
% distribution of statistics from those trials.
%
% Rotation about the centroid is chosen over free translation deliberately. In
% these cells contacts sit preferentially near the footprint rim, and CD58 may
% also vary with radius; a null that moved regions to the cell centre would be
% comparing rim pixels against middle pixels and would report that shared
% geometry as association. Rotation approximately preserves each region's
% distance from the edge, so a purely radial coincidence cannot survive it.
%
% Regions rotate INDEPENDENTLY of one another, not as one rigid constellation.
% See rotationNull below for why: the rigid version was tried first and could
% not place a single trial on a real cell.
%
% WHAT A POSITIVE RESULT DOES NOT YET PROVE
%
% The cell body sits ABOVE the bilayer and attenuates both the 485 excitation
% and the Atto488 emission. Where the membrane is closer there is more material
% in the light path, which would depress CD58 counts with no biology involved.
% That confound can only push the association NEGATIVE, so an enrichment
% (ratio > 1) cannot be explained by it, whereas a depletion (ratio < 1)
% probably should be. Both directions are reported, and the one-sided p-value
% is stated for the enrichment direction only.
%
% cfg fields
%   binning             'sliding4x4'
%   contactThresholdNm  30 - absolute height below which a pixel is "in contact"
%   minRegionAreaUm2    0.25 - about one diffraction spot; smaller blobs are noise
%   minMeasuredFraction 0.35 - reject regions built mostly from uninverted pixels
%   medianFilterPixels  5 - NaN-aware median filter applied before thresholding
%   clusterSigma        5 - CD58 cluster cut, in robust sigma above the median
%   minClusterAreaUm2   0.10
%   nTrials             2000 - randomisation trials
%   jitterFraction      0.10 - translational jitter, as a fraction of footprint radius
%   minInsideFraction   0.95 - a rotated region is accepted only if this much of
%                       it lands inside the footprint
%   angleGridSize       720 - resolution of the one-off sweep for admissible
%                       rotation angles, i.e. half a degree
%   randomSeed          7 - fixed so a rerun reproduces the p-value exactly
%   outputDir, namePrefix, fontSize, showFigures
%
% See also IMMUNE_CELL_MIET_PROFILE_FIGURE, IMMUNE_CELL_MIET_TOPOGRAPHY_FIGURE.

    if nargin < 2 || isempty(cfg); cfg = struct(); end
    cfg = withDefaults(cfg, struct( ...
        'binning', 'sliding4x4', 'contactThresholdNm', 30, ...
        'minRegionAreaUm2', 0.25, 'minMeasuredFraction', 0.35, ...
        'medianFilterPixels', 5, 'clusterSigma', 5, 'minClusterAreaUm2', 0.10, ...
        'nTrials', 2000, 'jitterFraction', 0.10, 'minInsideFraction', 0.95, ...
        'angleGridSize', 720, ...
        'randomSeed', 7, 'displayClip', [0.02 0.98], ...
        'colormap', 'turbo', 'cd58Colormap', 'bone', ...
        'outputDir', '', 'namePrefix', '', 'fontSize', 12, 'showFigures', false));

    heightMapsMat = char(heightMapsMat);
    loaded = load(heightMapsMat, 'heightMaps');
    heightMaps = loaded.heightMaps;
    if ~isfield(heightMaps.stages, cfg.binning)
        error('immune_cell_MIET_contact_regions:NoBinning', ...
            'Binning %s is not in %s', cfg.binning, heightMapsMat);
    end
    stage = heightMaps.stages.(cfg.binning);
    if ~stage.available
        error('immune_cell_MIET_contact_regions:StageUnavailable', ...
            'Binning %s is unavailable in %s', cfg.binning, heightMapsMat);
    end

    acquisitionDir = acquisitionFolder(heightMapsMat);
    [~, acquisition] = fileparts(acquisitionDir);
    if isempty(cfg.outputDir); cfg.outputDir = fileparts(heightMapsMat); end
    if ~isfolder(cfg.outputDir); mkdir(cfg.outputDir); end
    if isempty(cfg.namePrefix)
        cfg.namePrefix = sprintf('immune_cell_MIET_contacts_%s_%s', ...
            regexprep(acquisition, '^_', ''), cfg.binning);
    end

    heightNm = double(stage.heightNm);
    heightNm(~stage.heightMask) = NaN;
    footprint = logical(stage.displayMask);
    if ~any(footprint(:)); footprint = isfinite(heightNm); end
    cd58Full = loadCd58(acquisitionDir, size(heightNm));
    pxUm = pixelSizeFromHeader(acquisitionDir);
    if isempty(pxUm) || ~isfinite(pxUm) || pxUm <= 0
        pxUm = 1; lateralLabel = 'pixels';
    else
        lateralLabel = '\mum';
    end

    [rows, cols] = croppedRange(footprint, 6);
    height = heightNm(rows, cols);
    mask = footprint(rows, cols);
    cd58 = cd58Full(rows, cols);
    measured = isfinite(height) & mask;
    if nnz(measured) < 200
        error('immune_cell_MIET_contact_regions:TooFewPixels', ...
            'Only %d usable height pixels in %s', nnz(measured), heightMapsMat);
    end
    [nRow, nCol] = size(height);

    % ---- contact regions -------------------------------------------------
    % The height map is median-filtered first. Near tau_0 the calibration curve
    % is nearly flat, so dz/dtau reaches several hundred nm/ns and single pixels
    % scatter wildly; without the filter the threshold would carve those single
    % pixels out as "contacts". The filter is NaN-aware, so it never invents
    % height where none was measured.
    smoothed = nanMedianFilter(height, cfg.medianFilterPixels);
    smoothed(~mask) = NaN;

    pxAreaUm2 = pxUm^2;
    lowMask = smoothed < cfg.contactThresholdNm & mask;
    [labels, nRaw] = connectedComponents(lowMask);

    minAreaPx = max(1, ceil(cfg.minRegionAreaUm2 / pxAreaUm2));
    minClusterPx = max(1, ceil(cfg.minClusterAreaUm2 / pxAreaUm2));
    edgeDistance = distanceToBoundary(mask);

    regions = struct('label', {}, 'areaPx', {}, 'areaUm2', {}, ...
        'medianHeightNm', {}, 'minHeightNm', {}, 'measuredFraction', {}, ...
        'meanCd58', {}, 'medianCd58', {}, 'edgeDistanceUm', {}, ...
        'centroidXUm', {}, 'centroidYUm', {});
    contactMask = false(nRow, nCol);
    kept = 0;
    for label = 1:nRaw
        pixels = find(labels == label);
        if numel(pixels) < minAreaPx; continue; end
        frac = nnz(measured(pixels)) / numel(pixels);
        if frac < cfg.minMeasuredFraction; continue; end
        kept = kept + 1;
        contactMask(pixels) = true;
        [r, c] = ind2sub([nRow nCol], pixels);
        regions(kept).label = kept;
        regions(kept).areaPx = numel(pixels);
        regions(kept).areaUm2 = numel(pixels) * pxAreaUm2;
        regions(kept).medianHeightNm = median(height(pixels), 'omitnan');
        regions(kept).minHeightNm = min(height(pixels), [], 'omitnan');
        regions(kept).measuredFraction = frac;
        regions(kept).meanCd58 = mean(cd58(pixels), 'omitnan');
        regions(kept).medianCd58 = median(cd58(pixels), 'omitnan');
        regions(kept).edgeDistanceUm = mean(edgeDistance(pixels)) * pxUm;
        regions(kept).centroidXUm = (mean(c) - 1) * pxUm;
        regions(kept).centroidYUm = (mean(r) - 1) * pxUm;
    end

    % ---- CD58 clusters ---------------------------------------------------
    % Clusters are cut at a robust sigma above the CD58 median rather than at a
    % quantile: a quantile returns the same number of "clusters" from a flat
    % image as from a punctate one, which is exactly the failure mode the
    % absolute height threshold above was chosen to avoid.
    cd58Smooth = nanMedianFilter(cd58, 3);
    inFootprint = cd58Smooth(mask);
    inFootprint = inFootprint(isfinite(inFootprint));
    cd58Median = median(inFootprint);
    cd58Sigma = 1.4826 * median(abs(inFootprint - cd58Median));
    if ~(cd58Sigma > 0); cd58Sigma = max(eps, std(inFootprint)); end
    clusterCut = cd58Median + cfg.clusterSigma * cd58Sigma;
    [clusterLabels, nClusterRaw] = connectedComponents(cd58Smooth > clusterCut & mask);

    clusters = struct('areaUm2', {}, 'peakCd58', {}, 'centroidRow', {}, ...
        'centroidCol', {}, 'centroidXUm', {}, 'centroidYUm', {}, ...
        'inContact', {}, 'medianHeightNm', {});
    nCluster = 0;
    for label = 1:nClusterRaw
        pixels = find(clusterLabels == label);
        if numel(pixels) < minClusterPx; continue; end
        nCluster = nCluster + 1;
        [r, c] = ind2sub([nRow nCol], pixels);
        clusters(nCluster).areaUm2 = numel(pixels) * pxAreaUm2;
        clusters(nCluster).peakCd58 = max(cd58(pixels));
        clusters(nCluster).centroidRow = mean(r);
        clusters(nCluster).centroidCol = mean(c);
        clusters(nCluster).centroidXUm = (mean(c) - 1) * pxUm;
        clusters(nCluster).centroidYUm = (mean(r) - 1) * pxUm;
        clusters(nCluster).inContact = contactMask( ...
            clampIndex(mean(r), nRow), clampIndex(mean(c), nCol));
        clusters(nCluster).medianHeightNm = median(height(pixels), 'omitnan');
    end

    % ---- observed statistics and the randomisation -----------------------
    stats = contactStatistics(contactMask, mask, cd58, clusters, nRow, nCol);
    null = struct('trials', 0, 'attempts', 0, 'acceptanceRate', NaN, ...
        'medianRatio', [], 'meanRatio', [], 'clusterFraction', [], ...
        'insideFraction', [], 'meanInsideFraction', NaN, ...
        'angleAvailability', [], 'maskPixels', [], 'enrichment', [], ...
        'observedMaskPixels', NaN, 'maskPixelRatio', NaN);
    if kept > 0 && stats.outsidePixels > 0
        null = rotationNull(contactMask, mask, cd58, clusters, cfg);
    end

    % The MEAN ratio is the primary statistic and the median ratio only a coarse
    % cross-check. CD58 here is a photon count with a median of a few per pixel,
    % so a ratio of medians can only take values like 3/4, 4/4 or 5/4: on 151013
    % both the observed value and almost every null trial landed on exactly
    % 0.750, which returns p = 1 on both sides from pure quantisation rather
    % than from any absence of effect. The mean is continuous and is the
    % efficient estimator for counts, so it is what the figure and the p-value
    % report.
    p = struct();
    p.enrichedGreater = empiricalP(null.enrichment, stats.enrichment, 'greater');
    p.enrichedLess = empiricalP(null.enrichment, stats.enrichment, 'less');
    p.enrichedTwoSided = empiricalP(null.enrichment, stats.enrichment, 'two');
    p.meanRatioGreater = empiricalP(null.meanRatio, stats.meanRatio, 'greater');
    p.meanRatioLess = empiricalP(null.meanRatio, stats.meanRatio, 'less');
    p.medianRatioGreater = empiricalP(null.medianRatio, stats.medianRatio, 'greater');
    p.medianRatioLess = empiricalP(null.medianRatio, stats.medianRatio, 'less');
    p.clusterFractionGreater = empiricalP(null.clusterFraction, stats.clusterFraction, 'greater');

    out = struct();
    out.heightMapsMat = heightMapsMat;
    out.acquisition = acquisition;
    out.binning = cfg.binning;
    out.pixelSizeUm = pxUm;
    out.contactThresholdNm = cfg.contactThresholdNm;
    out.footprintAreaUm2 = nnz(mask) * pxAreaUm2;
    out.measuredFraction = nnz(measured) / max(1, nnz(mask));
    out.nRegions = kept;
    out.nRegionsBeforeFilter = nRaw;
    out.contactAreaUm2 = nnz(contactMask) * pxAreaUm2;
    out.contactAreaFraction = nnz(contactMask) / max(1, nnz(mask));
    out.regions = regions;
    out.nClusters = nCluster;
    out.clusters = clusters;
    out.clusterCut = clusterCut;
    out.cd58Median = cd58Median;
    out.cd58Sigma = cd58Sigma;
    out.stats = stats;
    out.null = null;
    out.p = p;
    out.config = cfg;
    out.figure = fullfile(cfg.outputDir, sprintf('%s.png', cfg.namePrefix));

    renderContacts(height, cd58, mask, contactMask, clusters, stats, null, p, ...
        pxUm, lateralLabel, cfg, heightMaps, acquisition, out, out.figure);
    fprintf(['  %-11s %d region(s), %.2f um2 (%.1f%% of footprint), CD58 enrichment ' ...
        '%.3f (enrich %s, deplete %s)\n'], cfg.binning, kept, out.contactAreaUm2, ...
        100 * out.contactAreaFraction, stats.enrichment, ...
        formatP(p.enrichedGreater), formatP(p.enrichedLess));
    fprintf('  %s\n', out.figure);
end

% -------------------------------------------------------------- statistics

function stats = contactStatistics(contactMask, footprint, cd58, clusters, nRow, nCol)
% CD58 inside the contact regions against CD58 elsewhere in the same footprint.
% Both a median ratio (robust, tests general enrichment) and a mean ratio
% (sensitive to bright puncta) are returned, because they answer different
% questions and can disagree: a few bright clusters move the mean and leave the
% median alone.
    inside = cd58(contactMask & footprint);
    outside = cd58(~contactMask & footprint);
    whole = cd58(footprint);
    inside = inside(isfinite(inside));
    outside = outside(isfinite(outside));
    whole = whole(isfinite(whole));

    stats = struct();
    stats.insidePixels = numel(inside);
    stats.outsidePixels = numel(outside);
    stats.medianInside = safeStat(@median, inside);
    stats.medianOutside = safeStat(@median, outside);
    stats.meanInside = safeStat(@mean, inside);
    stats.meanOutside = safeStat(@mean, outside);
    stats.footprintMean = safeStat(@mean, whole);
    stats.medianRatio = ratio(stats.medianInside, stats.medianOutside);
    stats.meanRatio = ratio(stats.meanInside, stats.meanOutside);
    % THE TESTED STATISTIC IS NORMALISED BY THE WHOLE FOOTPRINT, NOT BY THE
    % COMPLEMENT OF THE MASK.
    %
    % inside/outside looks like the natural quantity and is the one worth
    % reading, but it cannot be compared against this null. Independently
    % rotated regions may land on one another, so a null mask covers less
    % ground than the observed one - measured at 81% of it on 133830 - which
    % leaves part of the rim in the "outside" pool, lowers the null's outside
    % mean and pushes the null ratio up. Read against that, an unremarkable
    % observation looks strongly depleted: 133830 returned p = 0.0005 for
    % depletion purely from this. Dividing by the footprint mean instead gives
    % a denominator identical in every trial and in the observation, so only
    % the numerator - the CD58 actually covered by the regions - varies.
    stats.enrichment = ratio(stats.meanInside, stats.footprintMean);

    if isempty(clusters)
        stats.clusterFraction = NaN;
        stats.clustersInContact = 0;
    else
        rowsC = clampIndex([clusters.centroidRow], nRow);
        colsC = clampIndex([clusters.centroidCol], nCol);
        hit = contactMask(sub2ind([nRow nCol], rowsC(:), colsC(:)));
        stats.clustersInContact = nnz(hit);
        stats.clusterFraction = nnz(hit) / numel(hit);
    end
end

function null = rotationNull(contactMask, footprint, cd58, clusters, cfg)
% Each region is rotated about the footprint centroid through its OWN random
% angle, with a small translational jitter. The CD58 image, and the cluster
% centroids in it, stay exactly where they are; only the contact regions move.
%
% WHY EACH REGION ROTATES SEPARATELY AND NOT THE WHOLE PATTERN AS ONE BODY
%
% Rotating the entire configuration rigidly is the tidier null on paper, and it
% was tried first. It does not work on these cells. Contacts sit near the rim of
% a footprint that is nowhere near circular, so a single angle that keeps every
% region inside the cell almost never exists: on 151013 it accepted 0 of 80000
% attempts, which is no null distribution at all. Rotating region by region
% keeps everything that matters for the comparison - each region's shape, its
% area, and its distance from the centroid, and therefore roughly its distance
% from the edge - and gives up only the correlations BETWEEN regions, which are
% not what is being tested.
%
% VALID ANGLES ARE FOUND ONCE, NOT REDISCOVERED EVERY TRIAL
%
% Only a small part of the circle keeps a rim-hugging region inside the cell -
% on 151013 about 2% of angles - so drawing angles at random inside the trial
% loop spends almost all of its time on rejected draws, which cost roughly four
% minutes a cell and two hours a session. Each region's admissible angles are
% therefore swept once on a fixed grid, and every trial samples from that list.
% The null is the same; it is reached about two hundred times faster.
%
% A region with no admissible angle at all keeps its best-found placement rather
% than being dropped, so a trial always holds the full complement of regions.
% angleAvailability and meanInsideFraction are returned so a cell where rotation
% had little room to move anything is visible rather than silent.
    [nRow, nCol] = size(contactMask);
    [footRow, footCol] = find(footprint);
    centre = [mean(footRow), mean(footCol)];
    radius = sqrt(nnz(footprint) / pi);
    jitter = cfg.jitterFraction * radius;

    [labels, nRegion] = connectedComponents(contactMask);
    offsets = cell(nRegion, 1);
    for index = 1:nRegion
        [r, c] = ind2sub([nRow nCol], find(labels == index));
        offsets{index} = [r - centre(1), c - centre(2)];
    end

    % ---- sweep the circle once per region --------------------------------
    angleGrid = (0:cfg.angleGridSize - 1) * (2 * pi / cfg.angleGridSize);
    validAngles = cell(nRegion, 1);
    fallbackAngle = zeros(nRegion, 1);
    availability = zeros(nRegion, 1);
    for index = 1:nRegion
        fractions = zeros(numel(angleGrid), 1);
        for k = 1:numel(angleGrid)
            [~, fractions(k)] = placeRegion(offsets{index}, angleGrid(k), ...
                centre, nRow, nCol, footprint);
        end
        validAngles{index} = angleGrid(fractions >= cfg.minInsideFraction);
        availability(index) = numel(validAngles{index}) / numel(angleGrid);
        [~, best] = max(fractions);
        fallbackAngle(index) = angleGrid(best);
    end

    stream = RandStream('mt19937ar', 'Seed', cfg.randomSeed);
    maxJitterRetries = 8;

    null = struct('trials', 0, 'attempts', 0, 'acceptanceRate', NaN, ...
        'medianRatio', nan(cfg.nTrials, 1), 'meanRatio', nan(cfg.nTrials, 1), ...
        'clusterFraction', nan(cfg.nTrials, 1), 'insideFraction', nan(cfg.nTrials, 1), ...
        'maskPixels', nan(cfg.nTrials, 1), 'enrichment', nan(cfg.nTrials, 1));

    placements = 0;
    accepted = 0;
    for trial = 1:cfg.nTrials
        trialMask = false(nRow, nCol);
        insideTotal = 0;
        pixelTotal = 0;
        for index = 1:nRegion
            choices = validAngles{index};
            bestLinear = [];
            bestFraction = -1;
            for retry = 1:maxJitterRetries
                placements = placements + 1;
                if isempty(choices)
                    theta = fallbackAngle(index);
                else
                    theta = choices(randi(stream, numel(choices)));
                end
                shift = jitter * (2 * rand(stream, 1, 2) - 1);
                [linear, fraction] = placeRegion(offsets{index}, theta, ...
                    centre + shift, nRow, nCol, footprint);
                if fraction > bestFraction
                    bestFraction = fraction;
                    bestLinear = linear;
                end
                if fraction >= cfg.minInsideFraction
                    accepted = accepted + 1;
                    break;
                end
            end
            trialMask(bestLinear) = true;
            insideTotal = insideTotal + bestFraction * size(offsets{index}, 1);
            pixelTotal = pixelTotal + size(offsets{index}, 1);
        end

        trialStats = contactStatistics(trialMask, footprint, cd58, clusters, nRow, nCol);
        null.medianRatio(trial) = trialStats.medianRatio;
        null.meanRatio(trial) = trialStats.meanRatio;
        null.enrichment(trial) = trialStats.enrichment;
        null.clusterFraction(trial) = trialStats.clusterFraction;
        null.insideFraction(trial) = insideTotal / max(1, pixelTotal);
        % Independently rotated regions can land on top of one another, which
        % the observed configuration never does. A null mask that is
        % systematically smaller than the observed one would leave more of the
        % contact-like area in the "outside" pool and bias the null ratio, so
        % the covered area is recorded and compared rather than assumed equal.
        null.maskPixels(trial) = nnz(trialMask);
    end

    null.trials = cfg.nTrials;
    null.attempts = placements;
    null.acceptanceRate = accepted / max(1, placements);
    null.angleAvailability = availability;
    null.observedMaskPixels = nnz(contactMask);
    null.maskPixelRatio = median(null.maskPixels, 'omitnan') / max(1, nnz(contactMask));
    null.meanInsideFraction = mean(null.insideFraction, 'omitnan');
end

function [linear, fraction] = placeRegion(offsets, theta, centre, nRow, nCol, footprint)
% Rotate one region's pixel offsets by theta about centre and return the linear
% indices that land inside the footprint, plus the fraction that did.
    rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    moved = offsets * rotation.';
    r = round(moved(:, 1) + centre(1));
    c = round(moved(:, 2) + centre(2));
    onGrid = r >= 1 & r <= nRow & c >= 1 & c <= nCol;
    indices = zeros(numel(r), 1);
    indices(onGrid) = sub2ind([nRow nCol], r(onGrid), c(onGrid));
    good = onGrid;
    good(onGrid) = footprint(indices(onGrid));
    linear = indices(good);
    fraction = nnz(good) / numel(good);
end

function p = empiricalP(nullValues, observed, side)
% (1 + #{null at least as extreme}) / (1 + n). The added one keeps p away from
% zero, which a finite randomisation can never justify.
    nullValues = nullValues(isfinite(nullValues));
    if isempty(nullValues) || ~isfinite(observed); p = NaN; return; end
    n = numel(nullValues);
    switch side
        case 'greater'
            p = (1 + nnz(nullValues >= observed)) / (1 + n);
        case 'less'
            p = (1 + nnz(nullValues <= observed)) / (1 + n);
        otherwise
            centre = median(nullValues);
            p = (1 + nnz(abs(nullValues - centre) >= abs(observed - centre))) / (1 + n);
    end
end

% --------------------------------------------------------------- rendering

function renderContacts(height, cd58, mask, contactMask, clusters, stats, ...
        null, p, pxUm, lateralLabel, cfg, heightMaps, acquisition, out, outputFile)
    fs = cfg.fontSize;
    [nRow, nCol] = size(height);
    xUm = (0:nCol - 1) * pxUm;
    yUm = (0:nRow - 1) * pxUm;

    visibility = 'off';
    if cfg.showFigures; visibility = 'on'; end
    screen = get(groot, 'ScreenSize');
    wanted = [1400 900];
    figureSize = [min(wanted(1), screen(3) - 80), min(wanted(2), screen(4) - 120)];
    h = figure('Color', 'w', 'Visible', visibility, 'Position', [40 40 figureSize]);

    % One outline colour cannot serve both maps. White reads on turbo but
    % vanishes into the bright end of bone, so the CD58 panel gets a saturated
    % orange instead - far enough from the magenta cluster markers to stay
    % distinguishable from them.
    contactColour = [1.00 1.00 1.00];
    contactOnCd58 = [0.92 0.41 0.10];
    clusterColour = [1.00 0.25 0.85];

    ax1 = axes('Parent', h, 'Position', [0.045 0.600 0.26 0.31]);
    drawMap(ax1, xUm, yUm, height, perceptualColormapOrBuiltin(cfg.colormap), ...
        robustLimits(height, cfg.displayClip), lateralLabel, 'MIET height [nm]', fs);
    drawOutline(ax1, contactMask, xUm, yUm, contactColour, 1.4);
    title(ax1, sprintf('height, contacts < %g nm outlined', cfg.contactThresholdNm), ...
        'FontSize', fs, 'FontWeight', 'normal');

    ax2 = axes('Parent', h, 'Position', [0.375 0.600 0.26 0.31]);
    drawMap(ax2, xUm, yUm, cd58, perceptualColormapOrBuiltin(cfg.cd58Colormap), ...
        robustLimits(cd58, cfg.displayClip), lateralLabel, 'CD58 [photons]', fs);
    drawOutline(ax2, contactMask, xUm, yUm, contactOnCd58, 1.4);
    for index = 1:numel(clusters)
        plot(ax2, clusters(index).centroidXUm, clusters(index).centroidYUm, 'o', ...
            'MarkerEdgeColor', clusterColour, 'LineWidth', 1.4, 'MarkerSize', 9);
    end
    title(ax2, sprintf('CD58, %d cluster(s) above %.1f photons', ...
        numel(clusters), out.clusterCut), 'FontSize', fs, 'FontWeight', 'normal');

    % Where the two masks agree and where they do not. This panel answers the
    % question by eye; the histogram below answers it with a number.
    ax3 = axes('Parent', h, 'Position', [0.705 0.600 0.26 0.31]);
    grey = rescaleTo01(cd58, robustLimits(cd58, cfg.displayClip));
    overlay = repmat(0.55 * grey, 1, 1, 3);
    overlay = paint(overlay, contactMask & mask, [0.15 0.75 1.00], 0.55);
    overlay(repmat(~mask, 1, 1, 3)) = 1;
    image(ax3, [min(xUm) max(xUm)], [min(yUm) max(yUm)], overlay);
    axis(ax3, 'image');
    set(ax3, 'YDir', 'reverse', 'FontSize', fs - 1, 'Box', 'off');
    hold(ax3, 'on');
    for index = 1:numel(clusters)
        marker = 'o';
        if clusters(index).inContact; marker = 'p'; end
        plot(ax3, clusters(index).centroidXUm, clusters(index).centroidYUm, marker, ...
            'MarkerEdgeColor', clusterColour, 'LineWidth', 1.4, 'MarkerSize', 10);
    end
    xlabel(ax3, sprintf('x [%s]', lateralLabel), 'FontSize', fs);
    ylabel(ax3, sprintf('y [%s]', lateralLabel), 'FontSize', fs);
    title(ax3, sprintf('contacts (blue) on CD58; %d/%d cluster(s) inside', ...
        stats.clustersInContact, numel(clusters)), 'FontSize', fs, 'FontWeight', 'normal');

    % ---- null distribution ----------------------------------------------
    ax4 = axes('Parent', h, 'Position', [0.065 0.075 0.38 0.28]);
    if null.trials > 0 && any(isfinite(null.enrichment))
        histogram(ax4, null.enrichment, 40, 'FaceColor', [0.62 0.66 0.72], ...
            'EdgeColor', 'none');
        hold(ax4, 'on');
        yLimit = ylim(ax4);
        plot(ax4, stats.enrichment * [1 1], yLimit, '-', ...
            'Color', [0.85 0.15 0.15], 'LineWidth', 2.2);
        ylim(ax4, yLimit);
        legend(ax4, {sprintf('null, %d rotations', null.trials), ...
            sprintf('observed %.3f', stats.enrichment)}, ...
            'Location', 'northoutside', 'Orientation', 'horizontal', ...
            'FontSize', fs - 1, 'Box', 'off');
    else
        text(ax4, 0.5, 0.5, 'no contact regions - no null distribution', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', fs);
        set(ax4, 'XTick', [], 'YTick', []);
    end
    xlabel(ax4, 'CD58 mean in contacts / whole footprint', 'FontSize', fs);
    ylabel(ax4, 'trials', 'FontSize', fs);
    set(ax4, 'FontSize', fs - 1, 'Box', 'off');
    grid(ax4, 'on');

    % ---- per-region view -------------------------------------------------
    ax5 = axes('Parent', h, 'Position', [0.565 0.075 0.38 0.28]);
    if ~isempty(out.regions)
        areas = [out.regions.areaUm2];
        cd58Region = [out.regions.meanCd58];
        edge = [out.regions.edgeDistanceUm];
        sizes = 30 + 220 * rescaleTo01(areas(:).', [min(areas) max(areas) + eps]);
        scatter(ax5, edge, cd58Region, sizes, [0.15 0.45 0.80], 'filled', ...
            'MarkerFaceAlpha', 0.7);
        hold(ax5, 'on');
        plot(ax5, xlim(ax5), stats.footprintMean * [1 1], '--', ...
            'Color', [0.85 0.15 0.15], 'LineWidth', 1.6);
        legend(ax5, {'contact region (marker size = area)', ...
            'CD58 footprint mean'}, 'Location', 'northoutside', ...
            'Orientation', 'horizontal', 'FontSize', fs - 1, 'Box', 'off');
    else
        text(ax5, 0.5, 0.5, 'no contact regions passed the filters', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', fs);
        set(ax5, 'XTick', [], 'YTick', []);
    end
    xlabel(ax5, sprintf('mean distance from footprint edge [%s]', lateralLabel), ...
        'FontSize', fs);
    ylabel(ax5, 'region CD58 mean [photons]', 'FontSize', fs);
    set(ax5, 'FontSize', fs - 1, 'Box', 'off');
    grid(ax5, 'on');

    calib = heightMaps.calibration;
    annotation(h, 'textbox', [0.02 0.955 0.96 0.042], 'String', sprintf( ...
        ['%s   contact regions   |   %s   |   z < %g nm, area > %.2f um2   |   ' ...
         'tau_0 %.2f ns, qy %.3f'], regexprep(acquisition, '^_', ''), ...
        prettyBinning(cfg.binning), cfg.contactThresholdNm, cfg.minRegionAreaUm2, ...
        calibField(calib, 'tauFreeNs'), calibField(calib, 'quantumYield')), ...
        'EdgeColor', 'none', 'FontSize', fs + 3, 'FontWeight', 'bold', ...
        'Interpreter', 'none', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');

    annotation(h, 'textbox', [0.02 0.455 0.96 0.035], 'String', sprintf( ...
        ['%d contact region(s), %.2f um2 = %.1f%% of footprint      |      ' ...
         'CD58 mean %.2f in contacts / %.2f over footprint  =  %.3f'], ...
        out.nRegions, out.contactAreaUm2, 100 * out.contactAreaFraction, ...
        stats.meanInside, stats.footprintMean, stats.enrichment), ...
        'EdgeColor', 'none', 'FontSize', fs + 1, 'FontWeight', 'bold', ...
        'Interpreter', 'none', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');
    annotation(h, 'textbox', [0.02 0.420 0.96 0.035], 'String', sprintf( ...
        ['enriched p = %s,  depleted p = %s      |      clusters in contact ' ...
         '%d/%d (p = %s)      |      null %d rotations, %.0f%% of null mask ' ...
         'inside, null/observed area %.2f'], ...
        formatP(p.enrichedGreater), formatP(p.enrichedLess), ...
        stats.clustersInContact, numel(clusters), formatP(p.clusterFractionGreater), ...
        null.trials, 100 * nullField(null, 'meanInsideFraction'), ...
        nullField(null, 'maskPixelRatio')), ...
        'EdgeColor', 'none', 'FontSize', fs, 'Interpreter', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

    set(h, 'PaperPositionMode', 'auto');
    print(h, outputFile, '-dpng', '-r150');
    if ~cfg.showFigures; close(h); end
end

function drawMap(ax, xUm, yUm, data, map, limits, lateralLabel, barLabel, fs)
    handle = imagesc(ax, [min(xUm) max(xUm)], [min(yUm) max(yUm)], data);
    set(handle, 'AlphaData', isfinite(data));
    colormap(ax, map);
    if all(isfinite(limits)) && limits(2) > limits(1)
        caxis(ax, limits); %#ok<CAXIS>
    end
    axis(ax, 'image');
    set(ax, 'YDir', 'reverse', 'FontSize', fs - 1, 'Box', 'off', ...
        'Color', [0.93 0.93 0.93]);
    hold(ax, 'on');
    xlabel(ax, sprintf('x [%s]', lateralLabel), 'FontSize', fs);
    ylabel(ax, sprintf('y [%s]', lateralLabel), 'FontSize', fs);
    colourBar = colorbar(ax, 'FontSize', fs - 2);
    colourBar.Label.String = barLabel;
    colourBar.Label.FontSize = fs - 1;
end

function drawOutline(ax, binaryMask, xUm, yUm, colour, width)
% contourc at the half level traces the mask boundary without the Image
% Processing Toolbox, and returns every disconnected region separately.
    if ~any(binaryMask(:)); return; end
    matrix = contourc(xUm, yUm, double(binaryMask), [0.5 0.5]);
    index = 1;
    while index < size(matrix, 2)
        count = matrix(2, index);
        segment = matrix(:, index + 1:index + count);
        plot(ax, segment(1, :), segment(2, :), '-', 'Color', colour, ...
            'LineWidth', width);
        index = index + count + 1;
    end
end

function image3 = paint(image3, where, colour, strength)
    for channel = 1:3
        layer = image3(:, :, channel);
        layer(where) = (1 - strength) * layer(where) + strength * colour(channel);
        image3(:, :, channel) = layer;
    end
end

% ----------------------------------------------------------------- helpers

function out = nanMedianFilter(data, window)
% A true 2D median over a square window that ignores NaN, built by stacking the
% shifted copies. medfilt2 is not used: it needs the Image Processing Toolbox
% and it propagates NaN, which would eat the whole neighbourhood of every pixel
% that failed to invert - and on the sparser acquisitions that is most of them.
    window = max(1, round(window));
    if window <= 1; out = data; return; end
    radius = floor(window / 2);
    [nRow, nCol] = size(data);
    stack = nan(nRow, nCol, (2 * radius + 1)^2);
    slice = 0;
    for dRow = -radius:radius
        for dCol = -radius:radius
            slice = slice + 1;
            block = nan(nRow, nCol);
            rowsTo = max(1, 1 - dRow):min(nRow, nRow - dRow);
            colsTo = max(1, 1 - dCol):min(nCol, nCol - dCol);
            block(rowsTo, colsTo) = data(rowsTo + dRow, colsTo + dCol);
            stack(:, :, slice) = block;
        end
    end
    out = median(stack, 3, 'omitnan');
end

function [labels, nComponent] = connectedComponents(binaryMask)
% 8-connected labelling through the base-MATLAB graph object, so no Image
% Processing Toolbox licence is needed anywhere in this file. Each edge is added
% once, from a pixel to the four neighbours that follow it in scan order.
    [nRow, nCol] = size(binaryMask);
    labels = zeros(nRow, nCol);
    pixels = find(binaryMask);
    nComponent = 0;
    if isempty(pixels); return; end

    lookup = zeros(nRow, nCol);
    lookup(pixels) = 1:numel(pixels);
    [r, c] = ind2sub([nRow nCol], pixels);
    self = (1:numel(pixels)).';

    source = [];
    target = [];
    offsets = [0 1; 1 0; 1 1; 1 -1];
    for k = 1:size(offsets, 1)
        rr = r + offsets(k, 1);
        cc = c + offsets(k, 2);
        onGrid = rr >= 1 & rr <= nRow & cc >= 1 & cc <= nCol;
        neighbour = zeros(numel(rr), 1);
        neighbour(onGrid) = lookup(sub2ind([nRow nCol], rr(onGrid), cc(onGrid)));
        take = neighbour > 0;
        source = [source; self(take)];      %#ok<AGROW>
        target = [target; neighbour(take)]; %#ok<AGROW>
    end

    if isempty(source)
        labels(pixels) = 1:numel(pixels);
        nComponent = numel(pixels);
        return;
    end
    component = conncomp(graph(source, target, [], numel(pixels)));
    labels(pixels) = component;
    nComponent = max(component);
end

function distance = distanceToBoundary(footprint)
% Euclidean distance from every footprint pixel to the nearest boundary pixel,
% in pixels. bwdist would do this in one call but is Image Processing Toolbox;
% the boundary here is a few thousand pixels, so the direct computation is fast
% enough and carries no dependency.
    [nRow, nCol] = size(footprint);
    distance = zeros(nRow, nCol);
    if ~any(footprint(:)); return; end
    padded = false(nRow + 2, nCol + 2);
    padded(2:end - 1, 2:end - 1) = footprint;
    interior = padded(2:end - 1, 2:end - 1) & padded(1:end - 2, 2:end - 1) & ...
        padded(3:end, 2:end - 1) & padded(2:end - 1, 1:end - 2) & ...
        padded(2:end - 1, 3:end);
    boundary = footprint & ~interior;
    [boundaryRow, boundaryCol] = find(boundary);
    [pixelRow, pixelCol] = find(footprint);
    values = zeros(numel(pixelRow), 1);
    chunk = 20000;
    for start = 1:chunk:numel(pixelRow)
        stop = min(numel(pixelRow), start + chunk - 1);
        dr = pixelRow(start:stop) - boundaryRow.';
        dc = pixelCol(start:stop) - boundaryCol.';
        values(start:stop) = sqrt(min(dr.^2 + dc.^2, [], 2));
    end
    distance(sub2ind([nRow nCol], pixelRow, pixelCol)) = values;
end

function value = safeStat(fn, values)
    if isempty(values); value = NaN; else; value = fn(values); end
end

function value = ratio(inside, outside)
    if ~isfinite(inside) || ~isfinite(outside) || outside <= 0
        value = NaN;
    else
        value = inside / outside;
    end
end

function value = nullField(null, name)
    value = NaN;
    if isstruct(null) && isfield(null, name); value = null.(name); end
end

function text = formatP(p)
    if ~isfinite(p); text = 'n/a'; else; text = sprintf('%.4f', p); end
end

function scaled = rescaleTo01(data, limits)
    scaled = (double(data) - limits(1)) / max(eps, limits(2) - limits(1));
    scaled = max(0, min(1, scaled));
    scaled(~isfinite(scaled)) = 0;
end

function limits = robustLimits(data, clip)
    values = double(data(:));
    values = values(isfinite(values));
    if isempty(values); limits = [0 1]; return; end
    limits = [quantile(values, clip(1)), quantile(values, clip(2))];
    if ~(limits(2) > limits(1)); limits = [min(values), max(values) + eps]; end
end

function index = clampIndex(value, limit)
    index = max(1, min(limit, round(value)));
end

function cd58 = loadCd58(acquisitionDir, imageSize)
    ptu = fullfile(acquisitionDir, 'RawImage.ptu');
    if ~isfile(ptu)
        warning('immune_cell_MIET_contact_regions:NoPtu', ...
            'No RawImage.ptu in %s; CD58 statistics will be NaN.', acquisitionDir);
        cd58 = nan(imageSize);
        return;
    end
    try
        result = immune_cell_MIET_cd58_image(ptu, struct( ...
            'imageSize', imageSize, 'excitationNm', 485));
        cd58 = double(result.intensity);
    catch err
        warning('immune_cell_MIET_contact_regions:Cd58Failed', ...
            'CD58 unavailable (%s).', err.message);
        cd58 = nan(imageSize);
    end
end

function map = perceptualColormapOrBuiltin(name)
    switch lower(char(name))
        case {'viridis', 'magma', 'cividis', 'coolwarm', 'gray'}
            map = perceptualColormap(name, 256);
        otherwise
            map = feval(char(name), 256);
    end
end

function [rows, cols] = croppedRange(mask, margin)
    rowsAny = find(any(mask, 2));
    colsAny = find(any(mask, 1));
    if isempty(rowsAny); rowsAny = 1:size(mask, 1); end
    if isempty(colsAny); colsAny = 1:size(mask, 2); end
    rows = max(1, rowsAny(1) - margin):min(size(mask, 1), rowsAny(end) + margin);
    cols = max(1, colsAny(1) - margin):min(size(mask, 2), colsAny(end) + margin);
end

function folder = acquisitionFolder(heightMapsMat)
    folder = fileparts(heightMapsMat);
    for step = 1:3
        [parent, name] = fileparts(folder);
        if ~isempty(regexp(name, '^_\d{8}-\d{6}$', 'once')); return; end
        if isempty(parent) || strcmp(parent, folder); break; end
        folder = parent;
    end
end

function pxUm = pixelSizeFromHeader(acquisitionDir)
    pxUm = [];
    ptu = fullfile(acquisitionDir, 'RawImage.ptu');
    if ~isfile(ptu); return; end
    try
        head = PTU_Read_Head(ptu);
    catch
        return;
    end
    if ~isfield(head, 'ImgHdr_PixResol'); return; end
    value = double(head.ImgHdr_PixResol(1));
    if isfinite(value) && value > 0; pxUm = value; end
end

function value = calibField(calib, name)
    value = NaN;
    if isstruct(calib) && isfield(calib, 'params') && isfield(calib.params, name)
        value = double(calib.params.(name));
    elseif isstruct(calib) && isfield(calib, name)
        value = double(calib.(name));
    end
end

function text = prettyBinning(binning)
    switch binning
        case 'native';     text = 'native pixels';
        case 'sliding2x2'; text = '2x2 sliding TCSPC';
        case 'sliding4x4'; text = '4x4 sliding TCSPC';
        otherwise;         text = binning;
    end
end

function cfg = withDefaults(cfg, defaults)
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(cfg, names{index}) || isempty(cfg.(names{index}))
            cfg.(names{index}) = defaults.(names{index});
        end
    end
end
