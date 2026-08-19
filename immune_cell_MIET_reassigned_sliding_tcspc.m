function [pixelCube, pixelLinearIndices, intensity, meta] = ...
        immune_cell_MIET_reassigned_sliding_tcspc(ptu, ismResult, ...
        pixelMask, windowSize, photonBlockSize)
%IMMUNE_CELL_MIET_REASSIGNED_SLIDING_TCSPC Compact overlapping TCSPC sums.
%
% [pixelCube, pixelLinearIndices, intensity, meta] = ...
%     immune_cell_MIET_reassigned_sliding_tcspc(ptu, ismResult, ...
%     pixelMask, windowSize, photonBlockSize)
%
% Reassigns each photon with the same detector-shift convention used by
% immune_cell_MIET_reassigned_tcspc, then accumulates overlapping spatial
% windows with a one-pixel stride. pixelCube is [selectedPixel 1 time] and
% pixelLinearIndices gives the corresponding upper-left window anchors in
% the full [imageHeight imageWidth] grid. intensity is a full-size uint32
% image whose nonzero selected entries exactly equal sum(pixelCube, 3).
%
% Only complete windows are retained. Thus a 2-by-2 window produces an
% output grid with the original image dimensions, but anchors in the final
% row or column are invalid and remain zero in intensity. The physical
% centre of a window is displaced by (windowSize - 1) / 2 pixels from its
% stored upper-left anchor. Overlapping samples share photons and are not
% statistically independent.

    if nargin < 4 || isempty(windowSize)
        windowSize = [2 2];
    end
    if nargin < 5 || isempty(photonBlockSize)
        photonBlockSize = 1e6;
    end
    validateattributes(windowSize, {'numeric'}, ...
        {'real','finite','vector','numel',2,'integer','positive'}, ...
        mfilename, 'windowSize');
    validateattributes(photonBlockSize, {'numeric'}, ...
        {'real','finite','scalar','integer','positive'}, ...
        mfilename, 'photonBlockSize');
    if ~ismatrix(pixelMask) || ~(isnumeric(pixelMask) || islogical(pixelMask)) || ...
            ~isreal(pixelMask) || any(~isfinite(pixelMask(:)))
        error('immune_cell_MIET_reassigned_sliding_tcspc:PixelMask', ...
            'pixelMask must be a numeric or logical 2-D mask.');
    end
    pixelMask = logical(pixelMask);
    windowSize = double(windowSize(:)).';

    requiredPtu = {'head','im_col','im_line','im_tcspc','im_chan','Ngate'};
    for k = 1:numel(requiredPtu)
        if ~isfield(ptu, requiredPtu{k})
            error('immune_cell_MIET_reassigned_sliding_tcspc:PTUField', ...
                'Missing PTU field %s.', requiredPtu{k});
        end
    end
    requiredIsm = {'aprImage','channelIDs','shiftsToCenter'};
    for k = 1:numel(requiredIsm)
        if ~isfield(ismResult, requiredIsm{k})
            error('immune_cell_MIET_reassigned_sliding_tcspc:ISMField', ...
                'Missing ACO-ISM field %s.', requiredIsm{k});
        end
    end

    imageHeight = size(ismResult.aprImage, 1);
    imageWidth = size(ismResult.aprImage, 2);
    if ~isequal(size(pixelMask), [imageHeight imageWidth])
        error('immune_cell_MIET_reassigned_sliding_tcspc:MaskGeometry', ...
            'pixelMask must match size(ismResult.aprImage).');
    end
    if any(windowSize > [imageHeight imageWidth])
        error('immune_cell_MIET_reassigned_sliding_tcspc:WindowGeometry', ...
            'windowSize cannot exceed the reassigned image dimensions.');
    end

    nx = double(ptu.head.ImgHdr_PixX);
    ny = double(ptu.head.ImgHdr_PixY);
    % Detector stacks in this pipeline use [PixX, PixY, detector]. Test
    % that convention first because both alternatives match a square scan.
    if imageHeight == nx && imageWidth == ny
        rowField = 'im_col';
        columnField = 'im_line';
        coordinateMode = 'row=col_col=line';
    elseif imageHeight == ny && imageWidth == nx
        rowField = 'im_line';
        columnField = 'im_col';
        coordinateMode = 'row=line_col=col';
    else
        error('immune_cell_MIET_reassigned_sliding_tcspc:Geometry', ...
            'ACO-ISM image dimensions do not match the PTU header.');
    end

    windowHeight = windowSize(1);
    windowWidth = windowSize(2);
    maximumAnchorRow = imageHeight - windowHeight + 1;
    maximumAnchorColumn = imageWidth - windowWidth + 1;
    fullWindowSupport = false(imageHeight, imageWidth);
    fullWindowSupport(1:maximumAnchorRow, 1:maximumAnchorColumn) = true;
    selectedAnchorMask = pixelMask & fullWindowSupport;
    pixelLinearIndices = find(selectedAnchorMask);
    selectedPixelCount = numel(pixelLinearIndices);
    selectedLookup = zeros(imageHeight * imageWidth, 1, 'uint32');
    selectedLookup(pixelLinearIndices) = uint32(1:selectedPixelCount);

    gateLength = double(ptu.Ngate);
    if ~(isscalar(gateLength) && isfinite(gateLength) && gateLength >= 1)
        error('immune_cell_MIET_reassigned_sliding_tcspc:GateLength', ...
            'ptu.Ngate must specify a positive finite TCSPC-bin count.');
    end
    gateLength = round(gateLength);
    pixelCube = zeros(selectedPixelCount, 1, gateLength, 'uint16');
    channelLut = channelShiftLookup(ismResult.channelIDs);
    shifts = double(ismResult.shiftsToCenter);
    if size(shifts, 1) ~= numel(ismResult.channelIDs) || size(shifts, 2) < 2
        error('immune_cell_MIET_reassigned_sliding_tcspc:ShiftGeometry', ...
            'ismResult.shiftsToCenter must contain one two-axis shift per detector.');
    end

    oversampling = 2;
    alignedCount = min([numel(ptu.(rowField)), numel(ptu.(columnField)), ...
        numel(ptu.im_tcspc), numel(ptu.im_chan)]);
    windowPhotonContributionCount = 0;
    uniqueContributingPhotonCount = 0;
    % The MEX performs the same reduction as the MATLAB branch below with a
    % direct pointer increment per window hit. It is optional: if it is not
    % built, the MATLAB path runs and results are identical.
    useMex = exist('immune_cell_MIET_window_accumulate', 'file') == 3;
    accumulateDims = [imageHeight, maximumAnchorRow, maximumAnchorColumn, ...
        windowHeight, windowWidth, selectedPixelCount, gateLength];

    for first = 1:photonBlockSize:alignedCount
        last = min(first + photonBlockSize - 1, alignedCount);
        indices = first:last;
        row = double(ptu.(rowField)(indices));
        column = double(ptu.(columnField)(indices));
        timeBin = double(ptu.im_tcspc(indices));
        channel = double(ptu.im_chan(indices));
        if any(~isfinite(channel) | channel < 0 | channel > 255 | ...
                channel ~= round(channel))
            error('immune_cell_MIET_reassigned_sliding_tcspc:PhotonDetectorId', ...
                'Photon detector IDs must be finite integers in [0,255].');
        end
        shiftIndex = double(channelLut(channel + 1));
        if any(shiftIndex == 0)
            error('immune_cell_MIET_reassigned_sliding_tcspc:DetectorMapping', ...
                'Photon detector IDs are absent from the ACO-ISM result.');
        end
        rowOversampled = floor((row - 0.5 + ...
            shifts(shiftIndex, 1)) * oversampling + 0.5);
        columnOversampled = floor((column - 0.5 + ...
            shifts(shiftIndex, 2)) * oversampling + 0.5);
        valid = rowOversampled >= 1 & ...
            rowOversampled <= imageHeight * oversampling & ...
            columnOversampled >= 1 & ...
            columnOversampled <= imageWidth * oversampling & ...
            timeBin >= 1 & timeBin <= gateLength;
        if ~any(valid)
            continue;
        end

        rowBin = floor((rowOversampled(valid) - 1) / oversampling) + 1;
        columnBin = floor((columnOversampled(valid) - 1) / oversampling) + 1;
        timeBin = timeBin(valid);
        contributed = false(size(timeBin));

        % A source photon contributes to every selected upper-left anchor
        % whose spatial window contains that source pixel. For a 2-by-2
        % window this is at most four overlapping output samples.
        %
        % The destination indices for all window offsets are gathered first
        % and reduced in a single pass per photon block. Reducing per offset
        % instead cost one sort per offset (unique) for every block, which
        % profiled as the largest single self-time in the pipeline.
        if useMex
            [contributionCount, contributedMask, overflow] = ...
                immune_cell_MIET_window_accumulate(pixelCube, rowBin, ...
                columnBin, timeBin, selectedLookup, accumulateDims);
            if overflow
                error('immune_cell_MIET_reassigned_sliding_tcspc:Uint16Overflow', ...
                    ['A sliding-window TCSPC bin exceeds the uint16 count ' ...
                     'range. Reduce the spatial window or use a wider ' ...
                     'compact storage type.']);
            end
            windowPhotonContributionCount = ...
                windowPhotonContributionCount + contributionCount;
            uniqueContributingPhotonCount = ...
                uniqueContributingPhotonCount + nnz(contributedMask);
            continue;
        end

        offsetCount = windowHeight * windowWidth;
        linearParts = cell(offsetCount, 1);
        partCount = 0;
        for rowOffset = 0:windowHeight-1
            anchorRow = rowBin - rowOffset;
            validRow = anchorRow >= 1 & anchorRow <= maximumAnchorRow;
            for columnOffset = 0:windowWidth-1
                anchorColumn = columnBin - columnOffset;
                validAnchor = validRow & anchorColumn >= 1 & ...
                    anchorColumn <= maximumAnchorColumn;
                if ~any(validAnchor)
                    continue;
                end
                sourceIndex = find(validAnchor);
                fullAnchor = anchorRow(sourceIndex) + ...
                    (anchorColumn(sourceIndex) - 1) * imageHeight;
                compactPixel = selectedLookup(fullAnchor);
                selected = compactPixel > 0;
                if ~any(selected)
                    continue;
                end
                sourceIndex = sourceIndex(selected);
                compactPixel = compactPixel(selected);
                partCount = partCount + 1;
                linearParts{partCount} = double(compactPixel(:)) + ...
                    (timeBin(sourceIndex) - 1) * selectedPixelCount;
                contributed(sourceIndex) = true;
            end
        end
        if partCount > 0
            allLinear = vertcat(linearParts{1:partCount});
            clear linearParts
            % Dense scatter-add over the compact cube index space. accumarray
            % reduces duplicates directly, so no unique() sort is needed.
            delta = accumarray(allLinear, 1, [selectedPixelCount * gateLength, 1]);
            touched = find(delta);
            updated = double(pixelCube(touched)) + delta(touched);
            assertUint16Range(updated);
            pixelCube(touched) = uint16(updated);
            windowPhotonContributionCount = ...
                windowPhotonContributionCount + numel(allLinear);
        end
        uniqueContributingPhotonCount = uniqueContributingPhotonCount + ...
            nnz(contributed);
    end

    compactIntensity = zeros(selectedPixelCount, 1, 'double');
    for timeIndex = 1:gateLength
        compactIntensity = compactIntensity + ...
            double(pixelCube(:, 1, timeIndex));
    end
    if any(compactIntensity > double(intmax('uint32')))
        error('immune_cell_MIET_reassigned_sliding_tcspc:IntensityOverflow', ...
            'A sliding-window total exceeds the uint32 intensity range.');
    end
    intensity = zeros(imageHeight, imageWidth, 'uint32');
    intensity(pixelLinearIndices) = uint32(compactIntensity);

    dtNs = resolveTimeBinNs(ptu);
    meta = struct();
    meta.method = 'bounded overlapping reassigned spatial-window TCSPC sum';
    meta.coordinateMode = coordinateMode;
    meta.imageSize = [imageHeight imageWidth];
    meta.outputGridSize = [imageHeight imageWidth];
    meta.windowSize = windowSize;
    meta.windowArea = prod(windowSize);
    meta.stride = [1 1];
    meta.anchorConvention = 'upper-left array index';
    meta.windowCentreOffsetPixels = (windowSize - 1) / 2;
    meta.edgePolicy = ['only complete windows are fitted; incomplete final ' ...
        'row/column anchors remain absent and have zero intensity'];
    meta.fullWindowOnly = true;
    meta.selectedPixelCount = selectedPixelCount;
    meta.requestedMaskPixelCount = nnz(pixelMask);
    meta.edgeExcludedMaskPixelCount = nnz(pixelMask & ~fullWindowSupport);
    meta.windowPhotonContributionCount = windowPhotonContributionCount;
    meta.uniqueContributingPhotonCount = uniqueContributingPhotonCount;
    meta.overlappingWindowsSharePhotons = true;
    meta.compactCubeSize = size(pixelCube);
    meta.compactCubeClass = 'uint16';
    meta.intensitySize = size(intensity);
    meta.intensityClass = 'uint32';
    meta.gateLength = gateLength;
    meta.dtNs = dtNs;
    meta.tAxisNs = ((0:gateLength-1).' + 0.5) * dtNs;
    meta.photonBlockSize = photonBlockSize;
    meta.temporalDetectorAlignmentApplied = false;
end

function assertUint16Range(updated)
    if any(updated > double(intmax('uint16')))
        error('immune_cell_MIET_reassigned_sliding_tcspc:Uint16Overflow', ...
            ['A sliding-window TCSPC bin exceeds the uint16 count range. ' ...
             'Reduce the spatial window or use a wider compact storage type.']);
    end
end

function lut = channelShiftLookup(channelIds)
    channelIds = double(channelIds(:));
    if any(channelIds < 0 | channelIds > 255 | channelIds ~= round(channelIds))
        error('immune_cell_MIET_reassigned_sliding_tcspc:DetectorIds', ...
            'Detector IDs must be integers in [0,255].');
    end
    lut = zeros(256, 1, 'uint16');
    lut(channelIds + 1) = uint16(1:numel(channelIds));
end

function dtNs = resolveTimeBinNs(ptu)
    if isfield(ptu, 'Resolution_ns') && ~isempty(ptu.Resolution_ns)
        dtNs = double(ptu.Resolution_ns);
    else
        dtNs = double(ptu.head.MeasDesc_Resolution) * 1e9;
    end
    if ~(isscalar(dtNs) && isfinite(dtNs) && dtNs > 0)
        error('immune_cell_MIET_reassigned_sliding_tcspc:TimeResolution', ...
            'The PTU TCSPC time-bin width must be positive and finite.');
    end
end
