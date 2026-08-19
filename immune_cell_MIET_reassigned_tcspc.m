function [pixelCube, pixelLinearIndices, pooledCurve, meta, fullCube] = ...
        immune_cell_MIET_reassigned_tcspc(ptu, ismResult, pixelMask, ...
        pooledMask, buildFullCube, photonBlockSize)
%IMMUNE_CELL_MIET_REASSIGNED_TCSPC Bounded reassigned TCSPC accumulation.
% Histograms photons in bounded blocks into selected membrane pixels, a
% pooled outside-SLB curve, and optionally one detector-summed uint16 cube.

    if nargin < 5 || isempty(buildFullCube)
        buildFullCube = false;
    end
    if nargin < 6 || isempty(photonBlockSize)
        photonBlockSize = 1e6;
    end
    validateattributes(photonBlockSize, {'numeric'}, ...
        {'scalar','integer','positive'});
    buildFullCube = logical(buildFullCube);
    if ~isequal(size(pixelMask), size(pooledMask)) || ~ismatrix(pixelMask)
        error('immune_cell_MIET_reassigned_tcspc:MaskSize', ...
            'pixelMask and pooledMask must be equally sized 2-D masks.');
    end

    requiredPtu = {'head','im_col','im_line','im_tcspc','im_chan','Ngate'};
    for k = 1:numel(requiredPtu)
        if ~isfield(ptu, requiredPtu{k})
            error('immune_cell_MIET_reassigned_tcspc:PTUField', ...
                'Missing PTU field %s.', requiredPtu{k});
        end
    end
    requiredIsm = {'aprImage','channelIDs','shiftsToCenter'};
    for k = 1:numel(requiredIsm)
        if ~isfield(ismResult, requiredIsm{k})
            error('immune_cell_MIET_reassigned_tcspc:ISMField', ...
                'Missing ACO-ISM field %s.', requiredIsm{k});
        end
    end

    imageHeight = size(ismResult.aprImage, 1);
    imageWidth = size(ismResult.aprImage, 2);
    if ~isequal(size(pixelMask), [imageHeight imageWidth])
        error('immune_cell_MIET_reassigned_tcspc:MaskGeometry', ...
            'Masks must match size(ismResult.aprImage).');
    end
    nx = double(ptu.head.ImgHdr_PixX);
    ny = double(ptu.head.ImgHdr_PixY);
    % This pipeline creates detector stacks as [PixX, PixY, detector].
    % Prefer that convention for square scans, where shape is ambiguous.
    if imageHeight == nx && imageWidth == ny
        rowField = 'im_col';
        columnField = 'im_line';
        coordinateMode = 'row=col_col=line';
    elseif imageHeight == ny && imageWidth == nx
        rowField = 'im_line';
        columnField = 'im_col';
        coordinateMode = 'row=line_col=col';
    else
        error('immune_cell_MIET_reassigned_tcspc:Geometry', ...
            'ACO-ISM image dimensions do not match the PTU header.');
    end

    gateLength = round(double(ptu.Ngate));
    pixelLinearIndices = find(logical(pixelMask));
    selectedPixelCount = numel(pixelLinearIndices);
    selectedLookup = zeros(imageHeight * imageWidth, 1, 'uint32');
    selectedLookup(pixelLinearIndices) = uint32(1:selectedPixelCount);
    pixelCube = zeros(selectedPixelCount, 1, gateLength, 'uint16');
    if buildFullCube
        fullCube = zeros(imageHeight, imageWidth, gateLength, 'uint16');
    else
        fullCube = zeros(0, 0, 0, 'uint16');
    end
    pooledCurve = zeros(gateLength, 1, 'double');
    globalCurve = zeros(gateLength, 1, 'double');
    selectedPhotonCount = 0;
    channelLut = channelShiftLookup(ismResult.channelIDs);
    shifts = double(ismResult.shiftsToCenter);
    oversampling = 2;
    alignedCount = min([numel(ptu.(rowField)), numel(ptu.(columnField)), ...
        numel(ptu.im_tcspc), numel(ptu.im_chan)]);

    for first = 1:photonBlockSize:alignedCount
        last = min(first + photonBlockSize - 1, alignedCount);
        indices = first:last;
        row = double(ptu.(rowField)(indices));
        column = double(ptu.(columnField)(indices));
        timeBin = double(ptu.im_tcspc(indices));
        channel = double(ptu.im_chan(indices));
        shiftIndex = double(channelLut(channel + 1));
        if any(shiftIndex == 0)
            error('immune_cell_MIET_reassigned_tcspc:DetectorMapping', ...
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
        fullPixel = rowBin + (columnBin - 1) * imageHeight;
        globalCurve = globalCurve + accumarray(timeBin, 1, ...
            [gateLength 1], @sum, 0);

        pooledPhoton = logical(pooledMask(fullPixel));
        if any(pooledPhoton)
            pooledCurve = pooledCurve + accumarray(timeBin(pooledPhoton), 1, ...
                [gateLength 1], @sum, 0);
        end
        compactPixel = selectedLookup(fullPixel);
        selectedPhoton = compactPixel > 0;
        if any(selectedPhoton)
            compactLinear = double(compactPixel(selectedPhoton)) + ...
                (timeBin(selectedPhoton) - 1) * selectedPixelCount;
            [uniqueIndex, additions] = groupedCounts(compactLinear);
            updated = double(pixelCube(uniqueIndex)) + additions;
            assertUint16Range(updated, 'membrane TCSPC cube');
            pixelCube(uniqueIndex) = uint16(updated);
            selectedPhotonCount = selectedPhotonCount + nnz(selectedPhoton);
        end
        if buildFullCube
            fullLinear = fullPixel + (timeBin - 1) * imageHeight * imageWidth;
            [uniqueIndex, additions] = groupedCounts(fullLinear);
            updated = double(fullCube(uniqueIndex)) + additions;
            assertUint16Range(updated, 'full reassigned TCSPC cube');
            fullCube(uniqueIndex) = uint16(updated);
        end
    end

    dtNs = resolveTimeBinNs(ptu);
    meta = struct('coordinateMode', coordinateMode, ...
        'imageSize', [imageHeight imageWidth], ...
        'gateLength', gateLength, 'dtNs', dtNs, ...
        'tAxisNs', ((0:gateLength-1).' + 0.5) * dtNs, ...
        'selectedPixelCount', selectedPixelCount, ...
        'selectedPhotonCount', selectedPhotonCount, ...
        'pooledPhotonCount', sum(pooledCurve), ...
        'globalCurve', globalCurve, ...
        'fullCubeStored', buildFullCube, ...
        'fullCubeClass', 'uint16', ...
        'fullCubeSize', size(fullCube), ...
        'photonBlockSize', photonBlockSize, ...
        'temporalDetectorAlignmentApplied', false);
end

function [uniqueIndex, additions] = groupedCounts(linearIndex)
    linearIndex = double(linearIndex(:));
    [uniqueIndex, ~, group] = unique(linearIndex);
    additions = accumarray(group, 1);
end

function assertUint16Range(updated, label)
    if any(updated > double(intmax('uint16')))
        error('immune_cell_MIET_reassigned_tcspc:Uint16Overflow', ...
            '%s exceeds uint16 count range.', label);
    end
end

function lut = channelShiftLookup(channelIds)
    channelIds = double(channelIds(:));
    if any(channelIds < 0 | channelIds > 255 | channelIds ~= round(channelIds))
        error('immune_cell_MIET_reassigned_tcspc:DetectorIds', ...
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
end
