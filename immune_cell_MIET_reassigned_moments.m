function flim = immune_cell_MIET_reassigned_moments(ptu, ismResult, params)
%IMMUNE_CELL_MIET_REASSIGNED_MOMENTS Memory-bounded reassigned mean FLIM.
% Accumulates only the reassigned count and first two TCSPC-bin moments in
% bounded photon blocks. No full-photon double arrays or TCSPC cube are made.

    if nargin < 3 || isempty(params)
        params = struct();
    end
    params = setDefault(params, 'photonBlockSize', 1e6);
    params = setDefault(params, 'minCounts', 5);
    params = setDefault(params, 'oversampleXY', 2);
    validateattributes(params.photonBlockSize, {'numeric'}, ...
        {'scalar','integer','positive'});
    validateattributes(params.minCounts, {'numeric'}, ...
        {'scalar','nonnegative'});

    requiredPtu = {'head','Ngate','im_col','im_line','im_tcspc','im_chan'};
    for k = 1:numel(requiredPtu)
        if ~isfield(ptu, requiredPtu{k})
            error('immune_cell_MIET_reassigned_moments:PTUField', ...
                'Missing PTU field %s.', requiredPtu{k});
        end
    end
    requiredIsm = {'aprImage','channelIDs','shiftsToCenter'};
    for k = 1:numel(requiredIsm)
        if ~isfield(ismResult, requiredIsm{k})
            error('immune_cell_MIET_reassigned_moments:ISMField', ...
                'Missing ACO-ISM field %s.', requiredIsm{k});
        end
    end

    imageHeight = size(ismResult.aprImage, 1);
    imageWidth = size(ismResult.aprImage, 2);
    nx = double(ptu.head.ImgHdr_PixX);
    ny = double(ptu.head.ImgHdr_PixY);
    % This pipeline creates detector stacks as [PixX, PixY, detector].
    % Test that convention first because square scans satisfy both shape
    % checks and would otherwise swap the two APR shift components.
    if imageHeight == nx && imageWidth == ny
        coordinateMode = 'row=col_col=line';
        rowField = 'im_col';
        columnField = 'im_line';
    elseif imageHeight == ny && imageWidth == nx
        coordinateMode = 'row=line_col=col';
        rowField = 'im_line';
        columnField = 'im_col';
    else
        error('immune_cell_MIET_reassigned_moments:Geometry', ...
            'ACO-ISM image dimensions do not match the PTU header.');
    end

    channelLut = channelShiftLookup(ismResult.channelIDs);
    shifts = double(ismResult.shiftsToCenter);
    gateLength = round(double(ptu.Ngate));
    pixelCount = imageHeight * imageWidth;
    count = zeros(pixelCount, 1, 'double');
    sumBin = zeros(pixelCount, 1, 'double');
    sumBinSquared = zeros(pixelCount, 1, 'double');
    globalDecay = zeros(gateLength, 1, 'double');
    alignedCount = min([numel(ptu.(rowField)), numel(ptu.(columnField)), ...
        numel(ptu.im_tcspc), numel(ptu.im_chan)]);
    blockSize = double(params.photonBlockSize);
    oversampling = max(1, round(double(params.oversampleXY)));

    for first = 1:blockSize:alignedCount
        last = min(first + blockSize - 1, alignedCount);
        indices = first:last;
        row = double(ptu.(rowField)(indices));
        column = double(ptu.(columnField)(indices));
        timeBin = double(ptu.im_tcspc(indices));
        channel = double(ptu.im_chan(indices));
        shiftIndex = double(channelLut(channel + 1));
        if any(shiftIndex == 0)
            error('immune_cell_MIET_reassigned_moments:DetectorMapping', ...
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
        pixel = rowBin + (columnBin - 1) * imageHeight;
        count = count + accumarray(pixel, 1, [pixelCount 1], @sum, 0);
        sumBin = sumBin + accumarray(pixel, timeBin, ...
            [pixelCount 1], @sum, 0);
        sumBinSquared = sumBinSquared + accumarray(pixel, timeBin.^2, ...
            [pixelCount 1], @sum, 0);
        globalDecay = globalDecay + accumarray(timeBin, 1, ...
            [gateLength 1], @sum, 0);
    end

    [~, t0Bin] = max(globalDecay);
    dtNs = resolveTimeBinNs(ptu);
    meanArrival = zeros(pixelCount, 1, 'double');
    tauRms = zeros(pixelCount, 1, 'double');
    validPixel = count >= params.minCounts;
    meanBin = zeros(pixelCount, 1, 'double');
    meanBin(validPixel) = sumBin(validPixel) ./ count(validPixel);
    meanArrival(validPixel) = (meanBin(validPixel) - t0Bin) * dtNs;
    varianceBins = zeros(pixelCount, 1, 'double');
    varianceBins(validPixel) = sumBinSquared(validPixel) ./ ...
        count(validPixel) - meanBin(validPixel).^2;
    tauRms(validPixel) = sqrt(max(varianceBins(validPixel), 0)) * dtNs;

    stats = struct();
    stats.tag = reshape(count, [imageHeight imageWidth]);
    stats.meanArrival = reshape(meanArrival, [imageHeight imageWidth]);
    stats.tauMean = max(stats.meanArrival, 0);
    stats.tauRMS = reshape(tauRms, [imageHeight imageWidth]);
    stats.globalDecay = globalDecay;
    stats.tAxisNs = ((1:gateLength).' - 0.5) * dtNs;
    stats.t0Bin = double(t0Bin);
    stats.xyT = [];
    stats.rcT = [];

    flim = struct();
    flim.coordMode = coordinateMode;
    flim.tAxisNs = stats.tAxisNs;
    flim.reassigned = struct('total', stats, 'frames', []);
    flim.total = stats;
    flim.frames = [];
    flim.processing = struct('method', 'bounded reassigned photon moments', ...
        'photonBlockSize', blockSize, 'alignedPhotonCount', alignedCount, ...
        'storesPhotonCube', false);
end

function lut = channelShiftLookup(channelIds)
    channelIds = double(channelIds(:));
    if any(channelIds < 0 | channelIds > 255 | channelIds ~= round(channelIds))
        error('immune_cell_MIET_reassigned_moments:DetectorIds', ...
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

function s = setDefault(s, name, value)
    if ~isfield(s, name) || isempty(s.(name))
        s.(name) = value;
    end
end
