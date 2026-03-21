function flim = reassigned_flim(ptuOut, ismRes, params)
% REASSIGNED_FLIM
%
% Optimized FLIM reconstruction from PTU photons using APR shifts.
%
% Key optimization:
% If keepSameSize = true, we NEVER build the oversampled 3D cube.
% We compute oversampled photon bins, map them directly to native bins,
% and accumulate only the native-size uint16 cube.
%
% Storage policy:
%   - x-y-t cubes are stored as uint16
%   - tag / tauMean / tauRMS / meanArrival / globalDecay are stored as double
%
% Outputs:
%   flim.unassigned.total
%   flim.unassigned.frames
%   flim.reassigned.total
%   flim.reassigned.frames
%
% Compatibility:
%   flim.total  = flim.reassigned.total
%   flim.frames = flim.reassigned.frames

    if nargin < 3
        params = struct();
    end

    params = setDefault(params, 'oversampleXY', 2);
    params = setDefault(params, 'keepSameSize', true);
    params = setDefault(params, 'storeTotalCubes', true);
    params = setDefault(params, 'storeFrameCubes', false);
    params = setDefault(params, 'frameIndices', []);
    params = setDefault(params, 'minCounts', 20);
    params = setDefault(params, 'useBackground', true);
    params = setDefault(params, 'bgBins', []);
    params = setDefault(params, 't0Mode', 'auto');
    params = setDefault(params, 't0Bin', []);
    params = setDefault(params, 'overflowAction', 'error'); % 'error' or 'clip'

    % ------------------------------------------------------------
    % Basic dimensions
    % ------------------------------------------------------------
    nx = double(ptuOut.head.ImgHdr_PixX);
    ny = double(ptuOut.head.ImgHdr_PixY);

    if isfield(ptuOut, 'tcspc_pix') && ~isempty(ptuOut.tcspc_pix)
        Ng = size(ptuOut.tcspc_pix, 3);
    else
        Ng = double(max(ptuOut.im_tcspc));
    end

    imgH = size(ismRes.aprImage, 1);
    imgW = size(ismRes.aprImage, 2);

    dt_ns = ptuOut.head.MeasDesc_Resolution * 1e9;
    tAxisNs = ((1:Ng) - 0.5) * dt_ns;

    OS = max(1, round(params.oversampleXY));

    % ------------------------------------------------------------
    % Determine PTU coordinate -> image row/col convention
    % ------------------------------------------------------------
    if imgH == ny && imgW == nx
        coordMode = 'row=line_col=col';
    elseif imgH == nx && imgW == ny
        coordMode = 'row=col_col=line';
    else
        error(['Could not reconcile image size with PTU dimensions. ' ...
               'size(ismRes.aprImage) = [%d %d], but [PixY PixX] = [%d %d] and [PixX PixY] = [%d %d].'], ...
               imgH, imgW, ny, nx, nx, ny);
    end

    % ------------------------------------------------------------
    % Photon fields
    % ------------------------------------------------------------
    col0 = double(ptuOut.im_col);
    line0 = double(ptuOut.im_line);
    t0 = double(ptuOut.im_tcspc);
    ch = double(ptuOut.im_chan);

    if isfield(ptuOut, 'im_frame') && ~isempty(ptuOut.im_frame)
        fr = double(ptuOut.im_frame);
    else
        fr = ones(size(col0));
    end

    if isempty(params.frameIndices)
        frameKeep = unique(fr(:)).';
    else
        frameKeep = params.frameIndices(:).';
    end

    % ------------------------------------------------------------
    % Channel -> APR shift map
    % ------------------------------------------------------------
    [tf, loc] = ismember(ch, double(ismRes.channelIDs));
    assert(all(tf), 'Some photon channels in ptuOut are not present in ismRes.channelIDs.');

    dRow = ismRes.shiftsToCenter(loc, 1);
    dCol = ismRes.shiftsToCenter(loc, 2);

    % ------------------------------------------------------------
    % PTU coordinates -> image matrix row/col coordinates
    % ------------------------------------------------------------
    switch coordMode
        case 'row=line_col=col'
            row0 = line0;
            col1 = col0;
        case 'row=col_col=line'
            row0 = col0;
            col1 = line0;
        otherwise
            error('Unknown coordMode.');
    end

    % ------------------------------------------------------------
    % Keep selected frames and valid tcspc bins
    % ------------------------------------------------------------
    okBase = t0 >= 1 & t0 <= Ng & ismember(fr, frameKeep);

    row0 = row0(okBase);
    col1 = col1(okBase);
    t0   = t0(okBase);
    fr   = fr(okBase);
    dRow = dRow(okBase);
    dCol = dCol(okBase);

    % ------------------------------------------------------------
    % Convert photons to native bins directly
    % ------------------------------------------------------------
    [rb_un, cb_un, tb_un, fr_un] = map_photons_to_bins( ...
        row0, col1, zeros(size(dRow)), zeros(size(dCol)), t0, fr, ...
        imgH, imgW, Ng, OS, params.keepSameSize);

    [rb_re, cb_re, tb_re, fr_re] = map_photons_to_bins( ...
        row0, col1, dRow, dCol, t0, fr, ...
        imgH, imgW, Ng, OS, params.keepSameSize);

    if params.keepSameSize
        outH = imgH;
        outW = imgW;
    else
        outH = imgH * OS;
        outW = imgW * OS;
    end

    % ------------------------------------------------------------
    % TOTAL cubes
    % ------------------------------------------------------------
    cubeUn = build_cube_from_bins_uint16(rb_un, cb_un, tb_un, outH, outW, Ng, params.overflowAction);
    cubeRe = build_cube_from_bins_uint16(rb_re, cb_re, tb_re, outH, outW, Ng, params.overflowAction);

    statsUn = tcspc_cube_to_lifetime_double(cubeUn, tAxisNs, params);
    statsRe = tcspc_cube_to_lifetime_double(cubeRe, tAxisNs, params);

    % ------------------------------------------------------------
    % FRAMEWISE cubes and maps
    % ------------------------------------------------------------
    frameList = unique(fr(:)).';
    nFrames = numel(frameList);

    framesUn = init_frame_struct(outH, outW, Ng, nFrames, params.storeFrameCubes);
    framesRe = init_frame_struct(outH, outW, Ng, nFrames, params.storeFrameCubes);

    for k = 1:nFrames
        fk = frameList(k);

        indUn = (fr_un == fk);
        indRe = (fr_re == fk);

        if any(indUn)
            cubeF_un = build_cube_from_bins_uint16( ...
                rb_un(indUn), cb_un(indUn), tb_un(indUn), outH, outW, Ng, params.overflowAction);

            statsF_un = tcspc_cube_to_lifetime_double(cubeF_un, tAxisNs, params);
            framesUn = fill_frame_struct(framesUn, statsF_un, cubeF_un, k, params.storeFrameCubes);
        end

        if any(indRe)
            cubeF_re = build_cube_from_bins_uint16( ...
                rb_re(indRe), cb_re(indRe), tb_re(indRe), outH, outW, Ng, params.overflowAction);

            statsF_re = tcspc_cube_to_lifetime_double(cubeF_re, tAxisNs, params);
            framesRe = fill_frame_struct(framesRe, statsF_re, cubeF_re, k, params.storeFrameCubes);
        end
    end

    framesUn.frameNumbers = frameList;
    framesRe.frameNumbers = frameList;

    % ------------------------------------------------------------
    % Output
    % ------------------------------------------------------------
    flim = struct();
    flim.coordMode = coordMode;
    flim.tAxisNs = tAxisNs(:);

    flim.unassigned = struct();
    flim.unassigned.total = statsUn;
    if params.storeTotalCubes
        flim.unassigned.total.xyT = cubeUn;
        flim.unassigned.total.rcT = cubeUn;
    else
        flim.unassigned.total.xyT = [];
        flim.unassigned.total.rcT = [];
    end
    flim.unassigned.frames = framesUn;

    flim.reassigned = struct();
    flim.reassigned.total = statsRe;
    if params.storeTotalCubes
        flim.reassigned.total.xyT = cubeRe;
        flim.reassigned.total.rcT = cubeRe;
    else
        flim.reassigned.total.xyT = [];
        flim.reassigned.total.rcT = [];
    end
    flim.reassigned.frames = framesRe;

    flim.total = flim.reassigned.total;
    flim.frames = flim.reassigned.frames;
end


function [rb, cb, tb, frOut] = map_photons_to_bins(row0, col0, dRow, dCol, t0, fr, H, W, Ng, OS, keepSameSize)
% Map photons to final histogram bins.
%
% If keepSameSize = true:
%   compute oversampled bin, then map directly to its native downsampled bin
%   without ever building the oversampled cube.
%
% If keepSameSize = false:
%   return oversampled bins directly.

    row_os = floor((row0 - 0.5 + dRow) * OS + 0.5);
    col_os = floor((col0 - 0.5 + dCol) * OS + 0.5);
    tb = floor(t0);

    if keepSameSize
        ok = row_os >= 1 & row_os <= H*OS & col_os >= 1 & col_os <= W*OS & tb >= 1 & tb <= Ng;

        row_os = row_os(ok);
        col_os = col_os(ok);
        tb = tb(ok);
        frOut = fr(ok);

        rb = floor((row_os - 1) / OS) + 1;
        cb = floor((col_os - 1) / OS) + 1;
    else
        ok = row_os >= 1 & row_os <= H*OS & col_os >= 1 & col_os <= W*OS & tb >= 1 & tb <= Ng;

        rb = row_os(ok);
        cb = col_os(ok);
        tb = tb(ok);
        frOut = fr(ok);
    end
end


function cube = build_cube_from_bins_uint16(rb, cb, tb, H, W, Ng, overflowAction)
% Build [row,col,t] cube and store as uint16

    if isempty(rb)
        cube = zeros(H, W, Ng, 'uint16');
        return;
    end

    linIdx = sub2ind([H, W, Ng], rb, cb, tb);
    counts = accumarray(linIdx, 1, [H * W * Ng, 1], @sum, 0);
    cube = cast_counts_to_uint16(counts, [H, W, Ng], overflowAction);
end


function cube16 = cast_counts_to_uint16(counts, outSize, overflowAction)
    mx = max(counts(:));

    if mx > double(intmax('uint16'))
        switch lower(overflowAction)
            case 'error'
                error(['TCSPC cube exceeds uint16 range (max count = %g). ' ...
                       'Use fewer frames / lower oversampling, or switch overflowAction to ''clip''.'], mx);
            case 'clip'
                warning('TCSPC cube exceeds uint16 range (max count = %g). Clipping to uint16 max.', mx);
                counts = min(counts, double(intmax('uint16')));
            otherwise
                error('overflowAction must be ''error'' or ''clip''.');
        end
    end

    cube16 = reshape(uint16(counts), outSize);
end


function stats = tcspc_cube_to_lifetime_double(cube16, tAxisNs, params)
% Lifetime/statistics from uint16 cube, returned as double

    cube = double(cube16);
    [H, W, Ng] = size(cube);

    globalDecay = squeeze(sum(sum(cube, 1), 2));

    switch lower(params.t0Mode)
        case 'auto'
            [~, t0Bin] = max(globalDecay);
        case 'manual'
            assert(~isempty(params.t0Bin), 'Manual t0Mode requires params.t0Bin.');
            t0Bin = params.t0Bin;
        otherwise
            error('t0Mode must be ''auto'' or ''manual''.');
    end

    if params.useBackground
        if isempty(params.bgBins)
            nTail = max(3, round(0.1 * Ng));
            bgBins = (Ng - nTail + 1):Ng;
        else
            bgBins = params.bgBins;
        end
    else
        bgBins = [];
    end

    cubeCorr = cube;

    if ~isempty(bgBins)
        bg = mean(cubeCorr(:,:,bgBins), 3);
        cubeCorr = cubeCorr - bg;
        cubeCorr = max(cubeCorr, 0);
    end

    tag = sum(cubeCorr, 3);

    tShift = tAxisNs - tAxisNs(t0Bin);
    tShift3 = reshape(tShift, 1, 1, []);

    tmp1 = sum(cubeCorr .* tShift3, 3);
    tmp2 = sum(cubeCorr .* (tShift3.^2), 3);

    meanArrival = zeros(H, W, 'double');
    tauMean = zeros(H, W, 'double');
    tauRMS = zeros(H, W, 'double');

    valid = tag >= params.minCounts;

    meanArrival(valid) = tmp1(valid) ./ tag(valid);
    tauMean(valid) = max(meanArrival(valid), 0);

    varT = zeros(H, W, 'double');
    varT(valid) = tmp2(valid) ./ tag(valid) - meanArrival(valid).^2;
    varT = max(varT, 0);
    tauRMS(valid) = sqrt(varT(valid));

    stats = struct();
    stats.tag = tag;
    stats.meanArrival = meanArrival;
    stats.tauMean = tauMean;
    stats.tauRMS = tauRMS;
    stats.globalDecay = globalDecay;
    stats.tAxisNs = tAxisNs(:);
    stats.t0Bin = t0Bin;
end


function frames = init_frame_struct(H, W, Ng, nFrames, storeCubes)
    frames = struct();
    frames.tag = zeros(H, W, nFrames, 'double');
    frames.tauMean = zeros(H, W, nFrames, 'double');
    frames.tauRMS = zeros(H, W, nFrames, 'double');
    frames.meanArrival = zeros(H, W, nFrames, 'double');
    frames.globalDecay = zeros(Ng, nFrames, 'double');
    if storeCubes
        frames.xyT = cell(nFrames, 1);
        frames.rcT = cell(nFrames, 1);
    end
end


function frames = fill_frame_struct(frames, statsF, cubeF, k, storeCubes)
    frames.tag(:,:,k) = statsF.tag;
    frames.tauMean(:,:,k) = statsF.tauMean;
    frames.tauRMS(:,:,k) = statsF.tauRMS;
    frames.meanArrival(:,:,k) = statsF.meanArrival;
    frames.globalDecay(:,k) = statsF.globalDecay(:);
    if storeCubes
        frames.xyT{k} = cubeF;
        frames.rcT{k} = cubeF;
    end
end


function s = setDefault(s, fieldName, defaultValue)
    if ~isfield(s, fieldName) || isempty(s.(fieldName))
        s.(fieldName) = defaultValue;
    end
end