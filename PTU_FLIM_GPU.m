function out = PTU_FLIM_GPU(name, opts)
% PTU_FLIM_GPU
% End-to-end streaming PTU MultiFrame FLIM reader + reducer.
%
% This version is designed to be FAST in plain MATLAB:
%   - streams the PTU in chunks using PTU_Read
%   - computes per-frame FLIM maps on the fly
%   - accumulates final/global FLIM maps online
%   - avoids storing all photons unless explicitly requested
%   - avoids building dense tcspc cubes unless explicitly requested
%
% IMPORTANT
% ---------
% 1) This code assumes a MultiFrame Scan PTU file with:
%       head.ImgHdr_Ident == 9
%       head.ImgHdr_BiDirect == 0
% 2) PTU parsing remains CPU-side because PTU_Read is CPU code.
% 3) In plain MATLAB, the true bottleneck for "full GPU" FLIM is scatter-add
%    binning of photons into pixels. MATLAB does not provide a great native
%    GPU primitive for that. This implementation therefore focuses on the
%    fastest robust MATLAB-side path. If you later want a real large GPU
%    jump, replace the reducers here with a MEXCUDA atomicAdd kernel.
%
% OUTPUT SHAPE
% ------------
% Arrays are stored as [nx, ny, nCh, nFrame] or [nx, ny, nCh] to remain
% compatible with your existing code. Internally, the linear index is fixed
% consistently as:
%       pixIdx = x + (y-1)*nx
% where x = column, y = line.
%
% REQUIRED EXTERNAL FUNCTIONS
% ---------------------------
%   PTU_Read.m
%   PTU_Read_Head.m
%
% Example
% -------
% opts = struct;
% opts.photonsPerChunk = 5e6;
% opts.computePerFrame = true;
% opts.storeTcspcPix   = false;
% opts.storePhotonLists = false;
% opts.showWaitbar = false;
% out = PTU_FLIM_EndToEnd_GPU('myfile.ptu', opts);
%
% imagesc(out.tags(:,:,1)'); axis image; title('Intensity, channel 1');
%
% Narain-oriented version with online/global accumulation.

    if nargin < 2
        opts = struct();
    end
    opts = setDefaultOpts(opts);

    if ~ischar(name) && ~isstring(name)
        error('Input "name" must be a PTU filename.');
    end
    name = char(name);

    if numel(name) < 4 || ~strcmpi(name(end-2:end), 'ptu')
        error('Input file must be a .ptu file.');
    end

    head = PTU_Read_Head(name);
    if isempty(head)
        error('Could not read PTU header.');
    end

    if ~(head.ImgHdr_Ident == 9)
        error('This reader supports only ImgHdr_Ident == 9 (MultiFrame Scan).');
    end

    if ~(head.ImgHdr_BiDirect == 0)
        error('This version currently supports only monodirectional scan.');
    end

    % -------------------------------------------------------------
    % Geometry / acquisition parameters
    % -------------------------------------------------------------
    nx = double(head.ImgHdr_PixX);
    ny = double(head.ImgHdr_PixY);
    nPix = nx * ny;

    if isfield(head, 'ImgHdr_MaxFrames')
        nzGuess = double(head.ImgHdr_MaxFrames);
    else
        tim_p_frame = head.ImgHdr_TimePerPixel * nx * ny; % ms
        tot_time = head.TTResult_StopAfter;               % ms
        nzGuess = ceil(tot_time / tim_p_frame);
    end

    if isfinite(opts.maxFrames)
        nzAlloc = min(nzGuess, double(opts.maxFrames));
    else
        nzAlloc = nzGuess;
    end

    anzch = 32;

    rawRes_s = head.MeasDesc_Resolution;
    Resolution = max(1e9 * rawRes_s, opts.minLifetimeBin_ns); % ns per coarse bin
    chDiv = max(1, round((Resolution * 1e-9) / rawRes_s));
    Ngate = min(opts.maxNgate, ceil(1e9 * head.MeasDesc_GlobalResolution / Resolution) + 1);

    % Preserve compatibility with downstream code expecting coarse resolution
    head.MeasDesc_Resolution = Resolution * 1e-9;

    % Marker definitions
    LineStart = 4;
    LineStop  = 2;
    Frame     = 3;

    if isfield(head, 'ImgHdr_LineStart')
        LineStart = 2^(head.ImgHdr_LineStart - 1);
    end
    if isfield(head, 'ImgHdr_LineStop')
        LineStop = 2^(head.ImgHdr_LineStop - 1);
    end
    if isfield(head, 'ImgHdr_Frame')
        Frame = 2^(head.ImgHdr_Frame - 1);
    end

    % -------------------------------------------------------------
    % Discover active channels from a short read
    % -------------------------------------------------------------
    [~, ~, tmpchan0, tmpmarkers0] = PTU_Read(name, [1 min(1e4, head.TTResult_NumberOfRecords)], head);
    dind = unique(tmpchan0(tmpmarkers0 == 0));
    dind = dind(:)';
    nCh = numel(dind);

    if nCh == 0
        error('No active photon channels were found in the PTU file.');
    end

    chMapLUT = zeros(1, 256, 'uint16');
    chMapLUT(double(dind) + 1) = uint16(1:nCh);

    % -------------------------------------------------------------
    % Outputs / online accumulators
    % -------------------------------------------------------------
    if opts.computePerFrame
        tag = zeros(nx, ny, nCh, nzAlloc, 'single');
        tau = zeros(nx, ny, nCh, nzAlloc, 'single');
    else
        tag = [];
        tau = [];
    end

    globalCnt  = zeros(nPix, nCh, 'double');
    globalSum1 = zeros(nPix, nCh, 'double');
    globalSum2 = zeros(nPix, nCh, 'double');

    if opts.storeTcspcPix
        estBytes = double(nPix) * double(Ngate) * double(nCh) * 2;
        fprintf('Allocating tcspc_pix as uint16: %.2f GB\n', estBytes / 2^30);
        tcspc_pix = zeros(nx, ny, Ngate, nCh, 'uint16');
    else
        tcspc_pix = [];
    end

    if opts.storePhotonLists
        syncCell  = cell(nzAlloc, 1);
        tcspcCell = cell(nzAlloc, 1);
        chanCell  = cell(nzAlloc, 1);
        lineCell  = cell(nzAlloc, 1);
        colCell   = cell(nzAlloc, 1);
        frameCell = cell(nzAlloc, 1);
    else
        syncCell  = [];
        tcspcCell = [];
        chanCell  = [];
        lineCell  = [];
        colCell   = [];
        frameCell = [];
    end

    if opts.storeTimeCell
        SyncRate = 1 / head.MeasDesc_GlobalResolution;
        timeCell = cell(nCh, 1);
    else
        SyncRate = [];
        timeCell = [];
    end

    dtSum = 0;
    dtNum = 0;

    % -------------------------------------------------------------
    % Carry buffers between chunks
    % -------------------------------------------------------------
    yCarry      = zeros(0,1);
    tcspcCarry  = zeros(0,1,'uint16');
    chanCarry   = zeros(0,1,'uint8');
    markerCarry = zeros(0,1,'uint8');

    turnsStart = zeros(0,1);
    turnsStop  = zeros(0,1);

    cnt = 0;
    tend = 0;
    num = 1;
    frameIdx = 1;

    if opts.showWaitbar
        h = waitbar(0, sprintf('Frame %d / %d', frameIdx, nzAlloc));
    else
        h = [];
    end

    % -------------------------------------------------------------
    % Main chunk loop
    % -------------------------------------------------------------
    while num > 0

        [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 opts.photonsPerChunk], head);
        cnt = cnt + num;

        if num == 0
            break;
        end

        if ~isempty(yCarry)
            tmpy = tmpy + tend;
        end

        % Keep markers + valid photons
        indValid = (tmpmarkers > 0) | ((tmpchan < anzch) & (tmptcspc < Ngate * chDiv));

        yNew      = tmpy(indValid);
        tcspcNew  = floor(double(tmptcspc(indValid)) ./ chDiv) + 1;
        tcspcNew  = min(Ngate, max(1, tcspcNew));
        chanNew   = uint8(tmpchan(indValid));
        markerNew = uint8(tmpmarkers(indValid));

        % Append carry
        yCarry      = [yCarry; yNew]; %#ok<AGROW>
        tcspcCarry  = [tcspcCarry; uint16(tcspcNew)]; %#ok<AGROW>
        chanCarry   = [chanCarry; chanNew]; %#ok<AGROW>
        markerCarry = [markerCarry; markerNew]; %#ok<AGROW>

        % Update line markers
        if LineStart == LineStop
            tmpturns = yCarry(markerCarry == LineStart);
            if numel(turnsStart) > numel(turnsStop)
                turnsStart = [turnsStart; tmpturns(2:2:end)]; %#ok<AGROW>
                turnsStop  = [turnsStop;  tmpturns(1:2:end)]; %#ok<AGROW>
            else
                turnsStart = [turnsStart; tmpturns(1:2:end)]; %#ok<AGROW>
                turnsStop  = [turnsStop;  tmpturns(2:2:end)]; %#ok<AGROW>
            end
        else
            turnsStart = [turnsStart; yCarry(markerCarry == LineStart)]; %#ok<AGROW>
            turnsStop  = [turnsStop;  yCarry(markerCarry == LineStop)]; %#ok<AGROW>
        end

        % Frame boundaries available in current carry
        frameChange = yCarry(markerCarry == Frame);

        % Remove markers from photon carry
        indMarker = markerCarry ~= 0;
        yCarry(indMarker)      = [];
        tcspcCarry(indMarker)  = [];
        chanCarry(indMarker)   = [];
        markerCarry(indMarker) = [];

        if ~isempty(yCarry)
            tend = yCarry(end) + loc;
        else
            tend = loc;
        end

        % ---------------------------------------------------------
        % Process all complete frames now available
        % ---------------------------------------------------------
        for kf = 1:numel(frameChange)
            if frameIdx > nzAlloc
                break;
            end

            frameEnd = frameChange(kf);

            % Photons belonging to this completed frame
            indFramePhot = (yCarry < frameEnd);
            yf  = yCarry(indFramePhot);
            tf  = tcspcCarry(indFramePhot);
            chf = chanCarry(indFramePhot);

            yCarry(indFramePhot)      = [];
            tcspcCarry(indFramePhot)  = [];
            chanCarry(indFramePhot)   = [];
            markerCarry(indFramePhot) = [];

            % Line markers belonging to this frame
            indStartFrame = turnsStart < frameEnd;
            indStopFrame  = turnsStop  < frameEnd;

            startsF = turnsStart(indStartFrame);
            stopsF  = turnsStop(indStopFrame);

            turnsStart(indStartFrame) = [];
            turnsStop(indStopFrame)   = [];

            [startsUse, stopsUse] = pairStartStop(startsF, stopsF);

            if isempty(startsUse) || isempty(yf)
                frameIdx = frameIdx + 1;
                updateWaitbar(h, frameIdx, nzAlloc);
                continue;
            end

            % Restrict photons to valid temporal span
            indSpan = (yf >= startsUse(1)) & (yf <= stopsUse(end));
            yf  = yf(indSpan);
            tf  = tf(indSpan);
            chf = chf(indSpan);

            if isempty(yf)
                frameIdx = frameIdx + 1;
                updateWaitbar(h, frameIdx, nzAlloc);
                continue;
            end

            % Assign photons to lines
            lineIdx = discretize(double(yf(:)), [double(startsUse(:)); inf]);
            lineIdx = lineIdx(:);
            yf  = yf(:);
            tf  = tf(:);
            chf = chf(:);

            stopAtLine = stopsUse(lineIdx);
            stopAtLine = stopAtLine(:);
            ok = ~isnan(lineIdx) & (lineIdx >= 1) & (lineIdx <= ny) & (yf <= stopAtLine);

            yf      = yf(ok);
            tf      = tf(ok);
            chf     = chf(ok);
            lineIdx = lineIdx(ok);

            if isempty(yf)
                frameIdx = frameIdx + 1;
                updateWaitbar(h, frameIdx, nzAlloc);
                continue;
            end

            % Pixel coordinate along the line
            t1v = startsUse(lineIdx);
            t2v = stopsUse(lineIdx);
            t1v = t1v(:);
            t2v = t2v(:);
            frac = (double(yf(:)) - double(t1v)) ./ max(double(t2v - t1v), 1);
            frac = frac(:);

            col = 1 + floor(nx .* frac);
            col = min(nx, max(1, col));
            col = col(:);

            % Store photons only if explicitly requested
            if opts.storePhotonLists
                syncCell{frameIdx}  = yf(:);
                tcspcCell{frameIdx} = uint16(tf(:));
                chanCell{frameIdx}  = uint8(chf(:));
                lineCell{frameIdx}  = uint16(lineIdx(:));
                colCell{frameIdx}   = uint16(col(:));
                frameCell{frameIdx} = uint16(frameIdx * ones(numel(yf),1));
            end

            % Pixel time estimate
            dtThis = double(stopsUse(:) - startsUse(:));
            dtSum = dtSum + sum(dtThis);
            dtNum = dtNum + numel(dtThis);

            % Per-frame reduction
            [frameCnt, frameSum1, frameSum2] = reduceMomentsAllChannels( ...
                uint16(col(:)), uint16(lineIdx(:)), uint8(chf(:)), uint16(tf(:)), ...
                nx, ny, nPix, nCh, chMapLUT);

            if opts.computePerFrame
                tag(:,:,:,frameIdx) = reshape(single(frameCnt), [nx, ny, nCh]);
                tau(:,:,:,frameIdx) = momentsToTau(frameCnt, frameSum1, frameSum2, Resolution, nx, ny, nCh);
            end

            % Online global accumulation
            globalCnt  = globalCnt  + frameCnt;
            globalSum1 = globalSum1 + frameSum1;
            globalSum2 = globalSum2 + frameSum2;

            % Optional per-channel time vectors
            if opts.storeTimeCell
                localCh = chMapLUT(double(chf(:)) + 1);
                tAbs = round(double(yf(:)) / SyncRate / (Resolution * 1e-9)) + double(tf(:));
                for kc = 1:nCh
                    indc = (localCh == kc);
                    if any(indc)
                        timeCell{kc} = [timeCell{kc}; tAbs(indc)]; %#ok<AGROW>
                    end
                end
            end

            % Optional dense tcspc cube, updated online
            if opts.storeTcspcPix
                tcspc_pix = updateTcspcCube(tcspc_pix, uint16(col(:)), uint16(lineIdx(:)), ...
                    uint8(chf(:)), uint16(tf(:)), nx, nPix, nCh, chMapLUT, Ngate);
            end

            frameIdx = frameIdx + 1;
            updateWaitbar(h, frameIdx, nzAlloc);
        end

        if frameIdx > nzAlloc
            break;
        end
    end

    if ~isempty(h) && ishandle(h)
        close(h);
    end

    nFrames = frameIdx - 1;

    if opts.computePerFrame && nFrames < nzAlloc
        tag(:,:,:,nFrames+1:end) = [];
        tau(:,:,:,nFrames+1:end) = [];
    end

    % -------------------------------------------------------------
    % Final/global maps from online accumulators
    % -------------------------------------------------------------
    tags = reshape(single(globalCnt), [nx, ny, nCh]);
    taus = momentsToTau(globalCnt, globalSum1, globalSum2, Resolution, nx, ny, nCh);

    if dtNum > 0
        head.ImgHdr_PixelTime = 1e9 .* (dtSum / dtNum) ./ nx ./ head.TTResult_SyncRate;
        head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
    else
        head.ImgHdr_PixelTime = [];
        head.ImgHdr_DwellTime = [];
    end

    % -------------------------------------------------------------
    % Optional photon list concatenation
    % -------------------------------------------------------------
    if opts.storePhotonLists && nFrames > 0
        usedFrames = 1:nFrames;
        im_sync  = vertcat(syncCell{usedFrames});
        im_tcspc = vertcat(tcspcCell{usedFrames});
        im_chan  = vertcat(chanCell{usedFrames});
        im_line  = vertcat(lineCell{usedFrames});
        im_col   = vertcat(colCell{usedFrames});
        im_frame = vertcat(frameCell{usedFrames});
    else
        im_sync  = zeros(0,1);
        im_tcspc = uint16([]);
        im_chan  = uint8([]);
        im_line  = uint16([]);
        im_col   = uint16([]);
        im_frame = uint16([]);
    end

    if ~opts.storeTimeCell
        timeCell = cell(nCh, 1);
    end

    % -------------------------------------------------------------
    % Pack output
    % -------------------------------------------------------------
    out = struct();
    out.head = head;
    out.tag = tag;
    out.tau = tau;
    out.tags = tags;
    out.taus = taus;
    out.time = timeCell;
    out.tcspc_pix = tcspc_pix;
    out.im_sync = im_sync;
    out.im_tcspc = im_tcspc;
    out.im_chan = im_chan;
    out.im_line = im_line;
    out.im_col = im_col;
    out.im_frame = im_frame;
    out.dind = dind;
    out.nFrames = nFrames;
    out.Resolution_ns = Resolution;
    out.Ngate = Ngate;
    out.options = opts;
end

% =====================================================================
% Helpers
% =====================================================================

function opts = setDefaultOpts(opts)
    defaults.photonsPerChunk  = 5e6;
    defaults.storeTcspcPix    = false;
    defaults.computePerFrame  = true;
    defaults.storePhotonLists = false;
    defaults.storeTimeCell    = false;
    defaults.showWaitbar      = false;
    defaults.maxFrames        = inf;
    defaults.maxNgate         = 1024;
    defaults.minLifetimeBin_ns = 0.128;

    f = fieldnames(defaults);
    for k = 1:numel(f)
        if ~isfield(opts, f{k}) || isempty(opts.(f{k}))
            opts.(f{k}) = defaults.(f{k});
        end
    end
end

function [StartUse, StopUse] = pairStartStop(Turns1, Turns2)
    Turns1 = Turns1(:);
    Turns2 = Turns2(:);

    maxPairs = min(numel(Turns1), numel(Turns2));
    StartUse = zeros(maxPairs, 1, class(Turns1));
    StopUse  = zeros(maxPairs, 1, class(Turns2));

    i = 1;
    j = 1;
    k = 0;

    while i <= numel(Turns1) && j <= numel(Turns2)
        if Turns2(j) <= Turns1(i)
            j = j + 1;
        else
            k = k + 1;
            StartUse(k) = Turns1(i);
            StopUse(k)  = Turns2(j);
            i = i + 1;
            j = j + 1;
        end
    end

    StartUse = StartUse(1:k);
    StopUse  = StopUse(1:k);
end

function updateWaitbar(h, frameIdx, nzAlloc)
    if ~isempty(h) && ishandle(h)
        if mod(frameIdx, 10) == 0 || frameIdx == 2
            waitbar(min(frameIdx-1, nzAlloc) / max(nzAlloc,1), h, ...
                sprintf('Frame %d / %d', min(frameIdx, nzAlloc), nzAlloc));
            drawnow limitrate
        end
    end
end

function [cnt2, sum12, sum22] = reduceMomentsAllChannels(im_col, im_line, im_chan, im_tcspc, nx, ny, nPix, nCh, chMapLUT)
    % Consistent linear index for output shape [nx, ny, nCh]
    % x = column, y = line

    if isempty(im_col)
        cnt2  = zeros(nPix, nCh, 'double');
        sum12 = zeros(nPix, nCh, 'double');
        sum22 = zeros(nPix, nCh, 'double');
        return;
    end

    x = double(im_col(:));
    y = double(im_line(:));
    t = double(im_tcspc(:));

    chLocal = double(chMapLUT(double(im_chan(:)) + 1));
    chLocal = chLocal(:);
    x = x(:);
    y = y(:);
    t = t(:);

    valid = (chLocal >= 1) & (chLocal <= nCh) & (x >= 1) & (x <= nx) & (y >= 1) & (y <= ny);

    x = x(valid);
    y = y(valid);
    t = t(valid);
    chLocal = chLocal(valid);

    pixIdx = x + (y - 1) * nx;
    grpIdx = pixIdx + (chLocal - 1) * nPix;

    sz = [nPix * nCh, 1];
    cnt  = accumarray(grpIdx, 1,    sz, @sum, 0);
    sum1 = accumarray(grpIdx, t,    sz, @sum, 0);
    sum2 = accumarray(grpIdx, t.^2, sz, @sum, 0);

    cnt2  = reshape(cnt,  [nPix, nCh]);
    sum12 = reshape(sum1, [nPix, nCh]);
    sum22 = reshape(sum2, [nPix, nCh]);
end

function tauOut = momentsToTau(cnt2, sum12, sum22, Resolution, nx, ny, nCh)
    valid = cnt2 > 0;

    m1 = zeros(size(cnt2), 'double');
    m2 = zeros(size(cnt2), 'double');

    m1(valid) = (sum12(valid) ./ cnt2(valid)) * Resolution;
    m2(valid) = (sum22(valid) ./ cnt2(valid)) * (Resolution^2);

    varT = max(m2 - m1.^2, 0);
    tauOut = reshape(single(sqrt(varT)), [nx, ny, nCh]);
end

function tcspc_pix = updateTcspcCube(tcspc_pix, im_col, im_line, im_chan, im_tcspc, nx, nPix, nCh, chMapLUT, Ngate)
    if isempty(im_col)
        return;
    end

    x = double(im_col(:));
    y = double(im_line(:));
    g = double(im_tcspc(:));
    chLocal = double(chMapLUT(double(im_chan(:)) + 1));

    x = x(:);
    y = y(:);
    g = g(:);
    chLocal = chLocal(:);

    valid = (chLocal >= 1) & (chLocal <= nCh) & (x >= 1) & (x <= size(tcspc_pix,1)) & (y >= 1) & (y <= size(tcspc_pix,2)) & (g >= 1) & (g <= Ngate);
    x = x(valid);
    y = y(valid);
    g = g(valid);
    chLocal = chLocal(valid);

    pixIdx = x + (y - 1) * nx;
    grp3 = pixIdx + (g - 1) * nPix + (chLocal - 1) * (nPix * Ngate);

    counts = accumarray(grp3, 1, [nPix * Ngate * nCh, 1], @sum, 0);
    counts = reshape(counts, [nx, size(tcspc_pix,2), Ngate, nCh]);

    if max(counts(:)) > double(intmax('uint16'))
        warning('tcspc_pix exceeds uint16 range; clipping to uint16 max.');
        counts = min(counts, double(intmax('uint16')));
    end

    tcspc_pix = tcspc_pix + uint16(counts);
end

