function out = PTU_MultiFrameScanReadFast(name, photonsPerChunk, storeTcspcPix, useGPU)
% PTU_MultiFrameScanRead
%
% Fast reader for PTU MultiFrame Scan, monodirectional case:
%   head.ImgHdr_Ident == 9
%   head.ImgHdr_BiDirect == 0
%
% Memory strategy
% ---------------
% - Store all photons
% - Store framewise tag/tau as SINGLE
% - Store global tags/taus as SINGLE
% - Store tcspc_pix only if explicitly requested
%
% NOTE:
% tcspc_pix can still be extremely large even as uint16.
% Default is storeTcspcPix = false.
%
% GPU notes:
% - GPU acceleration is only used for tag/tau calculation when tcspc_pix
%   is being built (storeTcspcPix = true). Histogramming is still on CPU.
%
% For time series where each time point is stored in a separate PTU,
% call this reader per file. If no frame markers are present, remaining
% photons are treated as a single frame.

if nargin < 2 || isempty(photonsPerChunk)
    photonsPerChunk = 1e6;
end
if nargin < 3 || isempty(storeTcspcPix)
    storeTcspcPix = false;
end
if nargin < 4 || isempty(useGPU)
    useGPU = false;
end
useGPU = useGPU && gpuIsAvailable();
if useGPU && ~storeTcspcPix
    warning('PTU_MultiFrameScanReadFast:GPU', ...
        'useGPU requested but storeTcspcPix = false. GPU path requires tcspc_pix; proceeding on CPU.');
    useGPU = false;
end

out = struct();

if numel(name) < 3 || ~strcmpi(name(end-2:end), 'ptu')
    error('Input file must be a .ptu file.');
end

head = PTU_Read_Head(name);
if isempty(head)
    error('Could not read PTU header.');
end

if ~(head.ImgHdr_Ident == 9)
    error('This reader is intended for ImgHdr_Ident == 9 (MultiFrame Scan).');
end

if ~(head.ImgHdr_BiDirect == 0)
    error('This fast version currently supports only monodirectional scan.');
end

% ---------------------------------------------------------------------
% Dimensions
% ---------------------------------------------------------------------
nx = double(head.ImgHdr_PixX);
ny = double(head.ImgHdr_PixY);

if isfield(head, 'ImgHdr_MaxFrames')
    nz = double(head.ImgHdr_MaxFrames);
else
    tim_p_frame = head.ImgHdr_TimePerPixel * nx * ny; % ms
    tot_time = head.TTResult_StopAfter;               % ms
    nz = ceil(tot_time / tim_p_frame);
end

% ---------------------------------------------------------------------
% Channel and TCSPC settings
% ---------------------------------------------------------------------
anzch = 32;

rawRes_s = head.MeasDesc_Resolution;
Resolution = max(1e9 * rawRes_s, 0.128);             % ns
chDiv = max(1, round((Resolution * 1e-9) / rawRes_s));
Ngate = min(1024, ceil(1e9 * head.MeasDesc_GlobalResolution / Resolution) + 1);

% Preserve native timing metadata for downstream TCSPC/IRF reconstruction.
head.MeasDesc_Resolution_Original = rawRes_s;
head.TCSPC_chDiv = chDiv;

% Keep compatibility with downstream code that expects coarse resolution
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

% ---------------------------------------------------------------------
% First short read just to infer active channels
% ---------------------------------------------------------------------
[~, ~, tmpchan0, tmpmarkers0] = PTU_Read(name, [1 min(1e4, head.TTResult_NumberOfRecords)], head);
dind = unique(tmpchan0(tmpmarkers0 == 0));
maxch_n = numel(dind);

% ---------------------------------------------------------------------
% Preallocate outputs
% NOTE: kept in current [nx, ny, ...] shape for compatibility with your code
% even though standard MATLAB image layout would usually be [ny, nx, ...].
% ---------------------------------------------------------------------
tag = zeros(nx, ny, maxch_n, nz, 'single');
tau = zeros(nx, ny, maxch_n, nz, 'single');

% collect photon outputs in cells, concatenate once
syncCell  = cell(nz, 1);
tcspcCell = cell(nz, 1);
chanCell  = cell(nz, 1);
lineCell  = cell(nz, 1);
colCell   = cell(nz, 1);
frameCell = cell(nz, 1);

dtCell = cell(nz, 1);

% ---------------------------------------------------------------------
% Carry buffers between chunks
% ---------------------------------------------------------------------
yCarry      = [];
tcspcCarry  = [];
chanCarry   = [];
markerCarry = [];

turnsStart = [];
turnsStop  = [];

cnt = 0;
tend = 0;
num = 1;
frameIdx = 1;

h = waitbar(0, sprintf('Frame %d / %d', frameIdx, nz));

% ---------------------------------------------------------------------
% Main chunk loop
% ---------------------------------------------------------------------
while num > 0

    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photonsPerChunk], head);
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
    yCarry      = [yCarry; yNew];
    tcspcCarry  = [tcspcCarry; uint16(tcspcNew)];
    chanCarry   = [chanCarry; chanNew];
    markerCarry = [markerCarry; markerNew];

    % Update line markers
    if LineStart == LineStop
        tmpturns = yCarry(markerCarry == LineStart);
        if numel(turnsStart) > numel(turnsStop)
            turnsStart = [turnsStart; tmpturns(2:2:end)];
            turnsStop  = [turnsStop;  tmpturns(1:2:end)];
        else
            turnsStart = [turnsStart; tmpturns(1:2:end)];
            turnsStop  = [turnsStop;  tmpturns(2:2:end)];
        end
    else
        turnsStart = [turnsStart; yCarry(markerCarry == LineStart)];
        turnsStop  = [turnsStop;  yCarry(markerCarry == LineStop)];
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

    % -------------------------------------------------------------
    % Process all complete frames now available
    % -------------------------------------------------------------
    for kf = 1:numel(frameChange)
        if frameIdx > nz
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
            if mod(frameIdx, 10) == 0 || frameIdx == 2
                waitbar(min(frameIdx-1, nz)/nz, h, sprintf('Frame %d / %d', min(frameIdx, nz), nz));
                drawnow limitrate
            end
            continue;
        end

        % Restrict photons to valid temporal span
        indSpan = (yf >= startsUse(1)) & (yf <= stopsUse(end));
        yf  = yf(indSpan);
        tf  = tf(indSpan);
        chf = chf(indSpan);

        if isempty(yf)
            frameIdx = frameIdx + 1;
            if mod(frameIdx, 10) == 0 || frameIdx == 2
                waitbar(min(frameIdx-1, nz)/nz, h, sprintf('Frame %d / %d', min(frameIdx, nz), nz));
                drawnow limitrate
            end
            continue;
        end

        % Assign photons to lines
        lineIdx = discretize(double(yf), [double(startsUse(:)); inf]);
        ok = ~isnan(lineIdx) & (yf <= stopsUse(lineIdx)) & (lineIdx >= 1) & (lineIdx <= ny);

        yf      = yf(ok);
        tf      = tf(ok);
        chf     = chf(ok);
        lineIdx = lineIdx(ok);

        if isempty(yf)
            frameIdx = frameIdx + 1;
            if mod(frameIdx, 10) == 0 || frameIdx == 2
                waitbar(min(frameIdx-1, nz)/nz, h, sprintf('Frame %d / %d', min(frameIdx, nz), nz));
                drawnow limitrate
            end
            continue;
        end

        % Pixel coordinate along the line
        t1v = startsUse(lineIdx);
        t2v = stopsUse(lineIdx);
        frac = (double(yf) - double(t1v)) ./ double(t2v - t1v);

        col = 1 + floor(nx .* frac);
        col = min(nx, max(1, col));

        % Save photon-level outputs
        syncCell{frameIdx}  = yf(:);
        tcspcCell{frameIdx} = uint16(tf(:));
        chanCell{frameIdx}  = uint8(chf(:));
        lineCell{frameIdx}  = uint16(lineIdx(:));
        colCell{frameIdx}   = uint16(col(:));
        frameCell{frameIdx} = uint16(frameIdx * ones(numel(yf),1));

        % Dwell times
        dtCell{frameIdx} = double(stopsUse(:) - startsUse(:));

        % Fast per-frame tag/tau without building tcspc cube
        [tmptag, tmptau] = Process_Frame_fast( ...
            uint16(col(:)), uint16(lineIdx(:)), uint8(chf(:)), uint16(tf(:)), ...
            nx, ny, dind, Resolution);

        tag(:,:,:,frameIdx) = tmptag;
        tau(:,:,:,frameIdx) = tmptau;

        frameIdx = frameIdx + 1;
        if mod(frameIdx, 10) == 0 || frameIdx == 2
            waitbar(min(frameIdx-1, nz)/nz, h, sprintf('Frame %d / %d', min(frameIdx, nz), nz));
            drawnow limitrate
        end
    end
end

if ishandle(h)
    close(h);
end

% Trim unused trailing frames
if frameIdx <= nz
    tag(:,:,:,frameIdx:end) = [];
    tau(:,:,:,frameIdx:end) = [];
end

% ---------------------------------------------------------------------
% Concatenate photon outputs once
% ---------------------------------------------------------------------
usedFrames = 1:(frameIdx-1);

if isempty(usedFrames) || isempty(syncCell{usedFrames(1)})
    im_sync  = zeros(0,1);
    im_tcspc = uint16([]);
    im_chan  = uint8([]);
    im_line  = uint16([]);
    im_col   = uint16([]);
    im_frame = uint16([]);
else
    im_sync  = vertcat(syncCell{usedFrames});
    im_tcspc = vertcat(tcspcCell{usedFrames});
    im_chan  = vertcat(chanCell{usedFrames});
    im_line  = vertcat(lineCell{usedFrames});
    im_col   = vertcat(colCell{usedFrames});
    im_frame = vertcat(frameCell{usedFrames});
end

% ---------------------------------------------------------------------
% Pixel time
% ---------------------------------------------------------------------
dtCell = dtCell(~cellfun('isempty', dtCell));
if ~isempty(dtCell)
    dt_all = vertcat(dtCell{:});
    head.ImgHdr_PixelTime = 1e9 .* mean(dt_all) ./ nx ./ head.TTResult_SyncRate;
    head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
else
    head.ImgHdr_PixelTime = [];
    head.ImgHdr_DwellTime = [];
end

% ---------------------------------------------------------------------
% Build global tags, taus, time
% tcspc_pix is optional because it is often too large.
% ---------------------------------------------------------------------
SyncRate = 1 ./ head.MeasDesc_GlobalResolution;

tags = zeros(nx, ny, maxch_n, 'single');
taus = zeros(nx, ny, maxch_n, 'single');
time = cell(maxch_n, 1);

if storeTcspcPix
    estBytes = nx * ny * Ngate * maxch_n * 2; % uint16
    fprintf('Allocating tcspc_pix as uint16: %.2f GB\n', estBytes / 2^30);
    tcspc_pix = zeros(nx, ny, Ngate, maxch_n, 'uint16');
else
    tcspc_pix = [];
end

for ch = 1:maxch_n
    ind = (im_chan == dind(ch));

    if any(ind)
        % tag/tau from moments only
        [tag_ch, tau_ch, tcspc_ch] = compute_channel_stats( ...
            im_line(ind), im_col(ind), im_tcspc(ind), nx, ny, Ngate, Resolution, storeTcspcPix, useGPU);

        tags(:,:,ch) = tag_ch;
        taus(:,:,ch) = tau_ch;

        if storeTcspcPix
            tcspc_pix(:,:,:,ch) = tcspc_ch;
        end

        time{ch} = round(double(im_sync(ind)) / SyncRate / (Resolution * 1e-9)) + double(im_tcspc(ind));
    else
        time{ch} = [];
    end
end

out.head = head;
out.tag = tag;
out.tau = tau;
out.tags = tags;
out.taus = taus;
out.time = time;
out.tcspc_pix = tcspc_pix;
out.im_sync = im_sync;
out.im_tcspc = im_tcspc;
out.im_chan = im_chan;
out.im_line = im_line;
out.im_col = im_col;
out.im_frame = im_frame;
out.dind = dind;
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


function [tag, tau_rms] = Process_Frame_fast(im_col, im_line, im_chan, im_tcspc, nx, ny, dind, Resolution)
% Fast framewise tag/tau computation without dense tcspc cube
% Output kept as single to save RAM.

maxch_n = numel(dind);
tag = zeros(nx, ny, maxch_n, 'single');
tau_rms = zeros(nx, ny, maxch_n, 'single');

for k = 1:maxch_n
    ind = (im_chan == dind(k));
    if ~any(ind)
        continue;
    end

    pixIdx = sub2ind([nx, ny], double(im_line(ind)), double(im_col(ind)));
    tbin = double(im_tcspc(ind));

    nPix = nx * ny;
    cnt  = accumarray(pixIdx, 1,        [nPix, 1], @sum, 0);
    sum1 = accumarray(pixIdx, tbin,     [nPix, 1], @sum, 0);
    sum2 = accumarray(pixIdx, tbin.^2,  [nPix, 1], @sum, 0);

    tag2 = reshape(single(cnt), [nx, ny]);

    valid = cnt > 0;
    m1 = zeros(nPix, 1, 'double');
    m2 = zeros(nPix, 1, 'double');

    m1(valid) = (sum1(valid) ./ cnt(valid)) * Resolution;
    m2(valid) = (sum2(valid) ./ cnt(valid)) * (Resolution^2);

    varT = max(m2 - m1.^2, 0);
    tau2 = reshape(single(sqrt(varT)), [nx, ny]);

    tag(:,:,k) = tag2;
    tau_rms(:,:,k) = tau2;
end
end


function [tag_ch, tau_ch, tcspc_ch] = compute_channel_stats(im_line, im_col, im_tcspc, nx, ny, Ngate, Resolution, storeTcspc, useGPU)
% Global per-channel stats from photons.
% tag_ch, tau_ch are single
% tcspc_ch is uint16 if requested, else []

if storeTcspc && useGPU
    % Build tcspc_pix (CPU), then compute moments on GPU
    linIdx = sub2ind([nx, ny, Ngate], double(im_line), double(im_col), double(im_tcspc));
    counts3 = accumarray(linIdx, 1, [nx * ny * Ngate, 1], @sum, 0);

    if max(counts3) > double(intmax('uint16'))
        warning('tcspc_pix exceeds uint16 range; clipping to uint16 max.');
        counts3 = min(counts3, double(intmax('uint16')));
    end

    tcspc_ch = reshape(uint16(counts3), [nx, ny, Ngate]);

    tcspc_g = gpuArray(double(tcspc_ch));
    tag_g = sum(tcspc_g, 3);
    bins = gpuArray(reshape(1:Ngate, 1, 1, []));
    tmp1 = sum(tcspc_g .* bins, 3) * Resolution;
    tmp2 = sum(tcspc_g .* (bins.^2), 3) * (Resolution^2);

    denom = max(tag_g, 1);
    mean1 = tmp1 ./ denom;
    varT = max(tmp2 ./ denom - mean1.^2, 0);

    tag_ch = gather(single(tag_g));
    tau_ch = gather(single(sqrt(varT)));
else
    pixIdx = sub2ind([nx, ny], double(im_line), double(im_col));
    tbin = double(im_tcspc);

    nPix = nx * ny;
    cnt  = accumarray(pixIdx, 1,       [nPix, 1], @sum, 0);
    sum1 = accumarray(pixIdx, tbin,    [nPix, 1], @sum, 0);
    sum2 = accumarray(pixIdx, tbin.^2, [nPix, 1], @sum, 0);

    tag_ch = reshape(single(cnt), [nx, ny]);

    valid = cnt > 0;
    m1 = zeros(nPix, 1, 'double');
    m2 = zeros(nPix, 1, 'double');

    m1(valid) = (sum1(valid) ./ cnt(valid)) * Resolution;
    m2(valid) = (sum2(valid) ./ cnt(valid)) * (Resolution^2);

    varT = max(m2 - m1.^2, 0);
    tau_ch = reshape(single(sqrt(varT)), [nx, ny]);

    if storeTcspc
        linIdx = sub2ind([nx, ny, Ngate], double(im_line), double(im_col), double(im_tcspc));
        counts3 = accumarray(linIdx, 1, [nx * ny * Ngate, 1], @sum, 0);

        if max(counts3) > double(intmax('uint16'))
            warning('tcspc_pix exceeds uint16 range; clipping to uint16 max.');
            counts3 = min(counts3, double(intmax('uint16')));
        end

        tcspc_ch = reshape(uint16(counts3), [nx, ny, Ngate]);
    else
        tcspc_ch = [];
    end
end
end

function ok = gpuIsAvailable()
    ok = false;
    try
        ok = gpuDeviceCount > 0;
    catch
        ok = false;
    end
end
