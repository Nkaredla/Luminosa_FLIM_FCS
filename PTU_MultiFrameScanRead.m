function out = PTU_MultiFrameScanRead(name, photonsPerChunk)
% PTU_MultiFrameScanRead_Fast
%
% Fast reader for PTU MultiFrame Scan, monodirectional case:
%   head.ImgHdr_Ident == 9
%   head.ImgHdr_BiDirect == 0
%
% Output structure fields:
%   .head
%   .tag, .tau
%   .tags, .taus
%   .time
%   .tcspc_pix
%   .im_sync, .im_tcspc, .im_chan, .im_line, .im_col, .im_frame
%
% Notes
% -----
% This version is designed to preserve the logic of the original code,
% but make it much faster and more robust for chunked reading.

if nargin < 2 || isempty(photonsPerChunk)
    photonsPerChunk = 1e6;
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
nx = head.ImgHdr_PixX;
ny = head.ImgHdr_PixY;

if isfield(head, 'ImgHdr_MaxFrames')
    nz = head.ImgHdr_MaxFrames;
else
    tim_p_frame = head.ImgHdr_TimePerPixel*nx * ny; %ms
    tot_time = head.TTResult_StopAfter; % ms
    nz = ceil(tot_time / tim_p_frame);
end

% ---------------------------------------------------------------------
% Channel and TCSPC settings
% ---------------------------------------------------------------------
anzch = 32;

Resolution = max([1e9 * head.MeasDesc_Resolution, 0.128]); % ns
chDiv = 1e-9 * Resolution / head.MeasDesc_Resolution;
Ngate = ceil(1e9 * head.MeasDesc_GlobalResolution ./ Resolution) + 1;
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
[~, ~, tmpchan0, tmpmarkers0] = PTU_Read(name, [1 1e4], head);
dind = unique(tmpchan0(tmpmarkers0 == 0));
maxch_n = numel(dind);

% ---------------------------------------------------------------------
% Preallocate image outputs
% ---------------------------------------------------------------------
tag = zeros(nx, ny, maxch_n, nz, 'double');
tau = zeros(nx, ny, maxch_n, nz, 'double');

% We'll collect photon-resolved outputs in cells, concatenate once.
syncCell  = {};
tcspcCell = {};
chanCell  = {};
lineCell  = {};
colCell   = {};
frameCell = {};
timeCellAll = cell(maxch_n, 1);

% For dwell-time estimation
dt_all = [];

% ---------------------------------------------------------------------
% Carry buffers between chunks
% ---------------------------------------------------------------------
yCarry       = [];
tcspcCarry   = [];
chanCarry    = [];
markerCarry  = [];

turnsStart = [];
turnsStop  = [];

% Absolute offset for chunk stitching
cnt = 0;
tend = 0;
num = 1;

frameIdx = 1;

nphotTot = head.TTResult_NumberOfRecords;

h = waitbar(0, sprintf('Frame %d / %d', frameIdx, nz));

% ---------------------------------------------------------------------
% Main chunk loop
% ---------------------------------------------------------------------
while num

    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photonsPerChunk], head);

    cnt = cnt + num;

    if ~isempty(yCarry)
        tmpy = tmpy + tend;
    end

    % Keep markers + valid photons
    indValid = (tmpmarkers > 0) | ((tmpchan < anzch) & (tmptcspc < Ngate * chDiv));

    yNew       = tmpy(indValid);
    tcspcNew   = floor(tmptcspc(indValid) ./ chDiv) + 1;
    chanNew    = uint8(tmpchan(indValid));
    markerNew  = uint8(tmpmarkers(indValid));

    % Append to carry
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

        yf      = yCarry(indFramePhot);
        tf      = tcspcCarry(indFramePhot);
        chf     = chanCarry(indFramePhot);

        yCarry(indFramePhot)     = [];
        tcspcCarry(indFramePhot) = [];
        chanCarry(indFramePhot)  = [];
        markerCarry(indFramePhot)= [];

        % Line markers belonging to this frame
        indStartFrame = turnsStart < frameEnd;
        indStopFrame  = turnsStop  < frameEnd;

        startsF = turnsStart(indStartFrame);
        stopsF  = turnsStop(indStopFrame);

        turnsStart(indStartFrame) = [];
        turnsStop(indStopFrame)   = [];

        % Robust start/stop pairing
        [startsUse, stopsUse] = pairStartStop(startsF, stopsF);

        if isempty(startsUse) || isempty(yf)
            frameIdx = frameIdx + 1;
            waitbar(min(frameIdx-1, nz)/nz, h, sprintf('Frame %d / %d', min(frameIdx, nz), nz));
            drawnow limitrate
            continue;
        end

        % Restrict photons to spanned interval
        indSpan = (yf >= startsUse(1)) & (yf <= stopsUse(end));
        yf  = yf(indSpan);
        tf  = tf(indSpan);
        chf = chf(indSpan);

        if isempty(yf)
            frameIdx = frameIdx + 1;
            waitbar(min(frameIdx-1, nz)/nz, h, sprintf('Frame %d / %d', min(frameIdx, nz), nz));
            drawnow limitrate
            continue;
        end

        % Assign each photon to most recent line start
        lineIdx = discretize(double(yf), [double(startsUse(:)); inf]);
        ok = ~isnan(lineIdx) & (yf <= stopsUse(lineIdx));

        yf      = yf(ok);
        tf      = tf(ok);
        chf     = chf(ok);
        lineIdx = lineIdx(ok);

        if isempty(yf)
            frameIdx = frameIdx + 1;
            waitbar(min(frameIdx-1, nz)/nz, h, sprintf('Frame %d / %d', min(frameIdx, nz), nz));
            drawnow limitrate
            continue;
        end

        % Pixel coordinate
        t1v = startsUse(lineIdx);
        t2v = stopsUse(lineIdx);
        frac = (double(yf) - double(t1v)) ./ double(t2v - t1v);

        col = 1 + floor(nx .* frac);
        col = min(nx, max(1, col));

        % Save photon-level outputs
        nHere = numel(yf);
        syncCell{end+1,1}  = yf(:);
        tcspcCell{end+1,1} = uint16(tf(:));
        chanCell{end+1,1}  = uint8(chf(:));
        lineCell{end+1,1}  = uint16(lineIdx(:));
        colCell{end+1,1}   = uint16(col(:));
        frameCell{end+1,1} = uint16(frameIdx * ones(nHere,1));

        % Dwell times from actual valid lines only
        dt_all = [dt_all; double(stopsUse(:) - startsUse(:))]; %#ok<AGROW>

        % Frame-wise Process_Frame call for compatibility
        [tmptag, tmptau] = Process_Frame( ...
            yf(:), uint16(col(:)), uint16(lineIdx(:)), uint8(chf(:)), uint16(tf(:)), head);

        if ~isempty(tmptag)
            tag(:,:,:,frameIdx) = tmptag;
            tau(:,:,:,frameIdx) = tmptau;
        end

        frameIdx = frameIdx + 1;
        waitbar(min(frameIdx-1, nz)/nz, h, sprintf('Frame %d / %d', min(frameIdx, nz), nz));
        drawnow limitrate
    end
end

if nz>frameIdx
    tag(:,:,:,frameIdx+1:end) = [];
    tau(:,:,:,frameIdx+1:end) = [];
end

close(h);

% ---------------------------------------------------------------------
% Concatenate photon outputs once
% ---------------------------------------------------------------------
if isempty(syncCell)
    im_sync  = zeros(0,1);
    im_tcspc = uint16([]);
    im_chan  = uint8([]);
    im_line  = uint16([]);
    im_col   = uint16([]);
    im_frame = uint16([]);
else
    im_sync  = vertcat(syncCell{:});
    im_tcspc = vertcat(tcspcCell{:});
    im_chan  = vertcat(chanCell{:});
    im_line  = vertcat(lineCell{:});
    im_col   = vertcat(colCell{:});
    im_frame = vertcat(frameCell{:});
end

% ---------------------------------------------------------------------
% Pixel time
% ---------------------------------------------------------------------
if ~isempty(dt_all)
    head.ImgHdr_PixelTime = 1e9 .* mean(dt_all) ./ nx ./ head.TTResult_SyncRate;
    head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
else
    head.ImgHdr_PixelTime = [];
    head.ImgHdr_DwellTime = [];
end

% ---------------------------------------------------------------------
% Build tcspc_pix, tags, taus, time
% ---------------------------------------------------------------------
SyncRate = 1 ./ head.MeasDesc_GlobalResolution;

tcspc_pix = zeros(nx, ny, Ngate, maxch_n, 'uint16');
tags      = zeros(nx, ny, maxch_n, 'double');
taus      = zeros(nx, ny, maxch_n, 'double');
time      = cell(maxch_n, 1);

bin = permute(repmat((1:Ngate)', [1 nx, ny]), [2 3 1]) ;

for ch = 1:maxch_n
    ind = (im_chan == dind(ch));

    if any(ind)
        tcspc_pix(:,:,:,ch) = uint16(mHist3( ...
            double(im_line(ind)), ...
            double(im_col(ind)), ...
            double(im_tcspc(ind)), ...
            1:nx, 1:ny, 1:Ngate));

        time{ch} = round(im_sync(ind) * 1 / SyncRate / Resolution / 1e-9) + double(im_tcspc(ind));
        tags(:,:,ch) = double(sum(tcspc_pix(:,:,:,ch), 3));

        denom = tags(:,:,ch);
        denom(denom == 0) = NaN;

        m1 = double(sum(uint16(bin) .* tcspc_pix(:,:,:,ch), 3))*Resolution ./ denom;
        m2 = double(sum(uint16(bin).^2 .* tcspc_pix(:,:,:,ch), 3))*Resolution^2 ./ denom;
        taus(:,:,ch) = real(sqrt(max(m2 - m1.^2, 0)));
        tmpTau = taus(:,:,ch);
        tmpTau(isnan(tmpTau)) = 0;
        taus(:,:,ch) = tmpTau;
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
% Pair each start with the next stop that occurs after it.

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

function [tag, tau_rms, tcspc_pix, time] = Process_Frame(im_sync, im_col, im_line, im_chan, im_tcspc, head)

Resolution = max([head.MeasDesc_Resolution*1e9, 0.256]); % ns
chDiv      = max(1, round(1e-9 * Resolution / head.MeasDesc_Resolution));
SyncRate   = 1 / head.MeasDesc_GlobalResolution;

nx = head.ImgHdr_PixX;
ny = head.ImgHdr_PixY;

dind    = double(unique(im_chan(:)));
Ngate   = max(1, round(head.MeasDesc_GlobalResolution / head.MeasDesc_Resolution * ...
    (head.MeasDesc_Resolution / Resolution) * 1e9));
maxch_n = numel(dind);

tcspc_pix = zeros(nx, ny, Ngate, maxch_n, 'uint16');
time = cell(maxch_n, 1);
tag = zeros(nx, ny, maxch_n);
tau_rms = zeros(nx, ny, maxch_n);

bin = permute(repmat((1:Ngate)', [1 nx, ny]), [2,3,1]) ; % ns

for ch = 1:maxch_n
    ind = (im_chan == dind(ch));
    if ~any(ind)
        continue;
    end

    tcspc_bin = floor(double(im_tcspc(ind) - 1) / chDiv) + 1;
    tcspc_bin = min(Ngate, max(1, tcspc_bin));

    tcspc_pix(:,:,:,ch) = uint16(mHist3( ...
        double(im_line(ind)), ...
        double(im_col(ind)), ...
        double(tcspc_bin), ...
        1:nx, 1:ny, 1:Ngate));

    time{ch} = round(double(im_sync(ind)) / SyncRate / (Resolution * 1e-9)) + double(tcspc_bin);

    tag(:,:,ch) = double(sum(tcspc_pix(:,:,:,ch), 3));

    denom = tag(:,:,ch);
    valid = denom > 0;

    m1 = zeros(nx, ny);
    m2 = zeros(nx, ny);

    tmp1 = double(sum(uint16(bin) .* tcspc_pix(:,:,:,ch), 3))*Resolution;
    tmp2 = double(sum((uint16(bin.^2)) .* tcspc_pix(:,:,:,ch), 3))*Resolution^2;

    m1(valid) = tmp1(valid) ./ denom(valid);
    m2(valid) = tmp2(valid) ./ denom(valid);

    varT = m2 - m1.^2;
    varT = max(varT, 0);

    tau_rms(:,:,ch) = sqrt(varT);
end
end