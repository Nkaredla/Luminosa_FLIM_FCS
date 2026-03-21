function [head, im_sync, im_tcspc, im_chan, im_line, im_col, im_frame, num] = LPTU_LineScanRead(name, cnts, nx)
% LPTU_LineScanRead
% Reads photons from a PTU line-scan file.
%
% Input
% -----
% name : string
%     PTU filename
% cnts : scalar or 2-element vector
%     If scalar, number of records to read from the beginning.
%     If [start count], start record and number of records to read.
% nx : scalar
%     Number of pixels in x (default 50)
%
% Output
% ------
% head     : PTU header structure
% im_sync  : sync number of photon (global time)
% im_tcspc : TCSPC bin
% im_chan  : channel index
% im_line  : line index
% im_col   : pixel index along the line
% im_frame : empty (kept for compatibility)
% num      : number of records actually read

    if nargin < 3 || isempty(nx)
        nx = 50;
    end

    % Default outputs
    head     = [];
    im_sync  = [];
    im_tcspc = [];
    im_chan  = [];
    im_line  = [];
    im_col   = [];
    im_frame = [];
    num      = 0;

    % Check extension
    if numel(name) < 3 || ~strcmpi(name(end-2:end), 'ptu')
        return;
    end

    % Read header
    head = PTU_Read_Head(name);
    if isempty(head)
        return;
    end

    % Only monodirectional line scan
    if ~(head.ImgHdr_Ident == 9 && head.ImgHdr_BiDirect == 0)
        return;
    end

    % Marker definitions
    if ~isfield(head, 'ImgHdr_LineStart') || ~isfield(head, 'ImgHdr_LineStop')
        return;
    end

    LineStart = 2^(head.ImgHdr_LineStart - 1);
    LineStop  = 2^(head.ImgHdr_LineStop  - 1);

    % Original constants
    nz         = 1e3;
    anzch      = 32;
    Resolution = max([1e9 * head.MeasDesc_Resolution]);
    chDiv      = 1e-9 * Resolution / head.MeasDesc_Resolution;
    Ngate      = ceil(1e9 * head.MeasDesc_GlobalResolution ./ Resolution);

    % Read PTU records
    [tmpy, tmptcspc, tmpchan, tmpmarkers, num] = PTU_Read(name, cnts, head);

    % Keep markers + valid photons exactly as before
    ind = (tmpmarkers > 0) | ((tmpchan < anzch) & (tmptcspc < Ngate * chDiv));

    y       = tmpy(ind);
    tmpx    = floor(tmptcspc(ind) ./ chDiv) + 1;
    chan    = tmpchan(ind) + 1;
    markers = tmpmarkers(ind);

    % Extract raw line starts/stops
    if LineStart == LineStop
        tmpturns = y(markers == LineStart);
        Turns1   = tmpturns(1:2:end);
        Turns2   = tmpturns(2:2:end);
    else
        Turns1 = y(markers == LineStart);
        Turns2 = y(markers == LineStop);
    end

    % Remove markers once
    indMarker = (markers ~= 0);
    y(indMarker)    = [];
    tmpx(indMarker) = [];
    chan(indMarker) = [];

    % Make column vectors
    Turns1 = Turns1(:);
    Turns2 = Turns2(:);

    % No complete line markers
    if isempty(Turns1) || isempty(Turns2)
        head.ImgHdr_PixelTime = [];
        head.ImgHdr_DwellTime = [];
        return;
    end

    % ---------------------------------------------------------------
    % Robust temporal pairing:
    % pair each start with the next stop that occurs AFTER that start.
    %
    % This fixes chunked reads from the middle of a file, where the first
    % marker in the chunk may be a stop and/or the last marker may be a start.
    % ---------------------------------------------------------------
    maxPairs = min(numel(Turns1), numel(Turns2));
    StartUse = zeros(maxPairs, 1, class(Turns1));
    StopUse  = zeros(maxPairs, 1, class(Turns2));

    i = 1; % index into starts
    j = 1; % index into stops
    k = 0; % number of valid pairs

    while i <= numel(Turns1) && j <= numel(Turns2)
        if Turns2(j) <= Turns1(i)
            % stop occurs before current start -> stray stop in this chunk
            j = j + 1;
        else
            % valid pair: this start is closed by this stop
            k = k + 1;
            StartUse(k) = Turns1(i);
            StopUse(k)  = Turns2(j);
            i = i + 1;
            j = j + 1;
        end
    end

    StartUse = StartUse(1:k);
    StopUse  = StopUse(1:k);

    if isempty(StartUse)
        head.ImgHdr_PixelTime = [];
        head.ImgHdr_DwellTime = [];
        return;
    end

    nLines = numel(StartUse);

    % Restrict photons to overall valid time span
    ind = (y >= StartUse(1)) & (y <= StopUse(end));
    y    = y(ind);
    tmpx = tmpx(ind);
    chan = chan(ind);

    if isempty(y)
        dt = zeros(nz, 1);
        dt(1:nLines) = StopUse - StartUse;
        head.ImgHdr_PixelTime = 1e9 .* mean(dt) ./ nx ./ head.TTResult_SyncRate;
        head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
        return;
    end

    % Assign each photon to the interval beginning at the most recent start
    lineIdx = discretize(double(y), [double(StartUse); inf]);

    % Keep only photons that also lie before that line's stop
    ok = ~isnan(lineIdx) & (y <= StopUse(lineIdx));

    y       = y(ok);
    tmpx    = tmpx(ok);
    chan    = chan(ok);
    lineIdx = lineIdx(ok);

    % Build outputs
    im_sync  = y(:);
    im_tcspc = uint16(tmpx(:));
    im_chan  = uint8(chan(:));
    im_line  = uint16(lineIdx(:));

    t1v  = StartUse(lineIdx(:));
    t2v  = StopUse(lineIdx(:));
    frac = (double(y(:)) - double(t1v)) ./ double(t2v - t1v);

    % Keep your current safe clamping
    col = 1 + floor(nx .* frac);
    col = min(nx, max(1, col));
    im_col = uint16(col);

    % Reproduce original dt behavior
    dt = zeros(nz, 1);
    dt(1:nLines) = StopUse - StartUse;

    head.ImgHdr_PixelTime = 1e9 .* mean(dt) ./ nx ./ head.TTResult_SyncRate;
    head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
    im_frame = [];
end