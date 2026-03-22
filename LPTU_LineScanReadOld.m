function [head, im_sync, im_tcspc, im_chan, im_line, im_col, im_frame, num] = LPTU_LineScanReadOld(name, cnts, nx)
% reads the photons in a ptu file
%
% input:
%   name:     name of the ptu file
%   cnts:     photons to read
%             cnts(1): photon number to start from
%             cnts(2): number of photons to read
%
% output:
%   head:     header of the ptu file
%   im_sync:  sync number of photon (global time)
%   etc for the other outputs

if nargin < 3 || isempty(nx)
    nx = 50; % by default 50 pixels
end

if strcmp(name(end-2:end), 'ptu') % confirming that it is a .ptu file

    head = PTU_Read_Head(name); % Read PicoQuant Unified TTTR Files

    if ~isempty(head)

        if head.ImgHdr_Ident == 9 && head.ImgHdr_BiDirect == 0 % monodirectional line

            % [tmpsync, ~, tmpchan, tmpmarkers] = PTU_Read(name, [1 1e4], head);
            % dind = unique(tmpchan(~tmpmarkers));

            if isfield(head, 'ImgHdr_LineStart')
                LineStart = 2^(head.ImgHdr_LineStart - 1);
            end

            if isfield(head, 'ImgHdr_LineStop')
                LineStop = 2^(head.ImgHdr_LineStop - 1);
            end

            % SyncRate = head.TTResult_SyncRate; % laser frequency
            % TimeUnit = 1/SyncRate; % s
            % LineTimeSync = mode(diff(tmpsync(tmpmarkers == LineStart))); % time per line (double the time between start and stop)
            % LineFreq = 1/(LineTimeSync * TimeUnit); % scan frequency
            % TotMeasTime = head.TTResult_StopAfter; % total measurement time (ms)
            % PixTime = head.ImgHdr_TimePerPixel; % time per pixel (ms)

            nz         = 1e3;
            anzch      = 32;
            Resolution = max([1e9 * head.MeasDesc_Resolution]);
            chDiv      = 1e-9 * Resolution / head.MeasDesc_Resolution;
            Ngate      = ceil(1e9 * head.MeasDesc_GlobalResolution ./ Resolution);
            dt         = zeros(nz, 1);
            line       = 1;
            im_sync = [];
            im_tcspc = [];
            im_chan  = [];
            im_line  = [];
            im_col   = [];

            h = waitbar(0, 'Please wait ...');

            fprintf('\n\n');

            [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, cnts, head);

            % here while loop cut out
            ind = (tmpmarkers > 0) | ((tmpchan < anzch) & (tmptcspc < Ngate * chDiv));

            y       = [tmpy(ind)];                         %#ok<AGROW>
            tmpx    = [floor(tmptcspc(ind) ./ chDiv) + 1]; %#ok<AGROW>
            chan    = [tmpchan(ind) + 1];                 %#ok<AGROW>
            markers = [tmpmarkers(ind)];                  %#ok<AGROW>

            if LineStart == LineStop
                tmpturns = y(markers == LineStart);
                Turns1   = [tmpturns(1:2:end)]; %#ok<AGROW>
                Turns2   = [tmpturns(2:2:end)]; %#ok<AGROW>
            else
                Turns1 = [y(markers == LineStart)]; %#ok<AGROW>
                Turns2 = [y(markers == LineStop)];  %#ok<AGROW>
            end

            ind          = (markers ~= 0);
            y(ind)       = [];
            tmpx(ind)    = [];
            chan(ind)    = [];
            markers(ind) = [];

            if numel(Turns2) > 1
                for j = 1:numel(Turns2)-1

                    t1 = Turns1(1);
                    t2 = Turns2(1);

                    ind       = (y < t1);
                    y(ind)    = [];
                    tmpx(ind) = [];
                    chan(ind) = [];
                    markers(ind) = [];

                    ind = (y >= t1) & (y < t2);

                    im_sync  = [im_sync; y(ind)]; %#ok<*AGROW>
                    im_tcspc = [im_tcspc; uint16(tmpx(ind))];
                    im_chan  = [im_chan; uint8(chan(ind))];
                    im_line  = [im_line; uint16(line .* ones(sum(ind), 1))];
                    im_col   = [im_col; uint16(1 + floor(nx .* (y(ind) - t1) ./ (t2 - t1)))];

                    dt(line) = t2 - t1;
                    line     = line + 1;

                    waitbar(line / nz);
                    drawnow

                    Turns1(1) = [];
                    Turns2(1) = [];
                end
            end

            t1 = Turns1(end);
            t2 = Turns2(end);

            ind       = (y < t1);
            y(ind)    = [];
            tmpx(ind) = [];
            chan(ind) = [];

            ind = (y >= t1) & (y <= t2);

            im_sync  = [im_sync; y(ind)];
            im_tcspc = [im_tcspc; uint16(tmpx(ind))];
            im_chan  = [im_chan; uint8(chan(ind))];
            im_line  = [im_line; uint16(line .* ones(sum(ind), 1))];
            im_col   = [im_col; uint16(1 + floor(nx .* (y(ind) - t1) ./ (t2 - t1)))];

            dt(line) = t2 - t1;

            line = line + 1;
            waitbar(line / nz, h);
            drawnow

        end

        head.ImgHdr_PixelTime = 1e9 .* mean(dt) / nx / head.TTResult_SyncRate;
        head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
        im_frame = [];

        close(h);
    end
end
end