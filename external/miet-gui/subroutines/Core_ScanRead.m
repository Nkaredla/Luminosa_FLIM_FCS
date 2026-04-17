function [head, im_sync, im_tcspc, im_chan, im_line, im_col] = Core_ScanRead(name, deadtime)
% Core_ScanRead reads the .ht3 file of an image scan
% im_sync is a vector containing the sync number for each photon
% im_tcspc is a vector containing the time channel of each photon
% im_chan is a vector where the values represent the detector number of
% each photon
% im_line is a vector containing the line number of all the photons
% im_col is a vector containing the column number of all the photons

if nargin<2 || isempty(deadtime)
    deadtime = 100; % by default, deadtime is 100ns
end
nRemovedPhotons = 0;

if strcmp(name(end-2:end),'ht3')
    
    head = HT3_Read(name);
    
    if ~isempty(head)
        
        nx         = head.ImgHdr.PixX;                                              % x-pixels
        ny         = head.ImgHdr.PixY;                                              % y-pixels
        %             dt         = 1e-9*nx*head.ImgHdr.PixelTime*head.SyncRate;                   % number of pulses in one line
        xcord      = head.ImgHdr.X0+(1:nx)*head.ImgHdr.PixelSize;
        ycord      = head.ImgHdr.Y0+(1:ny)*head.ImgHdr.PixelSize;
        
        head.ImgHdr.Xcord = xcord;
        head.ImgHdr.Ycord = ycord;
        
        anzch      = 32;                                                    % max number of detectors
        Resolution = max([head.Resolution]);                                % Resolution of measurement
        chDiv      = Resolution/head.Resolution;                            % 1 in this case
        Ngate      = ceil(1e9/head.SyncRate./Resolution)+1;                 % TCPSC channels
        Timeunit   = 1./head.SyncRate;
        
        y        = [];
        tmpx     = [];
        chan     = [];
        markers  = [];
        line_chanind = [];
        im_sync  = zeros(round(head.Records),1);
        im_tcspc = [];
        im_chan  = [];
        im_line  = [];
        im_col   = [];
        Turns    = [];
        
        cnt      = 0;
        tend     = 0;
        line     = 1;
        
        h.waitbar = waitbar(0,'Please wait ...');
        
        fprintf('\n\n');
        
        c_cnt = 0;
        if head.ImgHdr.Pattern == 0 % Unidirectional scan
            
            [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);
            % Reading 5 million counts at a time
            while (num>0)
                
                cnt = cnt + num;
                if ~isempty(y)
                    tmpy = tmpy+tend;
                end;
                
                ind = (tmpmarkers>0)|((tmpchan<anzch)&(tmptcspc<Ngate*chDiv));
                
                y       = [y; tmpy(ind)];                         %#ok<AGROW>
                tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1;]; %#ok<AGROW>
                chan    = [chan; tmpchan(ind)+1];                 %#ok<AGROW>
                markers = [markers; tmpmarkers(ind)];             %#ok<AGROW>
                
                if isempty(line_chanind)
                    lmark = unique(chan(markers == 1),'legacy');
                    lmark(lmark == 9) = [];                  % The virtual photons where channel = 8 is triggered not by the piezo stage but by some external source
                    if numel(lmark)>1
                        ldiff   = zeros(numel(lmark),1);
                        numchan = lmark;
                        for idx = 1:numel(lmark)
                            numchan(idx) = sum(markers == 1 & chan == lmark(idx));
                            if numchan(idx)>1
                                ldiff(idx) = round(mean(diff(y(markers == 1 & chan == lmark(idx)))));
                            else
                                ldiff(idx) = 0;
                            end
                        end
                        chanratio = round((max(y)-min(y))./ldiff./numchan);
                        line_chanind = lmark(chanratio == 1);
                    else
                        line_chanind = lmark;
                    end
                end
                ind = markers==1 & chan==line_chanind;
                
                Turns = [Turns; y(ind)]; %#ok<AGROW>
                y(ind) = [];
                tmpx(ind) = [];
                chan(ind) = [];
                markers(ind) = [];
                
                tend  = y(end)+loc;
                %                 if deadtime>0
                %                     dind = unique(chan);
                %                     sync = y*Timeunit + tmpx*Resolution; % Note: this is a pseudo time, overcounts from previous reads are not taken into account.
                %                     flvp = false(numel(chan),numel(dind));
                %                     for ii = 1:numel(dind)
                %                         flvp(:,ii) = chan==dind(ii);
                %                     end
                %                     idx = any(flvp & circshift(flvp,[1 0]),2);
                %                     dsync=diff([0;sync]);
                %                     idx = idx & (dsync<deadtime);
                %                     tmpx(idx)       = [];
                %                     chan(idx)       = [];
                %                     y(idx)          = [];
                %                     markers(idx)    = [];
                %                     nRemovedPhotons = nRemovedPhotons + sum(idx);
                %                 end
                
                if numel(Turns)>1
                    for j=1:numel(Turns)-1
                        
                        dT = (Turns(2)-Turns(1));                           % Syncs in each line
                        t1 = Turns(1)+head.ImgHdr.TStartTo*dT;              % sync number for the begin of a to- scan
                        t2 = Turns(1)+head.ImgHdr.TStopTo*dT;               % sync number for the end of a to- scan
                        if ~isfield(head,'PixelDwellTime')
                            head.PixelDwellTime = (t2-t1)/nx/head.SyncRate;
                        end
                        ind = (markers~=0)|(y<t1);
                        y(ind)       = [];
                        tmpx(ind)    = [];
                        chan(ind)    = [];
                        markers(ind) = [];
                        
                        ind = (y>=t1)&(y<=t2);
                        
                        im_sync(c_cnt+1:c_cnt+sum(ind))   = y(ind);                        % % sync number of photons
                        im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW> % tcspc number of photons
                        im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW> % detector number of photons
                        im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW> % line number of photons
                        im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW> % separation of photons into pixels based on sync numbers
                        
                        c_cnt = c_cnt+sum(ind);
                        line = line +1;
                        waitbar(line/ny);
                        drawnow
                        
                        Turns(1) = [];
                        
                    end
                end
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);
                
            end;
            
            t1 = Turns(end)+head.ImgHdr.TStartTo*dT;
            t2 = Turns(end)+head.ImgHdr.TStopTo*dT;
            
            ind          = (y<t1);
            y(ind)       = [];
            tmpx(ind)    = [];
            chan(ind)    = [];
            
            ind = (y>=t1)&(y<=t2);
            
            %             im_sync   = [im_sync; y(ind)];
            im_sync(c_cnt+1:c_cnt+sum(ind))   = y(ind);
            im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
            im_chan   = [im_chan; uint8(chan(ind))];
            im_line   = [im_line; uint16(line.*ones(sum(ind),1))];
            im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];
            
            c_cnt = c_cnt+sum(ind);
            line = line +1;
            h.waitbar = waitbar(line/ny);
            drawnow
            
        else  % bidirectional scan
            
            [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);
            % tmpy            = sync value of each photon detected
            % tmptcspc        = tcspc channel of each photon detected
            % tmpchan         = detector channel for each detected photon
            % tmpmarkers      = marker indicating a virtual photon or
            %                   event. In this case, when the piezo shifts to the next line it is 1
            % num             = number of photons read actually.
            
            % tmpchan also has a different value when there is a line shift
            while (num>0)
                
                cnt = cnt + num;
                if ~isempty(y)
                    tmpy = tmpy+tend;
                end;
                
                ind = (tmpmarkers>0)|((tmpchan<anzch)&(tmptcspc<=Ngate*chDiv)); % indicator for photons
                
                y       = [y; tmpy(ind)];                         %#ok<AGROW>
                tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1]; %#ok<AGROW>
                chan    = [chan; tmpchan(ind)+1];                 %#ok<AGROW>
                markers = [markers; tmpmarkers(ind)];             %#ok<AGROW>
                if isempty(line_chanind)
                    lmark = unique(chan(markers == 1),'legacy');
                    lmark(lmark == 9) = [];                  % The virtual photons where channel = 8 is triggered not by the piezo stage but by some external source
                    if numel(lmark)>1
                        ldiff   = zeros(numel(lmark),1);
                        numchan = lmark;
                        for idx = 1:numel(lmark)
                            numchan(idx) = sum(markers == 1 & chan == lmark(idx));
                            if numchan(idx)>1
                                ldiff(idx) = round(mean(diff(y(markers == 1 & chan == lmark(idx)))));
                            else
                                ldiff(idx) = 0;
                            end
                        end
                        chanratio = round((max(y)-min(y))./ldiff./numchan);
                        line_chanind = lmark(chanratio == 1);
                    else
                        line_chanind = lmark;
                    end
                end
                
                ind = markers==1 & chan==line_chanind;
                
                Turns = [Turns; y(ind)]; %#ok<AGROW>
                y(ind) = [];
                tmpx(ind) = [];
                chan(ind) = [];
                markers(ind) = [];
                tend = y(end)+loc;
                if deadtime>0
                    dind = unique(chan,'legacy');
                    sync = y*Timeunit + tmpx*Resolution*1e-9; % Note: this is a pseudo time, overcounts from previous reads are not taken into account.
                    flvp = false(numel(chan),numel(dind));
                    for ii = 1:numel(dind)
                        flvp(:,ii) = chan==dind(ii);
                    end
                    idx = any(flvp & circshift(flvp,[1 0]),2);
                    dsync=diff([0;sync]);
                    idx = idx & (dsync<deadtime);
                    tmpx(idx)       = [];
                    chan(idx)       = [];
                    y(idx)          = [];
                    markers(idx)    = [];
                    nRemovedPhotons = nRemovedPhotons + sum(idx);
                end
                
                if numel(Turns)>2
                    for j=1:2:2*floor(numel(Turns)/2-1)
                        
                        dT = (Turns(2)-Turns(1));                           % Syncs in each line
                        t1 = Turns(1)+head.ImgHdr.TStartTo*dT;              % sync number for the begin of a to- scan
                        t2 = Turns(1)+head.ImgHdr.TStopTo*dT;               % sync number for the end of a to- scan
                        if ~isfield(head,'PixelDwellTime')
                            head.PixelDwellTime = (t2-t1)/nx/head.SyncRate;
                        end
                        ind = (markers~=0)|(y<t1);
                        y(ind)       = [];
                        tmpx(ind)    = [];
                        chan(ind)    = [];
                        markers(ind) = [];
                        
                        ind = (y>=t1)&(y<=t2);
                        
                        im_sync(c_cnt+1:c_cnt+sum(ind))   = y(ind);
                        %                         im_sync   = [im_sync; y(ind)];                        %#ok<AGROW> % sync number of photons
                        im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW> % sync number of photons
                        im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW> % detector number of photons
                        im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW> % line number of photons
                        im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW> % separation of photons into pixels based on sync numbers
                        
                        c_cnt = c_cnt+sum(ind);
                        line = line +1;
                        
                        t1 = Turns(1)+head.ImgHdr.TStartFro*dT;             % begin of fro- direction scan
                        t2 = Turns(1)+head.ImgHdr.TStopFro*dT;              % end of the fro- direction scan
                        
                        ind = (y<t1);
                        y(ind)       = [];
                        tmpx(ind)    = [];
                        chan(ind)    = [];
                        markers(ind) = [];
                        
                        ind = (y>=t1)&(y<=t2);
                        
                        %                         im_sync   = [im_sync; y(ind)];                        %#ok<AGROW> % sync number of photons
                        im_sync(c_cnt+1:c_cnt+sum(ind))   = y(ind);
                        im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW>
                        im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW>
                        im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                        im_col    = [im_col;  uint16(nx - floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                        
                        c_cnt = c_cnt+sum(ind);
                        line = line +1;
                        waitbar(line/ny);
                        drawnow
                        
                        Turns(1:2) = [];
                    end
                end
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);
                
            end;
            
            if ~isempty(Turns) % for the last one or two lines
                t1 = Turns(end-1)+head.ImgHdr.TStartTo*dT;
                t2 = Turns(end-1)+head.ImgHdr.TStopTo*dT;
                
                ind = (y<t1);
                y(ind)       = [];
                tmpx(ind)    = [];
                chan(ind)    = [];
                
                ind = (y>=t1)&(y<=t2);
                
                %                 im_sync   = [im_sync; y(ind)];
                im_sync(c_cnt+1:c_cnt+sum(ind))   = y(ind);
                im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                im_chan   = [im_chan; uint8(chan(ind))];
                im_line   = [im_line; uint16(line.*ones(sum(ind),1))];
                im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];
                
                c_cnt = c_cnt+sum(ind);
                line = line +1;
                
                t1 = Turns(end-1)+head.ImgHdr.TStartFro*dT;
                t2 = Turns(end-1)+head.ImgHdr.TStopFro*dT;
                
                ind = (y<t1);
                y(ind)       = [];
                tmpx(ind)    = [];
                chan(ind)    = [];
                
                ind = (y>=t1)&(y<=t2);
                
                im_sync(c_cnt+1:c_cnt+sum(ind))   = y(ind);
                %                 im_sync   = [im_sync; y(ind)];
                im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                im_chan   = [im_chan; uint8(chan(ind))];
                im_line   = [im_line; uint16(line.*ones(sum(ind),1))];
                im_col    = [im_col;  uint16(nx - floor(nx.*(y(ind)-t1)./(t2-t1)))];
                
                c_cnt = c_cnt+sum(ind);
                line = line +1;
                h.waitbar = waitbar(line/ny);
                drawnow
            end
        end;
        
        delete(h.waitbar);
        
    end
end
im_sync(im_sync==0)=[];
fprintf('Harp_tcspc removed %d of %d photons within dead time of %.2f ns\n',...
    nRemovedPhotons,numel(im_tcspc),deadtime*1e9);
if nargout<6
    %[head, im_sync, im_tcspc,  im_chan, im_line, im_col] = Core_ScanRead(name, deadtime)
    im_sync = im_tcspc;
    im_tcspc = im_chan;
    im_chan = im_line;
    im_line = im_col;
end
end
