function [head, im_sync, im_tcspc, im_chan, im_line, im_col, im_frame] = PTU_ScanRead(name, plt)
% head is the header of the file and it contains the experimental details
% such as pixel size, resolution, dwell time, etc.
% im_sync is a vector containing the sync number for each photon
% im_tcspc is a vector containing the time channel of each photon
% im_chan is a vector where the values represent the detector number of
% each photon
% im_line is a vector containing the line number of all the photons
% im_col is a vector containing the column number of all the photons
% im_frame is a vector containing the frame number of all the photons

if (nargin>1)&&~isempty(plt)
    plt = 1;
else
    plt = 0;
end

photons = 1e6;
if strcmp(name(end-2:end),'ptu')
    
    head = PTU_Read_Head(name);
    
    if ~isempty(head)
        nx = head.ImgHdr_PixX;
        ny = head.ImgHdr_PixY;
        
        if (head.ImgHdr_Ident == 1)||(head.ImgHdr_Ident == 6)
            
            anzch      = 32;
            Resolution = max([1e9*head.MeasDesc_Resolution,0.016]);
            chDiv      = 1e-9*Resolution/head.MeasDesc_Resolution;
            Ngate      = ceil(1e9*head.MeasDesc_GlobalResolution./Resolution)+1;
            head.MeasDesc_Resolution = Resolution*1e-9;
            
            LineStart = 4;
            LineStop  = 2;
            
            if isfield(head,'ImgHdr_LineStart')
                LineStart = 2^(head.ImgHdr_LineStart-1);
            end
            if isfield(head,'ImgHdr_LineStop')
                LineStop = 2^(head.ImgHdr_LineStop-1);
            end
            
            y        = [];
            tmpx     = [];
            chan     = [];
            markers  = [];
            dt       = zeros(ny,1);
            
            im_sync  = [];
            im_tcspc = [];
            im_chan  = [];
            im_line  = [];
            im_col   = [];
            im_frame = [];
            Turns1   = [];
            Turns2   = [];
            
            cnt      = 0;
            tend     = 0;
            line     = 1;
            
            h = waitbar(0,'Please wait ...');
            
            fprintf('\n\n');
            
            if head.ImgHdr_BiDirect == 0
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 5e6], head);
                
                while (num>0)
                    
                    cnt = cnt + num;
                    if ~isempty(y)
                        tmpy = tmpy+tend;
                    end
                    
                    ind = (tmpmarkers>0)|((tmpchan<anzch)&(tmptcspc<Ngate*chDiv));
                    
                    y       = [y; tmpy(ind)];                         %#ok<AGROW>
                    tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1]; %#ok<AGROW>
                    chan    = [chan; tmpchan(ind)+1];                 %#ok<AGROW>
                    markers = [markers; tmpmarkers(ind)];             %#ok<AGROW>
                    
                    if LineStart==LineStop
                        tmpturns = y(markers==LineStart);
                        if numel(Turns1)>numel(Turns2)           % first turn is a LineStop
                            Turns1 = [Turns1; tmpturns(2:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(1:2:end)];      %#ok<AGROW>
                        else
                            Turns1 = [Turns1; tmpturns(1:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(2:2:end)];      %#ok<AGROW>
                        end
                    else
                        Turns1 = [Turns1; y(markers==LineStart)]; %#ok<AGROW>
                        Turns2 = [Turns2; y(markers==LineStop)];  %#ok<AGROW>
                    end
                    
                    ind          = (markers~=0);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    markers(ind) = [];
                    
                    tend  = y(end)+loc;
                    
                    if numel(Turns2)>1
                        for j=1:numel(Turns2)-1
                            
                            t1 = Turns1(1);
                            t2 = Turns2(1);
                            
                            ind          = (y<t1);
                            y(ind)       = [];
                            tmpx(ind)    = [];
                            chan(ind)    = [];
                            markers(ind) = [];
                            
                            ind = (y>=t1)&(y<=t2);
                            
                            im_sync   = [im_sync; y(ind)]; %#ok<AGROW>
                            im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW>
                            im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW>
                            im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                            im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                            
                            dt(line)  = t2-t1;
                            line = line +1;
                            waitbar(line/ny);
                            drawnow
                            
                            Turns1(1) = [];
                            Turns2(1) = [];
                        end
                    end
                    
                    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 5e6], head);
                    
                end
                
                t1 = Turns1(end);
                t2 = Turns2(end);
                
                ind          = (y<t1);
                y(ind)       = [];
                tmpx(ind)    = [];
                chan(ind)    = [];
                
                ind = (y>=t1)&(y<=t2);
                
                im_sync   = [im_sync; y(ind)];
                im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                im_chan   = [im_chan; uint8(chan(ind))];
                im_line   = [im_line; uint16(line.*ones(sum(ind),1))];
                im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];
                dt(line)  = t2-t1;
                
                line = line +1;
                waitbar(line/ny,h);
                drawnow
                
            else  % bidirectional scan
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 5e6], head);
                
                while (num>0)
                    
                    cnt = cnt + num;
                    if ~isempty(y)
                        tmpy = tmpy+tend;
                    end
                    
                    ind = ((tmpchan<anzch)&(tmptcspc<=Ngate*chDiv));
                    
                    y       = [y; tmpy(ind)];                         %#ok<AGROW>
                    tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1;]; %#ok<AGROW>
                    chan    = [chan; tmpchan(ind)];                   %#ok<AGROW>
                    markers = [markers; tmpmarkers(ind)];             %#ok<AGROW>
                    
                    if LineStart==LineStop
                        tmpturns = y(markers==LineStart);
                        if numel(Turns1)>numel(Turns2)           % first turn is a LineStop
                            Turns1 = [Turns1; tmpturns(2:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(1:2:end)];      %#ok<AGROW>
                        else
                            Turns1 = [Turns1; tmpturns(1:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(2:2:end)];      %#ok<AGROW>
                        end
                    else
                        Turns1 = [Turns1; y(markers==LineStart)]; %#ok<AGROW>
                        Turns2 = [Turns2; y(markers==LineStop)];  %#ok<AGROW>
                    end
                    
                    ind          = (markers~=0);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    markers(ind) = [];
                    
                    tend = y(end)+loc;
                    
                    if numel(Turns2)>2
                        for j=1:2:2*floor(numel(Turns2)/2-1)
                            
                            t1 = Turns1(1);
                            t2 = Turns2(1);
                            
                            ind          = (y<t1);
                            y(ind)       = [];
                            tmpx(ind)    = [];
                            chan(ind)    = [];
                            markers(ind) = [];
                            
                            ind = (y>=t1)&(y<=t2);
                            
                            im_sync   = [im_sync; y(ind)];  %#ok<AGROW>
                            im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW>
                            im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW>
                            im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                            im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                            dt(line)  = t2-t1;
                            
                            line = line +1;
                            
                            t1 = Turns1(2);
                            t2 = Turns2(2);
                            
                            ind = (y<t1);
                            y(ind)       = [];
                            tmpx(ind)    = [];
                            chan(ind)    = [];
                            markers(ind) = [];
                            
                            ind = (y>=t1)&(y<=t2);
                            
                            im_sync   = [im_sync; y(ind)];  %#ok<AGROW>
                            im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW>
                            im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW>
                            im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                            im_col    = [im_col;  uint16(nx - floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                            dt(line)  = t2-t1;
                            
                            line = line +1;
                            waitbar(line/ny,h);
                            drawnow
                            
                            Turns1(1:2) = [];
                            Turns2(1:2) = [];
                        end
                    end
                    
                    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 5e6], head);
                    
                end
                
                if ~isempty(Turns2)
                    t1 = Turns1(end-1);
                    t2 = Turns2(end-1);
                    
                    ind = (y<t1);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    
                    ind = (y>=t1)&(y<=t2);
                    
                    im_sync   = [im_sync; y(ind)];
                    im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                    im_chan   = [im_chan; uint8(chan(ind))];
                    im_line   = [im_line; uint16(line.*ones(sum(ind),1))];
                    im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];
                    dt(line)  = t2-t1;
                    
                    line = line +1;
                    
                    t1 = Turns1(end);
                    t2 = Turns2(end);
                    
                    ind = (y<t1);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    
                    ind = (y>=t1)&(y<=t2);
                    
                    im_sync   = [im_sync; y(ind)];
                    im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                    im_chan   = [im_chan; uint8(chan(ind))];
                    im_line   = [im_line; uint16(line.*ones(sum(ind),1))];
                    im_col    = [im_col;  uint16(nx - floor(nx.*(y(ind)-t1)./(t2-t1)))];
                    dt(line)  = t2-t1;
                    
                    line = line +1;
                    waitbar(line/ny);
                    drawnow
                end
            end
            
            head.ImgHdr_PixelTime = 1e9.*mean(dt)/nx/head.TTResult_SyncRate;
            head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
            
            close(h);
            
        elseif head.ImgHdr_Ident == 3
            
            y        = [];
            tmpx     = [];
            chan     = [];
            marker   = [];
            
            dt       = zeros(ny,1);
            im_tcspc = [];
            im_chan  = [];
            im_line  = [];
            im_col   = [];
            im_frame = [];
            
            cnt      = 0;
            tend     = 0;
            line     = 1;
            n_frames = 0;
            f_times  = [];
            
            head.ImgHdr_X0       = 0;
            head.ImgHdr_Y0       = 0;
            head.ImgHdr_PixResol = 1;
            
            LineStart = 2^(head.ImgHdr_LineStart-1);
            LineStop  = 2^(head.ImgHdr_LineStop-1);
            Frame     = 2^(head.ImgHdr_Frame-1);
            
            h = waitbar(0, 'Please Wait...');
            
            if Frame < 1
                Frame = -1;
            end
            
            in_frame = false;
            
            if Frame < 1
                in_frame = true;
                n_frames = n_frames + 1;
            end
            
            [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 1e6], head);
            
            while (num>0)
                
                t_sync   = [];
                t_tcspc  = [];
                t_chan   = [];
                t_line   = [];
                t_col    = [];
                
                cnt = cnt + num;
                
                waitbar(cnt/head.TTResult_NumberOfRecords,h);
                
                tmpy = tmpy+tend;
                
                y       = [y; tmpy];                       %#ok<AGROW>
                tmpx    = [tmpx; tmptcspc];                %#ok<AGROW>
                chan    = [chan; tmpchan];                 %#ok<AGROW>
                marker  = [marker; tmpmarkers];            %#ok<AGROW>
                tend    = y(end)+loc;
                
                F  = y(bitand(marker,Frame)>0);
                
                while ~isempty(F)
                    
                    if ~in_frame
                        ind = (y<=F(1));
                        y(ind)       = [];
                        tmpx(ind)    = [];
                        chan(ind)    = [];
                        marker(ind)  = [];
                        line         = 1;
                        in_frame     = true;
                        n_frames     = n_frames + 1;
                        f_times      = [f_times; F(1)];
                        F(1)         = [];
                    end
                    
                    if ~isempty(F)
                        ind = y<F(1);
                        
                        f_y  = y(ind);
                        f_x  = tmpx(ind);
                        f_ch = chan(ind);
                        f_m  = marker(ind);
                        
                        y(ind)      = [];
                        tmpx(ind)   = [];
                        chan(ind)   = [];
                        marker(ind) = [];
                        
                    end
                    
                    L1 = f_y(bitand(f_m,LineStart)>0);
                    L2 = f_y(bitand(f_m,LineStop)>0);
                    
                    ll = line + numel(L2)-1; % this will be the last complete line in the data stack
                    
                    if ll > ny
                        L1 = L1(1:ny-line+1);
                        L2 = L2(1:ny-line+1);
                    end
                    
                    if numel(L1)>1
                        for j=1:numel(L2)
                            
                            ind = (f_y>L1(j))&(f_y<L2(j));
                            
                            t_sync   = [t_sync; y(ind)];
                            t_tcspc  = [t_tcspc; uint16(f_x(ind))];              %#ok<AGROW>
                            t_chan   = [t_chan; uint8(f_ch(ind))];                %#ok<AGROW>
                            t_line   = [t_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                            t_col    = [t_col;  uint16(1 + floor(nx.*(f_y(ind)-L1(j))./(L2(j)-L1(j))))];  %#ok<AGROW>
                            dt(line) = dt(line) + (L2(j)-L1(j));
                            line = line +1;
                        end
                    end
                    
                    if line>ny
                        in_frame = false;
                    end
                end
                
                im_sync   = [im_sync;  t_sync];
                im_tcspc  = [im_tcspc; t_tcspc];  %#ok<AGROW>
                im_chan   = [im_chan;  t_chan];   %#ok<AGROW>
                im_line   = [im_line;  t_line];   %#ok<AGROW>
                im_col    = [im_col;   t_col];    %#ok<AGROW>
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 1e6], head);
                
            end
            
            F  = y(bitand(marker,Frame)>0);
            
            t_sync   = [];
            t_tcspc  = [];
            t_chan   = [];
            t_line   = [];
            t_col    = [];
            
            if ~in_frame
                if isempty(F)
                    y       = [];
                    tmpx    = [];
                    chan    = [];
                    marker  = [];
                    line    = 1;
                else
                    ind = (y<=F(1));
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    marker(ind)  = [];
                    line         = 1;
                    n_frames     = n_frames + 1;
                    f_times      = [f_times; F(1)];
                end
            end
            
            f_y = y;
            f_x = tmpx;
            f_ch = chan;
            f_m  = marker;
            
            clear y tmpx chan;
            
            L1 = f_y(bitand(f_m,LineStart)>0);
            L2 = f_y(bitand(f_m,LineStop)>0);
            
            ll = line + numel(L2)-1; % this will be the last complete line in the data stack
            if ll > ny
                L1 = L1(1:ny-line+1);
                L2 = L2(1:ny-line+1);
            end
            
            if numel(L1)>1
                for j=1:numel(L2)
                    ind = (f_y>L1(j))&(f_y<L2(j));
                    
                    t_sync   = [t_sync; y(ind)];
                    t_tcspc  = [t_tcspc; uint16(f_x(ind))];              %#ok<AGROW>
                    t_chan   = [t_chan; uint8(f_ch(ind))];                %#ok<AGROW>
                    t_line   = [t_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                    t_col    = [t_col;  uint16(1 + floor(nx.*(f_y(ind)-L1(j))./(L2(j)-L1(j))))];  %#ok<AGROW>
                    dt(line) = dt(line) + (L2(j)-L1(j));
                    line = line +1;
                end
            end
            
            im_sync   = [im_sync; t_sync];
            im_tcspc  = [im_tcspc; t_tcspc];
            im_chan   = [im_chan;  t_chan];
            im_line   = [im_line;  t_line];
            im_col    = [im_col;   t_col];
            
            head.ImgHdr_FrameTime = 1e9.*mean(diff(f_times))/head.TTResult_SyncRate;
            head.ImgHdr_PixelTime = 1e9.*mean(dt)/nx/head.TTResult_SyncRate;
            head.ImgHdr_DwellTime = head.ImgHdr_PixelTime./n_frames;
            
            close(h);
            
            
        elseif head.ImgHdr_Ident == 9 % MultiFrame Scan
            
          if isfield(head,'ImgHdr_MaxFrames')
            nz = head.ImgHdr_MaxFrames;
          else
             tim_p_frame =  1./head.ImgHdr_LineFrequency*ny;
             tot_time = head.TTResult_StopAfter*1e-3;
             nz = ceil(tot_time/tim_p_frame);
              
          end
            
            [~, ~, tmpchan, tmpmarkers] = PTU_Read(name, [1 1e4], head);
            dind = unique(tmpchan(~tmpmarkers));
            tag = zeros(nx,ny,numel(dind),nz);
            tau = tag;
            
            anzch      = 32;
            Resolution = max([1e9*head.MeasDesc_Resolution,0.128]);
            chDiv      = 1e-9*Resolution/head.MeasDesc_Resolution;
            Ngate      = ceil(1e9*head.MeasDesc_GlobalResolution./Resolution)+1;
            head.MeasDesc_Resolution = Resolution*1e-9;
            
            LineStart = 4;
            LineStop  = 2;
            Frame     = 3;
            
            if isfield(head,'ImgHdr_LineStart')
                LineStart = 2^(head.ImgHdr_LineStart-1);
            end
            if isfield(head,'ImgHdr_LineStop')
                LineStop = 2^(head.ImgHdr_LineStop-1);
            end
            if isfield(head,'ImgHdr_Frame')
                Frame = 2^(head.ImgHdr_Frame-1);
            end
            y        = [];
            tmpx     = [];
            chan     = [];
            markers  = [];
            dt       = zeros(ny,1);
            
            nphot    = head.TTResult_NumberOfRecords;
            im_sync  = zeros(nphot,1);
            im_tcspc = uint16(im_sync);
            im_chan  = uint8(im_sync);
            im_line  = uint16(im_sync);
            im_col   = uint16(im_sync);
            im_frame = uint16(im_sync);
            
            cn_phot = 0;
            Turns1   = [];
            Turns2   = [];
            
            cnt      = 0;
            tend     = 0;
            line     = 1;
            frame    = 1;
            
            h = waitbar(0,['Frame ' num2str(frame) ' out of ' num2str(nz)]);
            %             h2 = waitbar(0,['Frame' num2str(frame)], 'Units','normalized', 'Position', [0.5 0.4 0.25 0.2]);
            
            fprintf('\n\n');
            
            if head.ImgHdr_BiDirect == 0
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photons], head);
                
                while (num>0)
                    
                    cnt = cnt + num;
                    if ~isempty(y)
                        tmpy = tmpy+tend;
                    end
                    
                    ind = (tmpmarkers>0)|((tmpchan<anzch)&(tmptcspc<Ngate*chDiv));
                    
                    y       = [y; tmpy(ind)];                         %#ok<AGROW>
                    tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1]; %#ok<AGROW>
                    chan    = uint8([chan; tmpchan(ind)+1]);                 %#ok<AGROW>
                    markers = uint8([markers; tmpmarkers(ind)]);             %#ok<AGROW>
                    
                    if LineStart==LineStop
                        tmpturns = y(markers==LineStart);
                        if numel(Turns1)>numel(Turns2)           % first turn is a LineStop
                            Turns1 = [Turns1; tmpturns(2:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(1:2:end)];      %#ok<AGROW>
                        else
                            Turns1 = [Turns1; tmpturns(1:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(2:2:end)];      %#ok<AGROW>
                        end
                    else
                        Turns1 = [Turns1; y(markers==LineStart)]; %#ok<AGROW>
                        Turns2 = [Turns2; y(markers==LineStop)];  %#ok<AGROW>
                    end
                    Framechange = y(markers==uint8(Frame));
                    
                    ind          = (markers~=0);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    markers(ind) = [];
                    
                    tend  = y(end)+loc;
                    
                    if numel(Framechange)>=1
                        for k = 1:numel(Framechange)
                            
                            line=1;
                            ind = y<Framechange(k);
                            
                            yf = y(ind);
                            tmpxf = tmpx(ind);
                            chanf = chan(ind);
                            %                             markersf = markers(ind);
                            
                            y(ind)       = [];
                            tmpx(ind)    = [];
                            chan(ind)    = [];
                            markers(ind) = [];
                            
                            
                            Turns2f = Turns2(Turns2<Framechange(k));
                            Turns1f = Turns1(Turns1<Framechange(k));
                            Turns2(Turns2<Framechange(k)) = [];
                            Turns1(Turns1<Framechange(k)) = [];
                            
                            %                             im_frame = [im_frame; uint16(frame*ones(sum(y<Framechange(1)),1))];  %#ok<AGROW>
                            
                            
                            if numel(Turns2f)>1
                                for j=1:numel(Turns2f)
                                    
                                    t1 = Turns1f(1);
                                    t2 = Turns2f(1);
                                    
                                    ind           = yf<=t1;
                                    yf(ind)       = [];
                                    tmpxf(ind)    = [];
                                    chanf(ind)    = [];
                                    %                                     markersf(ind) = [];
                                    
                                    ind = (yf>t1)&(yf<t2);
                                    
                                    im_frame(cn_phot+1:cn_phot+sum(ind))  = uint16(frame*ones(sum(ind),1));
                                    im_sync(cn_phot+1:cn_phot+sum(ind))   = yf(ind);
                                    im_tcspc(cn_phot+1:cn_phot+sum(ind))  = uint16(tmpxf(ind));
                                    im_chan(cn_phot+1:cn_phot+sum(ind))   = uint8(chanf(ind));
                                    im_line(cn_phot+1:cn_phot+sum(ind))   = uint16(line.*ones(sum(ind),1));
                                    im_col(cn_phot+1:cn_phot+sum(ind))    = uint16(1 + floor(nx.*(yf(ind)-t1)./(t2-t1)));
                                    
                                    cn_phot = cn_phot+sum(ind);
                                    
                                    dt(line)  = t2-t1;
                                    line = line +1;
                                    
                                    waitbar(frame/nz,h,sprintf(['Frame ' num2str(frame) ' out of ' num2str(nz)]));
                                    %                                     waitbar(line/ny,h2,sprintf(['Frame ' num2str(frame)]));
                                    drawnow
                                    
                                    Turns1f(1) = [];
                                    Turns2f(1) = [];
                                end
                            end
                            [tmptag, tmptau] = Process_Frame(im_sync(im_frame == frame),im_col(im_frame == frame),im_line(im_frame == frame),im_chan(im_frame == frame),im_tcspc(im_frame == frame), head);
                            if ~isempty(tmptag)
                                tag(:,:,:,frame) = tmptag;
                                tau(:,:,:,frame) = tmptau;
                            end
                            frame = frame+1;
                            
                            %                             Framechange(1)=[];
                        end
                    end
                    
                    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photons], head);
                    
                end
                head.ImgHdr_PixelTime = 1e9.*mean(dt)/nx/head.TTResult_SyncRate;
                head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
                
                
                close(h)
                im_frame(cn_phot+1:end) = [];
                im_sync(cn_phot+1:end) = [];
                im_tcspc(cn_phot+1:end) = [];
                im_chan(cn_phot+1:end) = [];
                im_line(cn_phot+1:end) = [];
                im_col(cn_phot+1:end) = [];
                %             close(h2)
                
            else  % bidirectional scan
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 5e6], head);
                
                while (num>0)
                    
                    cnt = cnt + num;
                    if ~isempty(y)
                        tmpy = tmpy+tend;
                    end
                    
                    ind = ((tmpchan<anzch)&(tmptcspc<=Ngate*chDiv));
                    
                    y       = [y; tmpy(ind)];                         %#ok<AGROW>
                    tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1;]; %#ok<AGROW>
                    chan    = [chan; tmpchan(ind)];                   %#ok<AGROW>
                    markers = [markers; tmpmarkers(ind)];             %#ok<AGROW>
                    
                    if LineStart==LineStop
                        tmpturns = y(markers==LineStart);
                        if numel(Turns1)>numel(Turns2)           % first turn is a LineStop
                            Turns1 = [Turns1; tmpturns(2:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(1:2:end)];      %#ok<AGROW>
                        else
                            Turns1 = [Turns1; tmpturns(1:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(2:2:end)];      %#ok<AGROW>
                        end
                    else
                        Turns1 = [Turns1; y(markers==LineStart)]; %#ok<AGROW>
                        Turns2 = [Turns2; y(markers==LineStop)];  %#ok<AGROW>
                    end
                    Framechange = y(markers==Frame);
                    
                    ind          = (markers~=0);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    markers(ind) = [];
                    
                    tend = y(end)+loc;
                    
                    
                    
                    
                    if numel(Framechange)>=1
                        for k = 1:numel(Framechange)
                            
                            line=1;
                            ind = y<Framechange(k);
                            
                            yf = y(ind);
                            tmpxf = tmpx(ind);
                            chanf = chan(ind);
                            markersf = markers(ind);
                            
                            y(ind)       = [];
                            tmpx(ind)    = [];
                            chan(ind)    = [];
                            markers(ind) = [];
                            
                            
                            Turns2f = Turns2(Turns2<=Framechange(k));
                            Turns1f = Turns1(Turns1<=Framechange(k));
                            Turns2(Turns2<=Framechange(k)) = [];
                            Turns1(Turns1<=Framechange(k)) = [];
                            
                            %                             im_frame = [im_frame; uint16(frame*ones(sum(y<Framechange(1)),1))];  %#ok<AGROW>
                            
                            
                            if numel(Turns2f)>2
                                for  j=1:2:2*floor(numel(Turns2f)/2-1)
                                    
                                    t1 = Turns1f(1);
                                    t2 = Turns2f(1);
                                    
                                    ind          = (yf<t1);
                                    yf(ind)       = [];
                                    tmpxf(ind)    = [];
                                    chanf(ind)    = [];
                                    markersf(ind) = [];
                                    
                                    ind = (yf>=t1)&(yf<=t2);
                                    
                                    im_frame(cn_phot+1:cn_phot+sum(ind))  = uint16(frame*ones(sum(ind),1));
                                    im_sync(cn_phot+1:cn_phot+sum(ind))   = yf(ind);
                                    im_tcspc(cn_phot+1:cn_phot+sum(ind))  = uint16(tmpxf(ind));
                                    im_chan(cn_phot+1:cn_phot+sum(ind))   = uint8(chanf(ind));
                                    im_line(cn_phot+1:cn_phot+sum(ind))   = uint16(line.*ones(sum(ind),1));
                                    im_col(cn_phot+1:cn_phot+sum(ind))    = uint16(1 + floor(nx.*(yf(ind)-t1)./(t2-t1)));
                                    
                                    cn_phot = cn_phot+sum(ind);
                                    dt(line)  = t2-t1;
                                    line = line +1;
                                    
                                    t1 = Turns1f(2);
                                    t2 = Turns2f(2);
                                    
                                    ind = (yf<t1);
                                    yf(ind)       = [];
                                    tmpxf(ind)    = [];
                                    chanf(ind)    = [];
                                    markersf(ind) = [];
                                    
                                    ind = (yf>=t1)&(yf<=t2);
                                    
                                    im_frame(cn_phot+1:cn_phot+sum(ind))  = uint16(frame*ones(sum(ind),1));
                                    im_sync(cn_phot+1:cn_phot+sum(ind))   = yf(ind);
                                    im_tcspc(cn_phot+1:cn_phot+sum(ind))  = uint16(tmpxf(ind));
                                    im_chan(cn_phot+1:cn_phot+sum(ind))   = uint8(chanf(ind));
                                    im_line(cn_phot+1:cn_phot+sum(ind))   = uint16(line.*ones(sum(ind),1));
                                    im_col(cn_phot+1:cn_phot+sum(ind))    = uint16(nx - floor(nx.*(yf(ind)-t1)./(t2-t1)));
                                    
                                    cn_phot = cn_phot+sum(ind);
                                    dt(line)  = t2-t1;
                                    
                                    line = line +1;
                                    waitbar(line/ny,h);
                                    drawnow
                                    
                                    Turns1f(1:2) = [];
                                    Turns2f(1:2) = [];
                                    
                                    
                                    waitbar(frame/head.ImgHdr_Frame,h,sprintf(['Frame ' num2str(frame) ' out of ' num2str(head.ImgHdr_Frame)]));
                                    %                                     waitbar(line/ny,h2,sprintf(['Frame ' num2str(frame)]));
                                    drawnow
                                    
%                                     Turns1f(1) = [];
%                                     Turns2f(1) = [];
                                end
                            end
                            [tmptag, tmptau] = Process_Frame(im_sync(im_frame == frame),im_col(im_frame == frame),im_line(im_frame == frame),im_chan(im_frame == frame),im_tcspc(im_frame == frame), head);
                            tag(:,:,:,frame) = tmptag;
                            tau(:,:,:,frame) = tmptau;
                            
                            frame = frame+1;
                            %                             Framechange(1)=[];
                        end
                    end
                    
                    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 5e6], head);
                    
                end
                
                if ~isempty(Turns2)
                    t1 = Turns1(end-1);
                    t2 = Turns2(end-1);
                    
                    ind = (y<t1);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    
                    ind = (y>=t1)&(y<=t2);
                    
                    im_frame(cn_phot+1:cn_phot+sum(ind))  = uint16(frame*ones(sum(ind),1));
                    im_sync(cn_phot+1:cn_phot+sum(ind))   = y(ind);
                    im_tcspc(cn_phot+1:cn_phot+sum(ind))  = uint16(tmpx(ind));
                    im_chan(cn_phot+1:cn_phot+sum(ind))   = uint8(chan(ind));
                    im_line(cn_phot+1:cn_phot+sum(ind))   = uint16(line.*ones(sum(ind),1));
                    im_col(cn_phot+1:cn_phot+sum(ind))    = uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)));
                    dt(line)  = t2-t1;
                    
                    line = line +1;
                    cn_phot = cn_phot+sum(ind);
                    
                    t1 = Turns1(end);
                    t2 = Turns2(end);
                    
                    ind = (y<t1);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    
                    ind = (y>=t1)&(y<=t2);
                    
                    im_frame(cn_phot+1:cn_phot+sum(ind))  = uint16(frame*ones(sum(ind),1));
                    im_sync(cn_phot+1:cn_phot+sum(ind))   = yf(ind);
                    im_tcspc(cn_phot+1:cn_phot+sum(ind))  = uint16(tmpxf(ind));
                    im_chan(cn_phot+1:cn_phot+sum(ind))   = uint8(chanf(ind));
                    im_line(cn_phot+1:cn_phot+sum(ind))   = uint16(line.*ones(sum(ind),1));
                    im_col(cn_phot+1:cn_phot+sum(ind))    = uint16(nx - floor(nx.*(y(ind)-t1)./(t2-t1)));
                    
                    cn_phot = cn_phot+sum(ind);
                    dt(line)  = t2-t1;
                    
                    line = line +1;
                    waitbar(line/ny);
                    drawnow
                end
                [tmptag, tmptau] = Process_Frame(im_sync(im_frame == frame),im_col(im_frame == frame),im_line(im_frame == frame),im_chan(im_frame == frame),im_tcspc(im_frame == frame), head);
               if frame<=nz&&~isempty(tmptag)
                tag(:,:,:,frame) = tmptag;
                tau(:,:,:,frame) = tmptau;
               end
                head.ImgHdr_PixelTime = 1e9.*mean(dt)/nx/head.TTResult_SyncRate;
                head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
                
                
                close(h);
                im_frame(cn_phot+1:end) = [];
                im_sync(cn_phot+1:end) = [];
                im_tcspc(cn_phot+1:end) = [];
                im_chan(cn_phot+1:end) = [];
                im_line(cn_phot+1:end) = [];
                im_col(cn_phot+1:end) = [];
                
            end
            
            
            SyncRate   = 1./head.MeasDesc_GlobalResolution;
            maxch_n = numel(dind);
            
            tcspc_pix = zeros(nx,ny,Ngate,maxch_n);
            %             time = zeros(numel(im_sync),maxch_n);
            time = {};
            tags = zeros(nx,ny,maxch_n);
            taus = tags;
            bin = permute(repmat((1:Ngate)',[1 nx,ny]),[2,3,1])*Resolution; % 3D time axis
            for ch = 1:maxch_n
                ind = im_chan==dind(ch)+1;
                tcspc_pix(:,:,:,ch) =  mHist3(double(im_line(ind)),double(im_col(ind)),double(im_tcspc(ind)),1:nx,1:ny,1:Ngate); % tcspc histograms for all the pixels at once!
                time{ch} = round(im_sync(ind)*1/SyncRate/Resolution/1e-9)+double(im_tcspc(ind)); % in tcspc bins
                tags(:,:,ch) = sum(tcspc_pix(:,:,:,ch),3);
                taus(:,:,ch) = real(sqrt((sum(bin.^2.*tcspc_pix(:,:,:,ch),3)./tags(:,:,ch))-(sum(bin.*tcspc_pix(:,:,:,ch),3)./tags(:,:,ch)).^2));
            end
            
            save([name(1:end-4) '_FLIM_data'],'tag','tau','tags','taus','time','tcspc_pix','head','im_sync','im_tcspc','im_line','im_col','im_chan','im_frame','-v7.3')
        end
        
    end
end


if ~isempty(head)&&(plt==1)
    x = head.ImgHdr_X0+(1:nx)*head.ImgHdr_PixResol;
    y = head.ImgHdr_Y0+(1:ny)*head.ImgHdr_PixResol;
    imagesc(x,y,sum(tags,3));
    set(gca,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse');
    xlabel('x / ???m');
    ylabel('y / ???m');
    title('Intensity')
    figure
    imagesc(x,y,sum(taus.*tags,3)./sum(tags,3));
    set(gca,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse');
    xlabel('x / ???m');
    ylabel('y / ???m');
    title('FLIM')
    
end
end

function [tag, tau, tcspc_pix, time] = Process_Frame(im_sync,im_col,im_line,im_chan,im_tcspc,head)

Resolution = max([head.MeasDesc_Resolution*1e9 0.256]); % resolution of 1 ns to calculate average lifetimes
chDiv      = ceil(1e-9*Resolution./head.MeasDesc_Resolution);
SyncRate   = 1./head.MeasDesc_GlobalResolution;
nx = head.ImgHdr_PixX;
ny = head.ImgHdr_PixY;
dind    = double(unique(im_chan));
Ngate   = round(head.MeasDesc_GlobalResolution./head.MeasDesc_Resolution*(head.MeasDesc_Resolution/Resolution)*1e9);
maxch_n = numel(dind);

tcspc_pix = zeros(nx,ny,Ngate,maxch_n);
% time = zeros(numel(im_sync),maxch_n);
time = {};
tag = zeros(nx,ny,maxch_n);
tau = tag;
bin = permute(repmat((1:Ngate)',[1 nx,ny]),[2,3,1])*Resolution; % 3D time axis
for ch = 1:maxch_n
    ind = im_chan==dind(ch);
    tcspc_pix(:,:,:,ch) =  mHist3(double(im_line(ind)),double(im_col(ind)),double(im_tcspc(ind)./chDiv),1:nx,1:ny,1:Ngate); % tcspc histograms for all the pixels at once!
    time{ch} = round(im_sync(ind)*1/SyncRate/Resolution/1e-9)+double(im_tcspc(ind)); % in tcspc bins
    tag(:,:,ch) = sum(tcspc_pix(:,:,:,ch),3);
    tau(:,:,ch) = real(sqrt((sum(bin.^2.*tcspc_pix(:,:,:,ch),3)./tag(:,:,ch))-(sum(bin.*tcspc_pix(:,:,:,ch),3)./tag(:,:,ch)).^2));
end
end
