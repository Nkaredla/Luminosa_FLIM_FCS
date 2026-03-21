function [bin,tcspcdata,head] = Harp_tcspc(name,resolution,deadtime,photons)
%[bin,tcspcdata,head] = Harp_tcspc(name,resolution,deadtime)
% This program gives the TCSPC Histogram for .ht3 files or .ptu files.
% name = filename
% resolution = binwidth of the TCSPC histogram
% deadtime [s] = filters photons within dead time of the detector. Typical value
% is 100e-9. By default, all photons will be read.
% This program needs HT3_Read.m, PTU_Read_Head.m, PTU_Read.m and mHist.m
% (c) Narain Karedla (2013)
% edited 31.07.2015 by Sebastian Isbaner

if nargin<4||isempty(photons)
    photons=1e6; % 1 million records at a time for processing.
end
if strcmp(name(end-2:end),'ht3')
    head=HT3_Read(name);
    Timeunit=1/(head.SyncRate);
    Resolution=head.Resolution*1e-9;
    NCounts=head.Records;
    
elseif strcmp(name(end-2:end),'ptu')
    head=PTU_Read_Head(name);
    Timeunit=1/(head.TTResult_SyncRate);
    Resolution=head.MeasDesc_Resolution;
    NCounts=head.TTResult_NumberOfRecords;
else
    disp('Not a valid file!')
    return
end

if nargin<2 || isempty(resolution)
    chdiv =1;
    resolution = Resolution;
else
    chdiv = ceil(resolution/Resolution);
    Resolution = max(Resolution,resolution);
end
if nargin<3 || isempty(deadtime)
    deadtime=0; % read all photons
end
Deadtime=deadtime;
NChannels=ceil(Timeunit/Resolution);
tcspcdata = [];
nRemovedPhotons=0;
tcspcfile = [name(1:end-4) '.ht3tcspc'];
if ( exist(tcspcfile, 'file') == 2) % check if tcspc data was saved
    load(tcspcfile, '-mat');
    if abs(resolution - Resolution)>1e-12 || deadtime ~= Deadtime
        delete(tcspcfile);
        [bin,tcspcdata,head] = Harp_tcspc(name,resolution,deadtime);
    end
else
    %     h = waitbar(0,'Please wait...');
    %
    num = 1;
    cnt = 0;
    bin = 0:NChannels-1;
    dind = [];
    while num>0
        if strcmp(name(end-2:end),'ht3')
            [sync, tcspc, chan, special, num ] = HT3_Read(name, [cnt+1, photons, head.length]);
        else
            if strcmp(name(end-2:end),'ptu')
                [sync, tcspc, chan, special, num] = PTU_Read(name, [cnt+1, photons],head);
            end
        end
        if isempty(dind)
            [dind,m]= unique(sort(chan(~special)),'legacy');
            occurence = diff([0;vertcat(m)]);
            dind(occurence<10)=[];
            dnum = length(dind);
            tcspcdata = zeros(NChannels,dnum);
        end
        tcspc(logical(special)) = [];
        chan(logical(special))  = [];
        sync(logical(special)) = [];
        tcspc = round(tcspc/chdiv);
        
        cnt = cnt+num;
        %         waitbar(cnt/NCounts,h)
        if (num>0) && (~isempty(chan))
            if deadtime>0
                %                 sync = sync*Timeunit + tcspc*Resolution; % Note: this is a pseudo time, overcounts from previous reads are not taken into account.
                %                 flvp = false(numel(chan),numel(dind));
                %                 for ii = 1:numel(dind)
                %                     flvp(:,ii) = chan==dind(ii);
                %                 end
                %                 idx = any(flvp & circshift(flvp,[1 0]),2);
                %                 dsync=diff([0;sync]);
                %                 idx = idx & (dsync<deadtime);
                %                 tcspc(idx)=[];
                %                 chan(idx)=[];
                %                 nRemovedPhotons = nRemovedPhotons + sum(idx);
                
                tttr = sync*Timeunit+tcspc*Resolution;
                difftttr = diff([0; tttr]);
                deadtime = 80e-9; % TCSPC deadtime + syncperiod
                ind = difftttr<deadtime; % photons arriving within the deatime + 1 sync period
                num = find(ind);
                t1 = sync(find(ind)-1); % first photon sync times
                t2 = sync(find(ind));   % second photon sync times
                tgate = ceil((tcspc(find(ind)-1)*resolution+80e-9)/Timeunit); % The nearest sync after the dead time
                xind = (t2-t1)<=(tgate);
                tmptcspc = tcspc(num(find(xind)));
                tmpchan  = chan(num(find(xind)));
                %                 tcspcdata(:,jj) = tcspcdata(:,jj)+ histc(tmptcspc,bin);
                
                for jj = 1:dnum
                    tcspcdata(:,jj) = tcspcdata(:,jj)+ mHist(tmptcspc(tmpchan == dind(jj)),bin); %#ok<AGROW>
                end
            else % No deadtime into consideration
                tmptcspc = tcspc;
                tmpchan = chan;
                for jj = 1:dnum
                    tcspcdata(:,jj) = tcspcdata(:,jj)+ histc(tmptcspc(tmpchan == dind(jj)),bin); %#ok<AGROW>
                end
            end
            semilogy(bin,tcspcdata)
            drawnow
        end
    end
    %     end
    %     fprintf('Harp_tcspc removed %d of %d photons within dead time of %.2f ns\n',...
    %         nRemovedPhotons,sum(tcspcdata(:))+nRemovedPhotons,deadtime*1e9);
        save(tcspcfile, 'bin', 'tcspcdata','head','Resolution','Deadtime','nRemovedPhotons');
%         close(h)
    semilogy(bin,tcspcdata)
    
end
