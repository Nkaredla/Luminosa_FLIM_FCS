function [data, tag, life_imm, life_imm_c]=Process_scan(name, num_PIE, timegate, flagint, cutoff, flagfit, pic)
% Input:
% name - filename of the .ht3-file (has to end on .ht3)
% num_PIE - number of pulses
% timegate - how to shift the bins of the tcspc curve so it start with peak
%            (see lines 69ff for a more detailed explanation)
% flagint - intensity image shows all photons (flagint=0 or []) or only
%           photons that were actually used for the calculations (flagint>1)
% cutoff - time in ns after the peak where the curve is purely exponential
% flagfit - if set to 1, lifetimes will be fitted, if set to 0, average
%           arrival time after cutoff is used as lifetime
% pic - set to true if you want plots of the intensity and lifetime images
%
% Output:
% data - structure containing timegates, tcspc-curve without the "bump" at
%        the end (tcspc), intensity image (tag) and lifetime image (life_imm),
%        binnumber of the peak (pos), cutoff in ns (cutoff) and cutoff in bins (shift)
% tag  - intensity image
% life_imm - lifetime image

if nargin<7 || isempty(pic)
    pic=0;
end

if nargin<6 || isempty(flagfit)
    flagfit = 1;         % By default, this program will fit the lifetime histograms for each pixel.
end

if nargin<5 || isempty(cutoff)
    cutoff=0.5; % ns
end

if nargin<4 || isempty(flagint) || flagint==0
    flagint=0; % if flagint ==1, The intensity image will represent photons after the cutoff
else
    flagint=1;
end

if nargin<3 || isempty(timegate)
    timegate=[];
end

if nargin<2 || isempty(num_PIE)
    num_PIE=1;
end

if exist([name(1:end-4),'_Core_Scan.mat'],'file')
    load([name(1:end-4),'_Core_Scan.mat'],'head','im_sync','im_tcspc','im_chan','im_line','im_col');
elseif strcmp(name(end-2:end),'ht3')
    disp('Read tcspc-data from ht3-file...');
    [head, im_sync, im_tcspc, im_chan, im_line, im_col] = Core_ScanRead(name);
    save([name(1:end-4),'_Core_Scan.mat'],'head','im_sync','im_tcspc','im_chan','im_line','im_col');
elseif strcmp(name(end-2:end),'ptu')
    disp('Read tcspc-data from ptu-file...');
    [head, im_sync, im_tcspc, im_chan, im_line, im_col] = PTU_ScanRead(name);
    head.ImgHdr.PixX = head.ImgHdr_PixX;
    head.ImgHdr.PixY = head.ImgHdr_PixY;
    head.ImgHdr.X0   = head.ImgHdr_X0;
    head.ImgHdr.Y0   = head.ImgHdr_Y0;
    head.ImgHdr.PixelSize = head.ImgHdr_PixResol;
    head.ImgHdr.PixelTime = head.ImgHdr_DwellTime; % pixel dwell time in ns
    head.Resolution  = 1e9*head.MeasDesc_Resolution;
    head.SyncRate    = 1./head.MeasDesc_GlobalResolution;
    save([name(1:end-4),'_Core_Scan.mat'],'head','im_sync','im_tcspc','im_chan','im_line','im_col');
else
    disp('You have to give an ht3- or ptu-file.');
    return
end

dind       = unique(im_chan);
maxch      = numel(dind);
maxres     = max([head.Resolution]);
Resolution = max([maxres 0.064]);
chDiv      = Resolution/maxres;
im_tcspc   = ceil(im_tcspc./chDiv);
Ngate      = double(max(im_tcspc));
pixel      = head.ImgHdr.PixelSize;
tcspcdata  = zeros(maxch, Ngate);
nx         = head.ImgHdr.PixX;
ny         = head.ImgHdr.PixY;
t          = (1:Ngate).*Resolution;
x0         = head.ImgHdr.X0;
y0         = head.ImgHdr.Y0;
% pxDwellTime= head.ImgHdr.PixelTime*1e-9; % pixel dwell time in seconds
% RepRate    = head.SyncRate;             % repetition rate of excitation pulses in Hz

deadtime   = 80; % ns
apd        = 33; % ns

for ch = 1:maxch
    tcspcdata(ch,:) = hist(double(im_tcspc(im_chan == dind(ch))),1:Ngate);
    indel(ch) = sum(tcspcdata(ch,:))<nx*ny; %#ok<*AGROW> deleting image planes containing less than nx*ny photons
end
tcspcdata(indel,:) = [];
dind(indel)    = [];
maxch_n        = numel(dind);
tcspcdata=tcspcdata.';

data.t=t;
data.tcspcdata=tcspcdata;
if isempty(timegate)
    disp('Determine time-gates...');
    [timegate, Ngate] = DetectTimeGates(tcspcdata, num_PIE, Resolution);
else
    Ngate = 1 + timegate(1,2)+timegate(1,4)-timegate(1,1);
end
%                 timegate(  1   ,:): time-window for pulse  1   in detection channel 1
%                 timegate(  2   ,:): time-window for pulse  1   in detection channel 2
%            ...  timegate(  n   ,:): time-window for pulse  1   in detection channel n
%                 timegate( n+1  ,:): time-window for pulse  2   in detection channel 1
%                 timegate( n+2  ,:): time-window for pulse  2   in detection channel 2
%            ...  timegate( 2*n  ,:): time-window for pulse  2   in detection channel n
%            ...  timegate(num_PIE*n,:): time-window for pulse num_PIE in detection channel n
%
%                 timegate(:, 1)    : begin of time-window
%                 timegate(:, 2)    : end of time-window
%                 timegate(:, 3)    : if > 0: continuation of time-window
%                 timegate(:, 4)    : if > 0: end of time-window

tcspc=zeros(Ngate,maxch_n,num_PIE); % Ngate timechannels, maxch_n detectors, num_PIE pulses in each.
ind=timegate(:,3)==0;
for j=1:num_PIE
    for i=1:maxch_n
        if ind((j-1)*maxch_n+i)
            tcspc(:,i,j)=tcspcdata(timegate((j-1)*maxch_n+i,1):timegate((j-1)*maxch_n+i,2),i);
        else
            tcspc(:,i,j)=[tcspcdata(timegate((j-1)*maxch_n+i,1):timegate((j-1)*maxch_n+i,2),i);...
                tcspcdata(timegate((j-1)*maxch_n+i,3):timegate((j-1)*maxch_n+i,4),i)];
        end
    end
end

% max_temp=zeros(size(tcspc,2),size(tcspc,3));
% for i=1:size(tcspc,2)
%     for j=1:size(tcspc,3)
%         [~,max_temp(i,j)]=max(tcspc(:,i,j));
%     end
% end
[~, ~, bgpix,~]= BG_val(tcspc,ceil(cutoff/Resolution),20,nx,ny);
% Ngate = size(tcspc,1);
% Ngate = min(minchan(:)); % cut away the "bump" at the end of the tcspc-curves
% tcspc = tcspc(1:Ngate,:,:);
% ind = false((numel(im_tcspc)),1);
% for j = 1:num_PIE
%     for i = 1:maxch_n
%         indn = mod(timegate((j-1)*maxch_n+i,1)+Ngate,NChannels)<timegate((j-1)*maxch_n+i,1); % true if the peak "hangs over" to beginning of tcspc curve
%         if indn
%             ind = ind|(im_tcspc>=timegate((j-1)*maxch_n+i,1) | (im_tcspc<= mod(timegate((j-1)*maxch_n+i,1)+Ngate,NChannels)) & im_chan==dind(i));
%         else
%             ind = ind|(im_tcspc>=timegate((j-1)*maxch_n+i,1) & (im_tcspc<=timegate((j-1)*maxch_n+i,1)+Ngate) & im_chan==dind(i));
%         end
%     end
% end
% delchan =~ind;
% im_tcspc(delchan) = [];
% im_col(delchan)   = [];
% im_line(delchan)  = [];
% im_chan(delchan)  = [];
% im_sync(delchan)  = [];

microtime=zeros(numel(im_tcspc),maxch_n,num_PIE);
im_pulse=zeros(numel(im_tcspc),1);
ind=timegate(:,3)==0;
for j=1:num_PIE
    for i=1:maxch_n % we want 1<=microtime<=Ngate
        if ind((j-1)*maxch_n+i)
            chan_ind=im_chan==dind(i);
            microtime(chan_ind,i,j)=im_tcspc(chan_ind)-timegate((j-1)*maxch_n+i,1)+1;
            im_pulse(chan_ind)=j;
        else
            chan_ind=im_chan==dind(i);
            timefirst=im_tcspc>timegate((j-1)*maxch_n+i,1)-1;
            timelast=im_tcspc<timegate((j-1)*maxch_n+i,end)+1;
            microtime(chan_ind & timefirst,i,j)=double(im_tcspc(chan_ind & timefirst))-timegate((j-1)*maxch_n+i,1)+1;
            microtime(chan_ind & timelast,i,j)=double(im_tcspc(chan_ind & timelast))+timegate((j-1)*maxch_n+i,2)-timegate((j-1)*maxch_n+i,1)+1;
            im_pulse((chan_ind & timefirst)|(chan_ind & timelast))=j;
        end
    end
end

data.timegate=timegate;
data.tcspc=tcspc;

c=0; pos=zeros(maxch_n*num_PIE,1);
for i=1:size(tcspc,2)
    for j=1:size(tcspc,3)
        c=c+1;
        [~,pos(c)]=max(tcspc(:,i,j));
    end
end

shift=ceil(cutoff/Resolution);
figure; plot((1:length(tcspc(:,1)))*Resolution,tcspc(:,1),'.b',([pos(1)+shift pos(1)+shift])*Resolution,[0 max(tcspc(:,1))],'-k');
title('cutoff in channel 1'); xlabel('time [ns]'); % plot where the cutoff is
drawnow

im_pixel=zeros(size(im_line));
for i=1:maxch_n
    ind = im_chan==dind(i);
    im_pixel(ind)=double(im_line(ind))+(double(im_col(ind))-1)*ny + nx*ny*(i-1)+ nx*ny*maxch_n*(double(im_pulse(ind))-1);
end

clearvars im_chan im_line im_col im_tcspc

tag=zeros(ny*nx*maxch_n*num_PIE,1); life_imm=NaN(numel(tag),1); life_imm_c=NaN(numel(tag),1); % numphot = tag;
if num_PIE>1
    microtime=sum(sum(microtime));
else
    microtime=sum(microtime,2);
end
if maxch~=0
    if flagint
        mic_ind=zeros(numel(im_sync),1);
        for j=1:num_PIE % find photons before the cutoff
            for ch=1:maxch_n
                mic_ind=mic_ind+double(microtime(:,ch,j)>pos((j-1)*maxch_n+ch)+shift);
            end
        end
        mic_ind=not(logical(mic_ind)); % index for photons before the cutoff
        tag=hist(im_pixel(~mic_ind),1:numel(tag)); % number of photons after cutoff per pixel
    else
        tag=hist(im_pixel,1:numel(tag));
    end
    indpix=unique((tag>100).*(1:numel(tag)).');
    if indpix(1)==0; indpix(1)=[]; end
    t0_pixel=zeros(numel(indpix),1);
    for j=1:num_PIE
        for ch=1:maxch_n
            ind=(nx*ny*(ch-1)+nx*ny*ch*(j-1))-1<indpix & indpix<(nx*ny*ch+nx*ny*ch*(j-1))+1;
            t0_pixel(ind)=pos((j-1)*maxch_n+ch)+shift+1;
        end
    end
    disp('Determine lifetimes:'); reverseStr = '';

    for i=1:numel(indpix)
        if flagfit
            [tmplife tmplife_c]= PixelTau(microtime(im_pixel==indpix(i)),im_sync(im_pixel==indpix(i)),t0_pixel(i),deadtime,apd,Resolution);
            
            
%             y = histc(microtime(im_pixel == indpix(i)),1:Ngate); % tcspc-histogram
            %             pidt = gradient(-log(1-cumsum(y)/sum(y)*(1-exp(-tag(indpix(i))/pxDwellTime/RepRate))));
            %             pidt = gradient(-log(1-cumsum(y)/sum(y)*(tag(indpix(i))/pxDwellTime/RepRate))); % 1-exp(-epsilon)=countrate/rep.rate
            %             [cx, k, ~, ~, ~, ~] = DistTailfit(pidt(t0_pixel(i):end), Resolution,[],[],[],[],1e9/RepRate);
            
            % Joerg's Pile-up correction
            
            %             eps = -log(1-(sum(y)./pxDwellTime./head.SyncRate));
            % %             eps = sum(y)./pxDwellTime./head.SyncRate;
            %             reco = y./sum(y);
            %             for jj=1:4
            %                 tmp = conv([reco; reco; reco],ones(deadtime,1));
            %                 tmp = tmp(size(reco,1):2*size(reco,1)-1);
            %                 reco = y.*exp(eps*tmp);
            %                 reco = reco/sum(reco);
            %             end
            
%             [cx, k, ~, ~, ~, ~] = DistTailfit(y(t0_pixel(i):end), Resolution);
            %              [cx, k, ~, ~, ~, ~] = DistTailfit(reco(t0_pixel(i):end), Resolution);
            life_imm(indpix(i)) = tmplife; life_imm_c(indpix(i))= tmplife_c;
%             life_imm(indpix(i)) = sum(cx)./sum(cx.*k);
        else
            life_imm(indpix(i))=(mean(microtime(im_pixel==indpix(i)))-t0_pixel(i))*Resolution;
        end
        msg = sprintf('Processed %d/%d pixels...\n', i, numel(indpix));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    tag=reshape(tag,[ny,nx,maxch_n,num_PIE]);
    life_imm=reshape(life_imm,[ny,nx,maxch_n,num_PIE]);
    life_imm_c=reshape(life_imm_c,[ny,nx,maxch_n,num_PIE]);
    
    if ~flagfit
        c=0;
        for i = 1:maxch_n
            for j = 1:num_PIE
                c = c+1;
                life_imm(:,:,i,j) = (life_imm(:,:,i,j).*tag(:,:,i,j)-round((Ngate-pos(c)-shift)/2)*Resolution*bgpix(i,j))...
                    ./(tag(:,:,i,j)-bgpix(i,j));
            end
        end
    end
end


data.tag=tag;
data.life_imm=life_imm;
data.life_imm_c = life_imm_c;
data.pos=pos;
data.shift=shift;
data.cutoff=cutoff;


if pic
    if maxch_n>0
        for j=1:num_PIE
            close all
            figure;
            set (gcf,'name','Intensity Image','NumberTitle','off')
            for ch=1:maxch_n
                subplot(1,maxch_n,ch)
                imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), data.tag(:,:,ch,j))
                colorbar
                axis equal
                axis tight
                colormap('hot')
            end
        end
        for j=1:num_PIE
            close all
            figure;
            set (gcf,'name','Lifetime Image','NumberTitle','off')
            for ch=1:maxch_n
                subplot(1,maxch_n,ch)
                imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), data.life_imm(:,:,ch,j))
                colorbar
                axis equal
                axis tight
                colormap('hot')
            end
        end
    end
end
save([name(1:end-4),'_PS'],'data','life_imm','tag');
end

