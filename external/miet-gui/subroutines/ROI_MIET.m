function res = ROI_MIET(name, num_PIE, timegate, flagint, cutoff, flagfit, min_cluster, pic)
% name = 'U:\Narain\from Anna\MIET_Actin_Vinculin_12h\Image_040.ht3';

if nargin<8 || isempty(pic)
    pic=0;
end

if nargin<7 || isempty(min_cluster)
    min_cluster = 10;
end

if nargin<6 || isempty(flagfit)
    flagfit = 1;         % By default, this program will fit the lifetime histograms for each pixel.
end

if nargin<5 || isempty(cutoff)
    cutoff=0.3; % ns
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
    load([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
elseif strcmp(name(end-2:end),'ht3')
    disp('Read tcspc-data from ht3-file...');
    [head, im_tcspc, im_chan, im_line, im_col] = Core_ScanRead(name);
    save([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
elseif strcmp(name(end-2:end),'ptu')
    disp('Read tcspc-data from ptu-file...');
    [head, im_tcspc, im_chan, im_line, im_col] = PTU_ScanRead(name);
    head.ImgHdr.PixX = head.ImgHdr_PixX;
    head.ImgHdr.PixY = head.ImgHdr_PixY;
    head.ImgHdr.X0   = head.ImgHdr_X0;
    head.ImgHdr.Y0   = head.ImgHdr_Y0;
    head.ImgHdr.PixelSize = head.ImgHdr_PixResol;
    head.ImgHdr.PixelTime = head.ImgHdr_DwellTime; % pixel dwell time in ns
    head.Resolution  = 1e9*head.MeasDesc_Resolution;
    head.SyncRate    = 1./head.MeasDesc_GlobalResolution;
    save([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
else
    disp('You have to give an ht3- or ptu-file.');
    return
end

dind       = unique(im_chan);
maxch      = numel(dind);
maxres     = max([head.Resolution]);
Resolution = max([maxres 0.032]);
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
pxDwellTime= head.ImgHdr.PixelTime*1e-9; % pixel dwell time in seconds
RepRate    = head.SyncRate;             % repetition rate of excitation pulses in Hz

for ch = 1:maxch
    tcspcdata(ch,:) = mHist(double(im_tcspc(im_chan == dind(ch))),1:Ngate);
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

max_temp=zeros(size(tcspc,2),size(tcspc,3));
for i=1:size(tcspc,2)
    for j=1:size(tcspc,3)
        [~,max_temp(i,j)]=max(tcspc(:,i,j));
    end
end
[~, minchan, bgpix,~]= BG_val(tcspc,ceil(cutoff/Resolution),20,nx,ny);
NChannels = size(tcspc,1);
Ngate = min(minchan(:)); % cut away the "bump" at the end of the tcspc-curves
tcspc = tcspc(1:Ngate,:,:);
ind = false((numel(im_tcspc)),1);
for j = 1:num_PIE
    for i = 1:maxch_n
        indn = mod(timegate((j-1)*maxch_n+i,1)+Ngate,NChannels)<timegate((j-1)*maxch_n+i,1); % true if the peak "hangs over" to beginning of tcspc curve
        if indn
            ind = ind|(im_tcspc>=timegate((j-1)*maxch_n+i,1) | (im_tcspc<= mod(timegate((j-1)*maxch_n+i,1)+Ngate,NChannels)) & im_chan==dind(i));
        else
            ind = ind|(im_tcspc>=timegate((j-1)*maxch_n+i,1) & (im_tcspc<=timegate((j-1)*maxch_n+i,1)+Ngate) & im_chan==dind(i));
        end
    end
end
delchan =~ind;
im_tcspc(delchan) = [];
im_col(delchan)   = [];
im_line(delchan)  = [];
im_chan(delchan)  = [];

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
im_pulse = uint8(im_pulse);

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

im_pixel=zeros(size(im_line));
for i=1:maxch_n
    ind = im_chan==dind(i);
    im_pixel(ind)=double(im_line(ind))+(double(im_col(ind))-1)*ny + nx*ny*(i-1)+ nx*ny*maxch_n*(double(im_pulse(ind))-1);
end

tag=zeros(ny*nx*maxch_n*num_PIE,1); life_imm=NaN(numel(tag),1); % numphot = tag;
if num_PIE>1
    microtime=sum(sum(microtime));
else
    microtime=sum(microtime,2);
end

if maxch~=0
    if flagint
        mic_ind=zeros(numel(im_tcspc),1);
        for j=1:num_PIE % find photons before the cutoff
            for ch=1:maxch_n
                mic_ind=mic_ind+double(microtime(:,ch,j)>pos((j-1)*maxch_n+ch)+shift);
            end
        end
        mic_ind=not(logical(mic_ind)); % index for photons before the cutoff
        tag=mHist(im_pixel(~mic_ind),1:numel(tag)); % number of photons after cutoff per pixel
    else
        tag=mHist(im_pixel,1:numel(tag));
    end
end
 tag=reshape(tag,[ny,nx,maxch_n,num_PIE]);
 impix1 = reshape(1:ny*nx,[ny nx]);
shift=ceil(cutoff/Resolution);

tag = tag(2:end,2:end);
impix1 = impix1(2:end,2:end);


% [grain, PixelList] = FRET_Grain_Image(tag,tag,1,0.9,min_cluster);
tsh = 0.9;

t = sort(tag(:));
baseline = t(ceil(tsh*numel(t))); % the top 15% pixels are considered only
d_gray = tag>baseline; % binary threshold image for donor only intensity

bw = d_gray; % binary file containing the region of interests in both the images
cc = bwconncomp(bw, 4); % get clusters from the binary file

for i = 1:cc.NumObjects
    num(i) = numel(cc.PixelIdxList{i});
end
ind = num<min_cluster;

PixelList = cell(cc.NumObjects-sum(ind),1);
grain = false(size(d_gray));
ind = find(~ind);
for i =  1:numel(ind)
    PixelList{i} = cc.PixelIdxList{ind(i)};
     grain(cc.PixelIdxList{ind(i)})=true; %
end

tcspc_da_roi = zeros(Ngate,numel(PixelList)); % tcspc of donor acceptor photons for each region of interest
tcspc_d_roi = zeros(Ngate,numel(PixelList)); % tcspc of donor only photons for each region of interest
tau_da = zeros(numel(PixelList),1); tau_d = tau_da; dtau_da = tau_da; dtau_d = tau_da;
[~,pos] = max(tcspc);

disp('Determine lifetimes:'); reverseStr = '';
for i = 1:numel(PixelList)
    list_da = impix1(PixelList{i}); % List of pixels for donor acceptor
%     list_d  = impix2(PixelList{i}); % List of pixels for acceptor
    tmpmicrotime= zeros(Ngate,numel(list_da));
    for j = 1:numel(list_da)
        tmpmicrotime(:,j) =  mHist(microtime(im_pixel == list_da(j)),1:Ngate);
    end
    [tcspc_da_roi(:,i) tau_da(i) dtau_da(i)] = Bootstrap_tau(bin,tmpmicrotime_da,shift,ceil(0.9*numel(tmpmicrotime_da)),100);
    [tcspc_d_roi(:,i) tau_d(i) dtau_d(i)] = Bootstrap_tau(bin,tmpmicrotime_d,shift,ceil(0.9*numel(tmpmicrotime_d)),100);
    msg = sprintf('Processed %d/%d pixels...\n', i, numel(PixelList));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
 