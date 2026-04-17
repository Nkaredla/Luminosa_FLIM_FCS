function [data, tag, tau1, tau_c]=Process_scanfast(name, num_PIE, timegate, flagint, cutoff, flagfit,task, pic)
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
% pic - set to true if you want plots of the intensity and lifetime
% imagesclo
%
% Output:
% data - structure containing timegates, tcspc-curve without the "bump" at
%        the end (tcspc), intensity image (tag) and lifetime image (life_imm),
%        binnumber of the peak (pos), cutoff in ns (cutoff) and cutoff in bins (shift)
% tag  - intensity image
% life_imm - lifetime image

if nargin<8 || isempty(pic)
    pic=0;
end

if nargin<7 || isempty(task)
    task = 0; % 0 if only intensity image is needed
    % 1 if lifetime image computation is also required
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

sync_pixel = head.ImgHdr.PixelTime*1e-9.*head.SyncRate; % syncs per pixel
treshold   = 0.01*sync_pixel; % treshold per pixel for dead-time correction

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


tcspc_pix = zeros(nx,ny,Ngate,maxch_n*num_PIE);
tcspc_tot = zeros(nx,ny,length(tcspcdata),maxch_n);
for ch = 1:maxch_n
    ind = im_chan==dind(ch);
    tcspc_tot(:,:,:,ch) =  mHist3(double(im_line(ind)),double(im_col(ind)),double(im_tcspc(ind)),1:nx,1:ny,1:size(tcspcdata,1)); % tcspc histograms for all the pixels at once!
    for pie = 1:num_PIE
        if timegate((pie-1)*maxch_n+ch,3)>0
        tcspc_pix(:,:,:,(pie-1)*maxch_n+ch) = tcspc_tot(:,:, [timegate((pie-1)*maxch_n+ch,1):timegate((pie-1)*maxch_n+ch,2) timegate((pie-1)*maxch_n+ch,3):timegate((pie-1)*maxch_n+ch,4)],ch);
        else
          tcspc_pix(:,:,:,(pie-1)*maxch_n+ch) = tcspc_tot(:,:, [timegate((pie-1)*maxch_n+ch,1):timegate((pie-1)*maxch_n+ch,2)],ch);   
        end
    end
end
% tcspc_pix contains the photons from (L1 D1), (L1 D2), (L2 D1) and (L2 D2)

if flagint
     [~,t0] = max(tcspcdata);
     t0 = t0+ceil(cutoff./Resolution);
     tag = sum(tcspc_pix(:,:,t0:end),3);
else
    tag = sum(tcspc_pix,3);
end
bin = permute(repmat((1:Ngate)',[1 nx,ny,maxch_n*num_PIE]),[3,2,1,4]).*Resolution;
tau1 = real(sqrt((sum(bin.^2.*tcspc_pix,3)./tag)-(sum(bin.*tcspc_pix,3)./tag).^2));
data.tag = tag;
data.tau1 = tau1;
data.bin = bin;
data.tcspc_pix = tcspc_pix;
data.head = head;

if task
for ch = 1:maxch_n
    ind = im_chan==dind(ch);
    time = im_sync(ind)*ceil(1/head.SyncRate/Resolution/1e-9)+double(im_tcspc(ind)); % in tcspc bins
    
    apd = ceil(apd./Resolution); deadtime = ceil(deadtime./Resolution); % in number of bins
    
    %     time_pix = mHist3(double(im_line(ind)),double(im_col(ind)),double(time),1:nx,1:ny,1:Ngate);
    
    % interphoton distance distribution of all the pixels!
    ipd = mHist3(double(im_line(ind)),double(im_col(ind)),diff([0; time(:)]),1:nx,1:ny,1:10*ceil(1/head.SyncRate/Resolution/1e-9)+apd+deadtime+1);
    ipd(:,:,end) = [];
    
    % mN for all the pixels
    N = 1:10; mN = zeros(nx,ny,numel(N));
    for i = 1:10
        range = (i-1)*ceil(1/head.SyncRate/Resolution/1e-9)+1:i*ceil(1/head.SyncRate/Resolution/1e-9); range = range+apd+deadtime;
        range = range(range<size(ipd,3));
        mN(:,:,i) = sum(ipd(:,:,range),3);
    end
    
    clearvars ipd % Just too much memory to keep holding it
    
    % Intensity Image
    tag = sum(tcspc_tot,3);
    
    % average Lifetime Image by variance calculation
    bin = permute(repmat((1:Ngate)',[1 nx,ny]),[2,3,1]).*Resolution; % 3D time axis
    tau1 = real(sqrt((sum(bin.^2.*tcspc_pix,3)./tag)-(sum(bin.*tcspc_pix,3)./tag).^2));
    
    % getting first estimate for hit rates
    nN = permute(repmat(N(:),[1 nx, ny]),[2,3,1]);
    epsilon1 = sum(mN.*nN,3)./sum(mN,3);
    epsilon1 = 1./epsilon1;
    % here comes the slowest part of correcting for deadtimes
    back = zeros(size(tcspc_pix));
    tic
    disp('Determining epsilon:'); reverseStr = '';
    epsilon = zeros(ny,nx);
    for x = 1:nx
        for y = 1:ny
            if epsilon1(y,x)>0.1 && tag(y,x)>treshold % minimum threshold
                epsilon(y,x) = Simplex('ExpFun0',epsilon1(y,x),[],[],[],[],N,squeeze(mN(y,x,:))./Ngate,[],[],1);
                epsilon(y,x) = 1./epsilon(y,x);
                back(y,x,:) = DeadTimeCorrection(squeeze(tcspc_pix(y,x,:)./sum(tcspc_pix(y,x,:))),deadtime,apd,epsilon(y,x));
                
            end
            msg = sprintf('Processed %d/%d pixels...\n', (x-1)*ny+y, nx*ny);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            
        end
    end
    toc
    tau2 = real(sqrt((sum(bin.^2.*back,3)./sum(back,3))-(sum(bin.*back,3)./sum(back,3)).^2)); % average lifetimes after correction
    
    if flagfit
        dt = Resolution;
        [~,t0] = max(back,[],3);
        t0 = t0+ceil(cutoff./dt);
        tau = zeros(size(tau1)); tau_c = tau;
        disp('Determining lifetimes:'); reverseStr = '';
        for x = 1:nx
            for y = 1:ny
                if tag(y,x)>100
                    [c, p, ~, ~, ~, ~] = DistTailfit(tcspc_pix(y,x,t0(y,x):end), dt);
                    % [~,c] = ExpFun(p,[t0:ceil(tsync./dt)].*dt,meas(t0:end),[],[],1);
                    tau(y,x) = sum(c)/sum(c.*p);
                    
                    % p = Simplex('ExpFun',p,[],[],[],[],[t0:ceil(tsync./dt)].*dt,back(t0:end),[],[],1);
                    if epsilon(y,x)>0
                        [c, p, ~, ~, ~, ~] = DistTailfit(back(y,x,t0(y,x):end), dt);
                        % [~,c] = ExpFun(p,[t0:ceil(tsync./dt)].*dt,back(t0:end),[],[],1);
                        tau_c(y,x) = sum(c)/sum(c.*p);
                    end
                end
                msg = sprintf('Processed %d/%d pixels...\n', (x-1)*ny+y, nx*ny);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
        
    end
    
    
    
end

nind = tau2==0|isnan(tau2);
% tau_c = tau2;
tau2(nind) = tau1(nind); % patch the zeros and nan lifetime values with the non corrected lifetime values.

nind = ~isfinite(epsilon)| isnan(epsilon)|epsilon>0.5 |epsilon==0;
epsilon_c = epsilon;
epsilon_c(nind) = 0;

nind = epsilon_c>0.1;
tag2 = epsilon_c.*nind.*sync_pixel;
tag2 = max(tag,tag2); % patched intensity image

data.tag=tag;
data.tau1=tau1;
data.tau2 = tau2;
data.tag_c = tag2;
data.epsilon1 = epsilon1;
data.epsilon  = epsilon;
data.tcspc_pix = tcspc_pix;
data.mN  = mN;
data.back = back;
if flagfit
    data.life_imm = tau;
    data.life_imm_c = tau_c;
end
% data.pos=pos;
% data.shift=shift;
% data.cutoff=cutoff;
end

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
save([name(1:end-4),'_PSDT'],'data');
end

