function res = lsFLCS(name,cnum,maxtime,timegates,para,irf,components,cutoff)
%
% Input parameters :
%
% name      : .ht3 or .ptu file which is to be correlated
% cnum      : number of pulses in one sync
% irf       : the instrument response function recorded (.ht3 or .ptu file)
% maxtime   : the maximum correlation time
% components: number of lifetime species present. The fitting takes place
%             according to the components specified. If components == 0,
%             the correlation takes place without any lifetime filters
%             (normal FCS)
% cutoff    : If complete fitting of the tcspc curve is desired, cutoff
%             must be []. For Time-gated-FLCS a positive cutoff given
%             declares the time gate after which the tail fitting is
%             performed.
% timegates : Indicating the begining and ending of each pulse region. If
%             given as [], they are automatically determined in the program
% resolution: The tcspc timechannel resolution which is used for lifetime
%             fitting and FLCS calculation.
%
% Output parameters :
%
% res.tau       : time axis for the tcspc;
% res.tcspc     : the tcspc histograms for each detector and each pulse.
% res.autotime  : the time axis for the correlations.
% res.cutoff    : the cutoff given by the user
% res.timegates : the timegates on the tcspc timescale for separating the
%                 tcspc histogram for multiple pulses.
% res.auto      : the total auto correlation matrix
% res.automean  : mean autocorrelation of all the bunches
% res.rate      : the rate of photons detected in each bunch
% res.head      : the header of the file
% res.time      : the time duration of each bunch

% (c) Sebastian Isbaner and Narain Karedla 2017

close all
nEvents=1e6; % 1 million records in one bunch
if nargin < 2 || isempty(cnum)
    cnum=1;
end

[bin, tcspcdata] = Harp_tcspc(name);
tcspcdata = tcspcdata(1:end-800);
% [head, im_sync, im_tcspc, im_chan, im_line, im_col, im_frame] = PTU_LineScanRead(name, [1 1e5]);
head = PTU_Read_Head(name);
head.SyncRate = head.TTResult_SyncRate;
head.Resolution = head.MeasDesc_Resolution;
[~, ~, im_chan, special] = PTU_Read(name,[1 1e3],head);
[dind,m]= unique(im_chan(special==0)+1,'legacy');
occurence = diff([0;vertcat(m)]);
dind(occurence<10)=[];
%dind=1; added by akshita
dnum = length(dind);
Resolution = head.MeasDesc_Resolution;
acqtime = head.TTResult_StopAfter;
macroresolution = 1/head.ImgHdr_LineFrequency;
Timeunit = head.MeasDesc_GlobalResolution;
NRecords = head.TTResult_NumberOfRecords;

if nargin<3 || isempty(maxtime)
    maxtime  = 10; % seconds
end

% Ncasc=1;
% Nsub=round(maxtime./macroresolution);
Nsub = 10;
Ncasc = ceil(log2(maxtime/macroresolution/Nsub));
autotime = cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
autotime = autotime*macroresolution;

if 3e4*maxtime*NRecords/acqtime>nEvents
    nEvents = round(3e4*maxtime*NRecords/acqtime);
end

% timegates
if nargin<4 || isempty(timegates)
    [t1, len] = AutodetectTimeGates(tcspcdata, cnum);
    timegates(:,1) = t1;
    timegates(:,2) = len;
else
    t1 = timegates(:,1);
    len = min(timegates(:,2));
end
%
if numel(t1)==cnum
    t1 = repmat(t1(:),dnum,1);
end
c=0;
for j=1:cnum
    for k=1:dnum
        c=c+1;
        tau(1:len,j,k) = bin(t1(c):t1(c)+len-1); %#ok<*AGROW>
        tcspc(1:len,j,k) = tcspcdata(t1(c):t1(c)+len-1,k);
    end
end

NChannels=ceil(Timeunit/Resolution);

tcspcdata2 = zeros(NChannels,dnum,ceil(NRecords/nEvents)); %NChannels: end of tcspc curves, equal to length(bin) or bin(end)+1, dnum: numer of detectors, last: number of record (also virtual) photons divided by bunch
semilogy(Resolution*bin,tcspcdata); drawnow
semilogy(Resolution*(1:numel(tau(:,:,1))),tcspc(:,:,1)); drawnow
res.tau = Resolution*tau;
res.tcspc = tcspc;
res.bin = bin;
res.tcspcdata = tcspcdata;

% cntlist = [1:nEvents:NRecords];   
% parfor jj = 1:ceil(NRecords/nEvents)
%
%     PTU_LineScanRead_cluster(name, [cntlist(jj) cntlist(jj)-1]);
% end

%%
if components>0 % Fluorescence Lifetime Correlation
    if length(components)<5 % assuming 5 will be the maximum number of lifetime components
        
        if (nargin<8 || isempty(cutoff)) % Fitting with IRF (FLCS)
            disp('Getting IRF data...')
            cutoff = [];
            if components ==1
                tcspcirf = [];
                decay = []; amp = []; idx = 0;
                for i = 1:cnum
                    for j = 1:dnum
                        param(:,idx*(components+1)+1:(idx+1)*(components+1)) = [ones(numel(tcspc(:,i,j)),1) tcspc(:,i,j)];
                        idx = idx+1;
                    end
                end
            else
                if nargin<6 || isempty(irf)
                    for i = 1:cnum
                        for j = 1:dnum
                            tcspcirf(:,i,j) = Calc_mIRF(head, tcspc(:,i,j).');
                        end
                    end
                    tcspcirf = max(tcspc(:))*tcspcirf./max(tcspcirf(:));
                    
                    irftcspc = tcspcirf;
                    clear tcspcirf;
                    for j=1:cnum
                        for k=1:dnum
                            tmpirf = irftcspc(:,j,k);
                            cutirf = 1e-3*max(tmpirf);
                            tmpirf(tmpirf<=1.5*cutirf)=1.5*cutirf;
                            tcspcirf(1:len,j,k) = (tmpirf-1.5*cutirf);
                        end
                    end
                else
                    if strcmpi(irf(end-3:end),'.ht3') || strcmpi(irf(end-3:end),'.ptu')
                        [~,irftcspc] = Harp_tcspc(irf,Resolution);
                        %         [t1,lenirf] = AutodetectTimeGates(irftcspc,cnum);
                        irftcspc = irftcspc(:,1);
                        %tcspcirf =
                        %irftcspc-max(irftcspc(end-200:end));%added by
                        %akshita from fluofit
                        %tcspcirf(tcspcirf<0)=0;
                        
                        for j=1:cnum
                            for k=1:dnum
                                %                                 tmpirf =
                                %                                 irftcspc(mod(t1(j):t1(j)+len-1,NChannels)+1,k);
                                %                                 this script was used only when the irf
                                %                                 was recorded at a different frequency
                                %                                 than the data
                                tmpirf = irftcspc(t1(j):t1(j)+len-1,k);
                                %                                 cutirf = 1.5*round(mean(irftcspc(round(mod(t1(j)+len/2:t1(j)+len-1,NChannels)+1),k)));
                                %                                 tmpirf(tmpirf<=cutirf)=cutirf;
                                %                                 tcspcirf(1:len,j,k) = tmpirf-cutirf;
                                %                                 tcspcirf(tcspcirf<0)=0;
                                bgirf = round(mean(irftcspc(t1(j)+len/2:t1(j)+len-1,k))); %backgrd is taken from irf from length/2 to end and mean of it
                                cutirf = bgirf + 6*round(sqrt(bgirf)); %since background is possion dist so 3times the sqrtof backgrd is cut still 3% probability is possible
                                tmpirf(tmpirf<=cutirf)=bgirf;
                                tcspcirf(1:len,j,k) = tmpirf-bgirf;
                                tcspcirf(tcspcirf<0)=0;
                            end
                        end
                        ind=1:size(tcspcirf,1);
                        tmp=ind>len;
                        tcspcirf(tmp==1,:,:)=[];
                        
                    elseif ischar(irf)
                        % load saved irf.mat
                        load(irf,'tcspcirf');
                        tmpirf = tcspcirf;
                        clear tcspcirf
                        c = 0;
                        for j=1:cnum
                            for k=1:dnum
                                c = c+1;
                                %                                 tcspcirf(1:len,j,k) = tmpirf(mod(t1(c):t1(c)+len-1,NChannels)+1,k);
                                tcspcirf(1:len,j,k) = tmpirf(t1(c):t1(c)+len-1,k);
                            end
                        end
                    else
                        tcspcirf = irf;
                        tcspcirf=circshift(tcspcirf,0);
                    end
                    tcspcirf=circshift(tcspcirf,1);
                end
                
                disp('Fitting Lifetimes...')
                param = zeros([length(tau),(components+1)*cnum*dnum]);
                idx = 0;
                for i=1:cnum
                    for j = 1:dnum
                        [~, ~, A, life, ~, ~, ~, zz] = ...
                            Fluofit(circshift(tcspcirf(1:end,i,j),0), tcspc(1:end,i,j), Timeunit*1e9/cnum,Resolution.*1e9,[0.03 1 2 3]);
                        
                        [~,ord] = sort(life);
                        %decay(:,i,j) = life(ord(end-components:end-1)); %selecting middle two lifetime components(when vesicles)
                        decay(:,i,j) = life(ord(end-components+1:end)); % selecting only last components
                        %decay(:,i,j) = life(ord(:));
                        
                        param(:,idx*(components+1)+1:(idx+1)*(components+1)) = zz(:,[1;ord(end-components+1:end)+1]);
                        %param(:,idx*(components+1)+1:(idx+1)*(components+1)) = zz(:,[1;ord(end-components:end-1)+1]);
                        %param(:,idx*(components+1)+1:(idx+1)*(components+1)) = zz(:,[1;ord(:)+1]);
                        % decay(:,i,j)=life;
                        %param(:,idx*(components+1)+1:(idx+1)*(components+1)) = zz;
                        
                        amp(:,i,j) = A(end-components+1:end)./sum(A);
                        %amp(:,i,j) = A(end-components:end-1)./sum(A);
                        %amp(:,i,j) = A(:)./sum(A);
                        
                        idx = idx+1;
                    end
                end
            end
            
        else % Fitting without IRF (time gated FLCS)
            tcspcirf = [];
            [~,pos] = max(tcspc);
            pos     = max(pos(:));
            shift   = ceil(cutoff*1e-9/Resolution);
            
            ind = tau(:,1,1) - tau(1,1,1) > pos+shift;
            tau = tau(ind,:,:);
            tcspc = tcspc(ind,:,:);
            
            param = zeros([length(tau),(components+1)*cnum*dnum]);
            idx = 0;
            if components ==1
                decay = []; amp = [];
                for i = 1:cnum
                    for j = 1:dnum
                        param(:,idx*(components+1)+1:(idx+1)*(components+1)) = [ones(numel(tcspc(:,i,j)),1) tcspc(:,i,j)];
                        idx = idx+1;
                    end
                end
            else
                %                       p = ones(components,1);
                p = [2 3];
                disp('Fitting Lifetimes...')
                for i = 1:cnum
                    for j = 1:dnum
                        decay(:,i,j) = Simplex('ExpFunMLE',p,[],[],1e-30,inf,(tau(:,i,j)-tau(1,i,j)).*resolution*1e9,tcspc(:,i,j),1,1,1);
                        [~, a, zz,z] = ExpFunMLE(decay(:,i,j),(tau(:,i,j)-tau(1,i,j)).*resolution*1e9,tcspc(:,i,j),1,1,1);
                        [decay(:,i,j),ord] = sort(decay(:,i,j));
                        param(:,idx*(components+1)+1:(idx+1)*(components+1)) = zz(:,[1;ord+1]);
                        amp(:,i,j) = a(2:end)./sum(a(2:end));
                        idx = idx +1;
                    end
                end
            end
        end
    else
        %         para = para(:);
        %         decay = repmat(para,[1 dnum cnum]); amp = [];
        %         para = 1e9./para.';
        %         tcspcirf = [];
        %         [~,pos] = max(tcspc);
        %         pos     = max(pos(:));
        %         shift   = ceil(cutoff*1e-9/Resolution);
        %
        %         ind = tau(:,1,1) - tau(1,1,1) > pos+shift;
        %         tau = tau(ind,:,:);
        %         tcspc = tcspc(ind,:,:);
        %         param = repmat([ones(length(tau),1) exp(-((tau(:,1,1)-tau(1,1,1)).*Resolution)*para)],[1 dnum*cnum]);
        param = components;
    end
    
    for k = 1:size(param,2)
        param(:,k) = param(:,k)./sum(param(:,k));
    end
    
    %save([name(1:end-4),'_tcspc'],'res','head','t','y','ymean','y_filtered','ymean_filtered');
    % Calculating the FLCS filters
    disp('Calculating Filters...')
    c = 0;
    master = param.*0;
    for i=1:cnum
        for j = 1:dnum
            master(:,c*(components+1)+1:(c+1)*(components+1)) = ...
                ((param(:,c*(components+1)+1:(c+1)*(components+1))'*diag(1./tcspc(:,i,j))*param(:,c*(components+1)+1:(c+1)*(components+1)))\(diag(1./tcspc(:,i,j))*param(:,c*(components+1)+1:(c+1)*(components+1)))')';
            subplot(cnum,dnum,c+1)
            plot(master(:,c*(components+1)+1:(c+1)*(components+1)))
            title(['Detector ',num2str(j),' Pulse ',num2str(i)])
            axis tight
            drawnow
            c = c+1;
        end
    end
else % Correlation without Lifetime filters (FCS)
    param     = [];
    tcspcirf  = [];
    decay     = [];
    amp       = [];
    if (nargin<8 || isempty(cutoff))
        cutoff = [];
    else % Time Gated FCS
        [~,pos] = max(tcspc);
        pos     = max(pos(:));
        shift   = ceil(cutoff*1e-9/Resolution);
        ind = tau(:,1,1) - tau(1,1,1) > pos+shift;
        tau = tau(ind,:,:);
        tcspc = tcspc(ind,:,:);
    end
    master = ones(size(tcspc));
end



%%
nPixels = head.ImgHdr_PixX; %nx spatial
% nPixels = head.nx*head.ny;
if nargin<5||isempty(para)
    pixel_n = [1:nPixels].';
    %     pixel_n(:,2) = pixel_n+1;
else
    pixelsub = para;
    pixelcasc = ceil(log2(nPixels/pixelsub));
    pixel_n = cumsum(reshape(repmat(2.^(0:pixelcasc-1),pixelsub,1),pixelcasc*pixelsub,1));
    ind = pixel_n>nPixels;
    pixel_n(ind) = [];
    %     pixel_n(end+1)=nPixel;
    %     pixel_n(:,2) = [pixel_n(2:end,1); nPixels+1];
end



X = repmat(1:nPixels,[length(pixel_n) 1]);
tmp = floor(X./repmat(pixel_n,[1 nPixels]));

for i = 1:nPixels
    line_idx(i) = find(tmp(:,i),1,'last');
end



master = master(:,[1 3 4]);
components = size(master,2)-1;
if components<5
    master = reshape(repmat(master,[1 1 length(pixel_n)]),[size(master,1) (components+1)*length(pixel_n)]);
end


auto = zeros(numel(autotime),dnum*cnum*length(pixel_n)*(components+1),dnum*cnum*length(pixel_n)*(components+1));
rate = [];
time = 0;
cnt = 0;
num = 1;
iBunch = 1;

while num>0
    
    [head, im_sync, im_tcspc, im_chan, ~, im_pixel, ~, num] = PTU_LineScanRead(name, [1 + cnt, nEvents]);
    
    cnt = cnt + num;
    if (num>0) && (~isempty(im_sync))
        dt = im_sync(end)-im_sync(1);
        time(iBunch) = dt*Timeunit;
        im_pixel = uint16(line_idx(im_pixel)).';
        
        flvp = zeros(numel(im_tcspc),cnum*dnum*length(pixel_n)*(components+1));
        
        % weighing each photon by the weights calculated in master
        idx = 0;
        for iPixel=1:length(pixel_n)
            for j=1:cnum
                for k = 1:dnum
                    ind =  im_pixel==iPixel & im_chan==dind(k) & ...
                        im_tcspc>=tau(1,j,k) & tau(end,j,k)>=im_tcspc;
                    %                 ind = (tmpx>=tau(1,i,j) & tmpx<=tau(end,i,j) & flv == dind(j));
                    flvp(ind,idx*(components+1)+1:(idx+1)*(components+1)) = master(im_tcspc(ind)-tau(1,j,k)+1,idx*(components+1)+1:(idx+1)*(components+1));
                    rate(j,k,iBunch) = sum(ind)/time(iBunch);
                    idx = idx+1;
                end
            end
        end
        
        
        if (num>0) && (~isempty(im_sync))
            for dd = 1:dnum
                tcspcdata2(:,dd,iBunch) = tcspcdata2(:,dd,iBunch)+mHist(double(im_tcspc(im_chan==dind(dd))),0:NChannels-1);
            end
        end
        %         rate(iBunch,:) = sum(flvp)/time(iBunch);
        
        tPhoton=round((im_sync.*Timeunit)./macroresolution);
        
        auto(:,:,:,iBunch) = tttr2xfcs(tPhoton, flvp, Ncasc, Nsub); % main function that produces the G(k,k',tau) 3D matrix. with k the number of pixels. these are auto and cross correlations.
        
    end
    tmpres.line_idx = line_idx;
    tmpres.autotime = autotime;
    tmpres.auto = auto;
    
    if components>1
        [G Gcross t xxi] = lsCrossRead(tmpres);
        Gnorm = G./repmat(mean(G(end-7:end,:,:),1),size(G,1),1);
        Gcnorm = Gcross./repmat(mean(Gcross(end-7:end,:,:),1),size(Gcross,1),1);
        Gnorm(isnan(Gnorm))=1;
        Gcnorm(isnan(Gcnorm))=1;
        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,(cnt/NRecords));
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        CombineImages(cat(3,Gnorm(:,:,end-1:end),Gcnorm(:,:,end-1:end)),2,2)
        colormap('jet')
        drawnow
    else
        %         automean=mean(reshape(auto,size(auto,1),[]),2);
        %         automean=0;
        %         for ii = 1:size(auto,2)
        %             automean=automean + auto(:,ii,ii,:);
        %         end
        %         automean = mean(squeeze(automean),2);
        %         semilogx(autotime,automean)
        
        G  = lsCrossRead(tmpres);
        Gnorm = G./repmat(mean(G(end-7:end,:,:),1),size(G,1),1);
        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,(cnt/NRecords));
        subplot('position', [0.1 0.1 0.7 0.8]);
        imagesc(Gnorm(:,:,2))
        axis off
        colormap('jet')
        drawnow;
    end
    
    iBunch = iBunch+1;
    
    
    
    res.line_idx = line_idx;
    res.bin2 = tau;
    res.tcspcdata2      = tcspcdata2(t1:t1+len);
    
    res.autotime = autotime;
    automean=0;
    for ii = 1:size(auto,2)
        automean=automean + auto(:,ii,ii,:);
    end
    automean = mean(squeeze(automean),2);
    res.automean        = automean;
    res.auto            = auto;
    res.rate            = rate;
    res.time            = time;
    head.NCounts        = cnt;
    res.head            = head;
    
    save([name(1:end-4),'_lsFLCS.mat'],'res','-v7.3')
    
    
    
end

end
