function res = lsFCS(name,cnum,maxtime,timegates,flagparallel)
%
% Input parameters :
%
% name          : .ht3 or .ptu file which is to be correlated
% cnum          : number of pulses in one sync
% irf           : the instrument response function recorded (.ht3 or .ptu file)
% maxtime       : the maximum correlation time
% components    : number of lifetime species present. The fitting takes place
%                 according to the components specified. If components == 0,
%                 the correlation takes place without any lifetime filters
%                 (normal FCS)
% cutoff        : If complete fitting of the tcspc curve is desired, cutoff
%                 must be []. For Time-gated-FLCS a positive cutoff given
%                 declares the time gate after which the tail fitting is
%                 performed.
% timegates     : Indicating the begining and ending of each pulse region. If
%                 given as [], they are automatically determined in the program
% flagparallel  : flag for parallel computing
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

% (c) Sebastian Isbaner & Narain Karedla 2018

close all
nEvents=2e6; % 1 million records in one bunch
if nargin < 2 || isempty(cnum)
    cnum=1;
end

if nargin<5||isempty(flagparallel)
    flagparallel = 0;
end

[bin, tcspcdata] = Harp_tcspc(name);
% [head, im_sync, im_tcspc, im_chan, im_line, im_col, im_frame] = PTU_LineScanRead(name, [1 1e6]);
% [dind,m]= unique(im_chan,'legacy');
% occurence = diff([0;vertcat(m)]);
head = PTU_Read_Head(name);
[~, ~, im_chan, special] = PTU_Read(name,[1 1e3],head);
[dind,m]= unique(sort(im_chan(special==0)+1),'legacy');
occurence = diff([0;vertcat(m)]);
dind(occurence<10)=[];
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

% if 1e5*maxtime*NRecords/acqtime>nEvents
%     nEvents = round(maxtime*NRecords/acqtime);
% end

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
for j=1:dnum
    for k=1:cnum
        c=c+1;
        tau(1:len,j,k) = bin(t1(c):t1(c)+len-1); %#ok<*AGROW>
        tcspc(1:len,j,k) = tcspcdata(t1(c):t1(c)+len-1,j);
    end
end

NChannels=ceil(Timeunit/Resolution);

tcspcdata2 = zeros(NChannels,dnum,ceil(NRecords/nEvents)); %NChannels: end of tcspc curves, equal to length(bin) or bin(end)+1, dnum: numer of detectors, last: number of record (also virtual) photons divided by bunch
semilogy(Resolution*bin,tcspcdata); drawnow
semilogy(Resolution*(1:length(tau(:,:,1))),reshape(tcspc,[length(tau(:,:,1)),cnum*dnum])); drawnow
res.tau = Resolution*tau;
res.tcspc = tcspc;
res.bin = bin;
res.tcspcdata = tcspcdata;
res.timegates = timegates;


if isfield(head,'CreatorSW_Name') && strcmpi(head.CreatorSW_Name,'Luminosa')
    nPixels = 50;
else
    nPixels = head.ImgHdr_PixX; %nx spatial
end

% nPixels = head.nx*head.ny;
% if nargin<5||isempty(para)
pixel_n = [1:nPixels].';

% else
%     pixelsub = 5;
%     pixelcasc = ceil(log2(nPixels/pixelsub));
%     pixel_n = cumsum(reshape(repmat(2.^(0:pixelcasc-1),pixelsub,1),pixelcasc*pixelsub,1));
%     ind = pixel_n>nPixels;
%     pixel_n(ind) = [];
%     %     pixel_n(end+1)=nPixel;
%     %     pixel_n(:,2) = [pixel_n(2:end,1); nPixels+1];
% end

% auto = zeros(numel(autotime),dnum*cnum*length(pixel_n),dnum*cnum*length(pixel_n));
rate = [];
time = 0;
cnt = 0;
num = 1;
iBunch = 1;

X = repmat(1:nPixels,[length(pixel_n) 1]);
tmp = floor(X./repmat(pixel_n,[1 nPixels]));

for i = 1:nPixels
    line_idx(i) = find(tmp(:,i),1,'last');
end

if flagparallel
    phot_read_idx = 0:nEvents:NRecords;
    phot_read_idx(end+1) = NRecords;

    npool = min([numel(phot_read_idx)-1,20]);

    MATLAB_2013b_or_newer = false;
    % Check MATLAB version
    MATLABversion = strsplit(version,'.');
    if(str2double(MATLABversion(1))>=8 && str2double(MATLABversion(2))>=2) % matlabpool -> parpool in MATLAB 2013b (8.2.x) and later
        MATLAB_2013b_or_newer = true;
    end
    if MATLAB_2013b_or_newer
        clust = parcluster;
        try
            p = gcp('nocreate');
            if isempty(p)
                nrRunningWorkers = 0;
            else
                nrRunningWorkers = p.NumWorkers;
            end

            if(nrRunningWorkers == 0)
                parpool(npool);
                p = gcp('nocreate');
                nrRunningWorkers = p.NumWorkers;
            end
            parallelProcessingAvailable = true;
            fprintf('Parallel processing available (%i workers).\n', nrRunningWorkers)
        catch
            parallelProcessingAvailable = false;
            fprintf(' Parallel processing unavailable.\n')
        end


        %           random_name = round(rand * 10000000);
        %     temp_name = [name(1:end-4) sprintf('_%i-%i.mat', random_name, jj)];

        addAttachedFiles(p,{'PTU_LineScanCorr.m','PTU_LineScanRead.m','PTU_Read.m','PTU_Read_Head.m','tttr2xfcs.m','cIntersect.m','cIntersect.cpp','cIntersect.mexw64','mexUtil.h'})




        for jj = 1:length(phot_read_idx)-1
            jobs(jj) = parfeval(p,@PTU_LineScanCorr, 1, ['\\jesrv2\AG-Enderlein' name(3:end)],  [phot_read_idx(jj)+1,nEvents],dind,cnum,pixel_n,tau,Timeunit,macroresolution,Ncasc,Nsub);
        end



        for jj = 1:length(phot_read_idx)-1
            %         get(jobs(jj), 'State')
            %             disp(jobs(jj).State)
            wait(jobs(jj))
        end

        auto = zeros(Ncasc*Nsub,length(pixel_n),length(phot_read_idx)-1);
        for jj = 1:length(phot_read_idx)-1
            outargs = fetchOutputs(jobs(jj));
            if dnum>1
                auto(:,:,jj) = squeeze(outargs(:,1,2,:)+outargs(:,2,1,:));
            else
                auto(:,:,jj) = squeeze(outargs);
            end
        end

    else % older versions of Matlab
        try
            nrRunningWorkers = matlabpool('size'); %#ok<*DPOOL>
            if(nrRunningWorkers == 0)
                matlabpool('open',npool);
            end
            parallelProcessingAvailable = true;
            fprintf('Parallel processing available (%i workers).\n', matlabpool('size'))
        catch
            parallelProcessingAvailable = false;
            fprintf('Parallel processing unavailable.\n')
        end
        clust = findResource('scheduler','type','jobmanager','Name','hal9001','LookupURL','134.76.92.49');

        %     random_name = round(rand * 10000000);
        %     temp_name = [name(1:end-4) sprintf('_%i-%i.mat', random_name, jj)];



        path_str = pwd;
        filedeps = cell(1,3);
        filedeps{1} = [path_str '\' 'PTU_LineScanRead.m'];
        filedeps{2} = [path_str '\' 'PTU_Read.m'];
        filedeps{3} = [path_str '\' 'PTU_Read_Head.m'];
        %     filedeps{4} = temp_name;


        for jj = 1:length(phot_read_idx)-1
            jobs{jj} = createJob(clust);
            set(jobs{jj}, 'JobData', jj);
            %         set(jobs{jj}, 'AdditionalPaths', {[path_str '\']});
            %         set(jobs{jj}, 'AttachedFiles', filedeps);
            set(jobs{jj}, 'FileDependencies', filedeps);
            createTask(jobs{jj}, @PTU_LineScanRead, 7, {name,  [phot_read_idx(jj)+1,nEvents]}); %file_name is short filename

            if (length(findJob(clust, 'State', 'queued')) < 20)
                submit(jobs{jj});
            end
        end
    end
    res.line_idx = line_idx;
    res.bin2 = tau;
    res.autotime = autotime;
    res.auto            = auto;
    head.NCounts        = cnt;
    res.head            = head;

    save([name(1:end-4),'_lsFCS.mat'],'res')


else

    while num>0

        if isfield(head,'CreatorSW_Name') && strcmpi(head.CreatorSW_Name,'Luminosa')
            [head, im_sync, im_tcspc, im_chan, ~, im_pixel, ~, num] = LPTU_LineScanRead(name, [1 + cnt, nEvents]);
        else
            [head, im_sync, im_tcspc, im_chan, ~, im_pixel, ~, num] = PTU_LineScanRead(name, [1 + cnt, nEvents]);
        end

        cnt = cnt + num;
        if (num>0) && (~isempty(im_sync))
            dt = im_sync(end)-im_sync(1);
            time(iBunch) = dt*Timeunit;
            flvp = zeros(length(im_chan),dnum*cnum*length(pixel_n));
            im_pixel = uint16(line_idx(im_pixel)).';


            r = 1;
            for iPixel=1:length(pixel_n)
                for j=1:dnum
                    for k=1:cnum
                        %                         flvp(:,r) = im_pixel==iPixel & im_chan==dind(k) & ...
                        %                             im_tcspc>=tau(1,j,k) & tau(end,j,k)>=im_tcspc; %flvp is the matrix that puts a 1 in the pixel that the photon corresponds to. Zo you get a matrix with "the number of photons in the bunch" rows and the number of pixels columns
                        %
                        flvp(:,r) = (im_pixel == iPixel) & ...
                            (im_chan  == uint8(dind(j))) & ...
                            (im_tcspc >= tau(1,j,k)) & ...
                            (im_tcspc <  tau(end,j,k) + 1 );
                        [tau(1,j,k) tau(end,j,k)]
                        r = r+1;
                    end
                end
            end
            if (num>0) && (~isempty(im_sync))
                for dd = 1:dnum
                    tcspcdata2(:,dd,iBunch) = tcspcdata2(:,dd,iBunch)+mHist(double(im_tcspc(im_chan==dind(dd))),0:NChannels-1);
                end
            end
            rate(iBunch,:) = sum(flvp)/time(iBunch);

            tPhoton=round((im_sync.*Timeunit)./macroresolution);

            auto(:,:,:,iBunch) = tttr2xfcsSym(tPhoton, flvp, Ncasc, Nsub); % main function that produces the G(k,k',tau) 3D matrix. with k the number of pixels. these are auto and cross correlations.

        end
%         subplot('position', [0.925 0.2 0.025 0.6]);
%         bar(0,(cnt/NRecords));
%         axis([-0.4 0.4 0 1]);
%         set(gca,'xtick',[],'ytick',[]);
%         subplot('position', [0.1 0.1 0.7 0.8]);
%         %     automean=mean(reshape(auto,size(auto,1),[]),2);
%         automean=0;
%         for ii = 1:size(auto,2)
%             automean=automean + auto(:,ii,ii,:);
%         end
%         automean = mean(squeeze(automean),2);
%         semilogx(autotime,automean)
%         drawnow;

        iBunch = iBunch+1;
    end


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

    save([name(1:end-4),'_lsFCS.mat'],'res')

end