name = 'D:\Luminosa\Data\Ciaofan\RawImage_2.ptu';

head = PTU_Read_Head(name);

out = PTU_MultiFrameScanReadFast(name); % read the channels

dind = out.dind;
ind = dind>=9 & dind<=31; % these are the PDA-23 detector channels
tags = out.tags(:,:,ind);

outNew = out;
outNew.dind = dind(ind);
outNew.tags = tags;

%% Pixel reassignment
params = struct();
params.imageSource = 'tags';      % use outPTU.tags(:,:,channel)
params.smoothSigma = 1;
params.useWindow = true;
params.normalizeImages = true;
params.upsampleReg = 20;

params.nIter = 1000;
params.checkEvery = 25;
params.stopTol = 1e-7;
params.minIter = 50;

params.pixelSize = head.ImgHdr_PixResol*1e3;            % nm/pixel in object plane
params.lambda = 690;              % nm
params.NA = 1.45;

params.showPlots = true;
results = run_ism_reconstruction_from_ptu(outNew, params);

flim = reassigned_flim(outNew, results, params);

int = results.acoImage;
lb = 10; % lower bound, 30 for cells, 10 for tissue
int = int-lb; %
int(int<0) = 0;
int = int./max(int(:)); % normalized intensity

%% Plotting figures
cmap = jet(256);
cmap = cmap(30:end-30,:);
tau = flim.total.tauRMS; % FLIM image ns, 10 for cells, 11 for tissue
% cim(tau ,log10(0.0001 + int.^3),[1.0 2.5],'v',cmap)
cim(flim.unassigned.total.tauRMS  , flim.unassigned.total.tag,[1 7],'v',cmap) % gamma correction
figure
cim(flim.reassigned.total.tauRMS  , results.acoImage,[1 7],'v',cmap) % gamma correction


%% TCSPC IRF calculation

tcspc_pix = flim.reassigned.total.xyT;   % same as flim.total.xyT
dtNs = outNew.head.MeasDesc_Resolution * 1e9;
pulsePeriodNs = outNew.head.MeasDesc_GlobalResolution * 1e9;
tAxisNs = flim.tAxisNs;
outNew.head.Resolution = dtNs;
outNew.head.SyncRate = outNew.head.TTResult_SyncRate;

tcspcIRF = Calc_mIRF(outNew.head, squeeze(sum(sum(tcspc_pix,1),2))');
tcspcIRF = squeeze(tcspcIRF);
tcspcIRF = tcspcIRF(:).';

%% Triexponential fitting and pattern matching
tau0 = [0.3 1.5 4.8];      % tri-exp
opts = struct();
opts.useGPU = true;
opts.mode = 'PIRLS';
opts.batchSize = 4000;
opts.includeBackground = true;
opts.pieIndex = 1;              % only relevant if tcspc_pix is 4D
opts.normalizePatterns = true;
opts.sortLifetimes = true;

outFLIM = GlobalMultiExpPatternMatchFromTCSPC( ...
    tcspc_pix, tcspcIRF, pulsePeriodNs, dtNs, tau0, opts);
