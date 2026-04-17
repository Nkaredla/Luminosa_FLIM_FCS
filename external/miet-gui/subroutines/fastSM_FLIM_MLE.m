function[head, FLIM, INT, SM]= fastSM_FLIM_MLE(name, IRF, pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, flagint, pic)

if nargin<14 || isempty(pic)
    pic=1;
end
if nargin<13 || isempty(flagint)
    flagint=1;
end
if nargin<12 || isempty(focpos)
    focpos=0;
end

if nargin<3 || isempty(pattern)
    pattern=[];
end

if nargin<2 || isempty(IRF)
    IRF=[];
    tcspcIRF = [];
end
 if ischar(d)
        d0 = [d0 0];
        d = 0;
 end
    
tailfit = 0; flagfit = 1; route = 1; % DistTail fitting
cutoff = 0.5; % 0.5 ns cutoff for tailfitting or for intensity image filtering

disp('Reading image...')

if exist([name(1:end-4),'_Core_Scan.mat'],'file')
    load([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
else
    [head, im_sync, im_tcspc, im_chan, im_line, im_col] = PTU_ScanRead(name);
    save([name(1:end-4),'_PTU_Scan.mat'],'head','im_sync','im_tcspc','im_chan','im_line','im_col');
end
dind       = unique(im_chan);
maxch      = numel(dind);
maxres     = max([head.MeasDesc_Resolution]);
Resolution = max([maxres 0.032]);
chDiv      = Resolution/maxres;
im_tcspc   = ceil(im_tcspc./chDiv);
Ngate      = double(max(im_tcspc));
pixel      = head.ImgHdr.PixelSize;
tcspc      = zeros(maxch, Ngate);
nx         = head.ImgHdr.PixX;
ny         = head.ImgHdr.PixY;
t          = (1:Ngate).*Resolution;
x0         = head.ImgHdr.X0;
y0         = head.ImgHdr.Y0;

for ch = 1:maxch
    tcspcdata(ch,:) = mHist(double(im_tcspc(im_chan == dind(ch))),1:Ngate);
    indel(ch) = sum(tcspcdata(ch,:))<500; %deleting image planes containing less than 500 photons
end
tcspcdata(indel,:) = [];
dind(indel)    = [];
maxch_n        = numel(dind);


%     disp('Determine time-gates...');
%     [timegate, Ngate] = DetectTimeGates(tcspcdata, 1, Resolution);
%
Ngate = size(tcspcdata,2);
tcspc_pix = zeros(ny,nx,Ngate,maxch_n); tag = zeros(ny,nx,maxch_n);
for ch = 1:maxch_n
    ind = im_chan==dind(ch);
    tcspc_pix(:,:,:,ch) =  mHist3(double(im_line(ind)),double(im_col(ind)),double(im_tcspc(ind)),1:nx,1:ny,1:Ngate); % tcspc histograms for all the pixels at once!
    %     bin = permute(repmat((1:Ngate)',[1 nx,ny]),[2,3,1]).*Resolution; % 3D time axis
    if flagint
        [~,t0] = max(tcspcdata(ch,:));
        t0 = t0+ceil(cutoff./Resolution);
        tag(:,:,ch) = sum(tcspc_pix(:,:,t0:t0+ceil(10/Resolution),ch),3);
    else
        tag(:,:,ch) = sum(tcspc_pix(:,:,:,ch),3);
    end
end
save([name(1:end-4),'_intensity_image.mat'],'tag');

disp(['Fitting patterns...'])
[field,~,imm_c, xc, yc, mask, al_res, be_res, theta_mol, phi_mol,err]=pattern_detect_auto(pixel,sum(tag,3),pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, 'MIET', pic);
if err==1
    INT=[]; SM=[];FLIM=[];
    return
end
[tcspc_sm, photon_sm] = Collect_tcspc(tcspc_pix,field);

% here comes the MLE part

if flagfit
    dt = Resolution;
    [~,t0] = max(tcspc_sm);
    t0 = t0+ceil(cutoff./dt);
    tau = zeros(numel(xc));
    if route
        disp('Determining lifetimes')
        if exist([name(1:end-4),'_tcspcIRF.mat'],'file')
            load([name(1:end-4),'_tcspcIRF.mat'],'tcspcIRF');
        else
            if isempty(tcspcIRF)
                disp('Calculating IRF...')
                if nargin<2 || isempty(IRF)
                    tcspcIRF = Calc_mIRF(head, sum(tcspc_sm.'));
                elseif strcmp(IRF(end-2:end),'ht3')
                    [binIRF, tcspcIRF, ~] = ht3TCSPC(IRF);
                    if numel(binIRF) ~= size(tcspc,2)
                        binIRF = binIRF+1;
                        num1 = numel(binIRF); p=size(tcspc,2);
                        Div = round(num1./p);
                        bin = 1:p;
                        tcspcIRF_n=zeros(numel(bin),size(tcspcIRF,2));
                        for i=0:numel(bin)-1
                            for j=1:size(tcspcIRF,2)
                                tcspcIRF_n(i+1,j)=sum(tcspcIRF(i*Div+1:(i+1)*Div,j));
                            end
                        end
                        tcspcIRF=tcspcIRF_n;
                        clear tcspcIRF_n;
                        for i=1:size(tcspcIRF,2)
                            %                                tmp1 = max(tcspcIRF(:,i));
                            %                                tmp2 = min(tcspcIRF(:,i));
                            %                                ord = ceil(log10(tmp1/tmp2));
                            %                                indtc =tcspcIRF<max(tcspcIRF).*(1/(10.^(ord-1)));
                            %                                tcspcIRF(indtc) = max(tcspcIRF).*(1/(10.^(ord-1)));
                            tcspcIRF(:,i) = tcspcIRF(:,i)-min(tcspcIRF(:,i));
                        end
                    end
                elseif size(IRF,1)==maxch_n && size(IRF,2)==Ngate
                    tcspcIRF=IRF;
                else
                    disp('invalid IRF, proceeding to calculate IRF')
                    tcspcIRF = Calc_mIRF(head, tcspc(ch,:));
                end
            end
            save([name(1:end-4),'_tcspcIRF.mat'],'tcspcIRF');
        end
        comp = 0;
        life_mat=cell(size(field,3),1);
        life_matav=zeros(size(field,3),1);
        disp('Fitting Lifetimes...')
        if size(tcspcIRF,1)>size(tcspcIRF,2)
            tcspcIRF=tcspcIRF.';
        end
        
        tcspcIRF = sum(tcspcIRF,3);
        
        for j=1:size(life_mat,1)
            if nargin<2 || isempty(IRF)
                
                [A, tau,~, ~, z, ~,~] = DistFluofit_extension(tcspcIRF, tcspc_sm(:,j), 1./(head.SyncRate.*1e-9), Resolution,[],1); %intlife, [zeros(length(intlife),1);5.*ones(length(intlife),1)]);
            else
                
                [A, tau,~, ~, z, ~,~] = DistFluofit_extension(tcspcIRF, tcspc_sm(:,j), 1./(head.SyncRate.*1e-9), Resolution,[],1); %intlife, [zeros(length(intlife),1);5.*ones(length(intlife),1)]);
            end
            life_mat{j}=1./tau; % check this!
            Amp{j}=A;
            fitz{j}=z;
            clear A; clear tau; clear z;
            comp=max(comp, numel(life_mat{j}(:)));
            life_matav(j)=sum(life_mat{j}.*Amp{j});
        end
    else
        comp = 0;
        life_mat=cell(size(field,3),1);
        life_matav=zeros(size(field,3),1);
        disp('Fitting Lifetimes...')
        for j=1:size(life_mat,1)
            if nargin<2 || isempty(IRF)
                
                [A, tau,~, ~, z, ~] = DistTailfit(tcspc_sm(t0(j):end,j), Resolution,1); %intlife, [zeros(length(intlife),1);5.*ones(length(intlife),1)]);
            else
                
                [A, tau,~, ~, z, ~] = DistTailfit(tcspc_sm(t0(j):end,j), Resolution,1);
            end
            life_mat{j}=1./tau; % check this!
            Amp{j}=A;
            fitz{j}=z;
            life_matav(j)=sum(A)./sum(tau.*A);
            clear A; clear tau; clear z;
            comp=max(comp, numel(life_mat{j}(:)));
        end
    end
    disp('Creating Lifetime Image...')
    
    lifeimm=zeros(size(field));
    
    for i = 1:numel(xc)
        lifeimm(:,:,i)=field(:,:,i).*life_matav(i);
    end
    lifeimm1= sum(lifeimm,3)./sum(field,3);
   
    %         for x=1:ny
    %             for y=1:nx
    %                 if sum(ismember(fieldcell{x,y},1:size(field,3)))>0
    %                     zind = fieldcell{x,y}(:);
    %                     ind = zind==0;
    %                     zind(ind) = [];
    %                     for i=1:numel(zind)
    %                         lifeimmcell{x,y}=life_matav(zind(i));
    %                     end
    %                     lifeimm(x,y)=mean(lifeimmcell{x,y}(:));
    %                 end
    %             end
    %         end
    
    
    %         for x = 1:nx
    %             for y = 1:ny
    %                 if tag(y,x)>100
    %                     [c, p, ~, ~, ~, ~] = DistTailfit(tcspc_pix(y,x,t0(y,x):end), dt);
    %                     tau(y,x) = sum(c)/sum(c.*p);
    %
    %                 end
    %
    %                 msg = sprintf('Processed %d/%d pixels...\n', (x-1)*ny+y, nx*ny);
    %                 fprintf([reverseStr, msg]);
    %                 reverseStr = repmat(sprintf('\b'), 1, length(msg));
    %             end
    %         end
    %
else
    bin = permute(repmat((1:Ngate)',[1 nx,ny]),[2,3,1]).*Resolution; % 3D time axis
    
    lifeimm = real(sqrt((sum(bin.^2.*tcspc_pix,3)./tag)-(sum(bin.*tcspc_pix,3)./tag).^2));
    
end


FLIM.life_mat    = life_mat;
% FLIM.life_mat2   = life_mat2;
FLIM.lifeimm1    = lifeimm1;
% FLIM.lifeimm2    = lifeimm2;
FLIM.Amp         = Amp;
% FLIM.Amp2        = Amp2;
FLIM.tcspc_sm    = tcspc_sm;
FLIM.photon_sm   = photon_sm;
FLIM.tcspcIRF    = tcspcIRF;
FLIM.life_matav  = life_matav;
% FLIM.life_matav2 = life_matav2;
FLIM.fitz        = fitz;
FLIM.tcspc       = tcspc;
FLIM.t           = t;

INT.field      = field;
INT.imm_c      = imm_c;
INT.mask       = mask;
INT.tag=tag;

SM.theta     = theta_mol;
SM.phi       = phi_mol;
SM.posx      = xc;
SM.posy      = yc;
SM.al_res    = al_res;
SM.be_res    = be_res;

disp ('File Processing Complete!')
if isempty(SM)&& isempty(FLIM)&& isempty(INT)
    return
else
    save([name(1:end-4) '_sm_flim.mat'],'FLIM','INT','SM','head');
end
end
