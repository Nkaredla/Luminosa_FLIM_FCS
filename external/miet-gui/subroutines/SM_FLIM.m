function [head, FLIM, INT, SM] = SM_FLIM(name, IRF, pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, flag, flagauto, pic)
% This program identifies single molecule patterns and combines the pixels
% of the identified patterns to obtain single molecules lifetimes.
%
% Input parameters:
% name is the '.ht3' file which is to be analysed.
% IRF is the '.ht3' file which contains the Instrument Response. If not
% given, the program calculates an IRF using parametric equations
%(see: Walther, K.A. et al, Mol. BioSyst. (2011) doi:10.1039/c0mb00132e)
% pattern corresponds to the laser polarization used for scanning the
% molecules. By default, it takes a linear laser polarization for
% identification. pattern = 'azimuthal' or 'radial' for azimuthally or
% radially polarized beams respectively.
% NA is the numerical aperture of the objective used.
% n0 is the vector of refractive indices of layers below the sample from
% bottom to top.
% n1 is the vector of refractive indices of layers above the surface from
% bottom to top.
% n is the refractive index of the medium containing the sample molecules.
% d0 is the vector of the thicknesses of the layers below the sample layer.
% length (d0)= length(n0)-1.
% d1 is the vector of the thicknesses of the layers above the sample layer.
% length (d1)= length(n1)-1.
% d is the thickness of the layer containing the sample molecules.
% lamex is the wavelength used for excitation
% focpos is the amount of defocussing.
%
% Output paramters:
% INT.tag is the intensity image
% INT.imm_c is the identified pattern image
% INT.field_comb is a matrix which contains separated single molecule patterns
% in different planes. The 4th dimension corresponds to different detection
% channels.
% INT.mask is a matrix containing the computed patterns used for identifying
% the molecules.
%
% FLIM.t is the time axis for the tcspc histograms
% FLIM.tcspc_sm contains the TCSPC histograms of all identified single molecules
% FLIM.photon_sm is the photon count of each identified single molecule.
% FLIM.life_mat is a cell containing the lifetime(s) of the identified molecules
% FLIM.Amp is the cell containing the Amplitudes of the fitted lifetime
% components.
% FLIM.lifeimm is the lifetime image of the identified molecules. The lifetimes
% shown in this image are taken b weighing the lifetimes by their
% amplitudes.
% FLIM.tcspcIRF is the calculated IRF function
% FLIM.life_matav is the vector containing the average lifetimes of the single
% molecules identified.
% FLIM.fitz is the cell containing the fitted lifetime curves
%
% SM.theta is the vector containing the polar orientation of all the molecules identified
% SM.phi is the vector containing polar orientation of all the molecules identified
% SM.posx is the vector containing the x-coordinates of all the molecules
% identified.
% SM.posy is the vector containing the y-coordinates of all the molecules
% identified.
% SM.al_res is the resolution of the polar angle used for creating pattern
% masks for molecule identification
% SM.be_res is the resolution of the azimuthal angle used for creating pattern
% masks for molecule identification
%
% Note: The program works only when there is one pulse in one laser sync
% period.

% (c) Narain Karedla, 2014. Last modified -- 06/08/2014

% name='U:\Anna\27-01-2014-Microtubles\Image_019.ht3';

if nargin<15 || isempty(pic)
    pic=1;
end
if nargin<14 || isempty(flagauto)
    flagauto=0;
end
if nargin<13 || isempty(flag)
    flag='lifetime';
end
if nargin<12 || isempty(focpos)
    focpos=0;
end

if nargin<3 || isempty(pattern)
    pattern=[];
end

if nargin<2 || isempty(IRF)
    IRF=[];
end
tailfit = 0;
%% Intensity image only...

if strcmpi(flag,'intensity')
    disp('Reading image...')
    
    if exist([name(1:end-4),'_Core_Scan.mat'],'file')
        load([name(1:end-4),'_Core_Scan.mat'],'head','im_sync','im_tcspc','im_chan','im_line','im_col');
    else
        [head, im_sync, im_tcspc, im_chan, im_line, im_col] = Core_ScanRead(name);
        save([name(1:end-4),'_Core_Scan.mat'],'head','im_sync','im_tcspc','im_chan','im_line','im_col');
    end
    dind       = unique(im_chan);
    maxch      = numel(dind);
    maxres     = max([head.Resolution]);
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
        tcspc(ch,:) = mHist(double(im_tcspc(im_chan == dind(ch))),1:Ngate);
        indel(ch) = sum(tcspc(ch,:))<500; %deleting image planes containing less than 500 photons
    end
    tcspc(indel,:) = [];
    dind(indel)    = [];
    maxch_n        = numel(dind);
    
    tag = zeros(ny,nx,maxch_n); % intensity image
    for ch=1:maxch_n
        chan_ind=im_chan==dind(ch);
        line_temp=im_line(chan_ind);
        col_temp=im_col(chan_ind);
        for y = 1:ny
            ind = (line_temp == y);
            tmp1 = col_temp(ind);
            for x = 1:nx
                tag(y,x,ch)=sum(tmp1 == x);
            end
        end
        clear line_temp; clear col_temp; clear chan_ind;
    end
    save([name(1:end-4),'_intensity_image.mat'],'tag');
    INT.tag=tag;
    FLIM.tcspc=tcspc;
    FLIM.t=t;
    SM=[];
    
    if pic==1
        if maxch_n>0
            figure;
            set (gcf,'name','Intensity Image','NumberTitle','off')
            for ch=1:maxch_n
                subplot(1,maxch_n,ch)
                imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), INT.tag(:,:,ch)) % check this!!
                colorbar
                axis equal
                axis tight
                colormap('hot')
            end
        end
    end
end
%% Pattern matching only...

if strcmpi(flag,'pattern match')
    disp('Reading image...')
    
    if exist([name(1:end-4),'_Core_Scan.mat'],'file')
        load([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
    else
        [head, im_tcspc, im_chan, im_line, im_col] = Core_ScanRead(name);
        save([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
    end
    dind       = unique(im_chan);
    maxch      = numel(dind);
    maxres     = max([head.Resolution]);
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
        tcspc(ch,:) = mHist(double(im_tcspc(im_chan == dind(ch))),1:Ngate);
        indel(ch) = sum(tcspc(ch,:))<500; %#ok<*AGROW> % deleting image planes containing less than 500 photons
    end
    tcspc(indel,:) = [];
    dind(indel)    = [];
    maxch_n        = numel(dind);
    
    tag = zeros(nx,ny,maxch_n); % intensity image
    for ch=1:maxch_n
        chan_ind=im_chan==dind(ch);
        line_temp=im_line(chan_ind);
        col_temp=im_col(chan_ind);
        for y = 1:ny
            ind = (line_temp == y);
            tmp1 = col_temp(ind);
            for x = 1:nx
                tag(y,x,ch)=sum(tmp1 == x);
            end
        end
        clear line_temp; clear col_temp; clear chan_ind;
    end
    save([name(1:end-4),'_intensity_image.mat'],'tag');
    
    %            indchan=ones(maxch_n,1);
    
    for ch=1:maxch_n
        disp(['Analyzing Channel  ', num2str(ch),' Out of ',num2str(maxch_n),' Channel(s)...']);
        [field,~,imm_c, xc, yc, mask, al_res, be_res, theta_mol, phi_mol,err]=pattern_detect_auto(pixel,tag(:,:,ch),pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, flag, pic);
        if err==1
            INT=[]; SM=[];FLIM=[];
            return
        end
        SM.xc{ch}    = xc;
        SM.yc{ch}    = yc;
        SM.al_res{ch}  = al_res;
        SM.be_res{ch}  = be_res;
        SM.theta{ch} = theta_mol;
        SM.phi{ch}   = phi_mol;
        INT.field{ch}  = field;
        INT.mask{ch}   = mask;
        INT.imm_c{ch}  = imm_c;
    end
    
    INT.tag    = tag;
    FLIM.tcspc = tcspc;
    FLIM.t     = t;
    
    if pic==1
        if maxch_n>0
            close all
            figure
            set (gcf,'name','Intensity Image','NumberTitle','off')
            for ch=1:maxch_n
                subplot(1,maxch_n,ch)
                imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), INT.tag(:,:,ch)) % check this!!
                colorbar
                axis equal
                axis tight
                colormap('hot')
            end
            
            figure
            set (gcf,'name','Fitted Patterns','NumberTitle','off')
            for ch=1:maxch_n
                subplot(1,maxch_n,ch)
                %                         mim(INT.imm_c{ch})
                imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), INT.imm_c{ch})
                axis equal
                axis tight
                colormap('hot')
            end
        end
    end
    
end
%% Pattern matching and lifetime analysis...
if strcmpi(flag,'combined')
    disp('Reading image...')
    
    if exist([name(1:end-4),'_Core_Scan.mat'],'file')
        load([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
    else
        [head, im_tcspc, im_chan, im_line, im_col] = Core_ScanRead(name);
        save([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
    end
    dind       = unique(im_chan);
    maxch      = numel(dind);
    maxres     = max([head.Resolution]);
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
    
    
    %     for ch = 1:maxch
    %         tcspc(ch,:) = mHist(double(im_tcspc(im_chan == dind(ch))),1:Ngate);
    %         indel(ch) = sum(tcspc(ch,:))<500; %#ok<*AGROW> % deleting image planes containing less than 500 photons
    %     end
    %     tcspc(indel,:) = [];
    %     dind(indel)    = [];
    %     maxch_n        = numel(dind);
    
    % Here I change to get the intensity images after a cutoff from the user
    [tag tcspc dind] = Intensity_im(im_chan, im_line, im_col, im_tcspc, 0.3);
    maxch_n        = numel(dind);
    
    %     tag = zeros(nx,ny,maxch_n); % intensity image
    %     for ch=1:maxch_n
    %         chan_ind=im_chan==dind(ch);
    %         line_temp=im_line(chan_ind);
    %         col_temp=im_col(chan_ind);
    %         for y = 1:ny
    %             ind = (line_temp == y);
    %             tmp1 = col_temp(ind);
    %             for x = 1:nx
    %                 tag(y,x,ch)=sum(tmp1 == x);
    %             end
    %         end
    %         clear line_temp; clear col_temp; clear chan_ind;
    %     end
    save([name(1:end-4),'_intensity_image.mat'],'tag');
    
    %     indchan=ones(maxch_n,1);
    
    disp(['Fitting patterns...'])
    [field,~,imm_c, xc, yc, mask, al_res, be_res, theta_mol, phi_mol,err]=pattern_detect_auto(pixel,sum(tag,3),pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, flag, pic);
    if err==1
        INT=[]; SM=[];FLIM=[];
        return
    end
    if numel(xc)>0
        for j=1:size(field,3)
            field(:,:,j)=field(:,:,j)';
        end
        %               field_comb(:,:,:,i)=field;
        %               clear field
        if strcmpi(IRF,'internal')
            bcgind=imm_c==0;
            bcgind=bcgind.*1;
        else
            bcgind=[];
        end
        
        [~,fieldcell,tcspc_sm,tcspcIRF]=fastsm_tcspc(bcgind, field, Ngate, im_line, im_col, im_tcspc, uint8(ones(numel(im_tcspc,1))));
        
        photon_sm=sum(tcspc_sm(:,:),2);
        if exist([name(1:end-4),'_tcspcIRF.mat'],'file')
            load([name(1:end-4),'_tcspcIRF.mat'],'tcspcIRF');
        else
            if isempty(tcspcIRF)
                disp('Calculating IRF...')
                if nargin<2 || isempty(IRF)
                    tcspcIRF = Calc_mIRF(head, sum(tcspc_sm));
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
                
                [A, tau,~, ~, z, ~,~] = DistFluofit_extension(tcspcIRF, tcspc_sm(j,:), 1./(head.SyncRate.*1e-9), Resolution,[],1); %intlife, [zeros(length(intlife),1);5.*ones(length(intlife),1)]);
            else
                
                [A, tau,~, ~, z, ~,~] = DistFluofit_extension(tcspcIRF, tcspc_sm(j,:), 1./(head.SyncRate.*1e-9), Resolution,[],1); %intlife, [zeros(length(intlife),1);5.*ones(length(intlife),1)]);
            end
            life_mat{j}=1./tau; % check this!
            Amp{j}=A;
            fitz{j}=z;
            clear A; clear tau; clear z;
            comp=max(comp, numel(life_mat{j}(:)));
            life_matav(j)=sum(life_mat{j}.*Amp{j});
        end
        disp('Creating Lifetime Image...')
        lifeimmcell=cell(ny,nx);
        lifeimm=zeros(size(lifeimmcell));
        for x=1:ny
            for y=1:nx
                if sum(ismember(fieldcell{x,y},1:size(field,3)))>0
                    zind = fieldcell{x,y}(:);
                    ind = zind==0;
                    zind(ind) = [];
                    for i=1:numel(zind)
                        lifeimmcell{x,y}=life_matav(zind(i));
                    end
                    lifeimm(x,y)=mean(lifeimmcell{x,y}(:));
                end
            end
        end
        lifeimm=lifeimm';
        FLIM.life_mat   = life_mat;
        FLIM.lifeimm    = lifeimm;
        FLIM.Amp        = Amp;
        FLIM.tcspc_sm   = tcspc_sm;
        FLIM.photon_sm  = photon_sm;
        FLIM.tcspcIRF   = tcspcIRF;
        FLIM.life_matav = life_matav;
        FLIM.fitz       = fitz;
        
        INT.field      = field;
        INT.imm_c      = imm_c;
        INT.mask       = mask;
        
        SM.theta     = theta_mol;
        SM.phi       = phi_mol;
        SM.posx      = xc;
        SM.posy      = yc;
        SM.al_res    = al_res;
        SM.be_res    = be_res;
        clear life_mat; clear lifeimm; clear Amp; clear tcspc_sm; clear photon_sm; clear tcspcIRF; clear field; clear imm_c;clear mask;
        
    else
        SM = [];
        disp('No patterns found')
    end
    FLIM.tcspc = tcspc;
    FLIM.t     = t;
    INT.tag    = tag;
    disp ('File Processing Complete!')
    if isempty(SM)&& isempty(FLIM)&& isempty(INT)
        return
    else
        save([name(1:end-4) '_sm_flim.mat'],'FLIM','INT','SM','head');
    end
end


%% Pattern matching and lifetime analysis...
if strcmpi(flag,'lifetime')|| strcmpi(flag,'MIET') || strcmpi(flag,'random')
    disp('Reading image...')
    
    if exist([name(1:end-4),'_Core_Scan.mat'],'file')
        load([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
    else
        [head, im_tcspc, im_chan, im_line, im_col] = Core_ScanRead(name);
        save([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
    end
    dind       = unique(im_chan);
    maxch      = numel(dind);
    maxres     = max([head.Resolution]);
    Resolution = max([maxres 0.032]);
    chDiv      = Resolution/maxres;
    im_tcspc   = ceil(im_tcspc./chDiv);
    Ngate      = double(max(im_tcspc));
    pixel      = head.ImgHdr.PixelSize;
%     tcspc      = zeros(maxch, Ngate);
    nx         = head.ImgHdr.PixX;
    ny         = head.ImgHdr.PixY;
    t          = (1:Ngate).*Resolution;
    x0         = head.ImgHdr.X0;
    y0         = head.ImgHdr.Y0;
    
    [tag tcspc dind] = Intensity_im(im_chan, im_line, im_col, im_tcspc, 0.3);
    maxch_n        = numel(dind);
    %     for ch = 1:maxch
    %         tcspc(ch,:) = mHist(double(im_tcspc(im_chan == dind(ch))),1:Ngate);
    %         indel(ch) = sum(tcspc(ch,:))<500; %#ok<*AGROW> % deleting image planes containing less than 500 photons
    %     end
    %     tcspc(indel,:) = [];
    %     dind(indel)    = [];
    %     maxch_n        = numel(dind);
    %
    %     tag = zeros(nx,ny,maxch_n); % intensity image
    %     for ch=1:maxch_n
    %         chan_ind=im_chan==dind(ch);
    %         line_temp=im_line(chan_ind);
    %         col_temp=im_col(chan_ind);
    %         for y = 1:ny
    %             ind = (line_temp == y);
    %             tmp1 = col_temp(ind);
    %             for x = 1:nx
    %                 tag(y,x,ch)=sum(tmp1 == x);
    %             end
    %         end
    %         clear line_temp; clear col_temp; clear chan_ind;
    %     end
    save([name(1:end-4),'_intensity_image.mat'],'tag');
    
    indchan=ones(maxch_n,1);
    
    for ch=1:maxch_n
        disp(['Analyzing Channel  ', num2str(ch),' Out of ',num2str(maxch_n),' Channel(s)...']);
        [field,~,imm_c, xc, yc, mask, al_res, be_res, theta_mol, phi_mol,err]=pattern_detect_auto(pixel,tag(:,:,ch),pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, flag, pic);
        if err==1
            INT=[]; SM=[];FLIM=[];
            return
        end
        
        if numel(xc)>0
            for j=1:size(field,3)
                field(:,:,j)=field(:,:,j)';
            end
            %               field_comb(:,:,:,i)=field;
            %               clear field
            if strcmpi(IRF,'internal')
                bcgind=imm_c==0;
                bcgind=bcgind.*1;
            else
                bcgind=[];
            end
            
            
            [~,fieldcell,tcspc_sm,tcspcIRF]=fastsm_tcspc(bcgind,field, Ngate, im_line, im_col, im_tcspc, im_chan);
            photon_sm=sum(tcspc_sm(:,:),2);
            if tailfit
                comp = 0;
                life_mat=cell(size(field,3),1);
                life_matav=zeros(size(field,3),1);
                disp('Fitting Lifetimes...')
                %                 [timegate, len] = DetectTimeGates(tcspc, 1,Resolution);
                [~,posmax] = max(tcspc);
                Chanmin = posmax+round(0.3/Resolution);
                for j=1:size(life_mat,1)
                    [A, tau,~, z, ~,~] = DistTailfit(tcspc_sm(j,[Chanmin:end, 1:posmax-20]),Resolution,1,[],[],[],0);
                    life_mat{j}=1./tau; % check this!
                    Amp{j}=A;
                    fitz{j}=z;
                    life_matav(j)=sum(A./tau)./sum(A);
                    clear A; clear tau; clear z;
                    comp=max(comp, numel(life_mat{j}(:)));
                    
                end
            else
                if exist([name(1:end-4),'_tcspcIRF.mat'],'file')
                    load([name(1:end-4),'_tcspcIRF.mat'],'tcspcIRF');
                else
                    if isempty(tcspcIRF)
                        disp('Calculating IRF...')
                        if nargin<2 || isempty(IRF)
                            tcspcIRF = Calc_mIRF(head, tcspc(ch,:));
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
                for j=1:size(life_mat,1)
                    if nargin<2 || isempty(IRF)
                        
                        [A, tau,~, ~, z, ~,~] = DistFluofit_extension(abs(tcspcIRF-mean(tcspcIRF(end-100:end))), tcspc_sm(j,:), 1./(head.SyncRate.*1e-9), Resolution,[],1); %intlife, [zeros(length(intlife),1);5.*ones(length(intlife),1)]);
                    else
                        
                        [A, tau,~, ~, z, ~,~] = DistFluofit_extension(abs(tcspcIRF-mean(tcspcIRF(end-100:end))), tcspc_sm(j,:), 1./(head.SyncRate.*1e-9), Resolution,[],1); %intlife, [zeros(length(intlife),1);5.*ones(length(intlife),1)]);
                    end
                    life_mat{j}=1./tau; % check this!
                    Amp{j}=A;
                    fitz{j}=z;
                    clear A; clear tau; clear z;
                    comp=max(comp, numel(life_mat{j}(:)));
                    life_matav(j)=sum(life_mat{j}.*Amp{j});
                end
            end
            disp('Creating Lifetime Image...')
            lifeimmcell=cell(ny,nx);
            lifeimm=zeros(size(lifeimmcell));
            for x=1:ny
                for y=1:nx
                    if sum(ismember(fieldcell{x,y},1:size(field,3)))>0
                        zind = fieldcell{x,y}(:);
                        ind = zind==0;
                        zind(ind) = [];
                        for i=1:numel(zind)
                            lifeimmcell{x,y}=life_matav(zind(i));
                        end
                        lifeimm(x,y)=mean(lifeimmcell{x,y}(:));
                    end
                end
            end
            lifeimm=lifeimm';
            FLIM.life_mat{ch}  = life_mat;
            FLIM.lifeimm{ch}   = lifeimm;
            FLIM.Amp{ch}       = Amp;
            FLIM.tcspc_sm{ch}  = tcspc_sm;
            FLIM.photon_sm{ch} = photon_sm;
            FLIM.tcspcIRF{ch}  = tcspcIRF;
            FLIM.life_matav{ch}= life_matav;
            FLIM.fitz{ch}      = fitz;
            
            INT.field{ch}      = field;
            INT.imm_c{ch}      = imm_c;
            INT.mask{ch}       = mask;
            
            SM.theta{ch}      = theta_mol;
            SM.phi{ch}        = phi_mol;
            SM.posx{ch}       = xc;
            SM.posy{ch}       = yc;
            SM.al_res(ch)     = al_res;
            SM.be_res(ch)     = be_res;
            clear life_mat; clear lifeimm; clear Amp; clear tcspc_sm; clear photon_sm; clear tcspcIRF; clear field; clear imm_c;clear mask;
            
        else
            indchan(ch)=0;
            maxch_n=maxch_n-1;
        end
        
    end
    if maxch_n==0
        SM=[];
    end
    
    FLIM.tcspc = tcspc;
    FLIM.t     = t;
    INT.tag    = tag;
    
    if pic==1
        if maxch_n>0
            close all
            figure;
            set (gcf,'name','Intensity Image','NumberTitle','off')
            for ch=1:maxch_n
                subplot(1,maxch_n,ch)
                imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), INT.tag(:,:,ch)) % check this!!
                colorbar
                axis equal
                axis tight
                colormap('hot')
            end
            
            figure;
            set (gcf,'name','Lifetime Image','NumberTitle','off')
            for ch=1:maxch_n
                subplot(1,maxch_n,ch)
                %                         mim(FLIM.lifeimm{ch})
                imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), FLIM.lifeimm{ch})
                axis equal
                axis tight
                colormap('hot')
            end
            
            figure
            set (gcf,'name','Fitted Patterns','NumberTitle','off')
            for ch=1:maxch_n
                subplot(1,maxch_n,ch)
                %                         mim(INT.imm_c{ch})
                imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), INT.imm_c{ch})
                axis equal
                axis tight
                colormap('hot')
            end
            
            figure;
            set (gcf,'name','Photon Count vs Lifetimes','NumberTitle','off')
            for ch=1:maxch_n
                subplot(1,maxch_n,ch)
                plot(FLIM.life_matav{ch},FLIM.photon_sm{ch},'o','MarkerSize',2)
                xlabel('Lifetimes(ns)')
                ylabel('Photons per molecule')
            end
            
            if strcmpi(pattern,'radial')
                figure;
                set(gcf,'name','Polar Orientation vs Lifetimes','NumberTitle','off')
                for ch=1:maxch_n
                    subplot(1,maxch_n,ch)
                    plot(SM.theta{ch}.*180./pi, FLIM.life_matav{ch},'o','MarkerSize',2)
                    xlabel('Inclination (degrees)')
                    ylabel('Lifetime (ns)')
                end
            end
        end
    else
        disp('No patterns found')
    end
end
disp ('File Processing Complete!')
if isempty(SM)&& isempty(FLIM)&& isempty(INT)
    return
else
    save([name(1:end-4) '_sm_flim.mat'],'FLIM','INT','SM','head');
end
end