function res = PTU_ROIAnalysis(name,automat)

if nargin<2||isempty(automat)
    automat = 0; % asks the user for each roi if 0, else detects
end
if exist([name(1:end-4),'_FLIM_data.mat'],'file')
    load([name(1:end-4),'_FLIM_data.mat']);
elseif strcmp(name(end-2:end),'ptu')
    disp('Read tcspc-data from ptu-file...');
    PTU_ScanRead(name);
    % head = PTU_Read_Head(name);
    load([name(1:end-4),'_FLIM_data.mat']);
else
    disp('You have to give a ptu-file.');
    return
end

dind       = unique(im_chan);
maxch      = numel(dind);
maxres     = max(1e9*head.MeasDesc_Resolution);
Resolution = max([maxres 0.064]);
chDiv      = round(Resolution/maxres);
im_tcspc   = floor(im_tcspc./chDiv);
Ngate      = double(max(im_tcspc))-1;
pixel      = head.ImgHdr_PixResol;
tcspc      = zeros(maxch, Ngate);
nx         = head.ImgHdr_PixX;
ny         = head.ImgHdr_PixY;
t          = (1:Ngate).*Resolution;
x0         = head.ImgHdr_X0;
y0         = head.ImgHdr_Y0;
frames     = 29;

if chDiv>1
    for ch = 1:maxch
        sizebefore = size(tcspc_pix(:,:,:,ch)); % number of channels before
        num = chDiv*floor(sizebefore(3)/chDiv);
        sizemiddle = [sizebefore(1:2),chDiv,floor(sizebefore(3)/chDiv)];
        sizeafter = sizebefore;
        sizeafter(3) = floor(sizeafter(3)/chDiv);
        Tcspc_pix(:,:,:,ch) = reshape(sum(reshape(tcspc_pix(:,:,1:num,ch),sizemiddle),3),sizeafter);
    end
else
    Tcspc_pix = tcspc_pix;
end

for ch = 1:maxch
    tcspc(ch,:) = mHist(double(im_tcspc(im_chan == dind(ch))),1:Ngate);
    indel(ch) = sum(tcspc(ch,:))<500; %#ok<*AGROW> % deleting image planes containing less than 500 photons
end
tcspc(indel,:) = [];
dind(indel)    = [];
maxch_n        = numel(dind);


% lineared intensity, lifetime and tscpc matrices
tag_lin = reshape(tag,[nx*ny, maxch_n, frames]);
% tags_lin = reshape(tags,[nx*ny,maxch_n]);
tau_lin = reshape(tau,[nx*ny, maxch_n, frames]);
% taus_lin = reshape(taus,[nx*ny,maxch_n]);
% Tcspc_pix_lin = reshape(Tcspc_pix,[nx*ny,Ngate]);

% pixel index for all frames, channels and xy position
im_pixel = double(im_line)+double(im_col-1)*ny+double(im_chan-1)*nx*ny+double(im_frame-1)*nx*ny*maxch_n;

if automat
    [~,posmax] = max(tcspc,[],2);
    posmax = repmat(max(posmax),[1 maxch_n]);
    cutoff = 1; % 1 ns cutoff from the peak of the tcspc
    shift = ceil(cutoff./Resolution);
    
    % get a pattern from the image
    %     [~,max_pos]    = max(tags(:));
    %     [mposx, mposy] = ind2sub(size(tags),max_pos);
    %     [~,~,wx,wy]    = Gauss2D(tags(mposx-9:mposx+9,mposy-9:mposy+9));
    %     [x,y]          = meshgrid(-1.5*max(wx,wy):1.5*max(wx,wy));
    %     pixels         = round(mean([wx,wy]).^2);
    %     ww             = mean([wx,wy]);
    %     mask           = exp(-2*((x/ww).^2+(y/ww).^2));
    % guess a pattern
    ww   = 0.3/pixel;
    [x,y]          = meshgrid(-1.5*ww:1.5*ww);
    pixels         = round(ww.^2);
    mask = exp(-2*((x/ww).^2+(y/ww).^2));
    % find single molecules
    [~, ~, ~, ~, xc, yc, ~, ~, ~, ~, imm_sep,imm_c] = FindPatternRapid(sum(tags,3),mask,mask>1e-6,[],pixels,0.5);
    ind = xc<6| xc>size(tags,1)-6| yc<6| yc>size(tags,1)-6;
    xc(ind)  =[]; yc(ind) =[];
    sm_num = numel(xc);
    for i=1:sm_num
        imm_sep(:,:,i) = imm_sep(:,:,i)./sum(sum(imm_sep(:,:,i)));
    end
    ind_imm = imm_sep>1e-3;
    
    %     tcspc_pix_lin = reshape(Tcspc_pix,[numel(tags),Ngate]);
    disp('Processing Molecules:');  reverseStr = '';
    for i=1:sm_num
        lifeimm=0*tags;
        coord_list{i} = find(ind_imm(:,:,i));
        for ch = 1:maxch_n
            coord_frame = reshape(repmat(nx*ny*ch*(0:frames-1),[numel(coord_list{i}),1]),[numel(coord_list{i})*frames,1])+repmat(coord_list{i},[frames,1]);
            ind = ismember(im_pixel,coord_frame);
            tcspc_roi_frame(:,:,i,ch) = uint8(mHist2(double(im_tcspc(ind)),ceil(im_pixel(ind)./nx/ny),1:Ngate,1:frames));
            tcspc_roi(:,i,ch) = sum(tcspc_roi_frame(:,:,i,ch),1);
            tag_trace(i,ch,1:frames) = sum(tcspc_roi_frame(:,:,i,ch),2);
            tau_trace(i,ch,1:frames) = nansum(tau_lin(coord_list{i},ch,:).*tag_lin(coord_list{i},ch,:),1)./sum(tag_lin(coord_list{i},ch,:),1);
            COORD{i,ch}=coord_list{i};
        end
        
%         disp('Fitting lifetime curves...')
        for ch=1:maxch_n
            A=[]; decay= []; z=[];
            [Amp, decay,~,taufit] = DistTailfit(tcspc_roi(posmax(ch)+shift:end,i,ch),Resolution,1);
            life_mat=1./decay;
            life_matav = sum(Amp)./sum(Amp.*decay);
%             disp(['The average lifetime in channel ',num2str(ch),'/', num2str(maxch_n),' is ',num2str(life_matav),' ns']);
            LIFE_MAT{i,ch}=life_mat; AMP{i,ch}=Amp; LIFE_MATAV(i,ch)=life_matav;
            TAU_FIT(:,i,ch) = taufit;
            %         [posx posy posz] = ind2sub([ny,nx,maxch_n],coord+ny*nx*(ch-1));
            lifeimm(coord_list{i}+nx*ny*(ch-1))=life_matav;
        end
        FIELD{i}=ind_imm(:,:,i);
        LIFEIMM(:,:,:,i)=lifeimm;
        clear life_mat Amp life_matav lifeimm coord_list;
        
        msg = sprintf('Processed %d/%d molecules...\n', i, sm_num);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    ROI = [];
    photon_sm=squeeze(sum(tcspc_roi,1));
else
    
    % from here starts the roi analysis
    run='y';
    c=0; % Number of regions selected. Do not touch it!
    [~,posmax] = max(tcspc,[],2);
    cutoff = 1; % 1 ns cutoff from the peak of the tcspc
    shift = ceil(cutoff./Resolution);
    
    while strcmpi(run,'y')
        c=c+1;
        close all
        mim(sum(tags,3));
        hold on
        if exist('ROI','var')
            for i=1:c-1;
                plot(ROI{i}(1,:),ROI{i}(2,:),'b.')
            end
        end
        roi=[]; field=[];
        disp('Select region, you can drag across and resize. To complete selection, double click mouse left.')
        h= imellipse;
        wait(h);
        pos(c,:)=getPosition(h);
        clear h;
        pos_center(c,:)=[pos(c,1)+pos(c,3)/2, pos(c,2)+pos(c,4)/2];
        theta=0:pi/50:2*pi;
        roi=[pos_center(c,1)+(pos(c,3)/2.*cos(theta)); pos_center(c,2)+(pos(c,4)/2.*sin(theta))];
        roi=round(roi);
        plot(roi(1,:),roi(2,:),'r.')
        hold off
        xc=unique(roi(1,:));
        pixels=0;
        for i=1:numel(xc)
            field(i,:)=[xc(i),min(roi(2,roi(1,:)==xc(i))), max(roi(2,roi(1,:)==xc(i)))];
        end
        tmp=[];
        indd=field(:,3)>size(tags,1);
        field(indd,3)=size(tags,1);
        indu=field(:,2)<1;
        field(indu,2)=1;
        indl=field(:,1)<0;
        field(indl,:)=[];
        indr=field(:,1)>size(tags,2);
        field(indr,:)=[];
        for i=1:size(field,1)
            pixels=pixels+field(i,3)-field(i,2)+1;
        end
        field_val=zeros(pixels,maxch_n);
        
        disp('Collecting photons...')
        
        lifeimm=zeros(size(tags));
        for ch=1:maxch_n
            xval=[]; yval=[];
            for i=1:size(field,1)
                xval=[xval;repmat(field(i,1),numel(field(i,2):field(i,3)),1)];
                yval=[yval; (field(i,2):field(i,3)).'];
            end
            coord=sub2ind([ny nx],yval,xval);
            % all coordinates in each frame for the detector channel
            coord_frame = reshape(repmat(nx*ny*ch*(0:frames-1),[numel(coord),1]),[numel(coord)*frames,1])+repmat(coord,[frames,1]);
            
            ind = ismember(im_pixel,coord_frame);
            tcspc_roi_frame(:,:,c,ch) = uint8(mHist2(double(im_tcspc(ind)),ceil(im_pixel(ind)./nx/ny),1:Ngate,1:frames));
            
            tcspc_roi(:,c,ch) = sum(tcspc_roi_frame(:,:,c,ch),1);
            tag_trace(c,1:frames) = sum(tcspc_roi_frame(:,:,c,ch),2);
            tau_trace(c,1:frames) = nansum(tau_lin(coord,:).*tag_lin(coord,:),1)./sum(tag_lin(coord,:),1);
            COORD{c,ch}=coord;
            
        end
        
        disp('Fitting lifetime curves...')
        for ch=1:maxch_n
            A=[]; decay= []; z=[];
            [Amp, decay,~,taufit] = DistTailfit(tcspc_roi(posmax(ch)+shift:end,c,ch),Resolution,1);
            
            life_mat=1./decay;
            
            life_matav = sum(Amp)./sum(Amp.*decay);
            disp(['The average lifetime in channel ',num2str(ch),'/', num2str(maxch_n),' is ',num2str(life_matav),' ns']);
            LIFE_MAT{c,ch}=life_mat; AMP{c,ch}=Amp; LIFE_MATAV(c,ch)=life_matav;
            TAU_FIT(:,c,ch) = taufit;
            %         [posx posy posz] = ind2sub([ny,nx,maxch_n],coord+ny*nx*(ch-1));
            lifeimm(coord+nx*ny*(ch-1))=life_matav;
        end
        ROI{c}=roi; FIELD{c}=field;
        LIFEIMM(:,:,:,c)=lifeimm;
        clear roi; clear field; clear life_mat; clear Amp; clear life_matav; clear lifeimm; clear coord;
        prompt='Press (y) to select another region of interest or (n) to quit analysis...';
        run=input(prompt,'s');
    end
    photon_sm=sum(tcspc_roi,3);
end

res.tag = tag;
res.tags = tags;
res.tau = tau;
res.taus = taus;
res.tcspc_roi = tcspc_roi;
res.field = FIELD; % note, when you save a matrix into a cell, use mim(,fliplr(flipud(tmp.'))) to view the data; where tmp = FIELD{i}
res.coord = COORD;
res.tcspc_roi_frame = tcspc_roi_frame;
res.tag_trace = tag_trace;
res.tau_trace = tau_trace;
res.life_mat = LIFE_MAT;
res.life_matav = LIFE_MATAV;
res.tau_fit = TAU_FIT;
res.amp  = AMP;
res.roi = ROI;
res.lifeimm = LIFEIMM;
res.photon_sm = photon_sm;
res.Resolution = Resolution;
res.head = head;
res.t = t;
res.tcspc_pix = tcspc_pix;