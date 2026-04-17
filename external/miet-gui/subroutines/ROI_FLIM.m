function [head, FLIM, INT] = ROI_FLIM(name, IRF, pic)
% This program can be used to fit fluorescence lifetimes from multiple
% selected regions from the user.
%
% Input Parameters:
% name is the '.ht3' file which is to be analysed.
% IRF is the '.ht3' file which contains the Instrument Response. If not
% given, the program calculates an IRF using parametric equations 
% (see: Walther, K.A. et al, Mol. BioSyst. (2011) doi:10.1039/c0mb00132e)
% If pic=1 then the results are displayed.
%
% Output Paramters:    
% FLIM.t is the time axis for the tcspc histograms
% INT.tag is the intensity image
% INT.roi is the cell containing the coordinates of region of interest(s) selected
% INT.field is a matrix which contains separated single molecule patterns
% in different planes. The 4th dimension corresponds to different detection
% channels.
% FLIM.tcspc_sm contains the TCSPC histograms of all identified single molecules
% FLIM.photon_sm is the photon count of each identified single molecule.
% FLIM.life_mat is a cell containing the lifetime(s) of the identified molecules
% FLIM.Amp is the cell containing the Amplitudes of the fitted lifetime
% components.
% FLIM.lifeim is the lifetime image of the identified molecules. The lifetimes
% shown in this image are taken b weighing the lifetimes by their
% amplitudes.
% FLIM.tcspcIRF is the calculated IRF function
%
% Note: The program works only when there is one pulse in one laser sync
% period. 

% (c) Narain Karedla, 2014.

% name = 'U:\Narain\QDsAziModeFreespace\Image_004.ht3';
% IRF= 'U:\Narain\140305\Point_011.ht3';

    if nargin<3 || isempty(pic)
        pic=0;
    end

    if nargin<2 || isempty(IRF)
        IRF=[];
    end
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
    Resolution = max([maxres 0.064]);
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
    INT.tag=tag;
    FLIM.tcspc=tcspc;
    FLIM.t=t;
    
    disp('Computing IRF...') 
    tcspcIRF=[];
    if exist([name(1:end-4),'_tcspcIRF.mat'],'file')
        load([name(1:end-4),'_tcspcIRF.mat'],'tcspcIRF');
    else
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
                    tcspcIRF(:,i) = tcspcIRF(:,i)-min(tcspcIRF(:,i));
                end
            end
        else
            disp('invalid IRF, proceeding to calculate IRF')
            tcspcIRF = Calc_mIRF(head, tcspc(ch,:));
        end
        save([name(1:end-4),'_tcspcIRF.mat'],'tcspcIRF');
    end
    FLIM.tcspcIRF = tcspcIRF;

   
    run='y';
    c=0; % Number of regions selected. Do not touch it!
    while strcmpi(run,'y')
        c=c+1;
        close all
        wim(sum(tag,3));
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
        indd=field(:,3)>size(tag,1);
        field(indd,3)=size(tag,1);
        indu=field(:,2)<1;
        field(indu,2)=1;
        indl=field(:,1)<0;
        field(indl,:)=[];
        indr=field(:,1)>size(tag,2);
        field(indr,:)=[];
        for i=1:size(field,1)
            pixels=pixels+field(i,3)-field(i,2)+1;
        end
        field_val=zeros(pixels,maxch_n);  

        disp('Collecting photons...') 
        
        lifeimm=zeros(size(tag));
        tcspc_sm=zeros(Ngate,ch);
        for ch=1:maxch_n
           
            for i=1:size(field,1)
                tmp=[tmp; tag(field(i,2):field(i,3),xc(i))];
            end
            field_val(:,ch)=tmp; tmp=[];
            base = field_val(:,ch)>0.2.*(max(field_val(:,ch))-min(field_val(:,ch)))+min(field_val(:,ch)); % background filter
            %                  c=0;
            photon_sm(ch)=sum(field_val(base,ch));
            xval=[]; yval=[];
            for i=1:size(field,1)
             xval=[xval;repmat(field(i,1),numel(field(i,2):field(i,3)),1)];
             yval=[yval; (field(i,2):field(i,3)).'];
            end
            coord=[xval(base),yval(base)];

            for i=1:size(coord,1)
             ind1=im_line==coord(i,2);
             ind2=im_col==coord(i,1);
             tmp1 = double(im_tcspc).*(ind1.*ind2);
             tmp2 = double(im_chan).*(ind1.*ind2);
             tcspc_sm(:,ch)=tcspc_sm(:,ch)+mHist(double(tmp1(tmp2 == ch)),1:Ngate);
            end
            TCSPC_SM(:,c,ch)=tcspc_sm(:,ch);
            COORD{c,ch}=coord;
        end
        PHOTON_SM(:,c)=photon_sm;
        disp('Fitting lifetime curves...')
        for j=1
            %              [A, tau,~, ~, z, ~,~] = DistFluofit(tcspcIRF(ch,:), tcspc_sm(:,ch), 1./(head.SyncRate.*1e-9), Resolution,[],1); %intlife, [zeros(length(intlife),1);5.*ones(length(intlife),1)]);
            A=[]; tau= []; z=[];
            if nargin<2 || isempty(IRF)
                [A, tau,~, ~, taufit, ~,~] = DistFluofit_extension(tcspcIRF, tcspc_sm(:,j), 1./(head.SyncRate.*1e-9), Resolution,[],1); %intlife, [zeros(length(intlife),1);5.*ones(length(intlife),1)]);
            else
                if size(tcspcIRF,1)>size(tcspcIRF,2)
                    tcspcIRF=tcspcIRF.';
                end
                [A, tau,~, ~, taufit, ~,~] = DistFluofit_extension(abs(tcspcIRF(j,:)-mean(tcspcIRF(j,end-100:end))), tcspc_sm(:,j), 1./(head.SyncRate.*1e-9), Resolution,[],1); %intlife, [zeros(length(intlife),1);5.*ones(length(intlife),1)]);
            end
            
            life_mat=1./tau;
            Amp=A;
            life_matav=sum(life_mat.*Amp);
            disp('fitted lifetime(s) = '); 
            for i=1:numel(life_mat)
                disp([num2str(life_mat(i)),' ns'])
            end
            disp('Amplitude(s) = ');
            for i=1:numel(Amp)
                disp([num2str(Amp(i))])
            end
            disp(['The fitted average lifetime for the selected region of interest in channel ',num2str(ch),' out of ', num2str(maxch_n),' channel(s) is ',num2str(life_matav),' ns']);
            LIFE_MAT{c,j}=life_mat; AMP{c,j}=Amp;LIFE_MATAV(c,j)=life_matav;
            TAU_FIT(:,c,j) = taufit; 
            for i=1:size(coord,1)
                lifeimm(coord(i,2),coord(i,1),j)=life_matav(j);
            end
        end
        ROI{c}=roi; FIELD{c}=field; 
        LIFEIMM(:,:,:,c)=lifeimm;
        clear roi; clear field; clear tcspc_sm; clear life_mat; clear Amp; clear life_matav; clear lifeimm; clear coord;
        prompt='Press (y) to select another region of interest or (n) to quit analysis...';
        run=input(prompt,'s');           
    end
   
%    for ch=1:maxch_n
%        for x=1:size(LIFEIMM,1)
%            for y=1:size(LIFEIMM,2)
%                lifeim(x,y,ch)=pmean(LIFEIMM(x,y,ch,:));
%            end
%        end
%    end
     lifeim(:,:,:)=sum(LIFEIMM,4);


     INT.roi            = ROI;
     INT.field          = FIELD;
     INT.coord          = COORD;
     FLIM.tcspc_sm      = TCSPC_SM;
     FLIM.photon_sm     = PHOTON_SM;
     FLIM.life_mat      = LIFE_MAT;
     FLIM.amp           = AMP;
     FLIM.life_matav    = LIFE_MATAV;
     FLIM.lifeim        = lifeim;
     FLIM.taufit        = TAU_FIT;
  
     if pic==1
         close all
         figure
         set(gcf,'name','Regions Selected','NumberTitle','off')
         wim(sum(tag,3))
         hold on
         if exist('ROI','var')
             for i=1:c;
                 plot(ROI{i}(1,:),ROI{i}(2,:),'b.')
             end
         end
         hold off
         
         figure
         set (gcf,'name','Intensity Image','NumberTitle','off')
         for ch=1:maxch_n
             subplot(1,maxch_n,ch)
             imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), INT.tag(:,:,ch)) % check this!!
             axis equal
             axis tight
             colormap('hot')
         end
         
         figure
         set(gcf,'name','Lifetime Image','NumberTitle','off')
         for ch=1:maxch_n
             subplot(1,maxch_n,ch)
             imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny),FLIM.lifeim(:,:,ch))
             axis equal
             axis tight
             colormap('hot')
         end
     end
     save([name(1:end-4),'_Roi_Analysis.mat'],'head','INT','FLIM')
end
    
                      
     