function [head, FLIM, INT] = FLIM_Analysis(name, IRF, num_PIE, pic)
% This program analyzes FLIM image and calculates the lifetime image by
% fitting lifetime curves for each pixel which has > 30 photons.
% name is the '.ht3' file which is to be analysed.
% IRF is the '.ht3' file which contains the Instrument Response. If not
% given, the program calculates an IRF using parametric equations 
%(see: Walther, K.A. et al, Mol. BioSyst. (2011) doi:10.1039/c0mb00132e)
% num_PIE is the number of pulses in one TCSPC cycle  

% Output paramters:    
% INT.tag is the intensity image
% INT.tag_select is the intensity image for the ROI selected.
% FLIM.t is the time axis for the tcspc histograms
% FLIM.life_cell is a cell containing the lifetime(s) of each pixel
% FLIM.Amp is the cell containing the Amplitudes of the fitted lifetime
% components.
% FLIM.life_imm is the lifetime image. The lifetimes
% shown in this image are taken by weighing the lifetimes by their
% amplitudes.
% FLIM.tcspcIRF is the calculated IRF function

% (c) Narain Karedla, 2014.

% name='U:\Olaf\MicroTime\Olaf_2-11-2014\Area_28.ht3';

if nargin<4 || isempty(pic)
    pic=1;
end

if nargin<3 || isempty(num_PIE)
    num_PIE=1;
end

if nargin<2 || isempty(IRF)
    IRF=[];
end



disp('Reading image...')    
[head, im_tcspc, im_chan, im_line, im_col] = Core_ScanRead(name,num_PIE);
    maxch      = max(im_chan);
    maxres     = max([head.Resolution]);
    Resolution = max([maxres 0.032]);
    chDiv      = Resolution/maxres;
    im_tcspc   = ceil(im_tcspc./chDiv);
    Ngate      = double(max(im_tcspc));
%     pixel      = head.ImgHdr.PixelSize;
    tcspc      = zeros(maxch, Ngate);    
    nx         = head.ImgHdr.PixX;                                             
    ny         = head.ImgHdr.PixY;
    t          = (1:Ngate).*Resolution;    
    dind       = unique(im_chan);
    
    
    for ch = 1:maxch
        tcspc(ch,:) = mHist(double(im_tcspc(im_chan == dind(ch))),1:Ngate);
        indel(ch) = sum(tcspc(ch,:))<500; %#ok<*AGROW> % deleting channels containing less than 500 photons
    end
   tcspc(indel,:)=[];
   dind(indel)=[];
   maxch_n=size(tcspc,1);
   
   tag = zeros(nx,ny,maxch_n); % intensity image
   for i=1:maxch
       chan_ind=im_chan==dind(i);
       line_temp=im_line(chan_ind);
       col_temp=im_col(chan_ind);
       for y = 1:ny
           ind = (line_temp == y);
           tmp1 = col_temp(ind);
           for x = 1:nx
               tag(y,x,i)=sum(tmp1 == x);
           end
       end
       clear line_temp; clear col_temp; clear chan_ind;
   end
   tag(:,:,indel)=[];

%   tcspc_pix=zeros(size(tag,1),size(tag,2),maxch_n,Ngate);
  
  
  for ch=1:maxch_n
     disp(['Analyzing Channel  ', num2str(ch),' Out of ',num2str(maxch_n),' Channel(s)...']);
     disp('Computing IRF...')
         if nargin<2 || isempty(IRF)
             tcspcIRF = Calc_mIRF(head, tcspc(ch,:));                      
         elseif strcmp(IRF(end-2:end),'ht3')
                    [binIRF, tcspcIRF, headIRF] = ht3TCSPC(IRF); %#ok<*NASGU>
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
         else
             disp('invalid IRF, proceeding to calculate IRF')
             tcspcIRF = Calc_mIRF(head, tcspc(ch,:));
         end
         
         
         
         
       disp('Select a Rectangular ROI from the image');      
       [a,b]=PickPixel(tag);
       tag_select=tag(a(1):a(2),b(1):b(2));

       life_cell = cell(size(tag_select));
       Amp       = cell(size(tag_select));
       life_imm  = zeros(size(tag_select));     

       disp('Fitting lifetimes for each pixel...')
        for y=1:size(tag_select,1) 
            for x=1:size(tag_select,2)
                if tag_select(y,x,ch)>150  %Olaf: changed from 30 to 150                                    
                ind1 = (im_line == y+a(1)-1);
                ind2 = (im_col == x+b(1)-1);
                tmp1 = double(im_tcspc).*(ind1.*ind2);
                tmp2 = double(im_chan).*(ind1.*ind2);
%                 tcspc_pix(y,x,ch,:)=mHist(double(tmp1(tmp2 == ch)),1:Ngate);
                [A, tau,~, ~, ~, ~,~] = DistFluofit(tcspcIRF, mHist(double(tmp1(tmp2 == ch)),1:Ngate), 1./(head.SyncRate.*1e-9), Resolution,[],1);
                Amp{y,x}=A; life_cell{y,x}=tau; 
                life_imm(y,x,ch)=sum(tau.*A);
                close all
%                 clear A; clear tau; 
                end
            end
            
        end
  end
  INT.tag           = tag;
  INT.tag_select    = tag_select;
  FLIM.life_imm     = life_imm;
  FLIM.tcspc        = tcspc;
  FLIM.tcspcIRF     = tcspcIRF;
  FLIM.t            = t;
  FLIM.Amp          = Amp;
  FLIM.life_cell    = life_cell;
  

  if pic ==1
      figure;
      set (gcf,'name','Intensity Image','NumberTitle','off')
      for ch=1:maxch_n
          subplot(1,maxch_n,ch)
          mim(INT.tag_select(:,:,ch))
      end

      figure;
      set (gcf,'name','Lifetime Image','NumberTitle','off')
      for ch=1:maxch_n
          subplot(1,maxch_n,ch)
          mim(FLIM.life_imm(:,:,ch))
      end
      
      figure
      set(gcs,'name','Lifetime Intensity Plot','NumberTitle','off')
      for ch=1:maxch_n
          subplot(1,maxch_n,ch)
          inten=INT.tag_select(:,:,ch);
          life=FLIM.life_imm(:,:,ch);
          plot(life(:),inten(:),'o','markersize',2.5)
      end
  end
end

  
  
  
  
