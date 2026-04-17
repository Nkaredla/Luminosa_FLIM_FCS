function MIET_viewer(name, addr, format)
% MIET_viewer is for post viewing the analysed smMIET data. It displays all
% the figures and saves them to the specified address (addr) in the format
% assigned (format)

% (c) Narain Karedla, 2014

% name= 'U:\Anna\27-01-2014-Microtubles\smFLIM\Image_021.ht3';

if nargin<3 || isempty(format)
    format='-dpng';
end

if nargin<2 || isempty(addr)
    [addr,~,~]=fileparts(name);
elseif ~exist(addr,'file')
    mkdir(addr)
end


if strcmp(name(end-3:end),'.ht3')
    dname=name(1:end-4);
    [~,fname,~]=fileparts(name);
else
    disp('Invalid File!')
    return
end

if exist([dname,'_sm_flim.mat'],'file')||exist([dname,'_smflim_analysis.mat'],'file')
    %    load([dname,'_Core_Scan.mat'],'head')
    
    if exist([dname,'_sm_flim.mat'],'file')
        load([dname,'_intensity_image.mat'],'tag')
        load([dname,'_sm_flim.mat'],'head','SM','INT','FLIM');
    else
        load([dname,'_smflim_analysis.mat'],'head','INT','FLIM');
    end
    
    nx         = head.ImgHdr.PixX;
    ny         = head.ImgHdr.PixY;
    x0         = head.ImgHdr.X0;
    y0         = head.ImgHdr.Y0;
    pixel      = head.ImgHdr.PixelSize;
    maxch_n=size(INT.tag,3);
    
    close all
    figure
    set (gcf,'name','Intensity Image','NumberTitle','off')
    for ch=1:maxch_n
        subplot(maxch_n,1,ch)
        imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), INT.tag(:,:,ch))
        if maxch_n>1
            title(['Detector ',num2str(ch)],'FontSize',16)
        end
        set(gca,'FontSize',16)
        xlabel(['\mu','m'],'FontSize',20)
        ylabel(['\mu','m'],'FontSize',20)
        h = colorbar;
        ylabel(h,'Photon Counts','FontSize',16)
        set(gca,'FontSize',16)
        axis equal
        axis tight
        colormap('hot')
    end
    print(format,'-r300',[addr,'\',fname ,'_Intensity'])
    
    
    if isfield(INT,'imm_c')
        figure
        set (gcf,'name','Fitted Patterns','NumberTitle','off')
        for ch=1:maxch_n
            subplot(maxch_n,1,ch)
            if iscell(INT.imm_c)
            imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), INT.imm_c{ch})
            else
                imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), INT.imm_c(:,:,ch))
            end
            if maxch_n>1
                title(['Detector ',num2str(ch)],'FontSize',16)
            end
            set(gca,'FontSize',16)
            xlabel(['\mu','m'],'FontSize',20)
            ylabel(['\mu','m'],'FontSize',20)
            axis equal
            axis tight
            colormap('hot')
        end
        print(format,'-r300',[addr,'\',fname ,'_Patterns'])
    end
    
    if isfield(FLIM,'lifeimm')
        figure
        set (gcf,'name','Lifetime Image','NumberTitle','off')
        for ch=1:maxch_n
            subplot(maxch_n,1,ch)
             if iscell(FLIM.lifeimm)
            imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), FLIM.lifeimm{ch})
            else
                imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), FLIM.lifeimm(:,:,ch))
            end
            if maxch_n>1
                title(['Detector ',num2str(ch)],'FontSize',16)
            end
            set(gca,'FontSize',16)
            xlabel(['\mu','m'],'FontSize',20)
            ylabel(['\mu','m'],'FontSize',20)
            h = colorbar;
            ylabel(h,'Lifetime (ns)','FontSize',16)
            set(gca,'FontSize',16)
            axis equal
            axis tight
            colormap('hot')
        end
        print(format,'-r300',[addr,'\',fname ,'_Lifetime'])
        
        figure;
        set (gcf,'name','Photon Count vs Lifetimes','NumberTitle','off')
        for ch=1:maxch_n
            subplot(1,maxch_n,ch)
             if iscell(FLIM.photon_sm)
             plot(FLIM.life_matav{ch},FLIM.photon_sm{ch},'o','MarkerSize',5)
            else
                plot(FLIM.life_matav(:,ch),FLIM.photon_sm(:,ch),'o','MarkerSize',5)
            end
           if maxch_n>1
                title(['Detector ',num2str(ch)],'FontSize',12)
            end
            xlabel('Lifetimes(ns)','FontSize',12)
            ylabel('Photons per molecule','FontSize',12)
        end
        print(format,'-r300',[addr,'\',fname ,'_Photons_Lifetime'])
        if ~isempty(SM)
             % then it is radially polarized patterns
                figure;
                set(gcf,'name','Polar Orientation vs Lifetimes','NumberTitle','off')
                for ch=1:maxch_n
                    subplot(1,maxch_n,ch)
                    if iscell(FLIM.life_matav)
                         plot(SM.theta{ch}.*180./pi, FLIM.life_matav{ch},'o','MarkerSize',5)
                    else
                        plot(SM.theta(ch,:).*180./pi, FLIM.life_matav(:,ch),'o','MarkerSize',5)
                    end
                    if maxch_n>1
                        title(['Detector ',num2str(ch)],'FontSize',12)
                    end
                    xlabel('Inclination (degrees)','FontSize',12)
                    ylabel('Lifetime (ns)','FontSize',12)
                end
           
            print(format,'-r300',[addr,'\',fname ,'_Orientation_Lifetime'])
            
        end
    end
    if exist([dname,'_miet_analysis.mat'],'file')
        load([dname,'_miet_analysis.mat'],'posimm_c')
        figure
        set(gcf,'name','Height Image','NumberTitle','off')
        for ch=1:maxch_n
            subplot(ch,1,maxch_n)
            if maxch_n>1
                title(['Detector ',num2str(ch)],'FontSize',16)
            end
            imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), posimm_c(:,:,ch).')    %#ok<NODEF>
            set(gca,'FontSize',16)
            xlabel(['\mu','m'],'FontSize',20)
            ylabel(['\mu','m'],'FontSize',20)
            h = colorbar;
            ylabel(h,'Axial Height (nm)','FontSize',16)
            set(gca,'FontSize',16)
            axis equal
            axis tight
            colormap('copper');
        end
        print(format,'-r300',[addr,'\',fname ,'_Axial_Distances'])
    end
    display('Processing complete!')
else
    disp('Analysis for this file has not been performed. Please analyze the data!')
end