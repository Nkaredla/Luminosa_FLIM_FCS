
folderName = 'D:\Luminosa\Data\260319\';
filenames = dir(folderName);
for k = 5:numel(filenames)
    name = [folderName, filenames(k).name,filesep,'RawImage.ptu'];
    head = PTU_Read_Head(name);
    out = PTU_MultiFrameScanReadFast(name); % read the channels
    int = out.tags;
    tauc = out.taus;
    % intensity clipping
    lb = prctile(int(:),60); % lower bound, 30 for cells, 10 for tissue
    ub = prctile(int(:),99);
    int(int<lb) = lb; 
    int(int>ub) = ub;
    int = int-lb; %
    int = int./max(int(:)); % normalized intensity
    
    cmap = jet(256);
    cmap = cmap(30:end-30,:);
    
    cim(tauc ,log10(0.001 + int.^2),[3.5 6.5],'v',cmap)

    print( [name(1:end-4),'_FLIM_image.png'],'-dpng','-r300')
            
%     cim(tauc, int.^0.5,[3 6],'v',cmap) % gamma correction
end
