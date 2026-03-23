
folderName = 'D:\Luminosa\Data\260319\';
filenames = dir(folderName);
for k = 5:numel(filenames)
    name = [folderName, filenames(k).name,filesep,'RawImage.ptu'];
    head = PTU_Read_Head(name);
    out = PTU_MultiFrameScanReadFast(name); % read the channels
    int = out.tags;
    tauc = out.taus;
    % intensity clipping
    lb = prctile(int(:),30); % lower bound, 30 for cells, 10 for tissue
    ub = prctile(int(:),99.5);
    int(int<lb) = lb; 
    int(int>ub) = ub;
    int = int-lb; %
    int = int./max(int(:)); % normalized intensity
    
    cmap = jet(256);
    cmap = cmap(30:end-30,:);
    
    hIm = cim(tauc ,log10(0.001 + int.^2),[3.5 6.5],'v',cmap);
    addPTUScaleBar(hIm, head, size(tauc), 'bottomright', ...
    'BarLengthUm', 5, ...
    'Label', '5 \mum', ...
    'Color', 'w', ...
    'LineWidth', 3, ...
    'FontSize', 20);

    save([name(1:end-4),'_FLIM.mat'], 'out','head');

    print( [name(1:end-4),'_FLIM_image.png'],'-dpng','-r300')
            
%     cim(tauc, int.^0.5,[3 6],'v',cmap) % gamma correction
end
