function [tau_fit,tau_avg,tau_px_max,tau_px_mean] = FreeSpaceLifetime(name)
% FreeSpaceLifetime(name) finds the free space lifetime of a dye by
% evaluating a sample consisting of this dye in the same medium that will
% be later used in the MIET experiment.
% 'name' is the name of the ht3-file containing the data
% 'tau_fit' fitted lifetime (TCSPC curve of whole image)
% 'tau_avg' average lifetime (TCSPC curve of whole image)
% 'tau_px_max' most common lifetime in lifetime image
% 'tau_px_mean' average lifetime in lifetime image (averaged over photons 
%     of all pixels with >25 photons)

disp('Get free space lifetime from tcspc histogram of the whole image...');
% read photon counts of bins that contain at least 1 photon
if strcmp(name(end-2:end),'ht3') == 0
    errordlg('You have to specify an ht3-file.');
else
    [bin_original,tcspcdata_original,head] = ReadTCSPC(name);
    numChan = size(tcspcdata_original,2); % how many channels, i.e. how many detectors?
end
if isempty(head)
    errordlg('Could not read ht3-data.');
else
    area = head.ImgHdrSize~=0;  % a flag to see if the given ht3 file is an area file or a point measurement
end

% make tcspc-histogram (number of photons vs. bin number) by inserting
% rows for all bins that did not contain any photons and writing a zero
% as the corresponding photon count; either for IRF or for actual data
tcspcdata = zeros(max(bin_original),size(tcspcdata_original,2));
for i = 1:size(tcspcdata_original,2)
    tcspcdata(bin_original,i) = tcspcdata_original(:,i);
end
bin=(1:max(bin_original)).';
clear bin_original tcspcdata_original

% find time-gates, shift data accordingly (either IRF or actual data)
timegate = DetectTimeGates(tcspcdata, 1, head.Resolution);
for chanInd = 1:numChan
    if timegate(chanInd,3) == 0
        tcspcdata(:,chanInd) = tcspcdata(timegate(chanInd,1):timegate(chanInd,2),chanInd);
    else
        tcspcdata(:,chanInd) = tcspcdata([timegate(chanInd,1):timegate(chanInd,2) timegate(chanInd,3):timegate(chanInd,4)],chanInd);
    end
end
tcspcdata = sum(tcspcdata,2);  % combining all detector channels
[~, chan_maxcnt] = max(tcspcdata);
cutoff = chan_maxcnt + round(0.5/head.Resolution); % 0.5 ns shift.

% 1st method: tail fitting using a bi-exponential model
p = [1,3];
figure;
p = Simplex('ExpFun',p,[],[],[],[],(bin(cutoff:end)-bin(cutoff)+1)*head.Resolution,tcspcdata(cutoff:end),1,[],1) ;
[~,c,~,z]= ExpFun(p,(bin(cutoff:end)-bin(cutoff)+1)*head.Resolution,tcspcdata(cutoff:end),1,[],[]);
tau_fit = sum(p.*c(2:3))/sum(c(2:3)); % weighted average of the two exponentials

% 2nd method: average arrival time of the photons after the cutoff
[~,minchan]= BG_val(tcspcdata, 20);
tau_avg = sum((tcspcdata(cutoff:minchan)-c(1)).*(bin(cutoff:minchan)-bin(cutoff)+1)*head.Resolution)/sum(tcspcdata(cutoff:minchan)-c(1));

% plot results
close(gcf); figure; subplot(4,1,1:3)
plot(bin*head.Resolution, tcspcdata, 'ob', bin(cutoff:end)*head.Resolution, z, 'r');
axis tight; legend('data','fit'); title(sprintf('tau\\_fit = %.1f ns, tau\\_avg = %.1f ns',tau_fit,tau_avg));
subplot(4,1,4)
plot(bin(cutoff:end)*head.Resolution, (tcspcdata(cutoff:end)-z)./sqrt(abs(z)),bin*head.Resolution,0*bin)
axis tight

% 3rd method: bin lifetimes of single pixels, find maximum and mean of histogram
disp('Get free space lifetime from tcspc histograms of single pixels...');
if area % only possible if image has more than one pixel
    [data, tag, life_imm]=Process_scan(name, 1) ; %#ok<ASGLU>
    save([name(1:end-4),'_PS.mat'],'data','tag','life_imm');
    tag(tag<26)=0;
    %         nx = head.ImgHdr.PixX;
    %         ny = head.ImgHdr.PixY;
    if size(life_imm,4)>1
        for z = 1:size(life_imm,4)
            tmplife_imm = life_imm(:,:,:,z);
            tmptag      = tag(:,:,:,z);
            for x=1:size(tmplife_imm,1)
                for y=1:size(tmplife_imm,2)
                    tmplife_imm2(x,y,z)=sum(tmptag(x,y,:).*tmplife_imm(x,y,:))/sum(tmptag(x,y,:)); %#ok<*AGROW>
                end
            end
        end
        tmptag2 = sum(tag,3);
        for x=1:size(tmplife_imm2,1)
            for y=1:size(tmplife_imm2,2)
                life_imm(x,y,1)=sum(tmptag2(x,y,:).*tmplife_imm2(x,y,:))/sum(tmptag2(x,y,:));
            end
        end
        life_imm = life_imm(:,:,1); tag = sum(tmptag2,3);
    elseif size(life_imm,3)>1 % combine different detector channels
        for x=1:size(life_imm,1)
            for y=1:size(life_imm,2)
                life_imm(x,y,1)=sum(tag(x,y,:).*life_imm(x,y,:))/sum(tag(x,y,:));
            end
        end
        life_imm = life_imm(:,:,1); tag = sum(tag,3);
    end
    % sort pixel lifetimes into histogram with 30 bins:
    % 0                       <= tau < 1/29*max(max(life_imm))
    % 1/29*max(max(life_imm)) <= tau < 2/29*max(max(life_imm))
    %                       ...
    % 28/29*max(max(life_imm))<= tau < max(max(life_imm))
    %                      tau = max(max(life_imm))
    numberOfBins=30;
    histoRanges=min(min(life_imm))+(0:numberOfBins-1)/(numberOfBins-1)*(max(max(life_imm))-min(min(life_imm)));
    total_intensity=sum(tag(~isnan(life_imm))); % only use pixels where all channels were valid
    [~,binnumber]=histc(reshape(life_imm,1,size(life_imm,1)*size(life_imm,2)),histoRanges);
    histogram=zeros(numberOfBins,1);
    for i=1:numel(binnumber) % histogram: number of photons that arrived with a certain lifetime
        if binnumber(i)~=0
            histogram(binnumber(i))=histogram(binnumber(i))+tag(i)/total_intensity;
        end
    end
    histogram = histogram(1:end-1)./sum(histogram); % the last point contains only the maximum value (no range) -> throw away
    histoBinMidpoints = (histoRanges(1:end-1)+histoRanges(2:end))/2;
    [~,position_of_max] = max(histogram); % estimate for mean value
    %         stdDev = sqrt(sum((histoBinMidpoints-histoBinMidpoints(position_of_max)).^2.*histogram.')); % estimate for standard deviation
    %         fitParam = Simplex('Gauss',[histoBinMidpoints(position_of_max) stdDev],...
    %             [0 0],[10 10],[],[],histoBinMidpoints,histogram,[],[],1); % fit to find mean value and standard deviation
    %         [~,~,fitPoints,~] = Gauss(fitParam,histoBinMidpoints,histogram,[],[],1); % get fitted curve (for plotting)
    %         tau_pixels = fitParam(1);
    tau_px_max = histoBinMidpoints(position_of_max);
    tau_px_mean = histoBinMidpoints*histogram;
    figure; bar(histoBinMidpoints,histogram); %hold on; % plot
    %         plot(histoBinMidpoints,fitPoints,'ob'); hold off;
    title(sprintf('tau\\_px\\_max = %.1f ns, tau\\_px\\_mean = %.1f ns',tau_px_max,tau_px_mean));
    xlabel('time [ns]'); ylabel('relative frequency [-]');
else
    tau_px_max = NaN; tau_px_mean=NaN;
end
end