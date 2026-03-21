function [G, Gcross, Gcarp, GcarpCross, t, xxi] = lsCrossRead(res, flag)

if nargin<2||isempty(flag)
    flag = 25;
end

t = res.autotime;
pixels = max(res.line_idx);
components = ceil(size(res.auto,2)/pixels);

timel = res.time(end);
if timel<mean(res.time)/3
    res.auto = res.auto(:,:,:,1:end-1); % remove the last bunch
end

ind = ~isnan(res.auto(:,1,1,1)); % check the correlations in the end
t = t(ind);
res.auto = res.auto(ind,:,:,:);

% edof=1; %cutting pixels from the linescan beginning and end
% pixelst=pixels-edof; %Number of pixels in one single scan with the edges cut of
xiaxis = pixels*2-1; %199 is 100 pixel possitive and negative minus 1 of xi=0
srep=length(res.auto(:,1,1)); %number of scans

G=zeros(srep,xiaxis,components); %Initialization. Scan number and xi dimensions.
Gcross = zeros(srep,xiaxis,(components-1)*components);

counter = zeros(2*pixels-1,1);
counter_cross = counter;
c = 0;
for k = 1:components
    %     auto_comp = res.auto(:,k:components:pixels*components,k:components:pixels*components,:);
    for i=1:pixels; %pixel1
        for j=1:pixels; %pixel2
            xi=j-i; %difference xi in the scanning direction positive and flow direction negative
            G(:,xi+pixels,k)=G(:,xi+pixels,k)+mean(res.auto(:,(i-1)*components+k,(j-1)*components+k,:),4); %correlation of the xi of j and i combinations
            counter(xi+pixels) = counter(xi+pixels) + 1;
            
        end
    end
    
    for p = k+1:components
        c = c+1;
        for i=1:pixels; %pixel1
            for j=1:pixels; %pixel2
                xi=j-i; %difference xi in the scanning direction positive and flow direction negative
                Gcross(:,xi+pixels,c)=Gcross(:,xi+pixels,c)+mean(res.auto(:,(i-1)*components+k,(j-1)*components+p,:),4); %correlation of the xi of j and i combinations
                counter(xi+pixels) = counter(xi+pixels) + 1;
                counter_cross(xi+pixels) = counter_cross(xi+pixels) + 1;
            end
        end
        
        c =c+1;
        for i=1:pixels; %pixel1
            for j=1:pixels; %pixel2
                xi=j-i; %difference xi in the scanning direction positive and flow direction negative
                Gcross(:,xi+pixels,c)=Gcross(:,xi+pixels,c)+mean(res.auto(:,(i-1)*components+p,(j-1)*components+k,:),4); %correlation of the xi of j and i combinations
                counter(xi+pixels) = counter(xi+pixels) + 1;
                
            end
        end
    end
end

% --- Carpet correlations: same pixel, separated by component ---
nt = numel(t);

Gcarp = zeros(nt, pixels, components);

for k = 1:components
    for i = 1:pixels
        idx = (i-1)*components + k;
        Gcarp(:, i, k) = mean(res.auto(:, idx, idx, :), 4);
    end
end

% --- Optional: same-pixel cross-component carpets ---
ncross = components * (components - 1);
GcarpCross = zeros(nt, pixels, ncross);

c = 0;
for k = 1:components
    for p = k+1:components
        c = c + 1;
        for i = 1:pixels
            idxk = (i-1)*components + k;
            idxp = (i-1)*components + p;
            GcarpCross(:, i, c) = mean(res.auto(:, idxk, idxp, :), 4);
        end

        c = c + 1;
        for i = 1:pixels
            idxk = (i-1)*components + k;
            idxp = (i-1)*components + p;
            GcarpCross(:, i, c) = mean(res.auto(:, idxp, idxk, :), 4);
        end
    end
end


% for j = 1:components+1
%     y(:,j,:) = res.auto(ind,j,j+(components+1),:)+res.auto(ind,j+(components+1),j,:); % D1L1 x D2L1 + D2L1 x D1L1
%     %         y(:,2,j,:) = res.auto(:,j+(components+1),j+cnum*(components+1),:)+res.auto(:,j+cnum*(components+1),j+(components+1),:); % D1L2 x D2L1 + D2L2 x D1L1
%     for k = j+1:components+1
%         c = c+1;
%         ycross(:,c,:) = res.auto(ind,j,k+(components+1),:)+res.auto(ind,j+components+1,k,:);
%         c=c+1;
%         ycross(:,c,:) = res.auto(ind,k+(components+1),j,:)+res.auto(ind,k,j+components+1,:);
%     end
% end



xxi = zeros(pixels,1);
for i = 1:pixels
    ind = res.line_idx==i;
    xxi(i) = round(mean(find(ind)))-1;
end
xxi=[-flipud(xxi(2:end)); xxi];


ind = round(size(G,2)/2+(-flag:flag));
G = G(:,ind,:);
Gcross = Gcross(:,ind,:);
%Still have to do something with that there are the most versions of xi=0
%and as the outher case only one time xi=+and- 99
%
% normvectorend=ones(srep,xiaxis);
% xiplace=1:xiaxis;
% normvector=zeros(xiaxis,1);
% normvector(1:pixels)=1:pixels;
% normvector(101:xiaxis)=(pixels-1):-1:1;
normvectorend=repmat(counter(ind)',[size(G,1),1 components]);
G=G./normvectorend;

normvectorend=repmat(counter_cross(ind)',[size(Gcross,1),1 components*(components-1)]);
Gcross=Gcross./normvectorend;

%plotting G against tau (the number of scans)n for every xi
% figure(1)
%G curve number
% taux=[1:srep]; %x axis
% plot(taux,G(:,xxi+pixels))
xxi = xxi(ind);
