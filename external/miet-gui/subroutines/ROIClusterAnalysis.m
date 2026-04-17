function [val_mean,val_stdDev,field_val] = ROIClusterAnalysis(A,nBin,min_cluster,B)
% [val_mean,val_stdDev,field_val]=FinalAnalysisROI(A,nBin,B) accepts an 
% array A of real values, plots them using imagesc and lets the user draw 
% an arbitrary number of elliptical regions of interest. All pixels contained 
% in these ROIs are collected in a single histogram and displayed together 
% with a Gaussian curve with the same mean and standard deviation. The 
% latter two parameters are displayed above the histogram and are returned 
% by the function, the values of the pixels are returned as field_val. 
% The number of bins used for the histogram can be given as the optional 
% parameter nBin (default is nBin=10). Pixels with the value "NaN" are 
% ignored in the analysis.
% Furhtermore, it is possible to give a sencond array B of the same size as
% A which contains intensity information. This array is plotted using imagesc
% and the user is asked to choose a threshold: pixels whose values in B are 
% smaller than the threshold will be ignored in the analysis of A.

if nargin>3 && ~isempty(B) % do we want to reject dim pixels?
    flag=true;
    if size(B)~=size(A)
        errordlg(['The two arrays A and B (1st and 3rd input parameters)'...
            'have to have the same size.']);
    end
else
    flag=false;
end
if nargin<3 || isempty(min_cluster) % set minimum number of pixels for cluter identification
    min_cluster=10;
end

if nargin<2 || isempty(nBin) % set default for number of bins
    nBin=10;
end

if size(A,1)>1 && size(A,2)>1 % open image, let user select ROI, collect values from pixels
    run='y';
    c=0; % number of selected regions
    pixels=false(size(A)); % which pixels are in the ROIs?
    figure;
    while strcmpi(run,'y') % get as many ROIs as the user wants
        c=c+1;
        imagesc(A);
        hold on
        if exist('ROI','var') % show existing ROIs
           for i=1:c-1;
              plot(ROI{i}(1,:),ROI{i}(2,:),'w.','MarkerSize',5) 
           end
        end
        roi=[]; field=[];
        % let the user draw an elliptical ROI
        disp('Select a region, you can drag across and resize. To complete selection, double click mouse left.')
        h= imellipse;
        wait(h);
        % get the position of the ROI
        pos(c,:)=getPosition(h);
        clear h;
        pos_center(c,:)=[pos(c,1)+pos(c,3)/2, pos(c,2)+pos(c,4)/2]; % position of the center of the ellipse
        theta=0:pi/200:2*pi;
        % x- and y-coordinates of the pixels forming the border of the ROI
        roi=[pos_center(c,1)+(pos(c,3)/2.*cos(theta)); pos_center(c,2)+(pos(c,4)/2.*sin(theta))];
        roi=round(roi);
        plot(roi(1,:),roi(2,:),'r.','MarkerSize',5)
        hold off
        % in "field", save the x-coordinates that are covered by the ROI, and
        % for each x-coordinate the lowest and highest y-coordinate of the ROI
        xc=unique(roi(1,:));
        for i=1:numel(xc) % for each x-coordinate: lowest and highest y-coordinate
           field(i,:)=[xc(i),min(roi(2,roi(1,:)==xc(i))), max(roi(2,roi(1,:)==xc(i)))];  
        end
        % crop field to values that actually contain a pixel
        indd=field(:,3)>size(A,1);
        field(indd,3)=size(A,1);
        indu=field(:,2)<1;
        field(indu,2)=1;
        indl=field(:,1)<1;
        field(indl,:)=[];
        indr=field(:,1)>size(A,2);
        field(indr,:)=[];
        % how many pixels are contained in the ROI?
        for i=1:size(field,1) 
            pixels(field(i,2):field(i,3),field(i,1))=true;
        end
        % write ROI into cell arrays to use later
        ROI{c}=roi;
        clear roi field indd indu indl indr
        prompt='Press (y) to select another region of interest or (n) to start analysis...';
        run=input(prompt,'s');           
    end

    % make final nice image showing where ROIs are selected
    imagesc(A); hold on
    for i=1:c;
      plot(ROI{i}(1,:),ROI{i}(2,:),'w.','MarkerSize',10) 
    end
    % if a second matrix has been given: show it and ask for threshold
    if flag
        [pixels, PixelList] = Grain_Image(B.*pixels,min_cluster,0.9);
    end
    % write values of selected pixels into field_val
    field_val = A(pixels);
    figure; mim(A.*pixels);
    figure; scatter(B(pixels),field_val,3,'sk');
    xlabel('intensity [counts]'); ylabel('lifetime [ns]');
else % take already collected values and process them
    field_val = A;
end
field_val(isnan(field_val))=[]; % delete all NaN-entries so they are ignored in the evaluation
pixelInd = numel(field_val); % number of pixels actually used in evaluation
% histogram the values, determine mean and standard deviation, plot Gaussian curve
[N,X]=hist(field_val,nBin); figure; bar(X,N,1,'w'); hold on; x=linspace(min(field_val),max(field_val),100);
val_mean = mean(field_val);
val_stdDev = sqrt( sum((field_val-val_mean).^2)/(pixelInd-1) );
val_curve = exp(-(x-val_mean).^2/(2*val_stdDev^2))/sqrt(2*pi*val_stdDev^2)*numel(field_val)*(max(field_val)-min(field_val))/nBin;
plot(x,val_curve,'-r'); hold off;
title(sprintf('mean=%.3f, standard deviation=%.3f',val_mean,val_stdDev));
