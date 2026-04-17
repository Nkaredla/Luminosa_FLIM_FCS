function handle = mim(x,p1,p2,p3,clmp)
% handle=mim(...) offers several different possibilities for visualizing x:
%
% mim(x)
% - x 2-dimensional: x is plotted as a scaled image with colormap 'hot'
% - x 3- or 4-dimensional: the same but individually for each x-y-plane 
%   (i.e. they have different color ranges, they are scaled differently!); 
%   these images are plotted next to each other with the 3rd dimension from
%   left to right and the 4th dimension from top to bottom
%
% mim(x,p1)
% - x 2-dimensional, p1 has two entries: entries of p1 are used as lower
%   and upper bound of the color range with which x is plotted; colormap 'hot'
% - x 2-, 3-, 4-dimensional, p1='h': same as mim(x) but with a horizontal colorbar
% - x 2-, 3-, 4-dimensional, p1='v': same as mim(x) but with a vertical colorbar
% - x 2-, 3-, 4-dimensional, p1 a string: same as mim(x) but the text in p1 is
%   written in white on the image
% - x 2-dimensional, p1 has the same size as x: x is visualized "weighted with p1",
%   i.e. the entries in x encode hue, the entries in p1 encode brightness (this 
%   does not work if x is 3- or 4-dimensional); colormap 'jet' with black as lowest
%
% mim(x,p1,p2,p3)
% - x and p1 2-dimensional, same size, p2 has two entries: same as mim(x,p1) but 
%   the entries of p2 are used as lower and upper bound of the color range (values 
%   that are too small or too large will be set to p2(1) or p2(2), respectively),
%   colormap 'jet' with black as the lowest entry; brightness scales from 10-100%;
%   if p3='v', a vertical colorbar is drawn next to the image, if p3 is  empty or 
%   set to any other value, a horizontal colorbar is drawn underneath the image
% - x 2-dimensional, p1 has two entries, p2='h' or 'v': same as mim(x,p1) but 
%   with a horizontal or vertical colorbar, respectively; colormap 'hot'
% - x 2-dimensional, p1='h' or 'v', p2 has two entries: same as mim(x,p2) but 
%   with a horizontal or vertical colorbar, respectively; colormap 'hot'
%
% mim(x,p1,p2,p3,clmp)
% - x and p1 2-dimensional, same size, p2 has two entries, clmp has 3 columns: 
%   same as mim(x,p1,p2,p3) but the colormap is given by clmp
%
% The output 'handle' contains the handle of the axis where the image was drawn
% (important for example if you want to add a title).

if nargin==1
    nd = ndims(x);
    switch nd
        case 2
            imagesc(x);
            handle = gca;
            axis image
            colormap hot
            axis off
        case 3
            nj = size(x,3);
            j = 1; 
            tmp = (x(:,:,j)-min(min(x(:,:,j))))/(max(max(x(:,:,j)))-min(min(x(:,:,j))));
            for j=2:nj
                tmp = [tmp, (x(:,:,j)-min(min(x(:,:,j))))/(max(max(x(:,:,j)))-min(min(x(:,:,j))))];
            end
            handle = mim(tmp);
        case 4
            nj = size(x,3);
            nk = size(x,4);
            j = 1; k = 1;
            tmp = (x(:,:,j,k)-min(min(x(:,:,j,k))))/(max(max(x(:,:,j,k)))-min(min(x(:,:,j,k))));
            for j=2:nj
                tmp = [tmp, (x(:,:,j,k)-min(min(x(:,:,j,k))))/(max(max(x(:,:,j,k)))-min(min(x(:,:,j,k))))];
            end
            for k=2:nk
                j = 1;
                tmptmp = (x(:,:,j,k)-min(min(x(:,:,j,k))))/(max(max(x(:,:,j,k)))-min(min(x(:,:,j,k))));
                for j=2:nj
                    tmptmp = [tmptmp, (x(:,:,j,k)-min(min(x(:,:,j,k))))/(max(max(x(:,:,j,k)))-min(min(x(:,:,j,k))))];
                end
                tmp = [tmp; tmptmp];
            end
            handle = mim(tmp);
    end
end
if nargin==2
    if numel(p1)==2
        imagesc(x,p1);
        handle = gca;
        axis image
        colormap hot
        axis off
    elseif p1=='h'
        handle = mim(x);
        colorbar('h');
    elseif p1=='v'
        handle = mim(x);
        colorbar;
    elseif ischar(p1)
        handle = mim(x); 
        [a,b]=size(x);
        text(b*(1-0.025*length(p1)),0.06*a,p1,'FontName','Times','FontSize',16,'Color','w');
    else
        c = zeros(size(x,1),size(x,2),3);
        mm = min(p1(isfinite(p1)));
        tmp = max(p1(isfinite(p1)))-mm;
        if tmp>0
            p1 = (p1-mm)/tmp;
        else
            p1 = zeros(size(p1));
        end
        mm = min(x(isfinite(p1)));
        tmp = max(x(isfinite(p1)))-mm;
        if tmp>0
            x = 1 + (x-mm)/tmp*127;
        else
            x = 128*ones(size(x));
        end
        x(isnan(x)) = 1;
        x(x<1) = 1;
        x(x>128) = 128;
        jet_black=[zeros(1,3); jet(127)]; % colormap 'jet' with black as lowest color
        col=double(jet_black);
        for j=1:3
            c(:,:,j) = ((floor(x)+1-x).*reshape(col(floor(x),j),size(x)) + (x-floor(x)).*reshape(col(ceil(x),j),size(x))).*p1;
        end
        if tmp>0
            clf
            if size(x,1)>=size(x,2)
                imagesc(c);
                handle = gca;
                axis image
                axis off
                h = colorbar;
                tmp = mm + get(h,'ytick')*tmp;
                tmp = round(tmp/10.^(floor(log10(tmp(end)))-2))*10.^(floor(log10(tmp(end)))-2);
                set(h, 'yticklabel', num2str(tmp'))
            else
                imagesc(c);
                handle = gca;
                axis image
                axis off
                h = colorbar('h');
                tmp = mm + get(h,'xtick')*tmp;
                tmp = round(tmp/10.^(floor(log10(tmp(end)))-2))*10.^(floor(log10(tmp(end)))-2);
                set(h, 'xticklabel', num2str(tmp'))
            end
            shading flat
            axis image
            colormap(jet_black)            
        else
            imagesc(c);
            handle = gca;
            axis image
            axis off
        end
    end
end
if nargin>2
    if size(x,1)==size(p1,1) && size(x,2)==size(p1,2) % x=array that determines color, p1=array that determines brightness
        if numel(p2)==2 % range of color map
            c = zeros(size(x,1),size(x,2),3);
            mm = min(p1(isfinite(p1)));
            tmp = max(p1(isfinite(p1)))-mm;
            if tmp>0
                p1 = 0.1+0.9*(p1-mm)/tmp; % "brightness" now has values from 0.1 to 1
            else
                p1 = zeros(size(p1));
            end
            mm = p2(1);
            tmp = p2(2) - mm;
            if tmp>0
                x = 1 + (x-mm)/tmp*127; % if p2(1)==min(x) and p2(2)==max(x), x now has values from 1 to 128, otherwise ? to ?
            else
                x = 128*ones(size(x));
            end
            x(isnan(x)) = 1;
            x(x<1) = 1;         % cut x to range specified by p2
            x(x>128) = 128;     % cut x to range specified by p2
            if nargin>4
                col = clmp;
            else
                col = jet(127);         % rgb-values of colormap 'jet' with 127 different colors
                col = [zeros(1,3); col];% black as lowest color 
            end
            if length(col)~=128
                for j =1:size(col,2)
                    tmpc(:,j) = interp1(1:length(col),col(:,j),1:length(col)/129:length(col));
                end
                col = tmpc;
            end
            for j=1:3 % interpolate rgb-values for pixels weighted with brightness
                c(:,:,j) = ((floor(x)+1-x).*reshape(col(floor(x),j),size(x)) + (x-floor(x)).*reshape(col(ceil(x),j),size(x))).*p1;
            end
            if tmp>0
                clf   % clear current figure
                if size(x,1)>size(x,2) || (nargin>3 && ~isempty(p3) && p3=='v')
                    handle = subplot('position',[0.1 0.1 0.7 0.8]);
                    imagesc(c);
                    axis image
                    axis off
                    subplot('position',[0.8 0.1 0.1 0.8]);
                    pcolor([0 0.1]*tmp,mm+(1:size(x,1))/size(x,1)*tmp,(mm+(1:size(x,1))'/size(x,1)*tmp)*ones(1,2))
                    set(gca,'xtick',[],'yaxislocation','right')                    
                else
                    handle = subplot('position',[0.1 0.3 0.8 0.6]);
                    imagesc(c);
                    axis image
                    axis off
                    subplot('position',[0.1 0.1 0.8 0.1]);
                    pcolor(mm+(1:size(x,1))/size(x,1)*tmp,[0 0.1]*tmp,ones(2,1)*(mm+(1:size(x,1))/size(x,1)*tmp))
                    set(gca,'ytick',[],'xaxislocation','bottom')                    
                end
                shading flat
                axis image
                colormap(col);
            else
                imagesc(c);
                handle = gca;
                axis image
                axis off
                colormap(col);
            end
        end
    elseif numel(p1)==2
        handle = mim(x,p1);
        if p2=='h'
            colorbar('h');
        elseif p2=='v'
            colorbar;
        end
    elseif numel(p2)==2
        handle = mim(x,p2);
        if p1=='h'
            colorbar('h');
        elseif p1=='v'
            colorbar;
        end
    end
end
figure(gcf);
