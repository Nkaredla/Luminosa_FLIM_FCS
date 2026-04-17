function handle = cim(x,p1,p2,p3,clmp)

if nargin==1
    nd = ndims(x);
    switch nd
        case 2
            handle = imagesc(x);
            axis image
            if nargin>4
                colormap(clmp)
            else
                colormap hot
            end
            axis off
            %set(gca,'Position',[0 0 1 1])
        case 3
            nj = size(x,3);
            j = 1;
            tmp = (x(:,:,j)-min(min(x(:,:,j))))/(max(max(x(:,:,j)))-min(min(x(:,:,j))));
            for j=2:nj
                tmp = [tmp, (x(:,:,j)-min(min(x(:,:,j))))/(max(max(x(:,:,j)))-min(min(x(:,:,j))))];
            end
            mim(tmp)
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
            mim(tmp)
    end
end
if nargin==2
    if numel(p1)==2
        handle = imagesc(x,p1);
        axis image
        if nargin>4
            colormap(clmp)
        else
            colormap hot
        end
        axis off
    elseif p1=='h'
        mim(x);
        colorbar('h');
    elseif p1=='v'
        mim(x);
        colorbar;
    elseif ischar(p1)
        mim(x);
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
            x = 1 + (x-mm)/tmp*63;
        else
            x = 64*ones(size(x));
        end
        x(isnan(x)) = 1;
        x(x<1) = 1;
        x(x>64) = 64;
        if nargin>4
            col = clmp;
        else
            col = double(jet);  % rgb-values of colormap with 64 different colors
        end
        if length(col)~=64
            for j =1:size(col,2)
                tmpc(:,j) = interp1(1:length(col),col(:,j),1:length(col)/65:length(col));
            end
            col = tmpc;
        end
        %             col = double(jet);
        col = [[1 1 1]; col];
        
        %         jet_black=[zeros(1,3); jet];
        %         col=double(jet_black);
        for j=1:3
            c(:,:,j) = ((floor(x)+1-x).*reshape(col(floor(x),j),size(x)) + (x-floor(x)).*reshape(col(ceil(x),j),size(x))).*p1;
        end
        if tmp>0
            clf
            if size(x,1)>=size(x,2)
                handle = imagesc(c);
                axis image
                axis off
                h = colorbar;
                tmp = mm + get(h,'ytick')*tmp;
                tmp = round(tmp/10.^(floor(log10(tmp(end)))-2))*10.^(floor(log10(tmp(end)))-2);
                set(h, 'yticklabel', num2str(tmp'))
            else
                handle = imagesc(c);
                axis image
                axis off
                h = colorbar('h');
                tmp = mm + get(h,'xtick')*tmp;
                tmp = round(tmp/10.^(floor(log10(tmp(end)))-2))*10.^(floor(log10(tmp(end)))-2);
                set(h, 'xticklabel', num2str(tmp'))
            end
            shading flat
            axis image
            colormap(jet)
        else
            handle = imagesc(c);
            axis image
            axis off
        end
    end
end
if nargin==3
    clmp = p2;
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
            x = 1 + (x-mm)/tmp*63;
        else
            x = 64*ones(size(x));
        end
        x(isnan(x)) = 1;
        x(x<1) = 1;
        x(x>64) = 64;
        col = clmp;
       
        if length(col)~=64
            for j =1:size(col,2)
                tmpc(:,j) = interp1(1:length(col),col(:,j),1:length(col)/65:length(col));
            end
            col = tmpc;
        end
        %             col = double(jet);
        col = [[1 1 1]; col];
        
        %         jet_black=[zeros(1,3); jet];
        %         col=double(jet_black);
        for j=1:3
            c(:,:,j) = ((floor(x)+1-x).*reshape(col(floor(x),j),size(x)) + (x-floor(x)).*reshape(col(ceil(x),j),size(x))).*p1;
        end
        if tmp>0
            clf
            if size(x,1)>=size(x,2)
                handle = imagesc(c);
                axis image
                axis off
                h = colorbar;
                tmp = mm + get(h,'ytick')*tmp;
                tmp = round(tmp/10.^(floor(log10(tmp(end)))-2))*10.^(floor(log10(tmp(end)))-2);
                set(h, 'yticklabel', num2str(tmp'))
            else
                handle = imagesc(c);
                axis image
                axis off
                h = colorbar('h');
                tmp = mm + get(h,'xtick')*tmp;
                tmp = round(tmp/10.^(floor(log10(tmp(end)))-2))*10.^(floor(log10(tmp(end)))-2);
                set(h, 'xticklabel', num2str(tmp'))
            end
            shading flat
            axis image
            colormap(jet)
        else
            handle = imagesc(c);
            axis image
            axis off
        end
    
end
if nargin>3
    if size(x,1)==size(p1,1) && size(x,2)==size(p1,2) % x=array that determines color, p1=array that determines brightness
        if numel(p2)==2 % range of color map
            c = zeros(size(x,1),size(x,2),3);
            mm = min(p1(isfinite(p1)));
            tmp = max(p1(isfinite(p1)))-mm;
            if tmp>0
                p1 = (p1-mm)/tmp; % "brightness" now has values from 0 to 1
            else
                p1 = zeros(size(p1));
            end
            mm = p2(1);
            tmp = p2(2) - mm;
            if tmp>0
                x = 1 + (x-mm)/tmp*63; % if p2(1)==min(x) and p2(2)==max(x), x now has values from 1 to 64, otherwise ? to ?
            else
                x = 64*ones(size(x));
            end
            x(isnan(x)) = 1;
            x(x<1) = 1;         % cut x to range specified by p2
            x(x>64) = 64;       % cut x to range specified by p2
            if nargin>4
                col = clmp;
            else
                col = double(jet);  % rgb-values of colormap with 64 different colors
            end
            %             col = double(jet);
            if length(col)~=64
                for j =1:size(col,2)
                    tmpc(:,j) = interp1(1:length(col),col(:,j),1:length(col)/65:length(col));
                end
                col = tmpc;
            end
            col = [[1 1 1]; col];
            for j=1:3 % interpolate rgb-values for pixels weighted with brightness
                c(:,:,j) = ((floor(x)+1-x).*reshape(col(floor(x),j),size(x)) + (x-floor(x)).*reshape(col(ceil(x),j),size(x))).*p1;
            end
            if tmp>0
                clf   % clear current figure
                if size(x,1)>size(x,2) || nargin>3 && ~isempty(p3) && p3=='v'
                    subplot('position',[0.1 0.1 0.7 0.8]);
                    handle = imagesc(c);
                    axis image
                    axis off
                    h= subplot('position',[0.8 0.1 0.05 0.8]);
                    pcolor([0 0.05]*tmp,mm+(1:size(x,1))/size(x,1)*tmp,(mm+(1:size(x,1))'/size(x,1)*tmp)*ones(1,2))
                    set(gca,'xtick',[],'yaxislocation','right')
                    set(h,'linewidth',2)
                else
                    subplot('position',[0.1 0.3 0.8 0.6]);
                    handle = imagesc(c);
                    axis image
                    axis off
                    h = subplot('position',[0.1 0.1 0.8 0.05]);
                    pcolor(mm+(1:size(x,1))/size(x,1)*tmp,[0 0.05]*tmp,ones(2,1)*(mm+(1:size(x,1))/size(x,1)*tmp))
                    set(gca,'ytick',[],'xaxislocation','bottom')
                     set(h,'linewidth',2)
                end
                shading flat
                axis image
                if nargin>4
                    colormap(clmp)
                else
                    colormap(jet)
                end
                
            else
                handle = imagesc(c);
                axis image
                axis off
            end
        end
    elseif numel(p1)==2
        mim(x,p1);
        if p2=='h'
            colorbar('h');
        elseif p2=='v'
            colorbar;
        end
    elseif numel(p2)==2
        mim(x,p2);
        if p1=='h'
            colorbar('h');
        elseif p1=='v'
            colorbar;
        end
    end
end
figure(gcf);
end

