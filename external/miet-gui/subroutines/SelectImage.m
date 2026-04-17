function pix = SelectImage(image)
        if iscell(image)
            if numel(size(image))>2 
               error('Only one dimensional cells are executed!')
            elseif size(image,2)>1
               error('Only one dimensional cells are executed!')
            else
                for i=1:size(image,1)
                    [x(i),y(i)]=size(image{i});
                end
                maxx=max(x); maxy=max(y);
                im=zeros(maxx,maxy,i);
            end
            for i=1:size(image,1)
            im(1:x(i),1:y(i),i) = image{i};
            end
        end

pix=[];
si.size  = size(im,3);
si.frame = 1;
si.skip  = 0;
si.fig   = [];

    function KeyListener(~,event)
        switch event.Key
            case 'leftarrow'
                si.skip = -1;
                tmp = mod(si.frame+si.skip,si.size+1);
                if tmp == 0
                    si.frame = si.size;
                else
                    si.frame = tmp;
                end
                RefreshPlot();
            case 'rightarrow'
                si.skip = 1;
                tmp = mod(si.frame+si.skip,si.size+1);
                if tmp == 0
                    si.frame = 1;
                else
                    si.frame = tmp;
                end
                RefreshPlot();
            case 'return'
                Pickframe();
                close(si.fig)
        end
    end

    function RefreshPlot()
        imagesc(im(:,:,si.frame));
        title(['Frame # ',num2str(si.frame)]);
    end

    function Pickframe
        pix=si.frame;
    end

scrsz = get(0,'ScreenSize');
si.fig = figure('Position',[(scrsz(3)-1680/2)/2,(scrsz(4)-1050/2)/2,1050/2,1050/2]);
colormap hot;
axis equal;

set(si.fig,'KeyPressFcn',@KeyListener);

imagesc(im(:,:,si.frame)); 
waitfor(si.fig);
            
            
end        