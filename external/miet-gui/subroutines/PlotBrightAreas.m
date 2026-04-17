function PlotBrightAreas(interestingQuantity,intensity,flag)

if nargin>2 && ~isempty(flag) && flag~=0
    figure; imagesc(intensity); title('intensity'); colorbar; 
end

% escape all "_" in variable name, otherwise character after "_" would be lower index in title
protectedName=protectName(inputname(1)); 

% display the "interesting quantity"
cm = colormap(jet(128));    % colormap: jet with 128 hues
cm(1,:)=[0 0 0];            % black for lowest value (i.e. for NaN)
f=figure; imagesc(interestingQuantity); colorbar; colormap(cm);
title(sprintf('%s, adjust intensity threshold for ignoring dim pixels',protectedName));

% insert a text field where the user can enter the threshold
uicontrol('Style','edit','String','1','Callback',@ThresholdField);
set(f,'ToolBar','figure');


function ThresholdField(ObjectHandle,~) % callback function: do not show dim pixels in height image
    threshold = str2double(get(ObjectHandle,'String'));
    tmp = interestingQuantity; tmp(intensity<threshold)=NaN;
    imagesc(tmp); colorbar; set(f,'ToolBar','figure');
    title(sprintf('%s',protectedName));
end

function protectedName = protectName(name) % automatically escape all "_"
    C = strsplit(name,'_');
    protectedName = C{1};
    for index=2:numel(C)
       protectedName=[protectedName '\_' C{index}]; %#ok<AGROW>
    end
end

end

