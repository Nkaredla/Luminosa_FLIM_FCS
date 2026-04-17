function MIET_disp(height,tag,smooth,sh,mini)

if nargin<2||isempty(tag)
    tag=[];
end
if nargin<5||isempty(mini)
   if isempty(tag)
        tag = height;
        tmp  = sort(tag(:),'descend');
        num  = ceil(0.6*numel(tmp));
        mini = tmp(num);
    else
        tmp  = sort(tag(:),'descend');
        num  = ceil(0.4*numel(tmp));
        mini = tmp(num);
    end
end
if nargin<4|| isempty(sh)
    sh=27/64; % 27/64 will give nice colour
end
if nargin<3|| isempty(smooth)
    smooth =2; % smoothing on height profile
end

maxh = max(height(:));
minh = min(height(:));

height(isnan(height))=maxh;
tmpim = mConv2(height,Disk(smooth));
tmpim(tmpim<=minh)=NaN;
thres_index = tag<mini;
tmpim(thres_index)=NaN;
shadow=ones(size(height))*max(max(tmpim))*sh;
shadow(thres_index)=0;

makenima(shadow,tmpim);
c = colorbar;
ylabel(c,'height/nm')
end

function makenima(shadow, height)
% "MAKE Nice IMAges" -- creates one figure with two graphs in it:
% the top one is a 3D view of a height profile (data in "height") with its
% corresponding shadow (data in "shadow");
% the bottom one is a side view of the 3D height profile

% Create figure window and set nice colormap
figure1 = figure('Colormap',...
    [0.8 0.8 0.8;0.7805 0.7805 0.7805;0.7611 0.7611 0.7611;0.7416 0.7416 0.7416;0.7222 0.7222 0.7222;0.7027 0.7027 0.7027;0.6833 0.6833 0.6833;0.6638 0.6638 0.6638;0.6444 0.6444 0.6444;0.6249 0.6249 0.6249;0.6055 0.6055 0.6055;0.586 0.586 0.586;0.5666 0.5666 0.5666;0.5471 0.5471 0.5471;0.5277 0.5277 0.5277;0.5082 0.5082 0.5082;0.4888 0.4888 0.4888;0.4693 0.4693 0.4693;0.4499 0.4499 0.4499;0.4304 0.4304 0.4304;0.411 0.411 0.411;0.3915 0.3915 0.3915;0.3721 0.3721 0.3721;0.3526 0.3526 0.3526;0.3332 0.3332 0.3332;0.3137 0.3137 0.3137;0 0 0;0.07922 0.02275 0.01647;0.1584 0.04549 0.03294;0.2376 0.06824 0.04941;0.3169 0.09098 0.06588;0.3961 0.1137 0.08235;0.4484 0.1366 0.1039;0.5007 0.1595 0.1255;0.5529 0.1824 0.1471;0.6052 0.2052 0.1686;0.6575 0.2281 0.1902;0.7098 0.251 0.2118;0.7284 0.3103 0.226;0.7471 0.3696 0.2402;0.7657 0.4289 0.2544;0.7843 0.4882 0.2686;0.8029 0.5475 0.2828;0.8216 0.6069 0.2971;0.8402 0.6662 0.3113;0.8588 0.7255 0.3255;0.8729 0.7518 0.389;0.8871 0.778 0.4525;0.9012 0.8043 0.5161;0.9153 0.8306 0.5796;0.9294 0.8569 0.6431;0.9435 0.8831 0.7067;0.9576 0.9094 0.7702;0.9718 0.9357 0.8337;0.9859 0.962 0.8973;1 0.9882 0.9608;1 0.9897 0.9657;1 0.9912 0.9706;1 0.9926 0.9755;1 0.9941 0.9804;1 0.9956 0.9853;1 0.9971 0.9902;1 0.9985 0.9951;1 1 1],...
    'Color',[0 0 0]);

% Create axes of bottom plot (profile of 3D data)
axes1 = axes('Parent',figure1,'ZColor',[1 1 1],'YColor',[1 1 1],...
    'XColor',[1 1 1],...
    'Position',[0.06522 0.06634 0.8934 0.1388],...
    'Layer','top',...
    'Color',[0 0 0],...
    'CLim',[1 max(max(height))]);
view([0 0]);
box('on');
hold('all');

% Create surf of bottom plot (profile of 3D data)
surf(height,'Parent',axes1,'SpecularExponent',15,...
    'SpecularStrength',0.6,...
    'FaceLighting','phong',...
    'LineStyle','none',...
    'FaceColor','interp');

% Create axes of top plot (surface & shadow of 3D data)
axes2 = axes('Parent',figure1,'Position',[0.08799 0.305 0.7588 0.6424],...
    'PlotBoxAspectRatio',[3.5 3.5 1.156],...
    'Layer','top',...
    'DataAspectRatio',[0.35 0.35 1],...
    'Color',[0 0 0],...
    'CLim',[0 max(max(height))],...
    'CameraViewAngle',7.162);
view([-79 62]);
hold('all');

% Create "shadow" on the "bottom" of the top image
contourf(shadow,'Parent',axes2,'LineStyle','none');
%contour(shadow,'LineStyle','none','LineColor',[0.3137 0.3176 0.3137],...
%    'Fill','on',...
%    'Parent',axes2);

% Create light object (further specified later with 'SpecularStrength',
% 'DiffuseStrength', 'AmbientStrength' etc.)
light('Parent',axes2);

% Create height profile (surface object) on the top plot
surf(height,'Parent',axes2,'SpecularExponent',20,...
    'SpecularStrength',0.3,...
    'DiffuseStrength',0.3,...
    'AmbientStrength',0.7,...
    'BackFaceLighting','unlit',...
    'FaceLighting','phong',...
    'LineStyle','none',...            % get rid of black lines in contour
    'FaceColor','interp');            % smooth colours of surface
hold off
end

