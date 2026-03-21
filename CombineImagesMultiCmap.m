function imRGB = CombineImagesMultiCmap(imraw, n, m, cmaps, flag, labelx, labely, fsize, clims)
% CombineImagesMultiCmap
%
% Tile a stack of images into one mosaic, but map each subimage with its
% own colormap before combining.
%
% INPUTS
% ------
% imraw   : [a x b x c] stack of images
% n, m    : number of rows and columns in the tiled layout
% cmaps   : colormap specification
%           - cell array, one entry per image or per component
%           - each entry can be:
%               * a colormap matrix [Nc x 3]
%               * a function handle, e.g. @parula
%               * a char/string, e.g. 'hot', 'turbo', 'parula'
%           - if empty, defaults are used
%
% flag    : optional
%           'scale'  -> each image scaled independently to its own range
%           otherwise uses clims if provided, or global slice min/max
%
% labelx  : optional cell array of length m
% labely  : optional cell array of length n
% fsize   : optional font size
% clims   : optional display limits
%           - []                      : automatic
%           - [min max]               : same for all images
%           - [c x 2]                 : one [min max] pair per image
%
% OUTPUT
% ------
% imRGB   : RGB mosaic image
%
% EXAMPLE
% -------
% cmaps = {'red','green','blue'};
% CombineImagesMultiCmap(cat(3,A1,A2,A3),1,3,cmaps,'scale',{'comp1','comp2','comp3'},{},14);

[a,b,c] = size(imraw);

if nargin < 4 || isempty(cmaps)
    cmaps = default_component_cmaps_black(c);
end
if nargin < 5
    flag = [];
end
if nargin < 6
    labelx = [];
end
if nargin < 7
    labely = [];
end
if nargin < 8 || isempty(fsize)
    fsize = 12;
end
if nargin < 9
    clims = [];
end

imRGB = zeros(n*a, m*b, 3);

for j = 1:n
    for k = 1:m
        idx = (j-1)*m + k;
        if idx <= c
            img = double(imraw(:,:,idx));

            % choose colormap for this panel
            cmap = get_this_cmap(cmaps, idx);

            % choose limits
            if ~isempty(flag) && ischar(flag) && strcmpi(flag,'scale')
                clim = [min(img(:)) max(img(:))];
            else
                clim = get_this_clim(clims, img, idx);
            end

            rgb = map2rgb_local(img, cmap, clim);

            rr = (j-1)*a+1 : j*a;
            cc = (k-1)*b+1 : k*b;
            imRGB(rr, cc, :) = rgb;
        end
    end
end

image(imRGB);
axis image;
set(gca,'YDir','reverse');
hold on
set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');

% grid lines
for j = 1:n-1
    plot(0.5:size(imRGB,2)+0.5, (j*a+0.5)*ones(1,size(imRGB,2)+1), ...
        'w', 'LineWidth', 1);
end
for k = 1:m-1
    plot((k*b+0.5)*ones(1,size(imRGB,1)+1), 0.5:size(imRGB,1)+0.5, ...
        'w', 'LineWidth', 1);
end

% x labels
if ~isempty(labelx) && numel(labelx) == m
    for k = 1:m
        text((2*k-1)*b/2, -size(imRGB,1)/25, labelx{k}, ...
            'HorizontalAlignment', 'center', ...
            'FontSize', fsize, 'Color', 'w', 'FontWeight', 'bold');
    end
end

% y labels
if ~isempty(labely) && numel(labely) == n
    for j = 1:n
        text(-1, (2*j-1)*a/2, labely{j}, ...
            'HorizontalAlignment', 'right', ...
            'FontSize', fsize, 'Color', 'w', 'FontWeight', 'bold');
    end
end

hold off
set(gca,'Position',[0.1 0.1 0.88 0.88]);
axis off

if nargout == 0
    clear imRGB
end
end


% -------------------------------------------------------------------------
function rgb = map2rgb_local(img, cmap, clim)
img = double(img);

if isempty(clim) || numel(clim) ~= 2 || clim(2) <= clim(1)
    clim = [min(img(:)) max(img(:))];
    if clim(2) <= clim(1)
        clim = [clim(1) clim(1)+1];
    end
end

imgn = (img - clim(1)) / (clim(2) - clim(1));
imgn = min(max(imgn, 0), 1);

nC = size(cmap,1);
ind = round(1 + imgn*(nC-1));
ind = min(max(ind,1),nC);

rgb = zeros([size(img) 3]);
rgb(:,:,1) = reshape(cmap(ind,1), size(img));
rgb(:,:,2) = reshape(cmap(ind,2), size(img));
rgb(:,:,3) = reshape(cmap(ind,3), size(img));
end


% -------------------------------------------------------------------------
function cmap = get_this_cmap(cmaps, idx)

if isa(cmaps, 'function_handle')
    cmap = cmaps(256);
    return
end

if isnumeric(cmaps) && size(cmaps,2) == 3
    cmap = cmaps;
    return
end

if ischar(cmaps) || isstring(cmaps)
    cmap = resolve_cmap(cmaps);
    return
end

if iscell(cmaps)
    if numel(cmaps) == 1
        entry = cmaps{1};
    elseif idx <= numel(cmaps)
        entry = cmaps{idx};
    else
        entry = cmaps{mod(idx-1,numel(cmaps))+1};
    end

    if isa(entry, 'function_handle')
        cmap = entry(256);
    elseif isnumeric(entry) && size(entry,2) == 3
        cmap = entry;
    elseif ischar(entry) || isstring(entry)
        cmap = resolve_cmap_black(entry);
    else
        error('Unsupported colormap entry for image %d.', idx);
    end
else
    error('Unsupported cmaps input.');
end
end


% -------------------------------------------------------------------------
function clim = get_this_clim(clims, img, idx)

if isempty(clims)
    clim = [min(img(:)) max(img(:))];
    if clim(2) <= clim(1)
        clim = [clim(1) clim(1)+1];
    end
    return
end

if isnumeric(clims) && isequal(size(clims), [1 2])
    clim = clims;
    return
end

if isnumeric(clims) && size(clims,2) == 2
    if idx <= size(clims,1)
        clim = clims(idx,:);
    else
        clim = clims(end,:);
    end
    return
end

error('clims must be [] or [1x2] or [c x 2].');
end


% -------------------------------------------------------------------------
function cmaps = default_component_cmaps(c)
base = {
    white2color([1 0 0],256)      % red
    white2color([0 0.7 0],256)    % green
    white2color([0 0.35 1],256)   % blue
    white2color([1 0.5 0],256)    % orange
    white2color([0.75 0 0.75],256)% magenta
    white2color([0 0.75 0.75],256)% cyan
    white2color([0.6 0.4 0.1],256)% brown/gold
    turbo(256)
    hot(256)
    parula(256)
    };

cmaps = cell(1,c);
for i = 1:c
    cmaps{i} = base{mod(i-1,numel(base))+1};
end
end

function cmaps = default_component_cmaps_black(c)
base = {
    black2color([1 0 0],256)       % red
    black2color([0 0.9 0],256)     % green
    black2color([0 0.4 1],256)     % blue
    black2color([1 0.6 0],256)     % orange
    black2color([1 0 1],256)       % magenta
    black2color([0 1 1],256)       % cyan
    black2color([1 1 0],256)       % yellow
    black2color([1 1 1],256)       % white
    };

cmaps = cell(1,c);
for i = 1:c
    cmaps{i} = base{mod(i-1,numel(base))+1};
end
end

function cmap = resolve_cmap_black(name)
name = char(lower(string(name)));

switch name
    case 'red'
        cmap = black2color([1 0 0],256);
    case 'green'
        cmap = black2color([0 0.9 0],256);
    case 'blue'
        cmap = black2color([0 0.4 1],256);
    case 'orange'
        cmap = black2color([1 0.6 0],256);
    case 'magenta'
        cmap = black2color([1 0 1],256);
    case 'cyan'
        cmap = black2color([0 1 1],256);
    case 'yellow'
        cmap = black2color([1 1 0],256);
    case 'white'
        cmap = black2color([1 1 1],256);
    otherwise
        error('Unknown black-background colormap name: %s', name);
end
end


% -------------------------------------------------------------------------
function cmap = resolve_cmap(name)
name = char(lower(string(name)));

switch name
    case 'parula'
        cmap = parula(256);
    case 'turbo'
        cmap = turbo(256);
    case 'hot'
        cmap = hot(256);
    case 'gray'
        cmap = gray(256);
    case 'red'
        cmap = white2color([1 0 0],256);
    case 'green'
        cmap = white2color([0 0.7 0],256);
    case 'blue'
        cmap = white2color([0 0.35 1],256);
    case 'orange'
        cmap = white2color([1 0.5 0],256);
    case 'magenta'
        cmap = white2color([1 0 1],256);
    case 'cyan'
        cmap = white2color([0 1 1],256);
    case 'yellow'
        cmap = white2color([1 1 0],256);
    otherwise
        error('Unknown colormap name: %s', name);
end
end


% -------------------------------------------------------------------------
function cmap = white2color(rgb, N)
if nargin < 2
    N = 256;
end
rgb = rgb(:).';
w = [1 1 1];
t = linspace(0,1,N).';
cmap = (1-t).*w + t.*rgb;
cmap = min(max(cmap,0),1);
end

function cmap = black2color(rgb, N)
if nargin < 2
    N = 256;
end
rgb = rgb(:).';
k = [0 0 0];
t = linspace(0,1,N).';
cmap = (1-t).*k + t.*rgb;
cmap = min(max(cmap,0),1);
end