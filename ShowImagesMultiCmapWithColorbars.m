function ShowImagesMultiCmapWithColorbars(imraw, cmaps, clims, titlestr, blackbg)
[a,b,c] = size(imraw);

if nargin < 2 || isempty(cmaps)
    cmaps = {'red','green','blue','orange','magenta','cyan'};
end
if nargin < 3
    clims = [];
end
if nargin < 4
    titlestr = [];
end
if nargin < 5
    blackbg = true;
end

nr = ceil(sqrt(c));
nc = round(c/nr);
t = tiledlayout(nr,nc,'TileSpacing','compact','Padding','compact');

for k = 1:c
    ax = nexttile;
    imagesc(imraw(:,:,k));
    axis image off

    cmap = get_this_cmap_local(cmaps, k);
    colormap(ax, cmap);

    clim = get_this_clim_local(clims, imraw(:,:,k), k);
    caxis(ax, clim);

    cb = colorbar(ax);
    cb.Color = [1 1 1];

    if ~isempty(titlestr) && numel(titlestr) >= k
        title(titlestr{k}, 'Color', [1 1 1]);
    end

    if blackbg
        set(ax, 'Color', 'k');
    end
end

if blackbg
    set(gcf, 'Color', 'k');
end
end


function cmap = get_this_cmap_local(cmaps, idx)
if iscell(cmaps)
    entry = cmaps{min(idx,numel(cmaps))};
else
    entry = cmaps;
end

if isnumeric(entry) && size(entry,2)==3
    cmap = entry;
elseif isa(entry,'function_handle')
    cmap = entry(256);
else
    switch lower(char(string(entry)))
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
        otherwise
            cmap = parula(256);
    end
end
end


function clim = get_this_clim_local(clims, img, idx)
if isempty(clims)
    clim = [min(img(:)) max(img(:))];
elseif isequal(size(clims), [1 2])
    clim = clims;
elseif size(clims,2) == 2
    clim = clims(min(idx,size(clims,1)),:);
else
    error('Invalid clims format');
end

if clim(2) <= clim(1)
    clim = [clim(1) clim(1)+1];
end
end


function cmap = black2color(rgb, N)
if nargin < 2
    N = 256;
end
rgb = rgb(:).';
t = linspace(0,1,N).';
cmap = t .* rgb;
end