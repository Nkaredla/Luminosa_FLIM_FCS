function fig = PDA23im(values, layout)
%PDA23IM Plot 23-detector intensities on the Luminosa honeycomb layout.
%
%   PDA23im(values)
%   PDA23im(values, layout)
%
%   values can be:
%       23 x 1       one detector map
%       1 x 23       one detector map
%       23 x nPlane  one detector map per column
%       nPlane x 23  one detector map per row

if nargin < 2 || isempty(layout)
    layout = 'honeycomb23';
end

[detXY, ~] = detectorLayout(layout, 1);
nDet = size(detXY, 1);

values = double(values);
if isvector(values)
    values = values(:);
elseif size(values, 1) == nDet
    % Already [detector x map].
elseif size(values, 2) == nDet
    values = values.';
else
    error('PDA23im:BadInput', ...
        'values must be a 23-vector, [23 x nMap], or [nMap x 23]. Got [%s].', ...
        num2str(size(values)));
end

if size(values, 1) ~= nDet
    error('PDA23im:BadInput', ...
        'values must contain one value per detector. Expected %d detector rows, got %d.', ...
        nDet, size(values, 1));
end

nMap = size(values, 2);
fig = figure('Color','w');

if nMap == 1
    ax = axes('Parent', fig);
    plotDetectorHexMap(detXY, values(:,1), 'Parent', ax);
    colormap(ax, 'parula');
    colorbar(ax);
    title(ax, '23-detector intensity');
    return;
end

nCol = ceil(sqrt(nMap));
nRow = ceil(nMap / nCol);
tl = tiledlayout(fig, nRow, nCol, 'Padding', 'compact', 'TileSpacing', 'compact');
commonCLim = [min(values(:)), max(values(:))];
if commonCLim(1) == commonCLim(2)
    commonCLim = commonCLim + [-1 1]*max(abs(commonCLim(1)), 1)*eps;
end

for k = 1:nMap
    ax = nexttile(tl);
    plotDetectorHexMap(detXY, values(:,k), 'Parent', ax, 'CLim', commonCLim);
    colormap(ax, 'parula');
    title(ax, sprintf('map %d', k));
end

cb = colorbar(ax);
cb.Layout.Tile = 'east';
title(tl, sprintf('23-detector intensity (%d maps)', nMap));

