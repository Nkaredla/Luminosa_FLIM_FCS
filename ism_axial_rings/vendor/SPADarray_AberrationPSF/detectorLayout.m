function [detXY, idxGrid, info] = detectorLayout(layoutName, detPitch)
%DETECTORLAYOUT Return detector coordinates for supported SPAD layouts.
%
%   [detXY, idxGrid, info] = detectorLayout(layoutName, detPitch)
%
%   layoutName:
%       'honeycomb23'  default Luminosa 23-channel honeycomb layout
%       'honeycomb93'  93-channel honeycomb layout (same point-up hex
%                      packing rule as honeycomb23, scaled to 11 rows)
%       'square25'     legacy 5x5 square layout
%
%   detPitch is the detector pitch in the same sample-equivalent units used
%   by the forward model. detXY is [nDetector x 2]. idxGrid is a sparse
%   display grid; missing honeycomb slots are NaN.

    if nargin < 1 || isempty(layoutName)
        layoutName = 'honeycomb23';
    end
    if nargin < 2 || isempty(detPitch)
        detPitch = 1;
    end

    layoutKey = lower(strrep(char(layoutName), ' ', ''));
    layoutKey = strrep(layoutKey, '-', '');
    layoutKey = strrep(layoutKey, '_', '');

    switch layoutKey
        case {'honeycomb23','hex23','hexagonal23','luminosa23','23'}
            rowCounts = [5 4 5 4 5];
            [detXY, idxGrid] = hexRowLayout(rowCounts, detPitch);

            info.name = 'honeycomb23';
            info.description = '23-channel honeycomb detector layout';
            info.rowCounts = rowCounts;

        case {'honeycomb93','hex93','hexagonal93','luminosa93','93'}
            % Symmetric 11-row taper (flat-topped hexagon, widest row = 11)
            % reconstructed from the 93-pixel array photo; same point-up
            % hex-packing rule as honeycomb23. Replace rowCounts here (or
            % pass an explicit detXY to callers) if measured centers show
            % this doesn't match the real hardware.
            rowCounts = [7 7 8 9 10 11 10 9 8 7 7];
            [detXY, idxGrid] = hexRowLayout(rowCounts, detPitch);

            info.name = 'honeycomb93';
            info.description = '93-channel honeycomb detector layout';
            info.rowCounts = rowCounts;

        case {'square25','legacy25','square5x5','5x5','square','25'}
            coords = ((0:4) - 2) * detPitch;
            [gx, gy] = meshgrid(coords, coords);
            detXY = [gx(:) gy(:)];
            idxGrid = detectorIndexGridFromXY(detXY, false, false, false);

            info.name = 'square25';
            info.description = 'Legacy 25-channel square 5x5 detector layout';
            info.rowCounts = 5 * ones(1, 5);

        otherwise
            error('detectorLayout:UnknownLayout', ...
                'Unknown detector layout "%s". Use ''honeycomb23'', ''honeycomb93'', or ''square25''.', ...
                char(layoutName));
    end

    info.nDet = size(detXY, 1);
    info.detPitch = detPitch;
    info.gridSize = size(idxGrid);
    info.idxGrid = idxGrid;
end

function [detXY, idxGrid] = hexRowLayout(rowCounts, detPitch)
%HEXROWLAYOUT Point-up hexagonal packing from a per-row pixel-count vector.
%
%   Row r holds rowCounts(r) pixels centered on x=0, spaced by detPitch
%   (integer offsets for an odd row count, half-integer offsets for an
%   even one so adjacent rows nest). Rows are stacked top-to-bottom with
%   vertical pitch sqrt(3)/2*detPitch, the touching-hexagon spacing.

    yPitch = sqrt(3)/2 * detPitch;
    maxHalfWidth = (max(rowCounts)-1)/2;
    xGrid = (-maxHalfWidth:0.5:maxHalfWidth) * detPitch;
    yGrid = ((numel(rowCounts)-1)/2:-1:-(numel(rowCounts)-1)/2) * yPitch;

    detXY = zeros(sum(rowCounts), 2);
    idxGrid = nan(numel(rowCounts), numel(xGrid));
    k = 0;

    for r = 1:numel(rowCounts)
        nThisRow = rowCounts(r);
        if mod(nThisRow, 2) == 1
            xRow = (-(nThisRow-1)/2:(nThisRow-1)/2) * detPitch;
        else
            xRow = ((-nThisRow/2):(nThisRow/2-1)) * detPitch + detPitch/2;
        end

        for c = 1:nThisRow
            k = k + 1;
            detXY(k,:) = [xRow(c), yGrid(r)];

            [~, gridCol] = min(abs(xGrid - xRow(c)));
            idxGrid(r, gridCol) = k;
        end
    end
end
