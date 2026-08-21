function micro = detectorMicroimage(stack, sim, iy, ix, flipX, flipY, transposeGrid)
% Return detector microimage at scan location (iy,ix) using physical detector order.

    if nargin < 5 || isempty(flipX), flipX = false; end
    if nargin < 6 || isempty(flipY), flipY = false; end
    if nargin < 7 || isempty(transposeGrid), transposeGrid = false; end

    idxGrid = detectorIndexGridFromXY(sim.detXY, flipX, flipY, transposeGrid);
    vals = squeeze(stack(iy, ix, :));
    micro = nan(size(idxGrid));
    valid = isfinite(idxGrid) & idxGrid >= 1 & idxGrid <= numel(vals);
    micro(valid) = vals(idxGrid(valid));
end
