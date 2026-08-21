function canvas = combineDetectorImagesPhysical(stack, sim, scaleMode, flipX, flipY, transposeGrid, gapPx)
% Combine detector-resolved scanned images into a physical detector mosaic
% using the physical detector order from sim.detXY.
%
% scaleMode:
%   'none'   : keep raw scale
%   'each'   : normalize each detector image independently
%   'global' : normalize whole stack globally

    if nargin < 3 || isempty(scaleMode), scaleMode = 'none'; end
    if nargin < 4 || isempty(flipX), flipX = false; end
    if nargin < 5 || isempty(flipY), flipY = false; end
    if nargin < 6 || isempty(transposeGrid), transposeGrid = false; end
    if nargin < 7 || isempty(gapPx), gapPx = 0; end

    idxGrid = detectorIndexGridFromXY(sim.detXY, flipX, flipY, transposeGrid);

    [ny, nx, ~] = size(stack);
    nRows = size(idxGrid,1);
    nCols = size(idxGrid,2);

    hasSparseSlots = any(~isfinite(idxGrid(:)));
    if hasSparseSlots
        colStep = max(1, round((nx + gapPx) / 2));
        rowStep = ny + gapPx;
        canvas = zeros((nRows-1)*rowStep + ny, (nCols-1)*colStep + nx);
    else
        canvas = zeros(nRows*ny, nCols*nx);
    end

    if strcmpi(scaleMode, 'global')
        smax = max(stack(:));
        if smax > 0
            stack = stack / smax;
        end
    end

    for r = 1:nRows
        for c = 1:nCols
            k = idxGrid(r,c);
            if ~isfinite(k) || k < 1 || k > size(stack,3)
                continue;
            end

            img = stack(:,:,k);

            if strcmpi(scaleMode, 'each')
                m = max(img(:));
                if m > 0
                    img = img / m;
                end
            end

            if hasSparseSlots
                rr = (r-1)*rowStep + (1:ny);
                cc = (c-1)*colStep + (1:nx);
            else
                rr = (r-1)*ny + (1:ny);
                cc = (c-1)*nx + (1:nx);
            end
            canvas(rr,cc) = img;
        end
    end
end
