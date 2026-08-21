function idxGrid = detectorIndexGridFromXY(detXY, flipX, flipY, transposeGrid, coordTolerance)
% Build a detector index grid from detector coordinates.
%
% Default display convention:
%   columns: x increasing left -> right
%   rows   : y decreasing top -> bottom (so top row is +y)
%
% Optional flags let you test parity mismatches.

    if nargin < 2 || isempty(flipX), flipX = true; end
    if nargin < 3 || isempty(flipY), flipY = true; end
    if nargin < 4 || isempty(transposeGrid), transposeGrid = false; end
    if nargin < 5, coordTolerance = []; end

    if isempty(detXY) || size(detXY,2) ~= 2
        error('detectorIndexGridFromXY:BadDetectorCoordinates', ...
            'detXY must be an [nDetector x 2] coordinate array.');
    end

    [xVals, xBin] = clusteredCoordinateValues(detXY(:,1), coordTolerance);
    [yVals, yBin] = clusteredCoordinateValues(detXY(:,2), coordTolerance);

    % Default: x left->right increasing
    xOrder = 1:numel(xVals);
    if flipX
        xOrder = fliplr(xOrder);
    end

    % Default: top row = +y
    yOrder = fliplr(1:numel(yVals));
    if flipY
        yOrder = fliplr(yOrder);
    end

    idxGrid = nan(numel(yOrder), numel(xOrder));

    for k = 1:size(detXY,1)
        r = find(yOrder == yBin(k), 1, 'first');
        c = find(xOrder == xBin(k), 1, 'first');
        if ~isempty(r) && ~isempty(c) && ~isfinite(idxGrid(r,c))
            idxGrid(r,c) = k;
        end
    end

    if transposeGrid
        idxGrid = idxGrid.';
    end
end

function [centers, bins] = clusteredCoordinateValues(vals, coordTolerance)
    vals = vals(:);
    [sortedVals, order] = sort(vals);

    if isempty(coordTolerance)
        coordTolerance = estimateCoordinateTolerance(sortedVals);
    end

    binsSorted = zeros(size(sortedVals));
    centers = [];
    groupStart = 1;
    groupIdx = 0;

    for ii = 2:numel(sortedVals)
        if sortedVals(ii) - sortedVals(ii-1) > coordTolerance
            groupIdx = groupIdx + 1;
            groupMembers = groupStart:(ii-1);
            centers(groupIdx) = mean(sortedVals(groupMembers)); %#ok<AGROW>
            binsSorted(groupMembers) = groupIdx;
            groupStart = ii;
        end
    end

    groupIdx = groupIdx + 1;
    groupMembers = groupStart:numel(sortedVals);
    centers(groupIdx) = mean(sortedVals(groupMembers));
    binsSorted(groupMembers) = groupIdx;

    bins = zeros(size(vals));
    bins(order) = binsSorted;
end

function tol = estimateCoordinateTolerance(sortedVals)
    if numel(sortedVals) < 2
        tol = 1e-12;
        return;
    end

    diffs = diff(sortedVals);
    scale = max(max(sortedVals) - min(sortedVals), 1);
    positiveDiffs = diffs(diffs > eps(scale));

    if isempty(positiveDiffs)
        tol = 1e-12;
        return;
    end

    sortedDiffs = sort(positiveDiffs);
    if numel(sortedDiffs) < 2
        tol = max(sortedDiffs(1)/2, 1e-12);
        return;
    end

    ratios = sortedDiffs(2:end) ./ max(sortedDiffs(1:end-1), realmin);
    [maxRatio, splitIdx] = max(ratios);
    if maxRatio > 3
        tol = sqrt(sortedDiffs(splitIdx) * sortedDiffs(splitIdx+1));
    else
        tol = sortedDiffs(1) / 2;
    end

    tol = max(tol, 1e-12);
end
