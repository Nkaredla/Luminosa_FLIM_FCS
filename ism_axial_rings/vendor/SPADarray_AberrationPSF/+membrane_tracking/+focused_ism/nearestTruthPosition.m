function [nearestMol, truthError] = nearestTruthPosition(truthXY, xy, ...
        photonCounts, opts)
    import membrane_tracking.focused_ism.*

    candidateIds = find(photonCounts > 0);
    if isempty(candidateIds)
        truthError = NaN;
        nearestMol = NaN;
        return;
    end

    delta = bsxfun(@minus, truthXY(candidateIds,:), xy);
    delta = applyMinimumImage(delta, opts);
    d = hypot(delta(:,1), delta(:,2));
    [truthError, localIndex] = min(d);
    nearestMol = candidateIds(localIndex);
end
