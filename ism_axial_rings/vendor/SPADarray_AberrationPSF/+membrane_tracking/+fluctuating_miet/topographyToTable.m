function output = topographyToTable(topography)
    import membrane_tracking.fluctuating_miet.*

    inside = topography.insideMembrane;
    nInside = sum(inside(:));
    nSnapshots = numel(topography.snapshotFrame);
    frame = repelem(topography.snapshotFrame,nInside);
    timeS = repelem(topography.snapshotTimeS,nInside);
    xUm = repmat(topography.xGridUm(inside),nSnapshots,1);
    yUm = repmat(topography.yGridUm(inside),nSnapshots,1);
    estimatedHeightUm = zeros(nInside*nSnapshots,1);
    posteriorSigmaUm = zeros(nInside*nSnapshots,1);
    trueHeightUm = zeros(nInside*nSnapshots,1);
    for snapshot = 1:nSnapshots
        rows = (snapshot-1)*nInside+(1:nInside);
        layer = topography.posteriorHeightUm(:,:,snapshot);
        estimatedHeightUm(rows) = layer(inside);
        layer = topography.posteriorSigmaUm(:,:,snapshot);
        posteriorSigmaUm(rows) = layer(inside);
        layer = topography.trueHeightUm(:,:,snapshot);
        trueHeightUm(rows) = layer(inside);
    end
    output = table(frame,timeS,xUm,yUm,estimatedHeightUm, ...
        posteriorSigmaUm,trueHeightUm);
end
