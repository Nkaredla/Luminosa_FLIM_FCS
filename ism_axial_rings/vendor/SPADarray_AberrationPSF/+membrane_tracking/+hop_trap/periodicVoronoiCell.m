function cellId = periodicVoronoiCell(positionUm, mesh, boxSizeUm)
%PERIODICVORONOICELL Return nearest periodic seed for one or more positions.

    if ~mesh.isActive
        cellId = ones(size(positionUm, 1), 1);
        return;
    end

    nPosition = size(positionUm, 1);
    cellId = zeros(nPosition, 1);
    for k = 1:nPosition
        delta = positionUm(k,:) - mesh.seedPositionsUm;
        delta = delta - boxSizeUm .* round(delta ./ boxSizeUm);
        [~, cellId(k)] = min(sum(delta.^2, 2));
    end
end
