function plotPeriodicVoronoiMesh(ax, mesh, boxSizeUm)
%PLOTPERIODICVORONOIMESH Draw periodic Voronoi edges clipped to the box.

    if ~mesh.isActive
        return;
    end
    offsetsX = [-boxSizeUm(1), 0, boxSizeUm(1)];
    offsetsY = [-boxSizeUm(2), 0, boxSizeUm(2)];
    extended = zeros(0, 2);
    for xOffset = offsetsX
        for yOffset = offsetsY
            extended = [extended; mesh.seedPositionsUm + ...
                [xOffset yOffset]]; %#ok<AGROW>
        end
    end
    hold(ax, 'on');
    lines = voronoi(ax, extended(:,1), extended(:,2), '-');
    set(lines, 'Color', [0.78 0.78 0.78], 'LineWidth', 0.45);
    plot(ax, mesh.seedPositionsUm(:,1), mesh.seedPositionsUm(:,2), '.', ...
        'Color', [0.55 0.55 0.55], 'MarkerSize', 3);
    xlim(ax, 0.5 * boxSizeUm(1) * [-1 1]);
    ylim(ax, 0.5 * boxSizeUm(2) * [-1 1]);
end
