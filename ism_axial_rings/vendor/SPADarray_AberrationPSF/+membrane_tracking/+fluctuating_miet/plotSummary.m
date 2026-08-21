function fig = plotSummary(trajectories, frames, detector, ...
        localizationTable, trackTable, surfaceModel, diffusion, opts)
    import membrane_tracking.fluctuating_miet.*

    fig = figure('Color','w','Name','Fluctuating MIET membrane tracking', ...
        'Position',[50 50 1500 860]);
    layout = tiledlayout(fig, 2, 3, 'TileSpacing','compact', ...
        'Padding','compact');

    ax1 = nexttile(layout, 1);
    topography = surfaceModel.topography;
    topographyAvailable = isfield(topography,'isValid') && ...
        topography.isValid && isfield(topography,'snapshotFrame');
    if ~topographyAvailable
        gridVector = linspace(-opts.membraneRadiusUm, ...
            opts.membraneRadiusUm,opts.topographyGridSize);
        [gridX,gridY] = meshgrid(gridVector,gridVector);
        gridZ = opts.tipHeightUm+0.5*opts.curvaturePerUm* ...
            (gridX.^2+gridY.^2);
        gridZ(gridX.^2+gridY.^2 > opts.membraneRadiusUm^2) = NaN;
        topography = struct('xGridUm',gridX,'yGridUm',gridY, ...
            'posteriorHeightUm',reshape(gridZ,size(gridZ,1), ...
            size(gridZ,2),1),'snapshotFrame',numel(trajectories.timeS), ...
            'snapshotTimeS',trajectories.timeS(end),'rmseUm',NaN);
    end
    snapshot = max(1,numel(topography.snapshotFrame));
    surf(ax1, topography.xGridUm, topography.yGridUm, ...
        topography.posteriorHeightUm(:,:,snapshot), ...
        topography.posteriorHeightUm(:,:,snapshot), ...
        'EdgeColor','none','FaceAlpha',0.82);
    hold(ax1,'on');
    path = squeeze(trajectories.positionsXYUm(1,:,:)).';
    plot3(ax1, path(:,1), path(:,2), ...
        squeeze(trajectories.positionsZUm(1,:)).', '-', ...
        'Color',[0.1 0.1 0.1], 'LineWidth',1.1);
    if height(localizationTable) > 0
        scatter3(ax1, localizationTable.xUm, localizationTable.yUm, ...
            localizationTable.heightUm, 9, localizationTable.timeS, ...
            'filled','MarkerFaceAlpha',0.45);
    end
    hold(ax1,'off');
    axis(ax1,'tight'); grid(ax1,'on'); view(ax1, 38, 26);
    xlabel(ax1,'x [um]'); ylabel(ax1,'y [um]'); zlabel(ax1,'z [um]');
    if topographyAvailable
        title(ax1, sprintf('Posterior topography at %.3f s (RMSE %.1f nm)', ...
            topography.snapshotTimeS(snapshot),1000*topography.rmseUm));
    else
        title(ax1,'Topography not estimable; simulated mean shown');
    end

    ax2 = nexttile(layout, 2);
    [~, brightFrame] = max(sum(frames,1));
    plotDetectorHexMap(detector.detXY, frames(:,brightFrame), ...
        'Parent', ax2, 'CellRadius', detector.detectorHexRadiusUm);
    colorbar(ax2);
    title(ax2, sprintf('ISM microimage, frame %d', brightFrame));

    ax3 = nexttile(layout, 3);
    contourf(ax3, topography.xGridUm, topography.yGridUm, ...
        topography.posteriorHeightUm(:,:,snapshot), 16, ...
        'LineColor','none');
    hold(ax3,'on');
    plot(ax3,path(:,1),path(:,2),'-','Color',[0.1 0.1 0.1], ...
        'LineWidth',1.0);
    if height(trackTable) > 0
        trackIds = unique(trackTable.trackId).';
        colors = lines(max(numel(trackIds),1));
        for trackIndex = 1:numel(trackIds)
            rows = trackTable.trackId == trackIds(trackIndex);
            plot(ax3,trackTable.xUm(rows),trackTable.yUm(rows),'o-', ...
                'Color',colors(trackIndex,:), 'MarkerSize',3, ...
                'LineWidth',0.8);
        end
    end
    hold(ax3,'off'); axis(ax3,'equal'); axis(ax3,'tight');
    xlabel(ax3,'x [um]'); ylabel(ax3,'y [um]');
    title(ax3,'Particle positions on reconstructed membrane');
    colorbar(ax3);

    ax4 = nexttile(layout, 4);
    trueHeight = squeeze(trajectories.positionsZUm(1,:));
    plot(ax4,trajectories.timeS,trueHeight,'-','Color',[0.15 0.15 0.15], ...
        'LineWidth',1.1);
    hold(ax4,'on');
    if height(localizationTable) > 0
        errorbar(ax4,localizationTable.timeS,localizationTable.heightUm, ...
            localizationTable.heightSigmaUm,'.','Color',[0.15 0.48 0.78], ...
            'MarkerSize',6,'CapSize',0);
    end
    hold(ax4,'off'); grid(ax4,'on');
    xlabel(ax4,'time [s]'); ylabel(ax4,'height [um]');
    title(ax4,'MIET height along the trajectory');

    ax5 = nexttile(layout, 5);
    if isfield(surfaceModel,'identifiability') && ...
            isfield(surfaceModel.identifiability,'profileNegLogLikelihood')
        id = surfaceModel.identifiability;
        deltaNll = id.profileNegLogLikelihood - ...
            min(id.profileNegLogLikelihood);
        plot(ax5, id.curvatureGridPerUm, deltaNll, 'o-', ...
            'Color',[0.20 0.45 0.70],'LineWidth',1.4);
        hold(ax5,'on');
        yline(ax5, 1.92, 'k--');
        xline(ax5, opts.curvaturePerUm, 'r-','LineWidth',1.4);
        hold(ax5,'off'); grid(ax5,'on');
        xlabel(ax5,'curvature [1/um]');
        ylabel(ax5,'profile \Delta(-log L)');
        title(ax5, sprintf('Amplitude profile, confounding ratio %.2f', ...
            id.confoundingRatio));
    else
        radius = hypot(localizationTable.xUm,localizationTable.yUm);
        scatter(ax5,radius,localizationTable.heightUm,10,'filled');
        grid(ax5,'on');
        xlabel(ax5,'radial position [um]'); ylabel(ax5,'height [um]');
        title(ax5,'Mean membrane profile');
    end

    ax6 = nexttile(layout, 6);
    msd = diffusion.msdTable;
    plot(ax6, msd.lagTimeS, msd.surfaceMsdCorrectedUm2, 'o-', ...
        'Color',[0.80 0.25 0.18],'LineWidth',1.3);
    hold(ax6,'on');
    plot(ax6, msd.lagTimeS, 4*opts.diffusionUm2PerS*msd.lagTimeS, ...
        'k--','LineWidth',1.2);
    hold(ax6,'off'); grid(ax6,'on');
    xlabel(ax6,'lag time [s]'); ylabel(ax6,'surface MSD [um^2]');
    title(ax6,'Surface diffusion');

    colormap(fig,parula(256));
    sgtitle(layout, sprintf(['kappa %+.3g (true %+.3g) 1/um   |   ' ...
        'D %.3g (true %.3g) um^2/s'], surfaceModel.curvaturePerUm, ...
        opts.curvaturePerUm, diffusion.DFisherMLEUm2PerS, ...
        opts.diffusionUm2PerS), 'FontWeight','bold');
end
