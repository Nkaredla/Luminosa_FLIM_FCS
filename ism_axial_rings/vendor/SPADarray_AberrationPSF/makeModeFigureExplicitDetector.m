function makeModeFigureExplicitDetector(sim, outFile)
%--------------------------------------------------------------------------
% makeModeFigureExplicitDetector
%
% Like makeModeFigure, but uses normalizedStackExplicitDetector and saves
% detector-channel scan-image mosaics as well.
%--------------------------------------------------------------------------

    scenarios = {
        'None',      struct();
        'Defocus',   struct('defocus', 0.18);
        'Astig x',   struct('astig_x', 0.18);
        'Coma x',    struct('coma_x', 0.18);
        'Spherical', struct('spherical', 0.18)};

    flipX = true;
    flipY = false;
    transposeGrid = false;
    nSc = size(scenarios,1);

    fig1 = figure('Color','w','Position',[100 100 1300 480]);

    for j = 1:nSc
        stack = normalizedStackExplicitDetector(sim, scenarios{j,2});

        [~, kCentre] = min(sum(sim.detXY.^2, 2));
        central = stack(:,:,kCentre);

        s = sum(stack,3);
        [~, idx] = max(s(:));
        [iy, ix] = ind2sub(size(s), idx);
        microVals = squeeze(stack(iy, ix, :));

        subplot(2, nSc, j);
        imagesc(sim.x, sim.y, central);
        axis image;
        set(gca,'YDir','normal');
        title(scenarios{j,1});
        xlabel('x (\mum)');
        ylabel('y (\mum)');
        colorbar;

        ax = subplot(2, nSc, j+nSc);
        plotDetectorHexMap(sim.detXY, microVals, 'Parent', ax, ...
            'EdgeColor', [0.35 0.35 0.35], 'LineWidth', 0.35);
        xlabel('');
        ylabel('');
        colorbar;
    end

    sgtitle('Explicit detector-plane integration model');
    exportgraphics(fig1, outFile, 'Resolution', 180);
%     close(fig1);

    fig2 = figure('Color','w','Position',[100 100 1300 420]);

    for j = 1:nSc
        stack = normalizedStackExplicitDetector(sim, scenarios{j,2});
        canvas = combineDetectorImagesPhysical(stack, sim, 'each', ...
            flipX, flipY, transposeGrid);

        subplot(1, nSc, j);
        imagesc(canvas);
        axis image off;
        set(gca,'YDir','normal');
        title([scenarios{j,1} ' composite']);
        colormap(gca, hot);
        colorbar;
    end

    sgtitle(sprintf('Composite %d-channel detector scan images: explicit detector integration', size(sim.detXY,1)));
    [p, n, ~] = fileparts(outFile);
    outFileComposite = fullfile(p, [n '_composite.png']);
    exportgraphics(fig2, outFileComposite, 'Resolution', 180);
%     close(fig2);
end
