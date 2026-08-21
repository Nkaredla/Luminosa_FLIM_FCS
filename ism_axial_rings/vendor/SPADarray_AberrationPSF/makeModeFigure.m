function makeModeFigure(sim, outFile)
%--------------------------------------------------------------------------
% Generate figure showing:
%   top row    = central detector-channel image
%   bottom row = detector microimage at brightest scan location
%
% Also saves a second figure:
%   *_composite.png
% containing the full detector-resolved scanned images tiled in physical
% detector order.
%--------------------------------------------------------------------------

    scenarios = {
        'None',      struct();
        'Defocus',   struct('defocus', 0.18);
        'Astig x',   struct('astig_x', 0.18);
        'Coma x',    struct('coma_x', 0.18);
        'Spherical', struct('spherical', 0.18)};

    flipX = false;
    flipY = true;
    transposeGrid = false;
    nSc = size(scenarios,1);

    %% -------- Figure 1: central channel + microimage --------
    fig1 = figure('Color','w','Position',[100 100 1300 480]);

    for j = 1:nSc
        stack = normalizedStack(sim, scenarios{j,2});

        % Central detector channel
        [~, kCentre] = min(sum(sim.detXY.^2, 2));
        central = stack(:,:,kCentre);

        % Brightest location in detector-summed image
        s = sum(stack,3);
        [~, idx] = max(s(:));
        [iy, ix] = ind2sub(size(s), idx);

        % Detector microimage values at the brightest scan location
        microVals = squeeze(stack(iy, ix, :));

        % Top row: central detector image
        subplot(2, nSc, j);
        imagesc(sim.x, sim.y, central);
        axis image;
%         set(gca,'YDir','normal');
        title(scenarios{j,1});
        xlabel('x (\mum)');
        ylabel('y (\mum)');
        colorbar;

        % Bottom row: microimage
        ax = subplot(2, nSc, j + nSc);
        plotDetectorHexMap(sim.detXY, microVals, 'Parent', ax, ...
            'EdgeColor', [0.35 0.35 0.35], 'LineWidth', 0.35);
        xlabel('');
        ylabel('');
        colorbar;
    end

    sgtitle('Detector-resolved signatures for representative aberrations (0.18 waves RMS)');
    exportgraphics(fig1, outFile, 'Resolution', 180);
    close(fig1);

    %% -------- Figure 2: composite tiled detector scan images --------
    fig2 = figure('Color','w','Position',[100 100 1300 420]);

    for j = 1:nSc
        stack = normalizedStack(sim, scenarios{j,2});
        canvas = combineDetectorImagesPhysical(stack, sim, 'each', ...
            flipX, flipY, transposeGrid);

        subplot(1, nSc, j);
        imagesc(canvas);
        axis image off;
%         set(gca,'YDir','normal');
        title([scenarios{j,1} ' composite']);
        colormap(gca, hot);
        colorbar;
    end

    sgtitle(sprintf('Composite %d-channel detector scan images in physical detector order', size(sim.detXY,1)));

    % Build second output filename automatically
    [p, n, ~] = fileparts(outFile);
    outFileComposite = fullfile(p, [n '_composite.png']);

    exportgraphics(fig2, outFileComposite, 'Resolution', 180);
%     close(fig2);
end
