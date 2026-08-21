function res = demoPhaseRetrieval3DBead(useRealData, rawData, planeZ, outDir, opts)
%--------------------------------------------------------------------------
% demoPhaseRetrieval3DBead
%
% PURPOSE
%   Run the multi-plane SPAD bead phase-retrieval workflow and save
%   publication-style phase and PSF comparison figures. By default this
%   simulates a bead recorded at z = [-0.1, 0, +0.1] um.
%
% USAGE
%   res = demoPhaseRetrieval3DBead();                         % simulation
%   res = demoPhaseRetrieval3DBead(true, rawData, planeZ);     % real data
%   res = demoPhaseRetrieval3DBead(false, [], [-0.2 0 0.2]);
%
% INPUT DATA LAYOUTS
%   rawData can be [ny nx nDetector nPlane], [ny nx nPlane nDetector], or an nDetector-cell array
%   whose cells contain [ny nx nPlane] volumes.
%--------------------------------------------------------------------------

    addParentFolderToPath();

    if nargin < 1 || isempty(useRealData), useRealData = false; end
    if ischar(useRealData) || isstring(useRealData)
        useRealData = any(strcmpi(char(useRealData), {'real','data','measured','true'}));
    end
    if nargin < 2, rawData = []; end
    if nargin < 3 || isempty(planeZ), planeZ = [-0.1 0 0.1]; end
    if nargin < 4 || isempty(outDir), outDir = 'output_matlab'; end
    if nargin < 5 || isempty(opts), opts = struct(); end

    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    opts.verbose = getOptionLocal(opts, 'verbose', true);
    opts.maxIter = getOptionLocal(opts, 'maxIter', 6);

    res = phaseRetrieval3DBead(useRealData, rawData, planeZ, opts);

    summaryFile = fullfile(outDir, 'phase_retrieval_3d_bead_summary.png');
    phaseFile = fullfile(outDir, 'phase_retrieval_2d_phase.png');
    psfFile = fullfile(outDir, 'phase_retrieval_array_psf_comparison.png');

    makePhaseRetrievalSummaryFigure(res, summaryFile);
    makePhase2DFigure(res, phaseFile);
    makeArrayVolumeComparisonFigure(res, psfFile);

    res.outputFiles.summary = summaryFile;
    res.outputFiles.phase2D = phaseFile;
    res.outputFiles.arrayComparison = psfFile;

    fprintf('[demoPhaseRetrieval3DBead] Saved %s\n', summaryFile);
    fprintf('[demoPhaseRetrieval3DBead] Saved %s\n', phaseFile);
    fprintf('[demoPhaseRetrieval3DBead] Saved %s\n', psfFile);
end

function addParentFolderToPath()
    thisDir = fileparts(mfilename('fullpath'));
    parentDir = fileparts(thisDir);

    if exist(fullfile(parentDir, 'defaultParams.m'), 'file') == 2
        addpath(parentDir);
    end
end

function makePhaseRetrievalSummaryFigure(res, outFile)
    sim = res.sim;
    nModes = numel(sim.modeOrder);
    xs = 1:nModes;

    est = res.estCoeffVector;
    hasTruth = isfield(res, 'trueCoeffVector');
    if hasTruth
        truth = res.trueCoeffVector;
    else
        truth = zeros(size(est));
    end

    fig = figure('Color','w','Position',[80 80 1300 760]);

    subplot(2,3,1);
    if hasTruth
        bar(xs-0.18, truth, 0.36, 'FaceColor', [0.75 0.25 0.20]);
        hold on;
        bar(xs+0.18, est, 0.36, 'FaceColor', [0.20 0.35 0.75]);
        legend({'true','retrieved'}, 'Box','off', 'Location','northwest');
    else
        bar(xs, est, 0.55, 'FaceColor', [0.20 0.35 0.75]);
    end
    grid on;
    set(gca, 'XTick', xs, 'XTickLabel', sim.modeOrder);
    xtickangle(35);
    ylabel('waves RMS');
    title('Recovered Zernike coefficients');

    subplot(2,3,2);
    plot(1:size(res.history,1), res.history(:,1), '-o', 'LineWidth', 1.3);
    grid on;
    xlabel('iteration');
    ylabel('residual norm');
    title('Gauss-Newton convergence');

    subplot(2,3,3);
    midPlane = ceil(numel(res.planeZ)/2);
    [~, kCentre] = min(sum(sim.detXY.^2, 2));
    dataImg = res.dataN(:,:,kCentre,midPlane);
    fitImg = res.fitStack(:,:,kCentre,midPlane);
    imagesc(sim.x, sim.y, [dataImg fitImg dataImg-fitImg]);
    axis image;
    set(gca,'YDir','normal');
    colorbar;
    title('Central SPAD: data | fit | residual');
    xlabel('x (um)');
    ylabel('y (um)');

    subplot(2,3,4);
    plotPupilPhaseMap(res.phase, 'Retrieved phase', []);

    subplot(2,3,5);
    if hasTruth
        plotPupilPhaseMap(res.truePhase, 'Simulated phase', []);
    else
        axis off;
        text(0.05, 0.5, 'Real-data mode: no simulated truth', 'FontSize', 11);
    end

    subplot(2,3,6);
    planeEnergy = squeeze(sum(sum(sum(res.dataN, 1), 2), 3));
    fitEnergy = squeeze(sum(sum(sum(res.fitStack, 1), 2), 3));
    plot(res.planeZ, planeEnergy, 'o-', 'LineWidth', 1.3);
    hold on;
    plot(res.planeZ, fitEnergy, 's-', 'LineWidth', 1.3);
    hold off;
    grid on;
    xlabel('recorded z plane (um)');
    ylabel('fraction of total counts');
    legend({'data','fit'}, 'Box','off', 'Location','best');
    title(sprintf('Axial energy, z offset = %.3f um', res.estZOffset));

    exportgraphics(fig, outFile, 'Resolution', 180);
    close(fig);
end

function makePhase2DFigure(res, outFile)
    hasTruth = isfield(res, 'truePhase');

    if hasTruth
        cmax = max(abs([res.truePhase.radians(:); res.phase.radians(:)]));
        clim = symmetricLimits(cmax);
        fig = figure('Color','w','Position',[120 120 900 420]);

        subplot(1,2,1);
        plotPupilPhaseMap(res.truePhase, 'Simulated pupil phase', clim);

        subplot(1,2,2);
        plotPupilPhaseMap(res.phase, 'Retrieved pupil phase', clim);
    else
        cmax = max(abs(res.phase.radians(:)));
        clim = symmetricLimits(cmax);
        fig = figure('Color','w','Position',[180 120 520 460]);
        plotPupilPhaseMap(res.phase, 'Retrieved pupil phase', clim);
    end

    colormap(fig, redBlueMap(256));
    exportgraphics(fig, outFile, 'Resolution', 220);
    close(fig);
end

function makeArrayVolumeComparisonFigure(res, outFile)
    sim = res.sim;
    hasTruth = isfield(res, 'trueCoeffs');

    if hasTruth
        truePsf = spadEffectivePSFArray(sim, res.trueCoeffs, 'each');
        fitPsf = spadEffectivePSFArray(sim, res.estCoeffs, 'each');
        renderArrayIsoOverlay( ...
            sim, truePsf, fitPsf, sim.z, outFile, ...
            'Calculated 3D effective PSFs: simulated truth vs retrieved fit', ...
            {'simulated','retrieved'});
    else
        inputVol = permute(res.dataN, [1 2 4 3]);   % y,x,z,ch
        fitVol = permute(res.fitStack, [1 2 4 3]);  % y,x,z,ch

        if numel(res.planeZ) >= 4
            renderArrayIsoOverlay( ...
                sim, inputVol, fitVol, res.planeZ, outFile, ...
                'Measured detector input volumes vs fitted volumes', ...
                {'input','fitted'});
        else
            renderArrayProjectionComparison( ...
                sim, inputVol, fitVol, outFile, ...
                'Measured detector input volumes vs fitted volumes', ...
                {'input','fitted'});
        end
    end
end

function renderArrayIsoOverlay(sim, volA, volB, zVec, outFile, titleText, labels)
    idxGrid = detectorIndexGridFromXY(sim.detXY, false, true, false);
    nRows = size(idxGrid, 1);
    nCols = size(idxGrid, 2);
    [X3D, Y3D, Z3D] = meshgrid(sim.x, sim.y, zVec);

    fig = figure('Color','w','Position',[40 40 1550 1250]);
    sgtitle(titleText, 'FontSize', 13, 'FontWeight', 'bold');

    for r = 1:nRows
        for c = 1:nCols
            k = idxGrid(r,c);
            ax = subplot(nRows, nCols, (r-1)*nCols + c);
            if ~isfinite(k) || k < 1 || k > size(sim.detXY,1)
                axis(ax, 'off');
                continue;
            end
            hold(ax, 'on');

            A = volA(:,:,:,k);
            B = volB(:,:,:,k);
            m = max([A(:); B(:)]);
            if m > 0
                A = A / m;
                B = B / m;
            end

            addIsoPatch(X3D, Y3D, Z3D, A, 0.25, [0.80 0.20 0.15], 0.32);
            addIsoPatch(X3D, Y3D, Z3D, B, 0.25, [0.15 0.32 0.85], 0.32);

            hold(ax, 'off');
            axis(ax, 'image');
            axis(ax, 'vis3d');
            axis(ax, 'off');
            view(ax, 3);
            axes(ax);
            light('Style', 'infinite', 'Position', [0 1 1]);
            light('Style', 'infinite', 'Position', [1 -1 -0.5]);
            lighting gouraud;
            dk = sim.detXY(k,:);
            title(ax, sprintf('(%.2f, %.2f)', dk(1), dk(2)), ...
                'FontSize', 6.5, 'Interpreter', 'none');
        end
    end

    annotation(fig, 'textbox', [0.02 0.01 0.96 0.035], ...
        'String', sprintf('Overlay per detector pixel: red = %s, blue = %s; isosurface = 25%% of per-channel peak', labels{1}, labels{2}), ...
        'EdgeColor', 'none', 'FontSize', 9, 'HorizontalAlignment', 'center');

    exportgraphics(fig, outFile, 'Resolution', 170);
    close(fig);
end

function addIsoPatch(X3D, Y3D, Z3D, V, level, color, alphaValue)
    if max(V(:)) <= level
        return;
    end

    try
        p = patch(isosurface(X3D, Y3D, Z3D, V, level));
        set(p, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alphaValue);
    catch
        % Sparse measured z stacks can be too thin for isosurface; callers
        % use projection rendering for those cases.
    end
end

function renderArrayProjectionComparison(sim, volA, volB, outFile, titleText, labels)
    idxGrid = detectorIndexGridFromXY(sim.detXY, false, true, false);
    nRows = size(idxGrid, 1);
    nCols = size(idxGrid, 2);
    ny = size(volA, 1);
    nx = size(volA, 2);

    gap = 2;
    tileW = 2*nx + gap;
    canvas = zeros(nRows*ny, nCols*tileW);

    for r = 1:nRows
        for c = 1:nCols
            k = idxGrid(r,c);
            if ~isfinite(k) || k < 1 || k > size(sim.detXY,1)
                continue;
            end
            A = max(volA(:,:,:,k), [], 3);
            B = max(volB(:,:,:,k), [], 3);
            m = max([A(:); B(:)]);
            if m > 0
                A = A / m;
                B = B / m;
            end

            rr = (r-1)*ny + (1:ny);
            cc0 = (c-1)*tileW;
            canvas(rr, cc0 + (1:nx)) = A;
            canvas(rr, cc0 + nx + gap + (1:nx)) = B;
        end
    end

    fig = figure('Color','w','Position',[80 80 1500 920]);
    imagesc(canvas);
    axis image off;
    colormap(fig, hot);
    colorbar;
    title(sprintf('%s: XY maximum projections', titleText), ...
        'FontSize', 12, 'FontWeight', 'bold');
    annotation(fig, 'textbox', [0.02 0.01 0.96 0.035], ...
        'String', sprintf('Each detector tile: left = %s, right = %s; intensities normalized per detector channel', labels{1}, labels{2}), ...
        'EdgeColor', 'none', 'FontSize', 9, 'HorizontalAlignment', 'center');

    exportgraphics(fig, outFile, 'Resolution', 180);
    close(fig);
end

function plotPupilPhaseMap(phase, titleText, clim)
    x = phase.rho .* cos(phase.phi);
    y = phase.rho .* sin(phase.phi);
    q = linspace(-1, 1, 220);
    [Xq, Yq] = meshgrid(q, q);
    Vq = griddata(x(:), y(:), phase.radians(:), Xq, Yq, 'linear');
    Vq(sqrt(Xq.^2 + Yq.^2) > 1) = NaN;

    imagesc(q, q, Vq);
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel('pupil x');
    ylabel('pupil y');
    title(titleText);
    cb = colorbar;
    cb.Label.String = 'phase (rad)';
    if ~isempty(clim)
        caxis(clim);
    end
end

function clim = symmetricLimits(cmax)
    if isempty(cmax) || ~isfinite(cmax) || cmax <= 0
        cmax = 1;
    end
    clim = [-cmax cmax];
end

function cmap = redBlueMap(n)
    if nargin < 1 || isempty(n)
        n = 256;
    end

    x = linspace(0, 1, n).';
    cmap = zeros(n, 3);
    cmap(:,1) = min(1, 2*x);
    cmap(:,2) = 1 - abs(2*x - 1);
    cmap(:,3) = min(1, 2*(1-x));
end

function v = getOptionLocal(opts, name, defaultValue)
    if isfield(opts, name) && ~isempty(opts.(name))
        v = opts.(name);
    else
        v = defaultValue;
    end
end
