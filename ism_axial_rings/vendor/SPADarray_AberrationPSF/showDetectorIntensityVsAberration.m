function out = showDetectorIntensityVsAberration(varargin)
%--------------------------------------------------------------------------
% showDetectorIntensityVsAberration
%
% PURPOSE
%   Show how the light from a single on-axis point emitter is distributed
%   across the 23-channel SPAD detector array, in focus and at several
%   axial (defocus) positions, for a set of single Zernike aberrations.
%
%   This is a forward demo of the same physics used by
%   estimateFullStackISMWavefront: the (aberrated) detection PSF is formed
%   with scalarPSFBessel and then integrated over each finite detector
%   pixel with detectorCollectionEfficiencyExplicit. The per-channel
%   collected intensity is drawn on the honeycomb layout with
%   plotDetectorHexMap.
%
% USAGE
%   showDetectorIntensityVsAberration();
%   showDetectorIntensityVsAberration('amplitudeWaves',0.2);
%   out = showDetectorIntensityVsAberration('detXY', measuredDetXY, ...
%             'zDefocusUm', -0.8:0.4:0.8);
%
% NAME-VALUE OPTIONS (experimental parameters)
%   'NA'                 objective NA                       (default 1.2)
%   'lamEm'              emission wavelength [um]           (default 0.520)
%   'nMedium'            immersion/medium index for the     (default 1.33)
%                        homogeneous scalar PSF (water)
%   'nImmersion'         water immersion index (stored)     (default 1.33)
%   'nGlass'             coverslip index (stored)           (default 1.518)
%   'nSample'            sample index, e.g. air (stored)    (default 1.0003)
%   'detPitchUm'         detector pitch in sample-equiv.    (default 0.18)
%                        units (ignored if 'detXY' given)
%   'detFillRatio'       detector size / pitch              (default 1.0)
%   'detXY'              [nDet x 2] measured detector        (default [] ->
%                        centres ("shifts"); overrides       honeycomb23)
%                        the nominal honeycomb layout
%   'detectorSubsamples' quadrature points per pixel        (default 7)
%   'modes'              cell list of single modes to show  (default below)
%   'amplitudeWaves'     aberration amplitude [waves RMS]   (default 0.15)
%   'zDefocusUm'         axial positions of the emitter     (default
%                        [-0.6 -0.3 0 0.3 0.6])
%   'fovXY'              lateral PSF field of view [um]     (default 1.6)
%   'nx'                 lateral samples (odd)              (default 81)
%   'normalize'          'mode' (each aberration row scaled  (default 'mode')
%                        to its own max across defocus),
%                        'global' (one shared scale), or
%                        'panel' (each map to its own max)
%   'colormap'           colormap name                      (default 'parula')
%   'showContinuousPSF'  also draw the detector-plane PSF   (default true)
%   'saveFigures'        export PNGs                        (default false)
%   'outputDir'          where to save                      (default
%                        output_matlab/detector_intensity_vs_aberration)
%
% OUTPUT
%   out.detVals     [nAberration x nZ x nDet] collected intensity (relative)
%   out.detXY       detector centres used
%   out.zDefocusUm  axial positions
%   out.labels      aberration row labels
%   out.sim         simulation struct used
%   out.psf         cell array of detection PSF stacks (one per aberration)
%--------------------------------------------------------------------------

    addpath(fileparts(mfilename('fullpath')));
    opts = parseInputs(varargin{:});

    % ----- build the simulation struct from the experimental parameters ---
    sim = defaultParams();
    sim.NA          = opts.NA;
    sim.nMedium     = opts.nMedium;       % drives the homogeneous scalar PSF
    sim.nImmersion  = opts.nImmersion;
    sim.nGlass      = opts.nGlass;
    sim.nSample     = opts.nSample;
    sim.sampleGeometry = 'homogeneous';
    sim.lamEm       = opts.lamEm;
    sim.lamRef      = opts.lamEm;          % express coeffs in waves at emission
    sim.detectorSubsamples = opts.detectorSubsamples;
    sim.detectorImageInverted = true;
    sim.detectorPixelShape = 'hex';

    % Lateral grid (shared sample-equivalent coordinates with the detector)
    sim.fovXY = opts.fovXY;
    sim.nx = opts.nx;  sim.ny = opts.nx;
    sim.x  = linspace(-sim.fovXY/2, sim.fovXY/2, sim.nx);
    sim.y  = sim.x;
    sim.dx = abs(sim.x(2)-sim.x(1));

    % Detector layout: measured centres ("shifts") or nominal honeycomb
    if isempty(opts.detXY)
        sim.detXY = detectorLayout('honeycomb23', opts.detPitchUm);
        sim.detPitch = opts.detPitchUm;
    else
        sim.detXY = double(opts.detXY);
        sim.detPitch = medianNearestPitch(sim.detXY);
    end
    sim.nDet = size(sim.detXY,1);
    sim.detSize = opts.detFillRatio * sim.detPitch;
    sim.detectorHexRadius = sim.detSize / sqrt(3);

    % Axial (defocus) positions are the PSF z planes
    sim.z  = opts.zDefocusUm(:).';
    sim.nz = numel(sim.z);

    % On-axis emitter -> centre sample pixel
    [~, xc] = min(abs(sim.x));
    [~, yc] = min(abs(sim.y));

    % ----- aberration scenarios ------------------------------------------
    scen = buildScenarios(opts);     % {label, coeffStruct}
    nAb  = size(scen,1);
    nZ   = sim.nz;

    % Two-line row labels: aberration name on top, induced amount beneath
    rowLabels = cell(nAb,1);
    for j = 1:nAb
        if isempty(fieldnames(scen{j,2}))
            rowLabels{j} = scen{j,1};                                % no aberration
        else
            rowLabels{j} = {scen{j,1}, sprintf('%.2f\\lambda', opts.amplitudeWaves)};
        end
    end

    detVals = zeros(nAb, nZ, sim.nDet);
    psfStack = cell(nAb,1);
    for j = 1:nAb
        coeffs = scen{j,2};
        hDet = psfBessel(sim, coeffs, sim.lamEm);                % [ny nx nz]
        eta  = detectorCollectionEfficiencyExplicit(sim, hDet);  % [ny nx nz nDet]
        for iz = 1:nZ
            detVals(j,iz,:) = squeeze(eta(yc,xc,iz,:));
        end
        psfStack{j} = hDet;
    end

    % ----- normalisation for display -------------------------------------
    %   'mode'   : each aberration row scaled to its own max across defocus
    %   'global' : every panel shares one whole-figure scale
    %   'panel'  : each panel auto-scaled to its own max
    dispVals = detVals;
    switch lower(opts.normalize)
        case 'global'
            s = max(detVals(:));  if s <= 0, s = 1; end
            dispVals = detVals / s;
            sharedCLim = [0 1];
            cbLabel = 'intensity (relative to global max)';
        case {'mode','row'}
            for j = 1:nAb
                s = max(reshape(detVals(j,:,:),1,[]));  if s <= 0, s = 1; end
                dispVals(j,:,:) = detVals(j,:,:) / s;
            end
            sharedCLim = [0 1];
            cbLabel = 'intensity (relative to each mode''s max)';
        case 'panel'
            sharedCLim = [];
            cbLabel = '';
        otherwise
            error('showDetectorIntensityVsAberration:BadNormalize', ...
                'normalize must be ''mode'', ''global'', or ''panel''.');
    end

    figH = figure('Color','w', 'Position', figurePositionForGrid(nAb, nZ));
    tl = tiledlayout(figH, nAb, nZ, 'Padding','compact', 'TileSpacing','compact');
    lastAx = [];
    for j = 1:nAb
        for iz = 1:nZ
            ax = nexttile(tl);
            vals = squeeze(dispVals(j,iz,:));
            plotDetectorHexMap(sim.detXY, vals, 'Parent', ax, ...
                'EdgeColor',[0.3 0.3 0.3], 'LineWidth',0.3);
            colormap(ax, opts.colormap);
            if ~isempty(sharedCLim), caxis(ax, sharedCLim); end
            addPanelLabels(ax, j, iz, rowLabels{j}, sim.z(iz));
            lastAx = ax;
        end
    end

    if ~isempty(sharedCLim) && ~isempty(lastAx)
        try
            cb = colorbar(lastAx);
            cb.Layout.Tile = 'east';
            cb.Label.String = cbLabel;
        catch
        end
    end
    title(tl, { sprintf(['Detector-plane intensity vs single aberration and defocus' ...
        '   (NA %.2f, \\lambda_{em} %.0f nm, n %.3f, pitch %.0f nm)'], ...
        sim.NA, 1000*sim.lamEm, sim.nMedium, 1000*sim.detPitch), ...
        ' ' }, 'FontWeight','bold');

    if opts.saveFigures
        ensureDir(opts.outputDir);
        exportgraphics(figH, fullfile(opts.outputDir, ...
            'detector_intensity_hexmaps.png'), 'Resolution', 180);
    end

    % ----- optional: continuous detector-plane PSF with pixel outlines ----
    if opts.showContinuousPSF
        figP = figure('Color','w', 'Position', figurePositionForGrid(nAb, nZ));
        tlp = tiledlayout(figP, nAb, nZ, 'Padding','compact', 'TileSpacing','compact');
        cellR = sim.detSize/sqrt(3);
        for j = 1:nAb
            hDet = psfStack{j};
            pmax = max(hDet(:));
            for iz = 1:nZ
                ax = nexttile(tlp);
                imagesc(ax, sim.x, sim.y, hDet(:,:,iz), [0 pmax]);
                axis(ax,'image'); set(ax,'YDir','normal');
                colormap(ax, opts.colormap);
                hold(ax,'on');
                drawHexOutlines(ax, sim.detXY, cellR);
                hold(ax,'off');
                set(ax,'XTick',[],'YTick',[]);
                if j==1, title(ax, sprintf('z = %+.2f \\mum', sim.z(iz))); end
                if iz==1
                    ylabel(ax, rowLabels{j}, 'Rotation',0, 'FontWeight','bold', ...
                        'HorizontalAlignment','right', 'VerticalAlignment','middle');
                end
            end
        end
        title(tlp, {'Detection PSF at the detector plane (hex = SPAD pixels)', ' '}, ...
            'FontWeight','bold');
        if opts.saveFigures
            exportgraphics(figP, fullfile(opts.outputDir, ...
                'detector_plane_psf.png'), 'Resolution', 180);
        end
    end

    % ----- pack output ----------------------------------------------------
    out = struct();
    out.detVals   = detVals;
    out.detXY     = sim.detXY;
    out.zDefocusUm = sim.z;
    out.labels    = scen(:,1);
    out.sim       = sim;
    out.psf       = psfStack;
    out.options   = opts;
end

% ========================================================================
function opts = parseInputs(varargin)
    p = inputParser;
    p.FunctionName = 'showDetectorIntensityVsAberration';
    addParameter(p, 'NA', 1.2);
    addParameter(p, 'lamEm', 0.520);
    addParameter(p, 'nMedium', 1.33);
    addParameter(p, 'nImmersion', 1.33);
    addParameter(p, 'nGlass', 1.518);
    addParameter(p, 'nSample', 1.0003);
    addParameter(p, 'detPitchUm', 0.18);
    addParameter(p, 'detFillRatio', 1.0);
    addParameter(p, 'detXY', []);
    addParameter(p, 'detectorSubsamples', 7);
    addParameter(p, 'modes', {'tilt_x','tilt_y','defocus', ...
        'astig_x','astig_y','coma_x','coma_y','spherical'});
    addParameter(p, 'amplitudeWaves', 0.15);
    addParameter(p, 'zDefocusUm', [-0.6 -0.3 0 0.3 0.6]);
    addParameter(p, 'fovXY', 1.6);
    addParameter(p, 'nx', 81);
    addParameter(p, 'normalize', 'mode');
    addParameter(p, 'colormap', 'parula');
    addParameter(p, 'showContinuousPSF', true);
    addParameter(p, 'saveFigures', false);
    addParameter(p, 'outputDir', '');
    parse(p, varargin{:});
    opts = p.Results;

    if ischar(opts.modes) || isstring(opts.modes)
        opts.modes = cellstr(opts.modes);
    end
    if mod(opts.nx,2) == 0, opts.nx = opts.nx + 1; end   % keep an on-axis pixel
    if isempty(opts.outputDir)
        root = fileparts(fileparts(mfilename('fullpath')));
        opts.outputDir = fullfile(root, 'output_matlab', ...
            'detector_intensity_vs_aberration');
    end
end

% ------------------------------------------------------------------------
function scen = buildScenarios(opts)
    A = opts.amplitudeWaves;
    scen = {'No aberration', struct()};
    for k = 1:numel(opts.modes)
        mode = opts.modes{k};
        scen(end+1,:) = { prettyMode(mode), struct(mode, A) }; %#ok<AGROW>
    end
end

% ------------------------------------------------------------------------
function name = prettyMode(mode)
    switch lower(mode)
        case 'tilt_x',    name = 'Tilt X';
        case 'tilt_y',    name = 'Tilt Y';
        case 'defocus',   name = 'Defocus';
        case 'astig_x',   name = 'Astig 0\circ';
        case 'astig_y',   name = 'Astig 45\circ';
        case 'coma_x',    name = 'Coma X';
        case 'coma_y',    name = 'Coma Y';
        case 'spherical', name = 'Spherical';
        otherwise,        name = strrep(mode,'_','\_');
    end
end

% ------------------------------------------------------------------------
function addPanelLabels(ax, j, iz, rowLabel, zValue)
    % Use axis decorations (title / ylabel) so tiledlayout reserves margin
    % and nothing is clipped. Re-enable the axes with hidden rulers so it
    % keeps plotDetectorHexMap's clean look while still carrying the labels.
    set(ax, 'Visible','on', 'Box','off', 'Color','none', ...
        'XColor','none', 'YColor','none', 'XTick',[], 'YTick',[]);
    if j == 1
        title(ax, sprintf('z = %+.2f \\mum', zValue), 'FontWeight','normal');
    end
    if iz == 1
        ylabel(ax, rowLabel, 'Rotation',0, 'Color','k', 'FontWeight','bold', ...
            'HorizontalAlignment','right', 'VerticalAlignment','middle');
    end
end

% ------------------------------------------------------------------------
function drawHexOutlines(ax, detXY, cellRadius)
    theta = (0:5)*pi/3 + pi/6;
    hx = cellRadius*cos(theta);
    hy = cellRadius*sin(theta);
    for k = 1:size(detXY,1)
        plot(ax, detXY(k,1)+[hx hx(1)], detXY(k,2)+[hy hy(1)], ...
            '-', 'Color',[1 1 1 0.5], 'LineWidth',0.4);
    end
end

% ------------------------------------------------------------------------
function pos = figurePositionForGrid(nRows, nCols)
    % Compact per-panel size, clamped to the screen and centred so tall
    % grids (many aberrations) still fit.
    panelW = 150; panelH = 120;          % per-tile pixels (kept small)
    marginW = 130; marginH = 95;         % row labels + title + colorbar
    w = panelW*nCols + marginW;
    h = panelH*nRows + marginH;
    scr = get(groot, 'ScreenSize');      % [x y width height]
    w = min(w, round(0.92*scr(3)));
    h = min(h, round(0.86*scr(4)));
    x = scr(1) + round((scr(3)-w)/2);
    y = scr(2) + max(round((scr(4)-h)/2), 40);
    pos = [x y w h];
end

% ------------------------------------------------------------------------
function pitch = medianNearestPitch(xy)
    n = size(xy,1);
    nearest = nan(n,1);
    for k = 1:n
        d = hypot(xy(:,1)-xy(k,1), xy(:,2)-xy(k,2));
        d(k) = inf;
        nearest(k) = min(d);
    end
    nearest = nearest(isfinite(nearest) & nearest>0);
    if isempty(nearest), pitch = 1; else, pitch = median(nearest); end
end

% ------------------------------------------------------------------------
function ensureDir(d)
    if ~isempty(d) && exist(d,'dir') ~= 7
        mkdir(d);
    end
end
