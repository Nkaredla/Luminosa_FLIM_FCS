function out = reconstruct3DSumPSFFromFocalFit(result, varargin)
%RECONSTRUCT3DSUMPSFFROMFOCALFIT 3D detector-summed PSF from a focal-plane aberration fit.
%
%   out = reconstruct3DSumPSFFromFocalFit(result)
%   out = reconstruct3DSumPSFFromFocalFit(result, 'zRangeUm', 3, 'nz', 61)
%
%   Takes the in-plane aberrations estimated by
%   estimateFocalPlaneCenterPointISMPoisson (result.fit.estCoeffs, in waves
%   RMS at result.fit.sim.lamRef) and reconstructs the 3D intensity PSF over
%   a lateral/axial grid using the SAME optical model the fit used (vectorial
%   Richards-Wolf, glass/air interface by default), then draws the xy, yz,
%   and zx sum-projections.
%
%   "Sum intensity PSF" is the SPAD-array detector-summed system PSF:
%   for each source position the excitation PSF times the total collection
%   efficiency of all detector pixels,
%
%       h_sum(x,y,z) = h_exc(x,y,z) * sum_k eta_k(x,y,z),
%
%   i.e. the ISM open-array ("sum image") PSF. Use 'psfType' to switch to the
%   bare detection or excitation PSF instead.
%
%   INPUT
%     result : struct returned by estimateFocalPlaneCenterPointISMPoisson
%              (must contain result.fit.estCoeffs and result.fit.sim), or a
%              coefficient struct together with a 'sim' option.
%
%   NAME-VALUE OPTIONS
%     'coeffs'        override the aberration coeff struct (default
%                     result.fit.estCoeffs)
%     'sim'           override the base sim struct (default result.fit.sim)
%     'excitationPolarizationMode' override the excitation polarization
%                     ('circular'/'xyAverage'/'x'/'y'); default inherits the
%                     fit's sim. Use 'circular' for a circularly polarized
%                     excitation laser (rotationally symmetric focus).
%     'collectionPolarizationMode' likewise for the collection path
%     'psfType'       'effective' (h_exc .* sum_k eta_k, default),
%                     'detectionSum' (sum_k eta_k), 'detection' (h_det),
%                     or 'excitation' (h_exc)
%     'projection'    'sum' (default) or 'max' projection along each axis
%     'fovXYUm'       lateral field of view [um]        (default 2.0)
%     'nx'            lateral samples (odd)             (default 101)
%     'zRangeUm'      total axial range [um]            (default 3.0)
%     'nz'            axial planes (odd)                (default 61)
%     'emitterHeightUm' source height above the glass   (default from sim)
%     'normalizeDisplay' scale each projection to its   (default true)
%                     own max for display
%     'logScale'      log10 display                     (default false)
%     'colormap'      colormap name                     (default 'inferno'
%                     if available, else 'hot')
%     'makeFigure'    draw the projection figure        (default true)
%     'saveFigures'   export a PNG                       (default false)
%     'outputDir'     export folder                      (default under
%                     output_matlab/focal_fit_3d_psf)
%     'verbose'       print a short summary              (default true)
%
%   OUTPUT
%     out.volume        [ny nx nz] reconstructed sum-intensity PSF
%     out.x, out.y, out.z   coordinate axes [um]
%     out.projections   struct with .xy [ny nx], .yz [ny nz], .zx [nz nx]
%     out.coeffs        aberration coeffs used
%     out.psfType, out.projectionMode
%     out.sim           sim used for the reconstruction (grid overridden)
%     out.figure        figure handle (or [])

    addpath(fileparts(mfilename('fullpath')));
    opts = parseInputs(varargin{:});

    [coeffs, baseSim] = resolveFitInputs(result, opts);
    sim = configurePSFGrid(baseSim, opts);
    emitterHeight = resolveEmitterHeight(sim, opts);

    [volume, hExc, hDetSum, hDet] = build3DVolume(sim, coeffs, emitterHeight);
    volume = selectPsfType(opts.psfType, volume, hExc, hDetSum, hDet);

    proj = projectVolume(volume, opts.projection);

    figHandle = [];
    if opts.makeFigure
        figHandle = plotProjections(proj, sim, coeffs, opts);
    end

    if opts.saveFigures && ~isempty(figHandle)
        ensureDir(opts.outputDir);
        exportgraphics(figHandle, fullfile(opts.outputDir, ...
            'focal_fit_3d_sum_psf_projections.png'), 'Resolution', 180);
    end

    if opts.verbose
        printSummary(coeffs, sim, opts);
    end

    out = struct();
    out.volume         = volume;
    out.x              = sim.x(:).';
    out.y              = sim.y(:).';
    out.z              = sim.z(:).';
    out.projections    = proj;
    out.coeffs         = coeffs;
    out.psfType        = opts.psfType;
    out.projectionMode = opts.projection;
    out.emitterHeightUm = emitterHeight;
    out.sim            = sim;
    out.figure         = figHandle;
    out.options        = opts;
end

% ============================================================================
function opts = parseInputs(varargin)
    p = inputParser;
    p.FunctionName = 'reconstruct3DSumPSFFromFocalFit';
    addParameter(p, 'coeffs', []);
    addParameter(p, 'sim', []);
    addParameter(p, 'excitationPolarizationMode', []);
    addParameter(p, 'collectionPolarizationMode', []);
    addParameter(p, 'psfType', 'effective');
    addParameter(p, 'projection', 'sum');
    addParameter(p, 'fovXYUm', 2.0);
    addParameter(p, 'nx', 101);
    addParameter(p, 'zRangeUm', 3.0);
    addParameter(p, 'nz', 61);
    addParameter(p, 'emitterHeightUm', []);
    addParameter(p, 'normalizeDisplay', true);
    addParameter(p, 'logScale', false);
    addParameter(p, 'colormap', '');
    addParameter(p, 'makeFigure', true);
    addParameter(p, 'saveFigures', false);
    addParameter(p, 'outputDir', '');
    addParameter(p, 'verbose', true);
    parse(p, varargin{:});
    opts = p.Results;

    opts.psfType = lower(char(opts.psfType));
    opts.projection = lower(char(opts.projection));
    if ~ismember(opts.psfType, {'effective','detectionsum','detection','excitation'})
        error('reconstruct3DSumPSFFromFocalFit:BadPsfType', ...
            'psfType must be effective, detectionSum, detection, or excitation.');
    end
    if ~ismember(opts.projection, {'sum','max'})
        error('reconstruct3DSumPSFFromFocalFit:BadProjection', ...
            'projection must be sum or max.');
    end
    if mod(opts.nx, 2) == 0, opts.nx = opts.nx + 1; end
    if mod(opts.nz, 2) == 0, opts.nz = opts.nz + 1; end
    if isempty(opts.outputDir)
        opts.outputDir = fullfile(fileparts(mfilename('fullpath')), ...
            'output_matlab', 'focal_fit_3d_psf');
    end
end

% ----------------------------------------------------------------------------
function [coeffs, baseSim] = resolveFitInputs(result, opts)
    if ~isempty(opts.coeffs)
        coeffs = opts.coeffs;
    elseif isstruct(result) && isfield(result, 'fit') && ...
            isfield(result.fit, 'estCoeffs')
        coeffs = result.fit.estCoeffs;
    elseif isstruct(result) && ~isfield(result, 'fit')
        coeffs = result;   % treat result itself as a coeff struct
    else
        error('reconstruct3DSumPSFFromFocalFit:NoCoeffs', ...
            'Provide a fit result with result.fit.estCoeffs or a ''coeffs'' option.');
    end

    if ~isempty(opts.sim)
        baseSim = opts.sim;
    elseif isstruct(result) && isfield(result, 'fit') && isfield(result.fit, 'sim')
        baseSim = result.fit.sim;
    else
        error('reconstruct3DSumPSFFromFocalFit:NoSim', ...
            'Provide result.fit.sim or a ''sim'' option describing the optics.');
    end
end

% ----------------------------------------------------------------------------
function sim = configurePSFGrid(baseSim, opts)
    sim = baseSim;
    sim.x = linspace(-opts.fovXYUm/2, opts.fovXYUm/2, opts.nx);
    sim.y = sim.x;
    sim.nx = opts.nx;
    sim.ny = opts.nx;
    sim.dx = abs(sim.x(2) - sim.x(1));
    sim.z = linspace(-opts.zRangeUm/2, opts.zRangeUm/2, opts.nz);
    sim.nz = opts.nz;
    sim.fovXY = opts.fovXYUm;
    sim.nzRange = opts.zRangeUm;
    if ~isfield(sim, 'includesVectorialPolarization') || ...
            isempty(sim.includesVectorialPolarization)
        sim.includesVectorialPolarization = true;
    end
    if ~isempty(opts.excitationPolarizationMode)
        sim.excitationPolarizationMode = char(opts.excitationPolarizationMode);
    end
    if ~isempty(opts.collectionPolarizationMode)
        sim.collectionPolarizationMode = char(opts.collectionPolarizationMode);
    end
end

% ----------------------------------------------------------------------------
function h = resolveEmitterHeight(sim, opts)
    if ~isempty(opts.emitterHeightUm)
        h = double(opts.emitterHeightUm);
    elseif isfield(sim, 'interfaceEmitterHeightUm') && ...
            ~isempty(sim.interfaceEmitterHeightUm)
        h = double(sim.interfaceEmitterHeightUm);
    else
        h = 0;
    end
end

% ----------------------------------------------------------------------------
function [volume, hExc, hDetSum, hDet] = build3DVolume(sim, coeffs, emitterHeight)
    coeffs = coeffStruct(sim, coeffs);
    zVec = sim.z(:).';

    isInterface = isfield(sim, 'sampleGeometry') && ...
        strcmpi(sim.sampleGeometry, 'airOnGlass');

    if isInterface
        simExc = setVectorialPolarizationForPath(sim, 'excitation');
        simExc.airInterfaceDirection = 'excitation';
        hExc = psfBesselAirInterface(simExc, coeffs, sim.lamExc, zVec, emitterHeight);

        simDet = setVectorialPolarizationForPath(sim, 'collection');
        simDet.airInterfaceDirection = 'collection';
        hDet = psfBesselAirInterface(simDet, coeffs, sim.lamEm, zVec, emitterHeight);
    else
        simExc = setVectorialPolarizationForPath(sim, 'excitation');
        hExc = psfBessel(simExc, coeffs, sim.lamExc);
        simDet = setVectorialPolarizationForPath(sim, 'collection');
        hDet = psfBessel(simDet, coeffs, sim.lamEm);
    end

    eta = detectorCollectionEfficiencyExplicit(sim, hDet);  % [ny nx nz nDet]
    hDetSum = sum(eta, 4);                                   % [ny nx nz]
    volume = hExc .* hDetSum;                                % effective sum PSF
end

% ----------------------------------------------------------------------------
function volume = selectPsfType(psfType, effective, hExc, hDetSum, hDet)
    switch psfType
        case 'effective'
            volume = effective;
        case 'detectionsum'
            volume = hDetSum;
        case 'detection'
            volume = hDet;
        case 'excitation'
            volume = hExc;
    end
end

% ----------------------------------------------------------------------------
function proj = projectVolume(volume, mode)
    switch mode
        case 'sum'
            xy = sum(volume, 3);
            yz = squeeze(sum(volume, 2));
            zx = squeeze(sum(volume, 1)).';   % -> [nz nx]
        case 'max'
            xy = max(volume, [], 3);
            yz = squeeze(max(volume, [], 2));
            zx = squeeze(max(volume, [], 1)).';
    end
    proj = struct('xy', xy, 'yz', yz, 'zx', zx);
end

% ----------------------------------------------------------------------------
function fig = plotProjections(proj, sim, coeffs, opts)
    cmap = resolveColormap(opts.colormap);
    x = sim.x(:).'; y = sim.y(:).'; z = sim.z(:).';

    fig = figure('Color', 'w', 'Position', [80 80 1280 430]);
    tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    % xy (top view): x horizontal, y vertical
    ax1 = nexttile(tl);
    drawProjection(ax1, x, y, proj.xy, opts);
    axis(ax1, 'image'); set(ax1, 'YDir', 'normal');
    xlabel(ax1, 'x [\mum]'); ylabel(ax1, 'y [\mum]');
    title(ax1, 'xy projection');
    colormap(ax1, cmap);

    % yz (side view): z horizontal, y vertical
    ax2 = nexttile(tl);
    drawProjection(ax2, z, y, proj.yz, opts);
    axis(ax2, 'image'); set(ax2, 'YDir', 'normal');
    xlabel(ax2, 'z [\mum]'); ylabel(ax2, 'y [\mum]');
    title(ax2, 'yz projection');
    colormap(ax2, cmap);

    % zx (side view): x horizontal, z vertical
    ax3 = nexttile(tl);
    drawProjection(ax3, x, z, proj.zx, opts);
    axis(ax3, 'image'); set(ax3, 'YDir', 'normal');
    xlabel(ax3, 'x [\mum]'); ylabel(ax3, 'z [\mum]');
    title(ax3, 'zx projection');
    colormap(ax3, cmap);

    cb = colorbar(ax3);
    cb.Layout.Tile = 'east';
    if opts.normalizeDisplay
        cb.Label.String = 'intensity (per-projection max)';
    else
        cb.Label.String = 'summed intensity';
    end

    title(tl, { sprintf('3D %s PSF (%s projection) from focal-plane fit', ...
        prettyPsfType(opts.psfType), opts.projection), ...
        coeffSummaryText(coeffs, sim) }, 'FontWeight', 'bold', 'Interpreter', 'tex');
end

% ----------------------------------------------------------------------------
function drawProjection(ax, xAxis, yAxis, img, opts)
    img = double(img);
    if opts.logScale
        img = log10(max(img, max(img(:)) * 1e-4));
    end
    if opts.normalizeDisplay
        m = max(img(:));
        if m > 0, img = img / m; end
    end
    imagesc(ax, xAxis, yAxis, img);
end

% ----------------------------------------------------------------------------
function cmap = resolveColormap(name)
    if isempty(name)
        candidates = {'inferno', 'magma', 'viridis', 'hot'};
    else
        candidates = {char(name), 'hot'};
    end
    for k = 1:numel(candidates)
        if exist(candidates{k}, 'file') == 2 || any(strcmp(candidates{k}, ...
                {'parula','hot','jet','gray','bone','hsv','cool','turbo'}))
            try
                cmap = feval(candidates{k}, 256);
                return;
            catch
            end
        end
    end
    cmap = hot(256);
end

% ----------------------------------------------------------------------------
function name = prettyPsfType(psfType)
    switch psfType
        case 'effective',    name = 'detector-summed effective';
        case 'detectionsum', name = 'detector-summed detection';
        case 'detection',    name = 'detection';
        case 'excitation',   name = 'excitation';
        otherwise,           name = psfType;
    end
end

% ----------------------------------------------------------------------------
function txt = coeffSummaryText(coeffs, sim)
    coeffs = coeffStruct(sim, coeffs);
    fn = fieldnames(coeffs);
    parts = {};
    for k = 1:numel(fn)
        v = coeffs.(fn{k});
        if isnumeric(v) && isscalar(v) && abs(v) > 1e-4
            parts{end+1} = sprintf('%s %+.3f\\lambda', ...
                strrep(fn{k}, '_', '\_'), v); %#ok<AGROW>
        end
    end
    if isempty(parts)
        txt = 'no significant aberration';
    else
        txt = strjoin(parts, ',  ');
    end
end

% ----------------------------------------------------------------------------
function printSummary(coeffs, sim, opts)
    coeffs = coeffStruct(sim, coeffs);
    fprintf('[reconstruct3DSumPSFFromFocalFit] %s PSF, %s projections, grid %dx%dx%d\n', ...
        prettyPsfType(opts.psfType), opts.projection, opts.nx, opts.nx, opts.nz);
    fn = fieldnames(coeffs);
    for k = 1:numel(fn)
        v = coeffs.(fn{k});
        if isnumeric(v) && isscalar(v) && abs(v) > 1e-4
            fprintf('    %-14s %+8.4f waves RMS (%+7.1f nm)\n', ...
                fn{k}, v, v * sim.lamRef * 1000);
        end
    end
end

% ----------------------------------------------------------------------------
function ensureDir(d)
    if ~isempty(d) && exist(d, 'dir') ~= 7
        mkdir(d);
    end
end
