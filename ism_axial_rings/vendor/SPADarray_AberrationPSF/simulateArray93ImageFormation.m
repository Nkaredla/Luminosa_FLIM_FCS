function out = simulateArray93ImageFormation(varargin)
%SIMULATEARRAY93IMAGEFORMATION Vectorial Richards-Wolf image formation on the 93-pixel array.
%
%   out = simulateArray93ImageFormation()
%   out = simulateArray93ImageFormation('defocusRad', 1)
%
%   Forms the on-axis point-emitter detection PSF with the vectorial
%   Richards-Wolf/Bessel model on the 93-channel honeycomb array
%   (detectorLayout('honeycomb93', ...)), then integrates that PSF over
%   each finite hex pixel with detectorCollectionEfficiencyExplicit -- the
%   same forward-model path showDetectorIntensityVsAberration.m uses for
%   the 23-array. Two cases are produced: in focus (no aberration) and with
%   a defocus wavefront term of 'defocusRad' radians.
%
%   EMITTER GEOMETRY (default)
%   By default the emitter sits on the coverslip glass with air on the
%   sample side, imaged by a high-NA oil-immersion objective through the
%   glass (immersion -> glass -> air). This is the layered vectorial
%   air-interface model (vectorialPSFBesselAirInterface, collection
%   direction, emitter height 0). Because the propagating-wave model needs
%   NA < nImmersion, a 1.45 NA objective uses oil (nImmersion 1.51), not
%   air. Pass 'sampleGeometry','homogeneous' for a single-medium test.
%
%   ARRAY SIZE / MAGNIFICATION
%   By default the sample-equivalent pixel pitch is chosen so the array's
%   net edge-to-edge diameter is 'arrayDiameterAU' Airy units (default
%   1.7 AU), i.e. the magnification is set to shrink the 93-pixel footprint
%   to that size (see detectorPitchForDiameterAU). Higher magnification ->
%   smaller pitch -> smaller AU footprint. Pass an explicit 'detPitchUm' to
%   override and pin the pitch directly.
%
%   DEFOCUS CONVENTION
%   zernikeModes.m's Z.defocus = sqrt(3)*(2*rho^2-1) has unit RMS over the
%   pupil, so coeffs.defocus is directly the RMS wavefront error in waves
%   at sim.lamRef (see vectorialPSFBessel.m, zernikePhaseMap.m: 1 wave =
%   2*pi rad). 'defocusRad' here is that same RMS wavefront error expressed
%   in radians, so internally coeffs.defocus = defocusRad / (2*pi).
%
%   NAME-VALUE OPTIONS
%     'NA'                 objective NA                       (default 1.45)
%     'lamEm'               emission wavelength [um]          (default 0.520)
%     'sampleGeometry'      'airOnGlass' (emitter on glass)    (default
%                           or 'homogeneous'                   'airOnGlass')
%     'emitterHeightUm'     emitter height above glass [um]    (default 0)
%     'nImmersion'          immersion index (oil)              (default 1.51)
%     'nGlass'              coverslip index                    (default 1.518)
%     'nSample'             sample-side index (air)            (default 1.0)
%     'coverslipThicknessUm' coverslip thickness [um]          (default 170)
%     'nMedium'             index for the homogeneous fallback (default 1.33)
%     'arrayDiameterAU'     target net array diameter [AU]     (default 1.7)
%                           (ignored if 'detPitchUm' is given)
%     'detPitchUm'          detector pitch in sample-equiv.    (default [] ->
%                           units; overrides arrayDiameterAU    derived from
%                           when supplied                       arrayDiameterAU)
%     'detFillRatio'        detector size / pitch              (default 1.0)
%     'detectorSubsamples'  quadrature points per pixel        (default 7)
%     'defocusRad'          RMS wavefront defocus [rad]        (default 1)
%     'verbose'             print pitch/diameter summary       (default true)
%     'fovXY'               lateral PSF field of view [um]     (default 2.6)
%     'nx'                  lateral samples (odd)              (default 131)
%     'colormap'            colormap name                      (default 'parula')
%     'showContinuousPSF'   also draw the detector-plane PSF   (default true)
%     'saveFigures'         export PNGs                        (default false)
%     'outputDir'           where to save                      (default
%                           output_matlab/array93_image_formation)
%
%   OUTPUT
%     out.sim              simulation struct used
%     out.detXY             [93 x 2] detector centers
%     out.psfFocus          in-focus detection PSF, [ny nx]
%     out.psfDefocus        defocused detection PSF, [ny nx]
%     out.detValsFocus      [93 x 1] collected intensity, in focus
%     out.detValsDefocus    [93 x 1] collected intensity, defocused
%     out.defocusRad        defocus amplitude used [rad]
%     out.defocusWaves      same amplitude in waves (coeffs.defocus)
%     out.detPitchUm        sample-equivalent pixel pitch used [um]
%     out.arrayDiameterAU   target net array diameter [AU]
%     out.netDiameterAU     achieved net edge-to-edge diameter [AU]
%     out.netDiameterUm     achieved net edge-to-edge diameter [um]
%     out.airyUnitUm        1 AU in sample space [um]
%     out.figHexMaps        figure handle, hex intensity maps
%     out.figContinuousPSF  figure handle, continuous detector-plane PSF (or [])

    addpath(fileparts(mfilename('fullpath')));
    opts = parseInputs(varargin{:});

    sim = defaultParams();
    sim.NA      = opts.NA;
    sim.lamEm   = opts.lamEm;
    sim.lamRef  = opts.lamEm;
    sim.includesVectorialPolarization = true;
    sim.diffractionModel = 'vectorial Richards-Wolf/Bessel';
    sim.vectorialPolarizationMode  = 'xyAverage';
    sim.collectionPolarizationMode = 'xyAverage';
    sim.detectorSubsamples = opts.detectorSubsamples;
    sim.detectorImageInverted = true;
    sim.detectorPixelShape = 'hex';

    sim.sampleGeometry = normalizeSampleGeometry(opts.sampleGeometry);
    if usesInterfaceGeometry(sim)
        % Point emitter sitting on the coverslip glass with air on the
        % sample side, viewed by a high-NA oil-immersion objective through
        % the glass (immersion -> glass -> air). The propagating-wave model
        % needs NA < nImmersion, so 1.45 NA pairs with oil (~1.51), not air.
        sim.nImmersion   = opts.nImmersion;
        sim.nGlass       = opts.nGlass;
        sim.nSample      = opts.nSample;
        sim.nDesignGlass = opts.nGlass;
        sim.coverslipThicknessUm       = opts.coverslipThicknessUm;
        sim.designCoverslipThicknessUm = opts.coverslipThicknessUm;
        sim.nMedium = sim.nImmersion;
        sim.interfaceEmitterHeightUm = opts.emitterHeightUm;
    else
        sim.nMedium = opts.nMedium;
    end

    sim.fovXY = opts.fovXY;
    sim.nx = opts.nx;  sim.ny = opts.nx;
    sim.x = linspace(-sim.fovXY/2, sim.fovXY/2, sim.nx);
    sim.y = sim.x;
    sim.dx = abs(sim.x(2) - sim.x(1));

    sim.detectorLayout = 'honeycomb93';
    sim.detPitch = opts.detPitchUm;
    [sim.detXY, sim.detectorIndexGrid, sim.detectorLayoutInfo] = ...
        detectorLayout('honeycomb93', sim.detPitch);
    sim.nDet = size(sim.detXY, 1);
    sim.detSize = opts.detFillRatio * sim.detPitch;
    sim.detectorHexRadius = sim.detSize / sqrt(3);

    sim.z = 0;   % in-focus stage; defocus is applied as a pupil phase term
    sim.nz = 1;

    [~, xc] = min(abs(sim.x));
    [~, yc] = min(abs(sim.y));

    defocusWaves = opts.defocusRad / (2*pi);

    % Detection PSF for the emitter on the glass: in focus and with a
    % vectorial-defocus wavefront term. The stage stays in focus (stageZ 0);
    % the "defocus" is the RMS-wavefront Zernike term.
    psfFocus   = detectionPSF(sim, struct());
    psfDefocus = detectionPSF(sim, struct('defocus', defocusWaves));

    etaFocus   = detectorCollectionEfficiencyExplicit(sim, psfFocus);
    etaDefocus = detectorCollectionEfficiencyExplicit(sim, psfDefocus);

    detValsFocus   = squeeze(etaFocus(yc, xc, 1, :));
    detValsDefocus = squeeze(etaDefocus(yc, xc, 1, :));

    % Achieved net (edge-to-edge) array diameter, for reporting/plot title.
    airyUnitUm = 1.22 * sim.lamEm / sim.NA;
    spanX = (max(sim.detXY(:,1)) - min(sim.detXY(:,1))) + sim.detSize;
    spanY = (max(sim.detXY(:,2)) - min(sim.detXY(:,2))) + 2*sim.detectorHexRadius;
    netDiameterUm = max(spanX, spanY);
    netDiameterAU = netDiameterUm / airyUnitUm;

    if opts.verbose
        fprintf(['[simulateArray93ImageFormation] pitch %.1f nm -> net diameter ' ...
            '%.3f um = %.2f AU (1 AU = %.1f nm, NA %.2f, lam_em %.0f nm).\n'], ...
            1000*sim.detPitch, netDiameterUm, netDiameterAU, ...
            1000*airyUnitUm, sim.NA, 1000*sim.lamEm);
    end

    figHexMaps = plotHexMapsFigure(sim, detValsFocus, detValsDefocus, opts, netDiameterAU);
    figContinuousPSF = [];
    if opts.showContinuousPSF
        figContinuousPSF = plotContinuousPSFFigure(sim, ...
            psfFocus(:,:,1), psfDefocus(:,:,1), opts);
    end

    if opts.saveFigures
        ensureDir(opts.outputDir);
        exportgraphics(figHexMaps, fullfile(opts.outputDir, ...
            'array93_hex_intensity.png'), 'Resolution', 180);
        if ~isempty(figContinuousPSF)
            exportgraphics(figContinuousPSF, fullfile(opts.outputDir, ...
                'array93_detector_plane_psf.png'), 'Resolution', 180);
        end
    end

    out = struct();
    out.sim               = sim;
    out.detXY             = sim.detXY;
    out.psfFocus           = psfFocus(:,:,1);
    out.psfDefocus         = psfDefocus(:,:,1);
    out.detValsFocus       = detValsFocus;
    out.detValsDefocus     = detValsDefocus;
    out.defocusRad         = opts.defocusRad;
    out.defocusWaves       = defocusWaves;
    out.detPitchUm          = sim.detPitch;
    out.arrayDiameterAU     = opts.arrayDiameterAU;
    out.netDiameterAU       = netDiameterAU;
    out.netDiameterUm       = netDiameterUm;
    out.airyUnitUm          = airyUnitUm;
    out.figHexMaps          = figHexMaps;
    out.figContinuousPSF    = figContinuousPSF;
    out.options             = opts;
end

% ============================================================================
function opts = parseInputs(varargin)
    p = inputParser;
    p.FunctionName = 'simulateArray93ImageFormation';
    addParameter(p, 'NA', 1.45);
    addParameter(p, 'lamEm', 0.520);
    addParameter(p, 'sampleGeometry', 'airOnGlass');
    addParameter(p, 'emitterHeightUm', 0);        % 0 = emitter on the glass
    addParameter(p, 'nImmersion', 1.51);          % oil immersion
    addParameter(p, 'nGlass', 1.518);             % coverslip
    addParameter(p, 'nSample', 1.0);              % air above the glass
    addParameter(p, 'coverslipThicknessUm', 170); % matched -> no cover aberration
    addParameter(p, 'nMedium', 1.33);             % homogeneous fallback only
    addParameter(p, 'arrayDiameterAU', 1.7);
    addParameter(p, 'detPitchUm', []);
    addParameter(p, 'detFillRatio', 1.0);
    addParameter(p, 'detectorSubsamples', 7);
    addParameter(p, 'defocusRad', 1);
    addParameter(p, 'verbose', true);
    addParameter(p, 'fovXY', 2.6);
    addParameter(p, 'nx', 131);
    addParameter(p, 'colormap', 'parula');
    addParameter(p, 'showContinuousPSF', true);
    addParameter(p, 'saveFigures', false);
    addParameter(p, 'outputDir', '');
    parse(p, varargin{:});
    opts = p.Results;

    if mod(opts.nx, 2) == 0, opts.nx = opts.nx + 1; end
    if isempty(opts.detPitchUm)
        % Set the magnification so the 93-array footprint spans the
        % requested Airy-unit diameter.
        opts.detPitchUm = detectorPitchForDiameterAU('honeycomb93', ...
            opts.arrayDiameterAU, opts.NA, opts.lamEm, opts.detFillRatio);
    end
    if isempty(opts.outputDir)
        opts.outputDir = fullfile(fileparts(mfilename('fullpath')), ...
            'output_matlab', 'array93_image_formation');
    end
end

% ----------------------------------------------------------------------------
function h = detectionPSF(sim, coeffs)
%DETECTIONPSF Collection-path detection PSF for the on-axis emitter.
    simDet = setVectorialPolarizationForPath(sim, 'collection');
    if usesInterfaceGeometry(sim)
        simDet.airInterfaceDirection = 'collection';
        h = psfBesselAirInterface(simDet, coeffs, sim.lamEm, ...
            0, sim.interfaceEmitterHeightUm);   % stageZ 0, emitter on glass
    else
        h = psfBessel(simDet, coeffs, sim.lamEm);
    end
end

% ----------------------------------------------------------------------------
function tf = usesInterfaceGeometry(sim)
    tf = isfield(sim, 'sampleGeometry') && ...
        strcmpi(sim.sampleGeometry, 'airOnGlass');
end

% ----------------------------------------------------------------------------
function geometry = normalizeSampleGeometry(value)
    key = lower(regexprep(char(value), '[\s_\-]', ''));
    switch key
        case 'homogeneous'
            geometry = 'homogeneous';
        case {'aironglass','interface','layeredinterface','glassinterface', ...
                'emitteronglass','onglass','glasswaterinterface','wateronglass'}
            geometry = 'airOnGlass';
        otherwise
            error('simulateArray93ImageFormation:BadSampleGeometry', ...
                ['sampleGeometry must be ''airOnGlass'' (emitter on the ' ...
                'coverslip) or ''homogeneous''.']);
    end
end

% ----------------------------------------------------------------------------
function fig = plotHexMapsFigure(sim, detValsFocus, detValsDefocus, opts, netDiameterAU)
    fig = figure('Color', 'w', 'Position', [100 100 760 420]);
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    sharedCLim = [0, max([detValsFocus(:); detValsDefocus(:)])];
    if sharedCLim(2) <= 0, sharedCLim(2) = 1; end

    ax1 = nexttile(tl);
    plotDetectorHexMap(sim.detXY, detValsFocus, 'Parent', ax1, ...
        'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.3);
    colormap(ax1, opts.colormap);
    caxis(ax1, sharedCLim);
    title(ax1, 'In focus');

    ax2 = nexttile(tl);
    plotDetectorHexMap(sim.detXY, detValsDefocus, 'Parent', ax2, ...
        'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.3);
    colormap(ax2, opts.colormap);
    caxis(ax2, sharedCLim);
    title(ax2, sprintf('Defocus %.2f rad (%.3f\\lambda)', ...
        opts.defocusRad, opts.defocusRad/(2*pi)));

    cb = colorbar(ax2);
    cb.Layout.Tile = 'east';
    cb.Label.String = 'collected intensity (relative)';

    if usesInterfaceGeometry(sim)
        geomLabel = sprintf('emitter on glass, oil n%.2f/air n%.2f', ...
            sim.nImmersion, sim.nSample);
    else
        geomLabel = sprintf('homogeneous n%.3f', sim.nMedium);
    end
    title(tl, sprintf(['93-pixel array detection intensity, on-axis emitter ' ...
        '(%s, NA %.2f, \\lambda_{em} %.0f nm, pitch %.0f nm, \\O %.2f AU)'], ...
        geomLabel, sim.NA, 1000*sim.lamEm, 1000*sim.detPitch, netDiameterAU), ...
        'FontWeight', 'bold');
end

% ----------------------------------------------------------------------------
function fig = plotContinuousPSFFigure(sim, hFocus, hDefocus, opts)
    fig = figure('Color', 'w', 'Position', [100 560 760 420]);
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    cellR = sim.detSize / sqrt(3);
    pmax = max([hFocus(:); hDefocus(:)]);

    ax1 = nexttile(tl);
    imagesc(ax1, sim.x, sim.y, hFocus, [0 pmax]);
    axis(ax1, 'image'); set(ax1, 'YDir', 'normal');
    colormap(ax1, opts.colormap);
    hold(ax1, 'on'); drawHexOutlines(ax1, sim.detXY, cellR); hold(ax1, 'off');
    set(ax1, 'XTick', [], 'YTick', []);
    title(ax1, 'In focus');

    ax2 = nexttile(tl);
    imagesc(ax2, sim.x, sim.y, hDefocus, [0 pmax]);
    axis(ax2, 'image'); set(ax2, 'YDir', 'normal');
    colormap(ax2, opts.colormap);
    hold(ax2, 'on'); drawHexOutlines(ax2, sim.detXY, cellR); hold(ax2, 'off');
    set(ax2, 'XTick', [], 'YTick', []);
    title(ax2, sprintf('Defocus %.2f rad', opts.defocusRad));

    title(tl, 'Detection PSF at the detector plane (hex = SPAD pixels)', ...
        'FontWeight', 'bold');
end

% ----------------------------------------------------------------------------
function drawHexOutlines(ax, detXY, cellRadius)
    theta = (0:5)*pi/3 + pi/6;
    hx = cellRadius*cos(theta);
    hy = cellRadius*sin(theta);
    for k = 1:size(detXY,1)
        plot(ax, detXY(k,1)+[hx hx(1)], detXY(k,2)+[hy hy(1)], ...
            '-', 'Color', [1 1 1 0.5], 'LineWidth', 0.4);
    end
end

% ----------------------------------------------------------------------------
function ensureDir(d)
    if ~isempty(d) && exist(d, 'dir') ~= 7
        mkdir(d);
    end
end
