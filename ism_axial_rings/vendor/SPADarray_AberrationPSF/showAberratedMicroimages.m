function out = showAberratedMicroimages(varargin)
%SHOWABERRATEDMICROIMAGES Aberrated PSF, array microimage and fitted centre vs emitter shift.
%
%   out = showAberratedMicroimages()
%   out = showAberratedMicroimages('aberrationMode','coma_x','amplitudeWaves',0.15, ...
%             'defocusUm',0.5,'shiftsUm',[-0.1 -0.05 0 0.05 0.1])
%
%   Companion diagnostic to simulateAberrationLocalizationBias. That routine
%   reports the bias as a number; this one shows the picture behind it. For each
%   requested emitter position it renders
%
%     row 1  the aberrated detection PSF in the sample-equivalent detector
%            plane, with the honeycomb pixel outlines drawn on top, so you can
%            see how much of the microimage the array actually catches;
%     row 2  the per-channel probabilities p_k on the honeycomb, i.e. the only
%            thing the fitter ever sees;
%     row 3  the cross-entropy landscape sum_k p_k^true log p_k^ideal(s) over
%            trial positions s, which is what the fit maximizes.
%
%   Each panel marks the true emitter position (white circle) and the fitted
%   position (red cross), so the bias is the gap between them. The landscape
%   panel is the useful one when a bias looks wrong: if it shows two comparable
%   maxima, the reported bias depends on which basin the solver landed in, and
%   the apparent-diffusion ratio derived from it is not trustworthy.
%
%   The fit uses exactly the model of simulateAberrationLocalizationBias: the
%   amplitude is free, so only the normalized channel-allocation vector matters
%   and the result is photon-count independent. No counts are sampled.
%
%   KEY NAME-VALUE OPTIONS
%     'aberrationMode'   Zernike mode name, or 'none' (default 'coma_x')
%     'amplitudeWaves'   RMS waves at sim.lamRef (default 0.15)
%     'defocusUm'        stage plane of the TRUE model (default 0.5)
%     'shiftsUm'         emitter x positions to render (default [-0.1 0 0.1])
%     'shiftAxis'        'x', 'y' or 'diagonal' (default 'x')
%     'landscapeHalfWidthUm'  extent of the cross-entropy map (default 0.25)
%     Optics and detector options match simulateAberrationLocalizationBias:
%       'detectorLayout', 'NA', 'lamEm', 'sampleGeometry', 'emitterHeightUm',
%       'nImmersion', 'nGlass', 'nSample', 'coverslipThicknessUm', 'nMedium',
%       'arrayDiameterAU', 'detPitchUm', 'detFillRatio', 'detectorSubsamples',
%       'fovXY', 'nx'
%     Output: 'makeFigure' (true), 'saveFigures' (false), 'outputDir' (''),
%       'verbose' (true)
%
%   OUTPUT
%     out.shiftsUm        emitter positions rendered [nShift x 2]
%     out.fittedUm        fitted positions [nShift x 2]
%     out.biasUm          fittedUm - shiftsUm
%     out.probability     [nShift x nDet] true channel probabilities
%     out.idealProbability[nShift x nDet] ideal probabilities at the fit
%     out.landscape       cross-entropy maps and their axes
%     out.multimodal      true where a competing maximum is within tolerance
%     out.psf             aberrated detection PSF slices
%     out.sim, out.figure

    addpath(fileparts(mfilename('fullpath')));
    opts = parseLocalInputs(varargin{:});
    sim = configureLocalSim(opts);

    coeffs = struct();
    if ~strcmpi(opts.aberrationMode, 'none') && opts.amplitudeWaves ~= 0
        coeffs.(opts.aberrationMode) = opts.amplitudeWaves;
    end

    trueModel  = localChannelModel(sim, coeffs, opts.defocusUm, opts);
    idealModel = localChannelModel(sim, struct(), 0, opts);

    shifts = shiftPositions(opts);
    nShift = size(shifts, 1);

    probability = zeros(nShift, sim.nDet);
    idealProb   = zeros(nShift, sim.nDet);
    fitted      = zeros(nShift, 2);
    multimodal  = false(nShift, 1);
    landscape   = cell(nShift, 1);

    for s = 1:nShift
        probability(s,:) = channelProbabilityAt(trueModel, shifts(s,1), shifts(s,2));
        [land, axesUm] = crossEntropyMap(probability(s,:), idealModel, shifts(s,:), opts);
        landscape{s} = struct('value', land, 'x', axesUm.x, 'y', axesUm.y);

        [fitted(s,:), multimodal(s)] = fitFromLandscape(land, axesUm, ...
            probability(s,:), idealModel, opts);
        idealProb(s,:) = channelProbabilityAt(idealModel, fitted(s,1), fitted(s,2));

        if opts.verbose
            fprintf(['[microimage] r = (%+6.3f, %+6.3f) um -> fit (%+6.3f, %+6.3f) um' ...
                ' | bias %6.1f nm%s\n'], shifts(s,1), shifts(s,2), ...
                fitted(s,1), fitted(s,2), ...
                1000*hypot(fitted(s,1)-shifts(s,1), fitted(s,2)-shifts(s,2)), ...
                repmat('  [MULTIMODAL]', 1, double(multimodal(s))));
        end
    end

    fig = [];
    if opts.makeFigure
        fig = renderPanels(sim, opts, trueModel, shifts, fitted, ...
            probability, landscape, multimodal);
        if opts.saveFigures
            if ~isempty(opts.outputDir) && exist(opts.outputDir, 'dir') ~= 7
                mkdir(opts.outputDir);
            end
            exportgraphics(fig, fullfile(opts.outputDir, sprintf( ...
                '%s_microimages_%s_%0.2fwaves_z%0.2f.png', sim.detectorLayout, ...
                opts.aberrationMode, opts.amplitudeWaves, opts.defocusUm)), ...
                'Resolution', 180);
        end
    end

    out = struct();
    out.sim = sim;
    out.coeffs = coeffs;
    out.shiftsUm = shifts;
    out.fittedUm = fitted;
    out.biasUm = fitted - shifts;
    out.probability = probability;
    out.idealProbability = idealProb;
    out.landscape = landscape;
    out.multimodal = multimodal;
    out.psf = trueModel.psf;
    out.figure = fig;
    out.options = opts;
end

% ============================================================================
function opts = parseLocalInputs(varargin)
    p = inputParser;
    p.FunctionName = 'showAberratedMicroimages';
    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'NA', 1.45);
    addParameter(p, 'lamEm', 0.520);
    addParameter(p, 'sampleGeometry', 'airOnGlass');
    addParameter(p, 'emitterHeightUm', 0);
    addParameter(p, 'nImmersion', 1.51);
    addParameter(p, 'nGlass', 1.518);
    addParameter(p, 'nSample', 1.33);
    addParameter(p, 'coverslipThicknessUm', 170);
    addParameter(p, 'nMedium', 1.33);
    addParameter(p, 'arrayDiameterAU', 1.7);
    addParameter(p, 'detPitchUm', []);
    addParameter(p, 'detFillRatio', 1.0);
    addParameter(p, 'detectorSubsamples', 7);
    addParameter(p, 'fovXY', 2.6);
    addParameter(p, 'nx', 121);
    addParameter(p, 'maxZernikeOrder', 4);

    addParameter(p, 'aberrationMode', 'coma_x');
    addParameter(p, 'amplitudeWaves', 0.15);
    addParameter(p, 'defocusUm', 0.5);
    addParameter(p, 'shiftsUm', [-0.10 0 0.10]);
    addParameter(p, 'shiftAxis', 'x');
    addParameter(p, 'landscapeHalfWidthUm', 0.25);
    addParameter(p, 'landscapeSamples', 81);
    addParameter(p, 'multimodalToleranceFraction', 0.25);
    addParameter(p, 'minProbability', 1e-12);

    addParameter(p, 'makeFigure', true);
    addParameter(p, 'saveFigures', false);
    addParameter(p, 'outputDir', '');
    addParameter(p, 'verbose', true);
    parse(p, varargin{:});
    opts = p.Results;

    opts.detectorLayout = char(opts.detectorLayout);
    opts.aberrationMode = char(opts.aberrationMode);
    opts.shiftAxis = lower(char(opts.shiftAxis));
    if mod(opts.nx, 2) == 0, opts.nx = opts.nx + 1; end
    if isempty(opts.detPitchUm)
        opts.detPitchUm = detectorPitchForDiameterAU(opts.detectorLayout, ...
            opts.arrayDiameterAU, opts.NA, opts.lamEm, opts.detFillRatio);
    end
    if isempty(opts.outputDir)
        opts.outputDir = fullfile(fileparts(mfilename('fullpath')), ...
            'output_matlab', 'aberration_localization_bias');
    end
end

% ----------------------------------------------------------------------------
function sim = configureLocalSim(opts)
%CONFIGURELOCALSIM Same optical configuration as simulateAberrationLocalizationBias.
    sim = defaultParams();
    sim.NA     = opts.NA;
    sim.lamEm  = opts.lamEm;
    sim.lamRef = opts.lamEm;
    sim.includesVectorialPolarization = true;
    sim.diffractionModel = 'vectorial Richards-Wolf/Bessel';
    sim.vectorialPolarizationMode  = 'xyAverage';
    sim.collectionPolarizationMode = 'xyAverage';
    sim.detectorSubsamples = opts.detectorSubsamples;
    sim.detectorImageInverted = true;
    sim.detectorPixelShape = 'hex';
    sim.M = max(sim.M, opts.maxZernikeOrder + 1);
    sim.modeOrder = zernikeModeNames(max(6, opts.maxZernikeOrder));

    key = lower(regexprep(char(opts.sampleGeometry), '[\s_\-]', ''));
    if strcmp(key, 'homogeneous')
        sim.sampleGeometry = 'homogeneous';
        sim.nMedium = opts.nMedium;
    else
        sim.sampleGeometry = 'airOnGlass';
        sim.nImmersion   = opts.nImmersion;
        sim.nGlass       = opts.nGlass;
        sim.nSample      = opts.nSample;
        sim.nDesignGlass = opts.nGlass;
        sim.coverslipThicknessUm       = opts.coverslipThicknessUm;
        sim.designCoverslipThicknessUm = opts.coverslipThicknessUm;
        sim.nMedium = sim.nImmersion;
        sim.interfaceEmitterHeightUm = opts.emitterHeightUm;
    end

    sim.fovXY = opts.fovXY;
    sim.nx = opts.nx;  sim.ny = opts.nx;
    sim.x = linspace(-sim.fovXY/2, sim.fovXY/2, sim.nx);
    sim.y = sim.x;
    sim.dx = abs(sim.x(2) - sim.x(1));

    sim.detectorLayout = opts.detectorLayout;
    sim.detPitch = opts.detPitchUm;
    [sim.detXY, sim.detectorIndexGrid, sim.detectorLayoutInfo] = ...
        detectorLayout(sim.detectorLayout, sim.detPitch);
    sim.nDet = size(sim.detXY, 1);
    sim.detSize = opts.detFillRatio * sim.detPitch;
    sim.detectorHexRadius = sim.detSize / sqrt(3);
end

% ----------------------------------------------------------------------------
function model = localChannelModel(sim, coeffs, planeZ, opts)
%LOCALCHANNELMODEL Normalized microimage and the PSF slice behind it.
    coeffs = coeffStruct(sim, coeffs);
    simDet = setVectorialPolarizationForPath(sim, 'collection');
    if strcmpi(sim.sampleGeometry, 'airOnGlass')
        simDet.airInterfaceDirection = 'collection';
        psf = psfBesselAirInterface(simDet, coeffs, sim.lamEm, ...
            planeZ, sim.interfaceEmitterHeightUm);
    else
        simDet.z = planeZ;
        simDet.nz = 1;
        psf = psfBessel(simDet, coeffs, sim.lamEm);
    end
    psf = psf(:,:,1);

    eta = detectorCollectionEfficiencyExplicit(sim, psf);
    eta = max(reshape(eta, numel(sim.y), numel(sim.x), sim.nDet), 0);
    total = sum(eta, 3);
    total(total < opts.minProbability) = opts.minProbability;

    model = struct();
    model.probability = eta ./ total;
    model.psf = psf;
    model.x = sim.x;
    model.y = sim.y;
    model.nDet = sim.nDet;
    model.planeZ = planeZ;
    model.coeffs = coeffs;
end

% ----------------------------------------------------------------------------
function pk = channelProbabilityAt(model, qx, qy)
    pk = zeros(1, model.nDet);
    for k = 1:model.nDet
        pk(k) = interp2(model.x, model.y, model.probability(:,:,k), ...
            qx, qy, 'cubic', 0);
    end
    pk = max(pk, 0);
    pk = pk / max(sum(pk), eps);
end

% ----------------------------------------------------------------------------
function shifts = shiftPositions(opts)
    v = opts.shiftsUm(:);
    switch opts.shiftAxis
        case 'x',        shifts = [v, zeros(size(v))];
        case 'y',        shifts = [zeros(size(v)), v];
        case 'diagonal', shifts = [v, v] / sqrt(2);
        otherwise
            error('showAberratedMicroimages:BadShiftAxis', ...
                'shiftAxis must be x, y or diagonal.');
    end
end

% ----------------------------------------------------------------------------
function [land, axesUm] = crossEntropyMap(pTrue, idealModel, centre, opts)
%CROSSENTROPYMAP The objective the fit maximizes, on a dense trial-position grid.
%   Evaluated exhaustively rather than by Newton descent, so a second maximum
%   cannot be missed.
    n = opts.landscapeSamples;
    h = opts.landscapeHalfWidthUm;
    axesUm.x = linspace(centre(1) - h, centre(1) + h, n);
    axesUm.y = linspace(centre(2) - h, centre(2) + h, n);
    [QX, QY] = meshgrid(axesUm.x, axesUm.y);

    land = zeros(n, n);
    for k = 1:idealModel.nDet
        pk = interp2(idealModel.x, idealModel.y, idealModel.probability(:,:,k), ...
            QX, QY, 'cubic', 0);
        land = land + pTrue(k) * log(max(pk, opts.minProbability));
    end
end

% ----------------------------------------------------------------------------
function [best, multimodal] = fitFromLandscape(land, axesUm, pTrue, idealModel, opts)
%FITFROMLANDSCAPE Global maximum of the landscape, refined by local Newton steps.
    [~, linear] = max(land(:));
    [iy, ix] = ind2sub(size(land), linear);
    best = [axesUm.x(ix), axesUm.y(iy)];
    best = polishMaximum(best, pTrue, idealModel, opts);

    % A competing maximum is a peak that is separated from the global one and
    % nearly as high. If one exists the bias is branch dependent.
    peak = max(land(:));
    floorValue = min(land(:));
    threshold = peak - opts.multimodalToleranceFraction * (peak - floorValue);
    high = land >= threshold;
    exclude = false(size(land));
    radius = max(2, round(size(land,1) / 12));
    rows = max(1, iy-radius):min(size(land,1), iy+radius);
    cols = max(1, ix-radius):min(size(land,2), ix+radius);
    exclude(rows, cols) = true;
    multimodal = any(high(:) & ~exclude(:));
end

% ----------------------------------------------------------------------------
function pos = polishMaximum(pos, pTrue, idealModel, opts)
%POLISHMAXIMUM Few Gauss-Newton steps with a consistently differenced objective.
    delta = 0.25 * (idealModel.x(2) - idealModel.x(1));
    for iter = 1:40
        p  = channelProbabilityRaw(idealModel, pos(1), pos(2), opts);
        px = channelProbabilityRaw(idealModel, pos(1)+delta, pos(2), opts);
        mx = channelProbabilityRaw(idealModel, pos(1)-delta, pos(2), opts);
        py = channelProbabilityRaw(idealModel, pos(1), pos(2)+delta, opts);
        my = channelProbabilityRaw(idealModel, pos(1), pos(2)-delta, opts);

        dpx = (px - mx) / (2*delta);
        dpy = (py - my) / (2*delta);
        sx = dpx ./ p;  sy = dpy ./ p;

        g = [sum(pTrue .* sx), sum(pTrue .* sy)];
        H = [sum(pTrue .* sx.^2), sum(pTrue .* sx .* sy); ...
             sum(pTrue .* sx .* sy), sum(pTrue .* sy.^2)];
        if rcond(H) < 1e-12, break; end
        step = (H \ g(:)).';
        if ~all(isfinite(step)), break; end
        pos = pos + step;
        if hypot(step(1), step(2)) < 1e-7, break; end
    end
end

function p = channelProbabilityRaw(model, qx, qy, opts)
    p = zeros(1, model.nDet);
    for k = 1:model.nDet
        p(k) = interp2(model.x, model.y, model.probability(:,:,k), ...
            qx, qy, 'cubic', 0);
    end
    p = max(p, opts.minProbability);
end

% ----------------------------------------------------------------------------
function fig = renderPanels(sim, opts, trueModel, shifts, fitted, ...
        probability, landscape, multimodal)
    nShift = size(shifts, 1);
    fig = figure('Color', 'w', 'Position', [50 50 min(320*nShift+120, 1800) 860]);
    tl = tiledlayout(fig, 3, nShift, 'Padding', 'compact', 'TileSpacing', 'compact');

    span = max(abs(sim.detXY(:))) + sim.detSize;
    viewLimit = 1.35 * span;

    % Row 1: PSF in the detector plane with the array footprint drawn on it.
    for s = 1:nShift
        ax = nexttile(tl, s);
        % detectorCollectionEfficiencyExplicit samples the centred PSF at
        % (d_k + r), so in detector coordinates d the microimage is h(d + r):
        % an emitter at +r puts its image at -r. Plotting the PSF against
        % (sim.x - r) reproduces that, and both markers are negated to match.
        imagesc(ax, sim.x - shifts(s,1), sim.y - shifts(s,2), trueModel.psf);
        set(ax, 'YDir', 'normal'); axis(ax, 'image');
        hold(ax, 'on');
        drawHoneycomb(ax, sim, [0.9 0.9 0.9]);
        plot(ax, -shifts(s,1), -shifts(s,2), 'wo', 'MarkerSize', 8, 'LineWidth', 1.4);
        plot(ax, -fitted(s,1), -fitted(s,2), ...
            'rx', 'MarkerSize', 10, 'LineWidth', 1.6);
        xlim(ax, [-viewLimit viewLimit]); ylim(ax, [-viewLimit viewLimit]);
        xlabel(ax, 'u [\mum]');
        if s == 1, ylabel(ax, 'PSF on array   v [\mum]'); end
        title(ax, sprintf('r = (%+.3f, %+.3f) \\mum', shifts(s,1), shifts(s,2)), ...
            'FontWeight', 'normal', 'FontSize', 9);
    end

    % Row 2: what the fitter actually sees, the 23 channel probabilities.
    for s = 1:nShift
        ax = nexttile(tl, nShift + s);
        drawHoneycombFilled(ax, sim, probability(s,:));
        hold(ax, 'on');
        plot(ax, -shifts(s,1), -shifts(s,2), 'wo', 'MarkerSize', 8, 'LineWidth', 1.4);
        plot(ax, -fitted(s,1), -fitted(s,2), 'rx', 'MarkerSize', 10, 'LineWidth', 1.6);
        axis(ax, 'image');
        xlim(ax, [-viewLimit viewLimit]); ylim(ax, [-viewLimit viewLimit]);
        colorbar(ax);
        xlabel(ax, 'u [\mum]');
        if s == 1, ylabel(ax, 'channel p_k   v [\mum]'); end
        title(ax, sprintf('max p_k = %.3f', max(probability(s,:))), ...
            'FontWeight', 'normal', 'FontSize', 9);
    end

    % Row 3: the objective, with true and fitted positions on it.
    for s = 1:nShift
        ax = nexttile(tl, 2*nShift + s);
        L = landscape{s};
        imagesc(ax, 1000*L.x, 1000*L.y, L.value);
        set(ax, 'YDir', 'normal'); axis(ax, 'image');
        hold(ax, 'on');
        contour(ax, 1000*L.x, 1000*L.y, L.value, 12, 'w', 'LineWidth', 0.4);
        plot(ax, 1000*shifts(s,1), 1000*shifts(s,2), 'wo', ...
            'MarkerSize', 8, 'LineWidth', 1.4);
        plot(ax, 1000*fitted(s,1), 1000*fitted(s,2), 'rx', ...
            'MarkerSize', 10, 'LineWidth', 1.6);
        xlabel(ax, 'trial s_x [nm]');
        if s == 1, ylabel(ax, 'cross-entropy   s_y [nm]'); end
        biasNm = 1000*hypot(fitted(s,1)-shifts(s,1), fitted(s,2)-shifts(s,2));
        flag = '';
        if multimodal(s), flag = '  MULTIMODAL'; end
        title(ax, sprintf('bias %.1f nm%s', biasNm, flag), ...
            'FontWeight', 'normal', 'FontSize', 9);
    end

    modeLabel = opts.aberrationMode;
    if strcmpi(modeLabel, 'none')
        headline = sprintf('No aberration, z = %.2f \\mum', opts.defocusUm);
    else
        headline = sprintf('%s = %.2f\\lambda RMS, z = %.2f \\mum', ...
            strrep(modeLabel, '_', '\_'), opts.amplitudeWaves, opts.defocusUm);
    end
    title(tl, sprintf(['%s (%s).  White circle = true position, red cross = ' ...
        'fitted position'], headline, sim.detectorLayout), 'FontWeight', 'bold');
end

% ----------------------------------------------------------------------------
function drawHoneycomb(ax, sim, colour)
    theta = (0:5)' * pi/3 + pi/6;    % point-up hex vertices
    for k = 1:sim.nDet
        vx = sim.detXY(k,1) + sim.detectorHexRadius * cos(theta);
        vy = sim.detXY(k,2) + sim.detectorHexRadius * sin(theta);
        plot(ax, [vx; vx(1)], [vy; vy(1)], '-', 'Color', colour, 'LineWidth', 0.6);
    end
end

% ----------------------------------------------------------------------------
function drawHoneycombFilled(ax, sim, values)
    theta = (0:5)' * pi/3 + pi/6;
    hold(ax, 'on');
    for k = 1:sim.nDet
        vx = sim.detXY(k,1) + sim.detectorHexRadius * cos(theta);
        vy = sim.detXY(k,2) + sim.detectorHexRadius * sin(theta);
        patch(ax, 'XData', vx, 'YData', vy, 'CData', values(k), ...
            'FaceColor', 'flat', 'EdgeColor', [0.35 0.35 0.35], 'LineWidth', 0.4);
    end
    clim(ax, [0 max(max(values), eps)]);
    set(ax, 'Color', 'w');
end
