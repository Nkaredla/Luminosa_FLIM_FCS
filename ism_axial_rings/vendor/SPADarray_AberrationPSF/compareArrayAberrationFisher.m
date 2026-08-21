function out = compareArrayAberrationFisher(varargin)
%COMPAREARRAYABERRATIONFISHER Compare aberration-estimation CRB: focused/defocused x 23/93 arrays.
%
%   out = compareArrayAberrationFisher()
%   out = compareArrayAberrationFisher('referencePhotons', 1e4)
%   out = compareArrayAberrationFisher('order23', 5, 'order93', 7)
%
%   Runs the Poisson Fisher-information / Cramer-Rao analysis
%   (simulateArrayFisherInformation) for four configurations and compares,
%   for every Zernike aberration mode, the flux-marginal CRB on the fitted
%   coefficient:
%
%       focused 23   - 23-channel array, no defocus,   modes to order 'order23' (5)
%       defocused 23 - 23-channel array, defocus bias, modes to order 'order23'
%       focused 93   - 93-channel array, no defocus,   modes to order 'order93' (7)
%       defocused 93 - 93-channel array, defocus bias, modes to order 'order93'
%
%   "focused" evaluates the Fisher information at zero aberration; "defocused"
%   evaluates it with a known defocus wavefront bias of 'defocusRad' radians
%   RMS present (defocus is still jointly fitted). This is the same defocus
%   term used by simulateArray93ImageFormation. Both arrays are sized to the
%   same Airy-unit footprint ('arrayDiameterAU') and use the same emitter-on-
%   glass high-NA oil optics, so the comparison isolates detector geometry,
%   focus condition, and mode order.
%
%   Lower CRB = better estimability of that mode. Modes of radial order 6-7
%   are only fitted by the 93 configurations, so the 23 series are absent
%   there. The per-mode CRB is the marginal bound from each array's own joint
%   fit of all its modes.
%
%   COMPUTE NOTE: this evaluates a finite-difference Jacobian over every
%   fitted mode for four configurations (order 7 -> 35 modes), so it is
%   compute-heavy. Reduce 'nx' or 'order93' for a faster preview.
%
%   KEY OPTIONS
%     'order23' (5), 'order93' (7)      radial orders fitted per array
%     'defocusRad' (1)                  defocus bias for the defocused configs
%     'photonCounts' (10.^(3:5))        photon budgets swept
%     'referencePhotons' (1e4)          photon count for the per-mode bar chart
%     'photonSweepMode' ('spherical')   mode shown in the CRB-vs-photons panel
%     Optics/geometry pass-throughs: 'NA' (1.45), 'lamEm' (0.520),
%       'sampleGeometry' ('airOnGlass'), 'emitterHeightUm' (0),
%       'nImmersion' (1.51), 'nGlass' (1.518), 'nSample' (1.0),
%       'coverslipThicknessUm' (170), 'arrayDiameterAU' (1.7),
%       'detFillRatio' (1.0), 'detectorSubsamples' (7), 'fovXY' (2.6),
%       'nx' (81), 'marginalizeFlux' (true), 'darkCountsPerDetector' (0),
%       'fdCoeff' (0.005).
%     Output: 'makeFigure' (true), 'saveFigures' (false), 'writeOutputs'
%       (false), 'outputDir' (''), 'verbose' (true).
%
%   OUTPUT
%     out.res23, out.res93     full simulateArrayFisherInformation results
%     out.modes                master mode list (93 order, superset)
%     out.modeOrder            radial order of each master mode
%     out.referencePhotons     photon count used for the bar chart
%     out.crb                  struct of [nMode x 1] CRB (waves): focused23,
%                              defocused23, focused93, defocused93
%     out.table                long-form comparison table
%     out.figPerMode           per-mode CRB comparison figure
%     out.figPhotonSweep       CRB-vs-photons figure for photonSweepMode

    addpath(fileparts(mfilename('fullpath')));
    opts = parseInputs(varargin{:});

    defocusWaves = opts.defocusRad / (2*pi);

    % Shared configuration: scenario 1 = focused (no aberration),
    % scenario 2 = defocused (defocus wavefront bias present).
    common = { ...
        'NA', opts.NA, 'lamEm', opts.lamEm, ...
        'sampleGeometry', opts.sampleGeometry, ...
        'emitterHeightUm', opts.emitterHeightUm, ...
        'nImmersion', opts.nImmersion, 'nGlass', opts.nGlass, ...
        'nSample', opts.nSample, 'nMedium', opts.nMedium, ...
        'coverslipThicknessUm', opts.coverslipThicknessUm, ...
        'arrayDiameterAU', opts.arrayDiameterAU, ...
        'detFillRatio', opts.detFillRatio, ...
        'detectorSubsamples', opts.detectorSubsamples, ...
        'fovXY', opts.fovXY, 'nx', opts.nx, ...
        'photonCounts', opts.photonCounts, ...
        'marginalizeFlux', opts.marginalizeFlux, ...
        'darkCountsPerDetector', opts.darkCountsPerDetector, ...
        'fdCoeff', opts.fdCoeff, ...
        'includeNoAberration', true, ...
        'aberrationModes', {'defocus'}, ...
        'amplitudeWaves', defocusWaves, ...
        'makeFigure', false, 'writeOutputs', false, ...
        'verbose', opts.verbose };

    if opts.verbose
        fprintf('[compareArrayAberrationFisher] running 23-channel array (order %d)...\n', ...
            opts.order23);
    end
    res23 = simulateArrayFisherInformation(common{:}, ...
        'detectorLayout', 'honeycomb23', 'maxZernikeOrder', opts.order23);

    if opts.verbose
        fprintf('[compareArrayAberrationFisher] running 93-channel array (order %d)...\n', ...
            opts.order93);
    end
    res93 = simulateArrayFisherInformation(common{:}, ...
        'detectorLayout', 'honeycomb93', 'maxZernikeOrder', opts.order93);

    % Master mode list = the higher-order (93) set, a superset of the 23 set.
    modes = res93.fitModes(:);
    nModes = numel(modes);
    modeOrder = cellfun(@modeRadialOrder, modes);

    [~, refIdx] = min(abs(opts.photonCounts - opts.referencePhotons));
    referencePhotons = opts.photonCounts(refIdx);

    crb = struct();
    crb.focused23   = alignCRB(res23, 1, refIdx, modes);
    crb.defocused23 = alignCRB(res23, 2, refIdx, modes);
    crb.focused93   = alignCRB(res93, 1, refIdx, modes);
    crb.defocused93 = alignCRB(res93, 2, refIdx, modes);

    resultsTable = buildComparisonTable(modes, modeOrder, referencePhotons, crb);

    figPerMode = [];
    figPhotonSweep = [];
    if opts.makeFigure
        figPerMode = plotPerModeComparison(modes, modeOrder, crb, ...
            referencePhotons, res93, opts);
        figPhotonSweep = plotPhotonSweepComparison(res23, res93, ...
            opts.photonSweepMode, opts);
    end

    if opts.writeOutputs || opts.saveFigures
        ensureDir(opts.outputDir);
        if opts.writeOutputs
            writetable(resultsTable, fullfile(opts.outputDir, ...
                'array_aberration_fisher_comparison.csv'));
        end
        if opts.saveFigures
            if ~isempty(figPerMode)
                exportgraphics(figPerMode, fullfile(opts.outputDir, ...
                    'array_aberration_crb_per_mode.png'), 'Resolution', 180);
            end
            if ~isempty(figPhotonSweep)
                exportgraphics(figPhotonSweep, fullfile(opts.outputDir, ...
                    'array_aberration_crb_vs_photons.png'), 'Resolution', 180);
            end
        end
    end

    if opts.verbose
        reportComparison(modes, modeOrder, crb, referencePhotons);
    end

    out = struct();
    out.res23             = res23;
    out.res93             = res93;
    out.modes             = modes;
    out.modeOrder         = modeOrder;
    out.referencePhotons  = referencePhotons;
    out.photonCounts      = opts.photonCounts;
    out.defocusRad        = opts.defocusRad;
    out.defocusWaves      = defocusWaves;
    out.crb               = crb;
    out.table             = resultsTable;
    out.figPerMode        = figPerMode;
    out.figPhotonSweep    = figPhotonSweep;
    out.options           = opts;
end

% ============================================================================
function opts = parseInputs(varargin)
    p = inputParser;
    p.FunctionName = 'compareArrayAberrationFisher';

    addParameter(p, 'order23', 5);
    addParameter(p, 'order93', 7);
    addParameter(p, 'defocusRad', 1);
    addParameter(p, 'photonCounts', 10.^(3:5));
    addParameter(p, 'referencePhotons', 1e4);
    addParameter(p, 'photonSweepMode', 'spherical');

    % Optics / geometry pass-throughs (shared by both arrays)
    addParameter(p, 'NA', 1.45);
    addParameter(p, 'lamEm', 0.520);
    addParameter(p, 'sampleGeometry', 'airOnGlass');
    addParameter(p, 'emitterHeightUm', 0);
    addParameter(p, 'nImmersion', 1.51);
    addParameter(p, 'nGlass', 1.518);
    addParameter(p, 'nSample', 1.0);
    addParameter(p, 'coverslipThicknessUm', 170);
    addParameter(p, 'nMedium', 1.33);
    addParameter(p, 'arrayDiameterAU', 1.7);
    addParameter(p, 'detFillRatio', 1.0);
    addParameter(p, 'detectorSubsamples', 7);
    addParameter(p, 'fovXY', 2.6);
    addParameter(p, 'nx', 81);
    addParameter(p, 'marginalizeFlux', true);
    addParameter(p, 'darkCountsPerDetector', 0);
    addParameter(p, 'fdCoeff', 0.005);

    addParameter(p, 'makeFigure', true);
    addParameter(p, 'saveFigures', false);
    addParameter(p, 'writeOutputs', false);
    addParameter(p, 'outputDir', '');
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;

    opts.photonSweepMode = char(opts.photonSweepMode);
    opts.photonCounts = opts.photonCounts(:).';
    if isempty(opts.outputDir)
        opts.outputDir = fullfile(fileparts(mfilename('fullpath')), ...
            'output_matlab', 'array_aberration_fisher_comparison');
    end
end

% ----------------------------------------------------------------------------
function v = alignCRB(res, scenarioIdx, photonIdx, targetModes)
%ALIGNCRB Map a result's per-mode CRB onto the master mode list (NaN if absent).
    v = nan(numel(targetModes), 1);
    src = reshape(res.crbWaves(scenarioIdx, photonIdx, :), [], 1);
    for i = 1:numel(targetModes)
        j = find(strcmp(res.fitModes, targetModes{i}), 1);
        if ~isempty(j)
            v(i) = src(j);
        end
    end
end

% ----------------------------------------------------------------------------
function o = modeRadialOrder(mode)
    o = NaN;
    for oo = 1:7
        if any(strcmp(zernikeModeNames(oo), mode))
            o = oo;
            return;
        end
    end
end

% ----------------------------------------------------------------------------
function T = buildComparisonTable(modes, modeOrder, referencePhotons, crb)
    configNames = {'focused23','defocused23','focused93','defocused93'};
    rows = {};
    for c = 1:numel(configNames)
        vec = crb.(configNames{c});
        for m = 1:numel(modes)
            rows(end+1, :) = { configNames{c}, modes{m}, modeOrder(m), ...
                referencePhotons, vec(m), 1000*vec(m) }; %#ok<AGROW>
        end
    end
    T = cell2table(rows, 'VariableNames', ...
        {'config', 'mode', 'radialOrder', 'photons', 'crbWaves', 'crbMilliWaves'});
end

% ----------------------------------------------------------------------------
function fig = plotPerModeComparison(modes, modeOrder, crb, referencePhotons, res93, opts)
    nModes = numel(modes);
    series = { crb.focused23,   'focused 23',   [0.20 0.44 0.74], 'o'; ...
               crb.defocused23, 'defocused 23', [0.30 0.65 0.90], 's'; ...
               crb.focused93,   'focused 93',   [0.85 0.33 0.10], '^'; ...
               crb.defocused93, 'defocused 93', [0.95 0.60 0.20], 'd' };

    fig = figure('Color', 'w', 'Position', [60 80 min(120+42*nModes, 1700) 560]);
    ax = axes('Parent', fig);
    hold(ax, 'on');
    x = (1:nModes).';

    hLines = gobjects(size(series,1), 1);
    for k = 1:size(series,1)
        y = 1000 * series{k,1}(:);          % milli-waves
        y(~isfinite(y)) = NaN;
        hLines(k) = plot(ax, x, y, series{k,4}, 'Color', series{k,3}, ...
            'MarkerFaceColor', series{k,3}, 'MarkerSize', 6, ...
            'LineStyle', '-', 'LineWidth', 0.75);
    end

    set(ax, 'YScale', 'log');
    grid(ax, 'on');
    drawOrderSeparators(ax, modeOrder);
    set(ax, 'XTick', x, 'XTickLabel', cellfun(@(m) strrep(m,'_','\_'), modes, ...
        'UniformOutput', false));
    try, xtickangle(ax, 45); catch, end
    xlim(ax, [0.5 nModes + 0.5]);
    ylabel(ax, 'flux-marginal CRB [m\lambda]');
    legend(ax, hLines, series(:,2), 'Location', 'northwest');
    title(ax, sprintf(['Per-mode aberration-estimation CRB at %.0f photons ' ...
        '(NA %.2f, footprint %.2f AU, defocus %.2f rad)'], ...
        referencePhotons, res93.sim.NA, res93.netDiameterAU, opts.defocusRad), ...
        'FontWeight', 'bold');
    hold(ax, 'off');
end

% ----------------------------------------------------------------------------
function drawOrderSeparators(ax, modeOrder)
    yl = ylim(ax);
    changeIdx = find(diff(modeOrder(:)) ~= 0);
    for c = changeIdx(:).'
        xv = c + 0.5;
        plot(ax, [xv xv], yl, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    end
    % Order labels at the top of each order block.
    orders = unique(modeOrder(isfinite(modeOrder)), 'stable');
    for o = orders(:).'
        idx = find(modeOrder == o);
        if isempty(idx), continue; end
        xc = mean([min(idx), max(idx)]);
        text(ax, xc, yl(2), sprintf('order %d', o), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 8, 'Color', [0.4 0.4 0.4]);
    end
    ylim(ax, yl);
end

% ----------------------------------------------------------------------------
function fig = plotPhotonSweepComparison(res23, res93, sweepMode, opts)
    fig = figure('Color', 'w', 'Position', [60 80 640 480]);
    ax = axes('Parent', fig);
    hold(ax, 'on');

    series = { res23, 1, 'focused 23',   [0.20 0.44 0.74], 'o'; ...
               res23, 2, 'defocused 23', [0.30 0.65 0.90], 's'; ...
               res93, 1, 'focused 93',   [0.85 0.33 0.10], '^'; ...
               res93, 2, 'defocused 93', [0.95 0.60 0.20], 'd' };

    hLines = gobjects(size(series,1), 1);
    labels = cell(size(series,1), 1);
    for k = 1:size(series,1)
        res = series{k,1};
        scen = series{k,2};
        j = find(strcmp(res.fitModes, sweepMode), 1);
        if isempty(j)
            y = nan(1, numel(res.photonCounts));
        else
            y = 1000 * reshape(res.crbWaves(scen, :, j), 1, []);
        end
        y(~isfinite(y)) = NaN;
        hLines(k) = plot(ax, res.photonCounts, y, series{k,5}, ...
            'Color', series{k,4}, 'MarkerFaceColor', series{k,4}, ...
            'MarkerSize', 6, 'LineStyle', '-', 'LineWidth', 1);
        labels{k} = series{k,3};
    end

    set(ax, 'XScale', 'log', 'YScale', 'log');
    grid(ax, 'on');
    xlabel(ax, 'photons per plane');
    ylabel(ax, 'flux-marginal CRB [m\lambda]');
    legend(ax, hLines, labels, 'Location', 'northeast');
    title(ax, sprintf('CRB vs photons for %s', strrep(sweepMode,'_','\_')), ...
        'FontWeight', 'bold');
    hold(ax, 'off');
end

% ----------------------------------------------------------------------------
function reportComparison(modes, modeOrder, crb, referencePhotons)
    fprintf('\n[compareArrayAberrationFisher] CRB (m-lambda) at %.0f photons\n', ...
        referencePhotons);
    fprintf('  %-22s %5s %10s %10s %10s %10s\n', 'mode', 'order', ...
        'foc23', 'def23', 'foc93', 'def93');
    for m = 1:numel(modes)
        fprintf('  %-22s %5d %10s %10s %10s %10s\n', modes{m}, modeOrder(m), ...
            fmtCrb(crb.focused23(m)), fmtCrb(crb.defocused23(m)), ...
            fmtCrb(crb.focused93(m)), fmtCrb(crb.defocused93(m)));
    end
end

function s = fmtCrb(v)
    if ~isfinite(v)
        s = '   -   ';
    else
        s = sprintf('%.3g', 1000*v);
    end
end

% ----------------------------------------------------------------------------
function ensureDir(d)
    if ~isempty(d) && exist(d, 'dir') ~= 7
        mkdir(d);
    end
end
