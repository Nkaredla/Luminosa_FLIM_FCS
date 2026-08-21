function out = simulateArrayFisherInformation(varargin)
%SIMULATEARRAYFISHERINFORMATION Poisson Fisher/CRB on a SPAD array vs photons and aberration.
%
%   out = simulateArrayFisherInformation('detectorLayout','honeycomb23','maxZernikeOrder',5)
%   out = simulateArrayFisherInformation('detectorLayout','honeycomb93','maxZernikeOrder',7)
%
%   Layout-parametrized Poisson Fisher-information / Cramer-Rao analysis for
%   the detector "microimage" of a point emitter, using the same conventions
%   as simulateRetrieveCenterMicroimageAberrations.m. It works on any
%   supported honeycomb layout ('honeycomb23' or 'honeycomb93') and any
%   Zernike order, and is the shared engine behind
%   simulateArray23FisherInformation, simulateArray93FisherInformation, and
%   compareArrayAberrationFisher.
%
%   The default optics place the emitter on the coverslip glass (air on the
%   sample side) imaged by a high-NA oil objective (layered vectorial
%   Richards-Wolf air-interface model). The array pitch is set from a target
%   Airy-unit footprint (arrayDiameterAU) so different layouts can be sized
%   to the same field.
%
%   MODEL / CONVENTIONS (identical to simulateRetrieveCenterMicroimageAberrations)
%   - The microimage is the nDet-vector of collected intensities for the
%     on-axis emitter, per-plane normalized to a probability vector p(theta).
%     Expected counts are mu_k = N*p_k + darkCounts.
%   - Poisson information I = J' diag(1/mu) J with J = d(mu)/d(theta), theta
%     the fitted Zernike coefficients in waves at sim.lamRef.
%   - With 'marginalizeFlux' true (default) the per-plane photon flux is a
%     nuisance and is projected out (flux-marginal CRB for an unknown-
%     brightness point source). Both marginal and conditional CRBs are kept.
%   - Poisson information is linear in the photon number, so the
%     probability-domain Jacobian is computed once per aberration operating
%     point and reused across all photon counts. CRB then falls as ~1/sqrt(N).
%
%   Detection-side model (no scanning excitation), matching
%   simulateArray93ImageFormation.
%
%   KEY NAME-VALUE OPTIONS
%     'detectorLayout'   'honeycomb23' or 'honeycomb93' (default 'honeycomb93')
%     'maxZernikeOrder'  highest radial order fitted; sets fitModes when
%                        'fitModes' is empty (default 4). The azimuthal
%                        harmonic count sim.M is raised to cover it.
%     'fitModes'         explicit fit-mode cell (overrides maxZernikeOrder)
%     'aberrationModes'  operating-point modes to evaluate Fisher at
%                        (default {'defocus','astig_x','coma_x','spherical'})
%     'includeNoAberration' add a zero-aberration baseline (default true)
%     'amplitudeWaves'   operating-point amplitude [waves] (default 0.15)
%     'photonCounts'     photon budgets per plane (default 10.^(2:6))
%     'planeZ'           stage plane(s) [um] (default 0, in focus)
%     'marginalizeFlux'  project out unknown flux (default true)
%     Optics/geometry: 'NA' (1.45), 'lamEm' (0.520), 'sampleGeometry'
%       ('airOnGlass'), 'emitterHeightUm' (0), 'nImmersion' (1.51),
%       'nGlass' (1.518), 'nSample' (1.0), 'coverslipThicknessUm' (170),
%       'nMedium' (1.33 fallback), 'arrayDiameterAU' (1.7), 'detPitchUm'
%       ([] -> derived), 'detFillRatio' (1.0), 'detectorSubsamples' (7),
%       'fovXY' (2.6), 'nx' (121).
%     Numerics: 'darkCountsPerDetector' (0), 'fdCoeff' (0.005),
%       'minExpectedCount' (1e-9), 'svdRelativeTolerance' (1e-4).
%     Output: 'colormap' ('lines'), 'makeFigure' (true), 'saveFigures'
%       (false), 'writeOutputs' (false), 'outputDir' (''), 'verbose' (true).
%
%   OUTPUT: see field comments at the end of the function.

    addpath(fileparts(mfilename('fullpath')));
    opts = parseInputs(varargin{:});

    sim = configureSim(opts);

    fitModes = opts.fitModes;
    nModes   = numel(fitModes);
    photons  = opts.photonCounts(:).';
    nPhoton  = numel(photons);
    planeZ   = opts.planeZ(:).';
    nPlane   = numel(planeZ);

    [~, xc] = min(abs(sim.x));
    [~, yc] = min(abs(sim.y));

    scenarios = buildScenarios(opts);          % {label, coeffStruct}
    nScen = size(scenarios, 1);

    darkPerDet = resolveDark(opts.darkCountsPerDetector, sim.nDet);
    darkStacked = repmat(darkPerDet, nPlane, 1);

    crbWaves     = nan(nScen, nPhoton, nModes);
    crbWavesCond = nan(nScen, nPhoton, nModes);
    totalInfo    = nan(nScen, nPhoton);
    condNumber   = nan(nScen, nPhoton);
    rankArr      = nan(nScen, nPhoton);
    infoAtMax    = cell(nScen, 1);

    for s = 1:nScen
        coeffs0 = scenarios{s, 2};

        probs0 = probMicroimage(sim, coeffs0, planeZ, yc, xc);   % [nDet x nPlane]
        probVec0 = probs0(:);
        dP = probabilityJacobian(sim, coeffs0, fitModes, planeZ, yc, xc, opts.fdCoeff);

        for j = 1:nPhoton
            N = photons(j);
            fish = fisherForPhotons(dP, probVec0, probs0, N, darkStacked, ...
                opts.marginalizeFlux, opts.minExpectedCount, ...
                opts.svdRelativeTolerance, nModes, nPlane, sim.nDet);

            crbWaves(s, j, :)     = fish.crbMarginal(:).';
            crbWavesCond(s, j, :) = fish.crbConditional(:).';
            totalInfo(s, j)       = trace(fish.informationMarginal);
            condNumber(s, j)      = fish.conditionNumber;
            rankArr(s, j)         = fish.rank;
            if j == nPhoton
                infoAtMax{s} = fish.informationMarginal;
            end
        end

        if opts.verbose
            reportScenario(scenarios{s,1}, fitModes, photons(end), ...
                squeeze(crbWaves(s, end, :)), rankArr(s, end), nModes, ...
                sim.detectorLayout);
        end
    end

    resultsTable = buildTable(scenarios(:,1), fitModes, photons, ...
        crbWaves, crbWavesCond, sim.detectorLayout);

    figCRB = [];
    if opts.makeFigure
        figCRB = plotCRBvsPhotons(scenarios(:,1), fitModes, photons, ...
            crbWaves, opts, sim.detectorLayout);
    end

    airyUnitUm = 1.22 * sim.lamEm / sim.NA;
    spanX = (max(sim.detXY(:,1)) - min(sim.detXY(:,1))) + sim.detSize;
    spanY = (max(sim.detXY(:,2)) - min(sim.detXY(:,2))) + 2*sim.detectorHexRadius;
    netDiameterAU = max(spanX, spanY) / airyUnitUm;

    if opts.writeOutputs || opts.saveFigures
        ensureDir(opts.outputDir);
        if opts.writeOutputs
            writetable(resultsTable, fullfile(opts.outputDir, ...
                sprintf('%s_fisher_crb.csv', sim.detectorLayout)));
        end
        if opts.saveFigures && ~isempty(figCRB)
            exportgraphics(figCRB, fullfile(opts.outputDir, ...
                sprintf('%s_fisher_crb_vs_photons.png', sim.detectorLayout)), ...
                'Resolution', 180);
        end
    end

    out = struct();
    out.sim                    = sim;               % configured array sim
    out.detectorLayout         = sim.detectorLayout;
    out.fitModes               = fitModes;          % estimated Zernike modes
    out.maxZernikeOrder        = opts.maxZernikeOrder;
    out.photonCounts           = photons;           % photon budgets used
    out.aberrationLabels       = scenarios(:,1);    % operating-point labels
    out.coeffScenarios         = scenarios(:,2);    % operating-point structs
    out.crbWaves               = crbWaves;          % [nScen x nPhoton x nMode] marginal CRB
    out.crbWavesConditional    = crbWavesCond;      % flux-conditional CRB
    out.totalInformation       = totalInfo;         % trace of marginal info
    out.conditionNumber        = condNumber;
    out.rank                   = rankArr;
    out.informationAtMaxPhotons = infoAtMax;
    out.table                  = resultsTable;
    out.figCRB                 = figCRB;
    out.detPitchUm             = sim.detPitch;
    out.airyUnitUm             = airyUnitUm;
    out.netDiameterAU          = netDiameterAU;
    out.options                = opts;
end

% ============================================================================
function opts = parseInputs(varargin)
    p = inputParser;
    p.FunctionName = 'simulateArrayFisherInformation';

    % Array / analysis
    addParameter(p, 'detectorLayout', 'honeycomb93');
    addParameter(p, 'maxZernikeOrder', 4);
    addParameter(p, 'fitModes', {});

    % Optics / geometry
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
    addParameter(p, 'detPitchUm', []);
    addParameter(p, 'detFillRatio', 1.0);
    addParameter(p, 'detectorSubsamples', 7);
    addParameter(p, 'fovXY', 2.6);
    addParameter(p, 'nx', 121);

    % Aberration operating points
    addParameter(p, 'aberrationModes', {'defocus','astig_x','coma_x','spherical'});
    addParameter(p, 'includeNoAberration', true);
    addParameter(p, 'amplitudeWaves', 0.15);

    % Photon / noise
    addParameter(p, 'photonCounts', 10.^(2:6));
    addParameter(p, 'planeZ', 0);
    addParameter(p, 'darkCountsPerDetector', 0);
    addParameter(p, 'marginalizeFlux', true);
    addParameter(p, 'fdCoeff', 0.005);
    addParameter(p, 'minExpectedCount', 1e-9);
    addParameter(p, 'svdRelativeTolerance', 1e-4);

    % Output / display
    addParameter(p, 'colormap', 'lines');
    addParameter(p, 'makeFigure', true);
    addParameter(p, 'saveFigures', false);
    addParameter(p, 'writeOutputs', false);
    addParameter(p, 'outputDir', '');
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;

    opts.detectorLayout = char(opts.detectorLayout);
    if mod(opts.nx, 2) == 0, opts.nx = opts.nx + 1; end
    if isempty(opts.fitModes)
        opts.fitModes = zernikeModeNames(opts.maxZernikeOrder);
    end
    if ischar(opts.fitModes) || isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    end
    if ischar(opts.aberrationModes) || isstring(opts.aberrationModes)
        opts.aberrationModes = cellstr(opts.aberrationModes);
    end
    if isempty(opts.detPitchUm)
        opts.detPitchUm = detectorPitchForDiameterAU(opts.detectorLayout, ...
            opts.arrayDiameterAU, opts.NA, opts.lamEm, opts.detFillRatio);
    end
    if isempty(opts.outputDir)
        opts.outputDir = fullfile(fileparts(mfilename('fullpath')), ...
            'output_matlab', 'array_fisher_information');
    end
end

% ----------------------------------------------------------------------------
function sim = configureSim(opts)
%CONFIGURESIM Build the array emitter-on-glass sim (mirrors simulateArray93ImageFormation).
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

    % Raise the azimuthal-harmonic bandwidth to cover the highest Zernike
    % order requested (m up to maxOrder, plus one for vectorial coupling).
    sim.M = max(sim.M, opts.maxZernikeOrder + 1);
    sim.modeOrder = zernikeModeNames(max(6, opts.maxZernikeOrder));

    sim.sampleGeometry = normalizeSampleGeometry(opts.sampleGeometry);
    if usesInterfaceGeometry(sim)
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

    sim.detectorLayout = opts.detectorLayout;
    sim.detPitch = opts.detPitchUm;
    [sim.detXY, sim.detectorIndexGrid, sim.detectorLayoutInfo] = ...
        detectorLayout(sim.detectorLayout, sim.detPitch);
    sim.nDet = size(sim.detXY, 1);
    sim.detSize = opts.detFillRatio * sim.detPitch;
    sim.detectorHexRadius = sim.detSize / sqrt(3);
end

% ----------------------------------------------------------------------------
function h = arrayDetectionPSF(sim, coeffs, planeZ)
%ARRAYDETECTIONPSF Collection-path detection PSF at the given stage plane(s).
    simDet = setVectorialPolarizationForPath(sim, 'collection');
    if usesInterfaceGeometry(sim)
        simDet.airInterfaceDirection = 'collection';
        h = psfBesselAirInterface(simDet, coeffs, sim.lamEm, ...
            planeZ, sim.interfaceEmitterHeightUm);
    else
        simDet.z = planeZ(:).';
        simDet.nz = numel(planeZ);
        h = psfBessel(simDet, coeffs, sim.lamEm);
    end
end

% ----------------------------------------------------------------------------
function probs = probMicroimage(sim, coeffs, planeZ, yc, xc)
%PROBMICROIMAGE Per-plane-normalized detector probability vector [nDet x nPlane].
    coeffs = coeffStruct(sim, coeffs);
    hDet = arrayDetectionPSF(sim, coeffs, planeZ);
    eta = detectorCollectionEfficiencyExplicit(sim, hDet);   % [ny nx nPlane nDet]
    nPlane = numel(planeZ);
    vals = zeros(sim.nDet, nPlane);
    for ip = 1:nPlane
        vals(:, ip) = squeeze(eta(yc, xc, ip, :));
    end
    probs = perPlaneNormalize(vals);
end

% ----------------------------------------------------------------------------
function dP = probabilityJacobian(sim, coeffs0, fitModes, planeZ, yc, xc, h)
%PROBABILITYJACOBIAN Central-difference d(prob)/d(coeff) over the fit modes.
    coeffs0 = coeffStruct(sim, coeffs0);
    p0 = probMicroimage(sim, coeffs0, planeZ, yc, xc);
    dP = zeros(numel(p0), numel(fitModes));
    for k = 1:numel(fitModes)
        cp = perturbCoeff(coeffs0, fitModes{k}, +h);
        cm = perturbCoeff(coeffs0, fitModes{k}, -h);
        pp = probMicroimage(sim, cp, planeZ, yc, xc);
        pm = probMicroimage(sim, cm, planeZ, yc, xc);
        dP(:, k) = (pp(:) - pm(:)) / (2 * h);
    end
end

function coeffs = perturbCoeff(coeffs0, modeName, delta)
    coeffs = coeffs0;
    if isfield(coeffs, modeName)
        coeffs.(modeName) = coeffs.(modeName) + delta;
    else
        coeffs.(modeName) = delta;
    end
end

% ----------------------------------------------------------------------------
function fish = fisherForPhotons(dP, probVec0, probs0, N, darkStacked, ...
        marginalizeFlux, minExpectedCount, svdTol, nModes, nPlane, nDet)
%FISHERFORPHOTONS Poisson Fisher/CRB at photon budget N (mirrors reference math).
    mu = max(N * probVec0 + darkStacked, minExpectedCount);
    Jcounts = N * dP;                                  % d(mu)/d(theta)
    weightedJcond = Jcounts ./ sqrt(mu);

    weightedJ = weightedJcond;
    if marginalizeFlux
        fluxSensitivity = perPlaneFluxSensitivity(probs0, nDet, nPlane);
        weightedFlux = fluxSensitivity ./ sqrt(mu);
        Q = orthonormalColumns(weightedFlux, svdTol);
        if ~isempty(Q)
            weightedJ = weightedJcond - Q * (Q.' * weightedJcond);
        end
    end

    marginal    = fisherStatsFromWeightedJ(weightedJ, nModes, svdTol);
    conditional = fisherStatsFromWeightedJ(weightedJcond, nModes, svdTol);

    fish = struct();
    fish.informationMarginal    = marginal.information;
    fish.informationConditional = conditional.information;
    fish.crbMarginal            = marginal.crbStd;
    fish.crbConditional         = conditional.crbStd;
    fish.conditionNumber        = marginal.conditionNumber;
    fish.rank                   = marginal.rank;
end

% ----------------------------------------------------------------------------
function sensitivity = perPlaneFluxSensitivity(probs, nDet, nPlane)
    sensitivity = zeros(nDet * nPlane, nPlane);
    for ip = 1:nPlane
        rows = (1:nDet) + (ip - 1) * nDet;
        sensitivity(rows, ip) = probs(:, ip);
    end
end

% ----------------------------------------------------------------------------
function Q = orthonormalColumns(A, relTol)
    if isempty(A)
        Q = zeros(size(A, 1), 0);
        return;
    end
    [Qfull, R] = qr(A, 0);
    d = abs(diag(R));
    if isempty(d)
        Q = zeros(size(A, 1), 0);
        return;
    end
    keep = d > relTol * max(max(d), eps);
    Q = Qfull(:, keep);
end

% ----------------------------------------------------------------------------
function stats = fisherStatsFromWeightedJ(weightedJ, nModes, relTol)
    I = weightedJ.' * weightedJ;
    I = 0.5 * (I + I.');

    [eigVecs, eigMat] = eig(I);
    eigVals = real(diag(eigMat));
    [eigVals, order] = sort(eigVals, 'descend');
    eigVecs = real(eigVecs(:, order));
    eigVals = max(eigVals, 0);
    relative = eigVals / max(eigVals(1), eps);
    rankValue = sum(relative > relTol);

    if rankValue >= nModes
        conditionNumber = eigVals(1) / max(eigVals(nModes), eps);
        covar = inv(I);
        crbStd = sqrt(max(diag(covar), 0));
    else
        conditionNumber = Inf;
        tolAbs = relTol * max(eigVals(1), eps);
        covar = pinv(I, tolAbs);
        crbStd = sqrt(max(diag(covar), 0));
        nullBasis = eigVecs(:, rankValue+1:end);
        for k = 1:numel(crbStd)
            if ~isempty(nullBasis) && norm(nullBasis(k, :)) > 1e-6
                crbStd(k) = Inf;
            end
        end
    end

    stats = struct();
    stats.information = I;
    stats.covariance = covar;
    stats.crbStd = crbStd(:);
    stats.eigenValues = eigVals(:);
    stats.rank = rankValue;
    stats.conditionNumber = conditionNumber;
    stats.isFullRank = rankValue >= nModes;
end

% ----------------------------------------------------------------------------
function y = perPlaneNormalize(values)
    values = max(double(values), 0);
    y = zeros(size(values));
    for ip = 1:size(values, 2)
        s = sum(values(:, ip));
        if s > eps
            y(:, ip) = values(:, ip) / s;
        end
    end
end

% ----------------------------------------------------------------------------
function dark = resolveDark(darkCounts, nDet)
    dark = double(darkCounts(:));
    if isscalar(dark)
        dark = repmat(dark, nDet, 1);
    elseif numel(dark) ~= nDet
        error('simulateArrayFisherInformation:BadDarkCounts', ...
            'darkCountsPerDetector must be scalar or contain %d values.', nDet);
    end
end

% ----------------------------------------------------------------------------
function scen = buildScenarios(opts)
    scen = {};
    if opts.includeNoAberration
        scen(end+1, :) = {'No aberration', struct()}; %#ok<AGROW>
    end
    A = opts.amplitudeWaves;
    for k = 1:numel(opts.aberrationModes)
        mode = opts.aberrationModes{k};
        scen(end+1, :) = { sprintf('%s = %.2f\\lambda', prettyMode(mode), A), ...
            struct(mode, A) }; %#ok<AGROW>
    end
    if isempty(scen)
        scen = {'No aberration', struct()};
    end
end

% ----------------------------------------------------------------------------
function name = prettyMode(mode)
    switch lower(mode)
        case 'tilt_x',    name = 'Tilt X';
        case 'tilt_y',    name = 'Tilt Y';
        case 'defocus',   name = 'Defocus';
        case 'astig_x',   name = 'Astig 0\circ';
        case 'astig_y',   name = 'Astig 45\circ';
        case 'coma_x',    name = 'Coma X';
        case 'coma_y',    name = 'Coma Y';
        case 'trefoil_x', name = 'Trefoil X';
        case 'trefoil_y', name = 'Trefoil Y';
        case 'spherical', name = 'Spherical';
        otherwise,        name = strrep(mode, '_', '\_');
    end
end

% ----------------------------------------------------------------------------
function T = buildTable(labels, fitModes, photons, crbWaves, crbWavesCond, layout)
    rows = {};
    for s = 1:numel(labels)
        rawLabel = regexprep(labels{s}, '\\[a-zA-Z]+|[{}]', '');
        for j = 1:numel(photons)
            for m = 1:numel(fitModes)
                rows(end+1, :) = { layout, rawLabel, photons(j), fitModes{m}, ...
                    crbWaves(s, j, m), crbWavesCond(s, j, m) }; %#ok<AGROW>
            end
        end
    end
    T = cell2table(rows, 'VariableNames', ...
        {'layout', 'aberration', 'photons', 'mode', 'crbWaves', 'crbWavesConditional'});
end

% ----------------------------------------------------------------------------
function reportScenario(label, fitModes, N, crb, rankValue, nModes, layout)
    plain = regexprep(label, '\\[a-zA-Z]+|[{}]', '');
    finite = crb(isfinite(crb));
    if isempty(finite)
        worst = Inf; worstName = 'n/a';
    else
        [worst, idx] = max(crb);
        worstName = fitModes{idx};
    end
    fprintf(['[%s Fisher] %-22s  N=%.0g: rank %d/%d, worst CRB %.4g ' ...
        'waves (%s)\n'], layout, plain, N, rankValue, nModes, worst, worstName);
end

% ----------------------------------------------------------------------------
function fig = plotCRBvsPhotons(labels, fitModes, photons, crbWaves, opts, layout)
    nScen = numel(labels);
    nCol = ceil(sqrt(nScen));
    nRow = ceil(nScen / nCol);
    fig = figure('Color', 'w', 'Position', [80 80 340*nCol 300*nRow]);
    tl = tiledlayout(fig, nRow, nCol, 'Padding', 'compact', 'TileSpacing', 'compact');
    cmap = feval(opts.colormap, numel(fitModes));

    lastAx = [];
    for s = 1:nScen
        ax = nexttile(tl);
        hold(ax, 'on');
        for m = 1:numel(fitModes)
            y = 1000 * squeeze(crbWaves(s, :, m));   % milli-waves
            y(~isfinite(y)) = NaN;
            plot(ax, photons, y, '-o', 'Color', cmap(m,:), ...
                'MarkerFaceColor', cmap(m,:), 'MarkerSize', 3, 'LineWidth', 1);
        end
        set(ax, 'XScale', 'log', 'YScale', 'log');
        grid(ax, 'on');
        xlabel(ax, 'photons per plane');
        ylabel(ax, 'CRB [m\lambda]');
        title(ax, labels{s}, 'FontWeight', 'normal');
        lastAx = ax;
    end

    if ~isempty(lastAx)
        lgd = legend(lastAx, cellfun(@(m) strrep(m,'_','\_'), fitModes, ...
            'UniformOutput', false));
        lgd.Layout.Tile = 'east';
    end
    title(tl, sprintf(['Flux-marginal Cramer-Rao bound on Zernike ' ...
        'coefficients (%s)'], layout), 'FontWeight', 'bold');
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
            error('simulateArrayFisherInformation:BadSampleGeometry', ...
                ['sampleGeometry must be ''airOnGlass'' (emitter on the ' ...
                'coverslip) or ''homogeneous''.']);
    end
end

% ----------------------------------------------------------------------------
function ensureDir(d)
    if ~isempty(d) && exist(d, 'dir') ~= 7
        mkdir(d);
    end
end
