function out = simulateAberrationLocalizationBias(varargin)
%SIMULATEABERRATIONLOCALIZATIONBIAS Localization bias from fitting an ideal PSF to aberrated data.
%
%   out = simulateAberrationLocalizationBias()
%   out = simulateAberrationLocalizationBias('aberrationModes',{'coma_x'}, ...
%             'amplitudeWaves',0.1,'defocusUm',[0 0.3])
%
%   A single molecule diffuses in the sample plane and is localized from the
%   SPAD-array microimage. The optics carry an unknown aberration, but the
%   fitter assumes an ideal (unaberrated, in-focus) PSF. This routine maps the
%   resulting systematic position error over the field and converts it into the
%   quantity that actually corrupts a diffusion measurement.
%
%   MODEL
%   Conditioned on a photon reaching the array, the microimage of an emitter at
%   r = (x,y) with stage defocus z is the channel-allocation vector
%
%       p_k(r,z) = eta_k(r,z) / sum_j eta_j(r,z)
%
%   with eta from detectorCollectionEfficiencyExplicit, i.e. the same
%   convention as +membrane_tracking/+focused_ism/ismDetectorChannelProbability.
%   Expected counts are mu_k = A*p_k, with A (captured signal) a free nuisance
%   parameter, since a diffusing molecule's brightness varies frame to frame
%   through the illumination profile, blinking and bleaching.
%
%   PSEUDO-TRUE POSITION
%   Because sum_k p_k = 1, the amplitude decouples exactly: profiling A out of
%   the Poisson log-likelihood leaves the cross-entropy, so the misspecified MLE
%   converges to
%
%       rHat(r) = argmax_s  sum_k p_k^true(r,z) * log p_k^ideal(s,0)
%
%   The bias field b(r) = rHat(r) - r is therefore independent of photon number.
%   It is a floor that no amount of signal removes, which is what a Cramer-Rao
%   analysis on its own does not show. Strictly this is the N -> infinity limit:
%   at finite N the MLE carries an additional O(1/N) bias on top of b, so b is
%   the asymptote, not the exact finite-N bias.
%
%   The argmax is GLOBAL. Once the true PSF is both aberrated and defocused the
%   cross-entropy becomes multimodal (a comatic microimage has a core and a
%   flare, and the ideal template can lock onto either), so the solver seeds
%   Newton with an exhaustive scan over the ideal-model nodes. Seeding at the
%   true position instead returns a local maximum whose basin changes across
%   the field, which puts a step discontinuity into b(r) and, after
%   differentiation, an arbitrarily large spurious apparent-D ratio.
%
%   WHAT THE BIAS DOES TO A DIFFUSION MEASUREMENT
%   A spatially uniform b is only an origin shift and is unobservable. The
%   damage comes from the structure of b:
%     - Jacobian J = I + db/dr maps true steps to observed ones, so an isotropic
%       D is measured as the tensor D*J*J'. Reported here as an apparent-D ratio
%       trace(J*J')/2 and an anisotropy sqrt(lmax/lmin) of J*J'. Note that the
%       reported anisotropy is the singular-value ratio of J, i.e. the ratio of
%       LENGTH scalings; the ratio of the apparent diffusion coefficients along
%       the two principal axes is its square.
%     - Spatial variation of J makes the apparent D position dependent, which
%       reads as heterogeneous or anomalous diffusion.
%     - Defocus dependence of b turns axial excursions into apparent lateral
%       motion.
%
%   ERROR BUDGET
%   With free amplitude and no background the per-photon Fisher information is
%   the multinomial form I1 = sum_k (grad p_k)(grad p_k)'/p_k, evaluated at the
%   true (aberrated) model. Precision is then sigma(N) = sigma1/sqrt(N), and the
%   crossover photon count where bias overtakes noise is N* = (sigma1/|b|)^2.
%   Past N* extra photons buy nothing until the aberration is calibrated.
%
%   BRIGHTNESS AND PHOTON STATISTICS
%   Nothing here is Monte Carlo: no counts are sampled and poissonSample is
%   never called. Both results are analytic limits of the Poisson model, so a
%   photophysics chain (laser power, saturation, quantum yield, blinking,
%   bleaching) would not change either one:
%     - the bias is the N -> infinity limit and is exactly N independent;
%     - the precision is exactly sigma1/sqrt(N), so photonCounts is applied as a
%       multiplier after the fact rather than simulated.
%   N is photons DETECTED ON THE ARRAY, since p_k is conditioned on the photon
%   landing there, and it is supplied by the caller rather than derived from a
%   brightness model. res.relativeCaptureEfficiency shows how the budget varies
%   across the field, so the achievable precision at r is
%   sigma1(r)/sqrt(N*relativeCapture(r)) rather than a uniform sigma1/sqrt(N).
%
%   This is also why only the collection path is modelled. For focused ISM the
%   counts are mu_k = N0*h_exc(r)*eta_k(r); the excitation term is common to all
%   channels, so it is absorbed entirely into the free amplitude A = N0*h_exc(r)
%   and cancels in p_k. An aberrated excitation focus therefore changes the
%   photon budget but CANNOT bias the fitted position, so long as amplitude is
%   free per frame -- the convention used by ismDetectorChannelProbability. If
%   brightness were instead held fixed across frames, the illumination gradient
%   would carry position information and excitation aberration would contribute
%   its own bias; that case is outside this model.
%
%   The 'No aberration' scenario at zero defocus is an exact self-test: true and
%   fitted models coincide, so the recovered bias must vanish to solver
%   tolerance. At nonzero defocus it is physics, not a control, and the reason
%   is worth stating because "a symmetric blur cannot shift a symmetric spot"
%   is the natural but wrong intuition. In a continuous, untruncated detector
%   the KL projection of a symmetric-about-r distribution onto a translation
%   family does sit at s = r, and there is no bias. The finite array breaks
%   that, and the mechanism is FLATNESS, not centroid shift:
%     - Defocus spreads the microimage, so the true channel distribution
%       p^true(r,z) is flatter across the 23 channels than the in-focus one.
%     - The ideal template's own flatness is controlled by |s|: the further the
%       emitter is from the array centre, the more of the microimage falls off
%       the array and the more even the surviving channels become. In the limit
%       |s| -> infinity every channel sees the same far tail and p^ideal -> 1/nDet.
%     - Matching a flatter pattern therefore pushes the fit OUTWARD: |s| > |r|.
%   The bias field is a radial EXPANSION, b ~ +(g-1) r with g > 1, so J ~ g*I
%   and the apparent D is inflated by roughly g^2 > 1. A perfect wavefront at
%   the wrong stage plane already reports a diffusion coefficient that is too
%   large, and it does so with zero bias at r = 0, which is why the centre
%   diagnostic reports bias 0 and D_app/D > 1 simultaneously.
%
%   Pushed far enough this stops being a bias and becomes a non-existence
%   result. Once the defocus exceeds roughly the depth of field lamEm/NA^2, the
%   true microimage is flatter than ANY finite-s ideal microimage, the
%   cross-entropy increases monotonically outward, and there is no interior
%   maximum: no pseudo-true position exists. The solver detects this (the
%   global argmax lands on the rim of the search window), marks those points
%   NaN and reports fractionValid, rather than returning the arbitrary value at
%   whichever radius the domain happened to be cut off.
%
%   KEY NAME-VALUE OPTIONS
%     'detectorLayout'   'honeycomb23' (default) or 'honeycomb93'
%     'aberrationModes'  modes to evaluate (default
%                        {'astig_x','coma_x','trefoil_x','spherical'}).
%                        tilt_x/tilt_y are excluded by default: they are a pure
%                        translation of the PSF and hence perfectly degenerate
%                        with position, giving a constant, unobservable shift.
%     'amplitudeWaves'   RMS-wave amplitudes at sim.lamRef (default [.05 .1 .15])
%     'includeNoAberration' add the zero-wavefront reference (default true)
%     'defocusUm'        stage planes [um] (default [0 0.25 0.5])
%     'evalHalfWidthUm'  half-width of the diffusion region (default 0.15)
%     'evalStride'       decimation of the fine grid for evaluation (default 1)
%     'photonCounts'     photon budgets for the error budget (default 10.^(2:5))
%     Optics: 'NA' (1.45), 'lamEm' (0.520), 'sampleGeometry' ('airOnGlass'),
%       'emitterHeightUm' (0), 'nImmersion' (1.51), 'nGlass' (1.518),
%       'nSample' (1.33), 'coverslipThicknessUm' (170), 'nMedium' (1.33),
%       'arrayDiameterAU' (1.7), 'detPitchUm' ([]), 'detFillRatio' (1.0),
%       'detectorSubsamples' (7), 'fovXY' (2.6), 'nx' (121)
%       NOTE the default nSample is 1.33, so 'airOnGlass' is really running as
%       a water-on-glass interface. Pass 'nSample',1.0 for an actual air
%       interface; that puts most of the NA 1.45 pupil beyond the critical
%       angle and makes the result strongly dependent on emitterHeightUm.
%     Numerics: 'maxNewtonIterations' (60), 'stepToleranceUm' (1e-7),
%       'minProbability' (1e-12), 'seedSearchHalfWidthUm' (0.40)
%     Output: 'makeFigure' (true), 'errorBudgetIncludeNoAberration' (false),
%       'saveFigures' (false), 'writeOutputs' (false), 'outputDir' (''),
%       'verbose' (true)
%
%   OUTPUT: see field comments at the end of the function.

    addpath(fileparts(mfilename('fullpath')));
    opts = parseInputs(varargin{:});

    sim = configureSim(opts);
    grid = buildEvaluationGrid(sim, opts);
    if isempty(opts.escapeRadiusUm)
        opts.escapeRadiusUm = grid.arrayRadiusUm;
    end
    depthOfFieldUm = sim.lamEm / sim.NA^2;
    if max(abs(opts.defocusUm)) > depthOfFieldUm
        warning('simulateAberrationLocalizationBias:DefocusBeyondDepthOfField', ...
            ['Requested defocus reaches %.2f um against a depth of field of ' ...
             '%.2f um. Past roughly one depth of field the microimage overfills ' ...
             'the array and its channel distribution flattens toward uniform. ' ...
             'The ideal template also tends to uniform, but only as |s| -> inf, ' ...
             'so the cross-entropy has no interior maximum and the fit runs to ' ...
             'the edge of the search domain. Those points are reported as ' ...
             'invalid (NaN), not as a large bias.'], ...
            max(abs(opts.defocusUm)), depthOfFieldUm);
    end

    idealModel = buildChannelModel(sim, struct(), 0, opts);

    scenarios = buildScenarios(opts);      % {label, coeffStruct, modeName, amplitude}
    nScen = size(scenarios, 1);
    defocus = opts.defocusUm(:).';
    nDefocus = numel(defocus);

    results = cell(nScen, nDefocus);
    for s = 1:nScen
        for d = 1:nDefocus
            trueModel = buildChannelModel(sim, scenarios{s,2}, defocus(d), opts);
            results{s,d} = analyzeScenario(trueModel, idealModel, grid, opts);
            if opts.verbose
                reportScenario(scenarios{s,1}, defocus(d), results{s,d}, grid);
            end
        end
    end

    resultsTable = buildTable(scenarios, defocus, results, opts, sim);

    figs = struct('biasField', [], 'errorBudget', [], ...
        'apparentDiffusion', [], 'radialProfile', []);
    if opts.makeFigure
        figs.biasField = plotBiasField(scenarios, defocus, results, grid, sim);
        figs.errorBudget = plotErrorBudget(scenarios, defocus, results, opts);
        figs.apparentDiffusion = plotApparentDiffusion(scenarios, defocus, results, opts);
        figs.radialProfile = plotRadialProfile(scenarios, defocus, results, grid);
    end

    if opts.writeOutputs || opts.saveFigures
        ensureOutputDir(opts.outputDir);
        if opts.writeOutputs
            writetable(resultsTable, fullfile(opts.outputDir, ...
                sprintf('%s_localization_bias.csv', sim.detectorLayout)));
        end
        if opts.saveFigures
            saveFigureSet(figs, opts.outputDir, sim.detectorLayout);
        end
    end

    out = struct();
    out.sim               = sim;             % configured array sim
    out.detectorLayout    = sim.detectorLayout;
    out.grid              = grid;            % evaluation grid (subset of sim grid)
    out.defocusUm         = defocus;         % stage planes evaluated
    out.scenarioLabels    = scenarios(:,1);  % aberration operating points
    out.coeffScenarios    = scenarios(:,2);
    out.scenarioModes     = scenarios(:,3);
    out.scenarioAmplitudes = cell2mat(scenarios(:,4));
    out.results           = results;         % [nScen x nDefocus] per-point maps
    out.photonCounts      = opts.photonCounts;
    out.table             = resultsTable;    % tidy summary, one row per scenario/defocus
    out.figures           = figs;
    out.idealModel        = idealModel;
    out.options           = opts;
end

% ============================================================================
function opts = parseInputs(varargin)
    p = inputParser;
    p.FunctionName = 'simulateAberrationLocalizationBias';

    addParameter(p, 'detectorLayout', 'honeycomb23');

    % Optics / geometry (matches simulateArrayFisherInformation defaults).
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

    % Aberration operating points.
    addParameter(p, 'aberrationModes', {'astig_x','coma_x','trefoil_x','spherical'});
    addParameter(p, 'amplitudeWaves', [0.05 0.10 0.15]);
    addParameter(p, 'includeNoAberration', true);
    % Defocus planes. The old default [0 0.25 0.5] is past the point where an
    % in-focus template can localize at all: at NA 1.45 the depth of field is
    % about lamEm/NA^2 = 0.25 um, and by 0.5 um the microimage overfills the
    % 1.7 AU array, flattens, and the misspecified MLE has no interior optimum.
    % Those planes produce divergent, not merely biased, estimates.
    addParameter(p, 'defocusUm', [0 0.06 0.12]);

    % Evaluation region. The array is only arrayDiameterAU Airy units across,
    % so the default half-width stays well inside the footprint; near the rim
    % the microimage falls off the array and localization degenerates.
    addParameter(p, 'evalHalfWidthUm', 0.15);
    addParameter(p, 'evalStride', 1);
    addParameter(p, 'radialBinUm', 0.05);

    % Error budget.
    addParameter(p, 'photonCounts', 10.^(2:5));

    % Numerics.
    addParameter(p, 'maxNewtonIterations', 60);
    addParameter(p, 'stepToleranceUm', 1e-7);
    addParameter(p, 'minProbability', 1e-12);
    % How far outside the evaluation region the coarse global search for the
    % cross-entropy maximum looks. Must comfortably exceed the largest bias
    % expected, otherwise the global maximiser can sit outside the search set.
    addParameter(p, 'seedSearchHalfWidthUm', 0.40);
    % A fit that lands further than this from the truth has not estimated a
    % position, it has run off. Empty means one array radius.
    addParameter(p, 'escapeRadiusUm', []);

    % Output.
    addParameter(p, 'makeFigure', true);
    addParameter(p, 'errorBudgetIncludeNoAberration', false);
    addParameter(p, 'saveFigures', false);
    addParameter(p, 'writeOutputs', false);
    addParameter(p, 'outputDir', '');
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;

    opts.detectorLayout = char(opts.detectorLayout);
    if mod(opts.nx, 2) == 0, opts.nx = opts.nx + 1; end
    if ischar(opts.aberrationModes) || isstring(opts.aberrationModes)
        opts.aberrationModes = cellstr(opts.aberrationModes);
    end
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
function sim = configureSim(opts)
%CONFIGURESIM Emitter-on-glass array sim, mirroring simulateArrayFisherInformation.
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
function grid = buildEvaluationGrid(sim, opts)
%BUILDEVALUATIONGRID Decimated subset of the fine grid covering the diffusion region.
%   Keeping the evaluation points on fine-grid nodes makes the true-model
%   probabilities and their gradients exact array lookups; only the ideal model,
%   queried at off-grid trial positions, needs interpolation.
    inside = find(abs(sim.x) <= opts.evalHalfWidthUm + 1e-12);
    if numel(inside) < 3
        error('simulateAberrationLocalizationBias:EvalRegionTooSmall', ...
            'evalHalfWidthUm is smaller than the grid pitch; increase it or nx.');
    end
    idx = inside(1:opts.evalStride:end);
    if numel(idx) < 3
        idx = inside;
    end

    grid = struct();
    grid.indexX = idx;
    grid.indexY = idx;
    grid.x = sim.x(idx);
    grid.y = sim.y(idx);
    grid.nX = numel(idx);
    grid.nY = numel(idx);
    [grid.X, grid.Y] = meshgrid(grid.x, grid.y);
    grid.R = hypot(grid.X, grid.Y);
    grid.nPoints = numel(grid.X);
    grid.spacingUm = abs(grid.x(2) - grid.x(1));
    [~, grid.centerIndexX] = min(abs(grid.x));
    [~, grid.centerIndexY] = min(abs(grid.y));

    grid.airyUnitUm = 1.22 * sim.lamEm / sim.NA;
    span = max( ...
        (max(sim.detXY(:,1)) - min(sim.detXY(:,1))) + sim.detSize, ...
        (max(sim.detXY(:,2)) - min(sim.detXY(:,2))) + 2*sim.detectorHexRadius);
    grid.arrayRadiusUm = span / 2;
    grid.fractionOfArrayRadius = max(grid.x) / grid.arrayRadiusUm;

    if grid.fractionOfArrayRadius > 0.7
        warning('simulateAberrationLocalizationBias:EvalRegionNearArrayEdge', ...
            ['Evaluation half-width %.3f um reaches %.0f%% of the array radius ' ...
             '(%.3f um). Near the rim the microimage spills off the array, so ' ...
             'bias and distortion there reflect vanishing position information ' ...
             'rather than the aberration.'], ...
            max(grid.x), 100*grid.fractionOfArrayRadius, grid.arrayRadiusUm);
    end

    grid.radialEdges = 0:opts.radialBinUm:(max(grid.R(:)) + opts.radialBinUm);
end

% ----------------------------------------------------------------------------
function model = buildChannelModel(sim, coeffs, planeZ, opts)
%BUILDCHANNELMODEL Normalized microimage p_k(x,y) and its lateral gradient.
%   p is indexed (y, x, channel) on the fine simulation grid.
    coeffs = coeffStruct(sim, coeffs);
    hDet = arrayDetectionPSF(sim, coeffs, planeZ);
    eta = detectorCollectionEfficiencyExplicit(sim, hDet);     % [ny nx 1 nDet]
    eta = max(reshape(eta, numel(sim.y), numel(sim.x), sim.nDet), 0);

    total = sum(eta, 3);
    total(total < opts.minProbability) = opts.minProbability;
    prob = eta ./ total;

    dPdx = zeros(size(prob));
    dPdy = zeros(size(prob));
    for k = 1:sim.nDet
        [dPdx(:,:,k), dPdy(:,:,k)] = gradient(prob(:,:,k), sim.dx, sim.dx);
    end

    model = struct();
    model.probability = prob;
    model.dPdx = dPdx;
    model.dPdy = dPdy;
    % Total array response, proportional to the flux captured by the whole
    % array. Normalizing it away is what makes the bias brightness-independent.
    % NOTE this is NOT an absolute fraction: detectorCollectionEfficiencyExplicit
    % returns a mean intensity per pixel, so the scale carries the PSF
    % normalization and omits the pixel-area factor. Use the centre-normalized
    % relativeCaptureEfficiency, which is free of both, to see how the photon
    % budget varies across the field.
    model.captureEfficiency = sum(eta, 3);
    model.x = sim.x;
    model.y = sim.y;
    model.nDet = sim.nDet;
    model.planeZ = planeZ;
    model.coeffs = coeffs;
    % Step used for central differences of the cubic interpolant when this
    % model is queried off-grid (see interpolateModel). A fraction of the grid
    % pitch keeps the difference inside one interpolation cell, so it
    % differentiates the interpolant rather than the sampled data.
    model.derivativeStepUm = 0.25 * sim.dx;
end

% ----------------------------------------------------------------------------
function h = arrayDetectionPSF(sim, coeffs, planeZ)
%ARRAYDETECTIONPSF Collection-path detection PSF at the given stage plane.
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
function res = analyzeScenario(trueModel, idealModel, grid, opts)
%ANALYZESCENARIO Bias field, distortion Jacobian and error budget for one case.
    pTrue = sampleGridPoints(trueModel.probability, grid);      % [nPts x nDet]

    fitted = solvePseudoTruePositions(pTrue, idealModel, grid, opts);

    biasX = reshape(fitted.x - grid.X(:), grid.nY, grid.nX);
    biasY = reshape(fitted.y - grid.Y(:), grid.nY, grid.nX);

    % NO INTERIOR OPTIMUM. As |s| grows the ideal microimage approaches the
    % uniform distribution: far from the array every channel sees the same far
    % tail of the PSF, so p_k^ideal -> 1/nDet. A strongly defocused true
    % microimage overfills the array and is ALSO nearly uniform. The
    % cross-entropy is then maximised by running s away to infinity, not by any
    % finite position, and the "fit" simply reports wherever the search domain
    % was cut off. Those points carry no position estimate at all, so they are
    % marked invalid rather than allowed into the bias field, where a single
    % escaped point would poison the Jacobian of its whole neighbourhood.
    escaped = reshape(fitted.escaped | ...
        hypot(fitted.x - grid.X(:), fitted.y - grid.Y(:)) > opts.escapeRadiusUm, ...
        grid.nY, grid.nX);
    biasX(escaped) = NaN;
    biasY(escaped) = NaN;
    biasMagnitude = hypot(biasX, biasY);
    valid = ~escaped;

    % A uniform bias is only an origin shift and is unobservable in tracking;
    % the residual after removing it is what distorts a trajectory.
    meanBias = [mean(biasX(:), 'omitnan'), mean(biasY(:), 'omitnan')];
    residualX = biasX - meanBias(1);
    residualY = biasY - meanBias(2);
    residualMagnitude = hypot(residualX, residualY);

    distortion = distortionMetrics(biasX, biasY, grid);
    fisher = positionFisher(trueModel, grid, opts);

    res = struct();
    res.biasX = biasX;                       % [um] systematic x error vs position
    res.biasY = biasY;
    res.biasMagnitude = biasMagnitude;
    res.escaped = escaped;                   % no interior optimum at these points
    res.fractionValid = mean(valid(:));
    res.meanBiasUm = meanBias;               % unobservable origin shift
    res.residualMagnitude = residualMagnitude;
    res.rmsBiasUm = sqrt(mean(biasMagnitude(:).^2, 'omitnan'));
    res.rmsResidualBiasUm = sqrt(mean(residualMagnitude(:).^2, 'omitnan'));
    res.maxBiasUm = max(biasMagnitude(:), [], 'omitnan');
    res.converged = fitted.converged;
    res.fractionConverged = mean(fitted.converged);
    finiteResidual = fitted.gradientNormUm(isfinite(fitted.gradientNormUm));
    if isempty(finiteResidual)
        res.maxNewtonResidualUm = NaN;
    else
        res.maxNewtonResidualUm = max(finiteResidual);
    end

    % A discontinuous bias field means the solver jumped between branches of a
    % multimodal cross-entropy; its Jacobian, and therefore the apparent
    % diffusion ratio, is then meaningless rather than merely inaccurate.
    stepX = diff(biasX, 1, 2); stepY = diff(biasY, 1, 2);
    stepXv = diff(biasX, 1, 1); stepYv = diff(biasY, 1, 1);
    jumps = [hypot(stepX(:), stepY(:)); hypot(stepXv(:), stepYv(:))];
    res.maxNeighbourBiasJumpUm = max([jumps(~isnan(jumps)); 0]);
    res.biasFieldSmooth = res.maxNeighbourBiasJumpUm < 0.5 * grid.spacingUm;

    res.jacobian = distortion.jacobian;              % [nY nX 2 2]
    res.apparentDiffusionRatio = distortion.ratio;   % D_app / D_true map
    res.anisotropy = distortion.anisotropy;
    res.meanApparentDiffusionRatio = mean(distortion.ratio(:), 'omitnan');
    res.stdApparentDiffusionRatio = std(distortion.ratio(:), 'omitnan');
    res.meanAnisotropy = mean(distortion.anisotropy(:), 'omitnan');

    res.sigmaOnePhotonUm = fisher.sigmaOnePhoton;    % [nY nX] CRB at N = 1
    res.meanSigmaOnePhotonUm = mean(fisher.sigmaOnePhoton(:), 'omitnan');
    res.precisionUm = res.meanSigmaOnePhotonUm ./ sqrt(opts.photonCounts(:).');
    res.crossoverPhotons = (res.meanSigmaOnePhotonUm / max(res.rmsBiasUm, eps))^2;
    res.crossoverPhotonsResidual = ...
        (res.meanSigmaOnePhotonUm / max(res.rmsResidualBiasUm, eps))^2;

    % Field centre, where position information is highest and the numbers are
    % least contaminated by the array rim.
    cy = grid.centerIndexY; cx = grid.centerIndexX;
    res.centerBiasUm = biasMagnitude(cy, cx);
    res.centerSigmaOnePhotonUm = fisher.sigmaOnePhoton(cy, cx);
    res.centerApparentDiffusionRatio = distortion.ratio(cy, cx);
    res.centerCrossoverPhotons = ...
        (res.centerSigmaOnePhotonUm / max(res.centerBiasUm, eps))^2;

    % Array capture fraction. photonCounts is applied as a uniform N, so this
    % map says how optimistic that is away from the centre: the achievable
    % precision at r scales as sigma1(r) / sqrt(N * relativeCapture(r)).
    capture = reshape(sampleGridPoints(trueModel.captureEfficiency, grid), ...
        grid.nY, grid.nX);
    res.captureEfficiency = capture;
    res.centerCaptureEfficiency = capture(cy, cx);
    res.relativeCaptureEfficiency = capture / max(capture(cy, cx), eps);

    res.radial = radialProfile(grid, biasMagnitude, fisher.sigmaOnePhoton, ...
        distortion.ratio);
end

% ----------------------------------------------------------------------------
function profile = radialProfile(grid, biasMagnitude, sigmaOnePhoton, ratio)
%RADIALPROFILE Bin the field metrics by radius so the usable field is visible.
    edges = grid.radialEdges;
    nBin = numel(edges) - 1;
    profile = struct('radiusUm', nan(1, nBin), 'biasUm', nan(1, nBin), ...
        'sigmaOnePhotonUm', nan(1, nBin), 'apparentDiffusionRatio', nan(1, nBin), ...
        'count', zeros(1, nBin));
    for b = 1:nBin
        m = grid.R >= edges(b) & grid.R < edges(b+1);
        profile.count(b) = nnz(m);
        if ~any(m(:)), continue; end
        profile.radiusUm(b) = mean(grid.R(m));
        profile.biasUm(b) = mean(biasMagnitude(m), 'omitnan');
        profile.sigmaOnePhotonUm(b) = mean(sigmaOnePhoton(m), 'omitnan');
        profile.apparentDiffusionRatio(b) = mean(ratio(m), 'omitnan');
    end
end

% ----------------------------------------------------------------------------
function values = sampleGridPoints(field3D, grid)
%SAMPLEGRIDPOINTS Extract [nPoints x nChannel] values at the evaluation nodes.
    sub = field3D(grid.indexY, grid.indexX, :);
    values = reshape(sub, grid.nPoints, size(field3D, 3));
end

% ----------------------------------------------------------------------------
function fitted = solvePseudoTruePositions(pTrue, idealModel, grid, opts)
%SOLVEPSEUDOTRUEPOSITIONS Maximize the cross-entropy sum_k pTrue_k log pIdeal_k(s).
%
%   Profiling the free amplitude out of the Poisson log-likelihood leaves this
%   cross-entropy, so the limit of the misspecified MLE is its maximizer. The
%   Gauss-Newton curvature -sum_k pTrue_k (grad p / p)(grad p / p)' is negative
%   semi-definite by construction, so the ascent direction is always valid.
%
%   The cross-entropy is NOT unimodal once the true PSF is strongly aberrated
%   and defocused: a comatic, defocused microimage has a core and a flare, and
%   the ideal template can lock onto either. Starting Newton at the true
%   position therefore returns whichever local maximum happens to sit in that
%   basin, and the basin boundary sweeps across the field as the emitter moves,
%   putting a step discontinuity in b(r). The theory calls for the GLOBAL
%   maximiser, so a coarse exhaustive scan over the ideal-model nodes seeds the
%   solver and Newton is used only to polish.
    nPts = grid.nPoints;
    active = true(nPts, 1);
    converged = false(nPts, 1);

    xLimit = [min(idealModel.x), max(idealModel.x)];
    yLimit = [min(idealModel.y), max(idealModel.y)];
    margin = 2 * (idealModel.x(2) - idealModel.x(1));

    [posX, posY, seedOnBoundary] = coarseGlobalSeed(pTrue, idealModel, grid, opts);
    objective = evaluateObjective(posX, posY, pTrue, idealModel, opts);
    gradientNorm = inf(nPts, 1);

    % The Newton step is built from a finite-differenced interpolant, so it
    % cannot be driven below the noise floor of that difference. Demanding
    % 1e-7 um from it just guarantees every point exhausts its iterations and
    % gets reported as unconverged. A thousandth of the grid pitch (about 20 pm
    % here) is many orders below any bias worth reporting.
    tolerance = max(opts.stepToleranceUm, 1e-3 * grid.spacingUm);

    for iter = 1:opts.maxNewtonIterations
        if ~any(active), break; end
        rows = find(active);

        [p, dpx, dpy] = interpolateModel(posX(rows), posY(rows), idealModel);
        p = max(p, opts.minProbability);
        w = pTrue(rows, :);

        gx = sum(w .* dpx ./ p, 2);
        gy = sum(w .* dpy ./ p, 2);
        hxx = sum(w .* (dpx ./ p).^2, 2);
        hyy = sum(w .* (dpy ./ p).^2, 2);
        hxy = sum(w .* (dpx ./ p) .* (dpy ./ p), 2);

        det2 = hxx .* hyy - hxy.^2;
        ridge = 1e-12 * max(hxx + hyy, eps);
        bad = ~(det2 > 0);
        hxx(bad) = hxx(bad) + ridge(bad);
        hyy(bad) = hyy(bad) + ridge(bad);
        det2 = hxx .* hyy - hxy.^2;
        det2(det2 <= 0) = eps;

        stepX = ( hyy .* gx - hxy .* gy) ./ det2;
        stepY = (-hxy .* gx + hxx .* gy) ./ det2;

        % Newton residual in position units: the length of the step the
        % quadratic model asks for. This is the honest stationarity test.
        gradientNorm(rows) = hypot(stepX, stepY);

        [posX, posY, objective, moved] = backtrackingAscent(rows, posX, posY, ...
            stepX, stepY, objective, pTrue, idealModel, opts, ...
            xLimit, yLimit, margin);

        % A point is converged only if the step the solver WANTED is small.
        % moved < tolerance alone also fires when the line search failed 12
        % times and the point simply stalled, which silently reported stalled
        % points as converged and let their (wrong) positions into the bias
        % field.
        settled = gradientNorm(rows) < tolerance;
        stalled = ~settled & (moved < opts.stepToleranceUm);
        converged(rows(settled)) = true;
        active(rows(settled | stalled)) = false;
    end

    fitted = struct('x', posX, 'y', posY, 'converged', converged, ...
        'escaped', seedOnBoundary, 'gradientNormUm', gradientNorm, ...
        'objective', objective);
end

% ----------------------------------------------------------------------------
function [seedX, seedY, onBoundary] = coarseGlobalSeed(pTrue, idealModel, grid, opts)
%COARSEGLOBALSEED Exhaustive argmax of the cross-entropy over ideal-model nodes.
%   Vectorised as one matrix product: the cross-entropy is linear in log p, so
%   the whole [nPoint x nNode] objective table is pTrue * logP'. Restricting the
%   node set to a window around the evaluation region keeps that table small.
%
%   onBoundary marks points whose best node lies on the rim of that window.
%   That is the signature of an objective with NO interior maximum: the ideal
%   microimage tends to the uniform distribution as |s| -> infinity, so if the
%   true microimage is flatter than anything the ideal model produces at finite
%   s, the cross-entropy increases monotonically outward and the reported
%   "position" is just wherever the window was cut. Widening the window does
%   not help; the estimator genuinely does not exist for those points.
    window = max(grid.x(end), grid.y(end)) + opts.seedSearchHalfWidthUm;
    ix = find(abs(idealModel.x) <= window);
    iy = find(abs(idealModel.y) <= window);
    if numel(ix) < 3 || numel(iy) < 3
        seedX = grid.X(:); seedY = grid.Y(:);
        onBoundary = false(grid.nPoints, 1);
        return;
    end

    logP = zeros(numel(iy)*numel(ix), idealModel.nDet);
    for k = 1:idealModel.nDet
        block = idealModel.probability(iy, ix, k);
        logP(:,k) = log(max(block(:), opts.minProbability));
    end

    score = pTrue * logP.';                 % [nPoint x nNode]
    [~, best] = max(score, [], 2);

    [NX, NY] = meshgrid(idealModel.x(ix), idealModel.y(iy));
    seedX = NX(best);  seedX = seedX(:);
    seedY = NY(best);  seedY = seedY(:);

    nRow = numel(iy);  nCol = numel(ix);
    [rowOfBest, colOfBest] = ind2sub([nRow nCol], best(:));
    onBoundary = rowOfBest <= 1 | rowOfBest >= nRow | ...
                 colOfBest <= 1 | colOfBest >= nCol;
end

% ----------------------------------------------------------------------------
function [posX, posY, objective, moved] = backtrackingAscent(rows, posX, posY, ...
        stepX, stepY, objective, pTrue, idealModel, opts, xLimit, yLimit, margin)
%BACKTRACKINGASCENT Halve the Newton step per point until the objective improves.
    scale = ones(numel(rows), 1);
    accepted = false(numel(rows), 1);
    bestX = posX(rows);
    bestY = posY(rows);
    bestObjective = objective(rows);

    for trial = 1:12
        pending = find(~accepted);
        if isempty(pending), break; end

        candX = clampTo(posX(rows(pending)) + scale(pending) .* stepX(pending), ...
            xLimit, margin);
        candY = clampTo(posY(rows(pending)) + scale(pending) .* stepY(pending), ...
            yLimit, margin);
        candObjective = evaluateObjective(candX, candY, pTrue(rows(pending), :), ...
            idealModel, opts);

        better = candObjective > bestObjective(pending);
        bestX(pending(better)) = candX(better);
        bestY(pending(better)) = candY(better);
        bestObjective(pending(better)) = candObjective(better);
        accepted(pending(better)) = true;
        scale(pending(~better)) = scale(pending(~better)) / 2;
    end

    moved = hypot(bestX - posX(rows), bestY - posY(rows));
    posX(rows) = bestX;
    posY(rows) = bestY;
    objective(rows) = bestObjective;
end

% ----------------------------------------------------------------------------
function value = evaluateObjective(qx, qy, pTrue, idealModel, opts)
    p = max(interpolateModel(qx, qy, idealModel), opts.minProbability);
    value = sum(pTrue .* log(p), 2);
end

% ----------------------------------------------------------------------------
function [p, dpx, dpy] = interpolateModel(qx, qy, model)
%INTERPOLATEMODEL Cubic interpolation of the ideal microimage at trial positions.
%
%   The derivatives are central differences OF THE SAME INTERPOLANT that
%   supplies p, not a separate interpolation of a pre-differentiated array.
%   That distinction matters: the objective the line search accepts and the
%   gradient the Newton step uses must belong to one function, otherwise the
%   solver stops where the (inconsistent) gradient vanishes rather than at a
%   stationary point of the objective. The mismatch is O(h^2) but it is locked
%   to the grid, so the bias field picks up a ripple with the grid period, and
%   differentiating that ripple in distortionMetrics multiplies it by 1/h and
%   corrupts the apparent-diffusion ratio.
    nQ = numel(qx);
    p = interpolateProbability(qx, qy, model);
    if nargout <= 1
        return;
    end
    delta = model.derivativeStepUm;
    pxPlus  = interpolateProbability(qx + delta, qy, model);
    pxMinus = interpolateProbability(qx - delta, qy, model);
    pyPlus  = interpolateProbability(qx, qy + delta, model);
    pyMinus = interpolateProbability(qx, qy - delta, model);
    dpx = (pxPlus - pxMinus) / (2*delta);
    dpy = (pyPlus - pyMinus) / (2*delta);
    if nQ == 0
        dpx = zeros(0, model.nDet);
        dpy = zeros(0, model.nDet);
    end
end

% ----------------------------------------------------------------------------
function p = interpolateProbability(qx, qy, model)
    p = zeros(numel(qx), model.nDet);
    for k = 1:model.nDet
        p(:,k) = interp2(model.x, model.y, model.probability(:,:,k), ...
            qx, qy, 'cubic', 0);
    end
end

% ----------------------------------------------------------------------------
function v = clampTo(v, limits, margin)
    v = min(max(v, limits(1) + margin), limits(2) - margin);
end

% ----------------------------------------------------------------------------
function distortion = distortionMetrics(biasX, biasY, grid)
%DISTORTIONMETRICS Jacobian of the apparent-position map and its diffusion effect.
%   Observed steps are dr_obs = J dr_true with J = I + db/dr, so an isotropic
%   true D is measured as the tensor D*J*J'.
    h = grid.spacingUm;
    [dbxdx, dbxdy] = gradient(biasX, h, h);
    [dbydx, dbydy] = gradient(biasY, h, h);

    jacobian = zeros(grid.nY, grid.nX, 2, 2);
    jacobian(:,:,1,1) = 1 + dbxdx;
    jacobian(:,:,1,2) = dbxdy;
    jacobian(:,:,2,1) = dbydx;
    jacobian(:,:,2,2) = 1 + dbydy;

    a = jacobian(:,:,1,1); b = jacobian(:,:,1,2);
    c = jacobian(:,:,2,1); d = jacobian(:,:,2,2);

    % Symmetric J*J', whose eigenvalues scale the apparent diffusion tensor.
    m11 = a.^2 + b.^2;
    m22 = c.^2 + d.^2;
    m12 = a.*c + b.*d;

    traceM = m11 + m22;
    detM = m11 .* m22 - m12.^2;
    disc = sqrt(max(traceM.^2 / 4 - detM, 0));
    lambdaMax = traceM / 2 + disc;
    lambdaMin = max(traceM / 2 - disc, 0);

    distortion = struct();
    distortion.jacobian = jacobian;
    distortion.ratio = traceM / 2;                                  % scalar D_app/D_true
    distortion.anisotropy = sqrt(lambdaMax ./ max(lambdaMin, eps));
end

% ----------------------------------------------------------------------------
function fisher = positionFisher(trueModel, grid, opts)
%POSITIONFISHER Per-photon multinomial information for (x,y) at each grid point.
%   With free amplitude and no background, sum_k p_k = 1 makes the amplitude
%   score orthogonal to the position scores, so the position information is
%   exactly I1 = sum_k (grad p_k)(grad p_k)'/p_k and sigma(N) = sigma1/sqrt(N).
    p = max(sampleGridPoints(trueModel.probability, grid), opts.minProbability);
    dpx = sampleGridPoints(trueModel.dPdx, grid);
    dpy = sampleGridPoints(trueModel.dPdy, grid);

    ixx = sum(dpx.^2 ./ p, 2);
    iyy = sum(dpy.^2 ./ p, 2);
    ixy = sum(dpx .* dpy ./ p, 2);

    det2 = ixx .* iyy - ixy.^2;
    varianceTrace = nan(size(det2));
    valid = det2 > 0;
    varianceTrace(valid) = (ixx(valid) + iyy(valid)) ./ det2(valid);

    fisher = struct();
    fisher.sigmaOnePhoton = reshape(sqrt(varianceTrace / 2), grid.nY, grid.nX);
    fisher.informationXX = reshape(ixx, grid.nY, grid.nX);
    fisher.informationYY = reshape(iyy, grid.nY, grid.nX);
    fisher.informationXY = reshape(ixy, grid.nY, grid.nX);
end

% ----------------------------------------------------------------------------
function scen = buildScenarios(opts)
    scen = {};
    if opts.includeNoAberration
        scen(end+1, :) = {'No aberration', struct(), 'none', 0}; %#ok<AGROW>
    end
    for k = 1:numel(opts.aberrationModes)
        mode = opts.aberrationModes{k};
        for a = 1:numel(opts.amplitudeWaves)
            amplitude = opts.amplitudeWaves(a);
            scen(end+1, :) = { ...
                sprintf('%s = %.2f\\lambda', prettyMode(mode), amplitude), ...
                struct(mode, amplitude), mode, amplitude}; %#ok<AGROW>
        end
    end
    if isempty(scen)
        scen = {'No aberration', struct(), 'none', 0};
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
function T = buildTable(scenarios, defocus, results, opts, sim)
    rows = {};
    for s = 1:size(scenarios, 1)
        plain = regexprep(scenarios{s,1}, '\\[a-zA-Z]+|[{}]', '');
        for d = 1:numel(defocus)
            r = results{s,d};
            rows(end+1, :) = { sim.detectorLayout, plain, scenarios{s,3}, ...
                scenarios{s,4}, defocus(d), ...
                1000 * r.centerBiasUm, 1000 * r.rmsBiasUm, 1000 * r.maxBiasUm, ...
                1000 * r.rmsResidualBiasUm, 1000 * norm(r.meanBiasUm), ...
                r.centerApparentDiffusionRatio, r.meanApparentDiffusionRatio, ...
                r.stdApparentDiffusionRatio, r.meanAnisotropy, ...
                1000 * r.centerSigmaOnePhotonUm, 1000 * r.meanSigmaOnePhotonUm, ...
                r.centerCrossoverPhotons, r.crossoverPhotons, ...
                r.fractionConverged, r.fractionValid }; %#ok<AGROW>
        end
    end
    T = cell2table(rows, 'VariableNames', {'layout', 'scenario', 'mode', ...
        'amplitudeWaves', 'defocusUm', 'centerBiasNm', 'rmsBiasNm', 'maxBiasNm', ...
        'rmsResidualBiasNm', 'meanBiasNm', 'centerDiffusionRatio', ...
        'apparentDiffusionRatio', 'apparentDiffusionRatioStd', 'anisotropy', ...
        'centerSigmaOnePhotonNm', 'sigmaOnePhotonNm', ...
        'centerCrossoverPhotons', 'crossoverPhotons', 'fractionConverged', ...
        'fractionValid'});
    T.Properties.UserData = opts;
end

% ----------------------------------------------------------------------------
function reportScenario(label, planeZ, res, grid)
    plain = regexprep(label, '\\[a-zA-Z]+|[{}]', '');
    fprintf(['[bias] %-22s z=%+5.2f um | centre: bias %6.1f nm, D_app/D %.3f, ' ...
        'sigma1 %5.1f nm, N* %.3g | field: RMS bias %6.1f nm, D_app/D %.3f\n'], ...
        plain, planeZ, 1000*res.centerBiasUm, res.centerApparentDiffusionRatio, ...
        1000*res.centerSigmaOnePhotonUm, res.centerCrossoverPhotons, ...
        1000*res.rmsBiasUm, res.meanApparentDiffusionRatio);
    if res.fractionValid < 1
        fprintf(['       ^ NO ESTIMATOR at %.0f%% of the field: the ' ...
            'cross-entropy has no interior maximum there, so those points have ' ...
            'no pseudo-true position at all. Reduce the defocus.\n'], ...
            100*(1 - res.fractionValid));
    end
    if res.fractionValid > 0 && ~res.biasFieldSmooth
        fprintf(['       ^ UNRELIABLE: largest neighbour bias jump %.1f nm ' ...
            'against a %.1f nm grid. The bias field is discontinuous, so the ' ...
            'D_app/D above is a numerical artefact.\n'], ...
            1000*res.maxNeighbourBiasJumpUm, 1000*grid.spacingUm);
    end
end

% ----------------------------------------------------------------------------
function fig = plotBiasField(scenarios, defocus, results, grid, sim)
%PLOTBIASFIELD Bias magnitude map with the displacement field overlaid.
    nScen = size(scenarios, 1);
    nDefocus = numel(defocus);
    fig = figure('Color', 'w', 'Position', [60 60 300*nDefocus+140 260*nScen]);
    tl = tiledlayout(fig, nScen, nDefocus, 'Padding', 'compact', ...
        'TileSpacing', 'compact');

    for s = 1:nScen
        for d = 1:nDefocus
            r = results{s,d};
            ax = nexttile(tl);
            imagesc(ax, grid.x, grid.y, 1000 * r.biasMagnitude);
            set(ax, 'YDir', 'normal');
            axis(ax, 'image');
            hold(ax, 'on');
            quiver(ax, grid.X, grid.Y, r.biasX, r.biasY, 0.9, 'k', 'LineWidth', 0.6);
            colorbar(ax);
            xlabel(ax, 'x [\mum]');
            ylabel(ax, 'y [\mum]');
            title(ax, sprintf('%s, z=%.2f \\mum', scenarios{s,1}, defocus(d)), ...
                'FontWeight', 'normal', 'FontSize', 9);
        end
    end
    title(tl, sprintf(['Localization bias from assuming an ideal PSF (%s), ' ...
        'colour = |bias| [nm]'], sim.detectorLayout), 'FontWeight', 'bold');
end

% ----------------------------------------------------------------------------
function fig = plotErrorBudget(scenarios, defocus, results, opts)
%PLOTERRORBUDGET Photon-independent bias floor against the shot-noise CRB.
%   The zero-wavefront reference is dropped by default. At z = 0 its bias floor
%   is a solver-tolerance number, which on a log axis stretches the y range by
%   many decades and squashes every curve that carries information. Set
%   'errorBudgetIncludeNoAberration' to true to put it back.
    show = 1:size(scenarios, 1);
    if ~opts.errorBudgetIncludeNoAberration
        show = show(~strcmp(scenarios(:,3), 'none'));
    end
    if isempty(show)
        show = 1:size(scenarios, 1);
    end

    nDefocus = numel(defocus);
    nCol = min(nDefocus, 3);
    nRow = ceil(nDefocus / nCol);
    fig = figure('Color', 'w', 'Position', [70 70 380*nCol 320*nRow]);
    tl = tiledlayout(fig, nRow, nCol, 'Padding', 'compact', 'TileSpacing', 'compact');
    cmap = lines(size(scenarios, 1));
    photons = opts.photonCounts(:).';

    lastAx = [];
    for d = 1:nDefocus
        ax = nexttile(tl);
        hold(ax, 'on');
        for s = show
            r = results{s,d};
            plot(ax, photons, 1000 * r.meanSigmaOnePhotonUm ./ sqrt(photons), ...
                '--', 'Color', cmap(s,:), 'LineWidth', 1);
            plot(ax, photons, repmat(1000 * r.rmsBiasUm, size(photons)), ...
                '-', 'Color', cmap(s,:), 'LineWidth', 1.4);
        end
        set(ax, 'XScale', 'log', 'YScale', 'log');
        grid(ax, 'on');
        xlabel(ax, 'detected photons N');
        ylabel(ax, 'localization error [nm]');
        title(ax, sprintf('z = %.2f \\mum', defocus(d)), 'FontWeight', 'normal');
        lastAx = ax;
    end

    if ~isempty(lastAx)
        lgd = legend(lastAx, repelem(scenarios(show,1), 2), 'FontSize', 7);
        lgd.Layout.Tile = 'east';
    end
    title(tl, ['Error budget: solid = aberration bias floor (photon ' ...
        'independent), dashed = shot-noise CRB'], 'FontWeight', 'bold');
end

% ----------------------------------------------------------------------------
function fig = plotApparentDiffusion(scenarios, defocus, results, opts)
%PLOTAPPARENTDIFFUSION Measured-to-true diffusion ratio versus aberration amplitude.
    modes = unique(scenarios(:,3), 'stable');
    modes = modes(~strcmp(modes, 'none'));
    fig = figure('Color', 'w', 'Position', [80 80 420 340]);
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on');
    cmap = lines(max(numel(modes), 1));
    markers = {'o', 's', '^', 'd', 'v', '>'};

    for m = 1:numel(modes)
        for d = 1:numel(defocus)
            sel = find(strcmp(scenarios(:,3), modes{m}));
            amplitude = cell2mat(scenarios(sel, 4));
            ratio = arrayfun(@(k) results{sel(k), d}.meanApparentDiffusionRatio, ...
                1:numel(sel));
            plot(ax, amplitude, ratio, ['-' markers{min(d, numel(markers))}], ...
                'Color', cmap(m,:), 'MarkerFaceColor', cmap(m,:), ...
                'MarkerSize', 4, 'LineWidth', 1, ...
                'DisplayName', sprintf('%s, z=%.2f', prettyMode(modes{m}), defocus(d)));
        end
    end
    yline(ax, 1, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    grid(ax, 'on');
    xlabel(ax, 'aberration amplitude [RMS waves]');
    ylabel(ax, 'apparent D / true D');
    title(ax, 'Diffusion bias from an unknown wavefront', 'FontWeight', 'normal');
    legend(ax, 'Location', 'best', 'FontSize', 7);
end

% ----------------------------------------------------------------------------
function fig = plotRadialProfile(scenarios, defocus, results, evalGrid)
%PLOTRADIALPROFILE Bias and precision versus radius, which sets the usable field.
    nScen = size(scenarios, 1);
    fig = figure('Color', 'w', 'Position', [90 90 900 300]);
    tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    cmap = lines(nScen);
    d = 1;   % first defocus plane

    fields = {'biasUm', 'sigmaOnePhotonUm', 'apparentDiffusionRatio'};
    labels = {'|bias| [nm]', '\sigma_1 [nm]', 'apparent D / true D'};
    scales = [1000, 1000, 1];

    lastAx = [];
    for f = 1:3
        ax = nexttile(tl);
        hold(ax, 'on');
        for s = 1:nScen
            profile = results{s,d}.radial;
            plot(ax, 1000 * profile.radiusUm, scales(f) * profile.(fields{f}), ...
                '-o', 'Color', cmap(s,:), 'MarkerFaceColor', cmap(s,:), ...
                'MarkerSize', 3, 'LineWidth', 1, 'DisplayName', scenarios{s,1});
        end
        grid(ax, 'on');
        xlabel(ax, 'radius from array centre [nm]');
        ylabel(ax, labels{f});
        if f == 3
            yline(ax, 1, 'k:', 'HandleVisibility', 'off');
        end
        lastAx = ax;
    end
    title(tl, sprintf(['Radial dependence at z = %.2f \\mum ' ...
        '(array radius %.0f nm, field covers %.0f%%)'], defocus(d), ...
        1000*evalGrid.arrayRadiusUm, 100*evalGrid.fractionOfArrayRadius), ...
        'FontWeight', 'bold');
    if ~isempty(lastAx)
        lgd = legend(lastAx, 'FontSize', 7);
        lgd.Layout.Tile = 'east';
    end
end

% ----------------------------------------------------------------------------
function saveFigureSet(figs, outputDir, layout)
    names = fieldnames(figs);
    for k = 1:numel(names)
        f = figs.(names{k});
        if isempty(f) || ~isgraphics(f), continue; end
        exportgraphics(f, fullfile(outputDir, ...
            sprintf('%s_%s.png', layout, names{k})), 'Resolution', 180);
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
            error('simulateAberrationLocalizationBias:BadSampleGeometry', ...
                ['sampleGeometry must be ''airOnGlass'' (emitter on the ' ...
                'coverslip) or ''homogeneous''.']);
    end
end

% ----------------------------------------------------------------------------
function ensureOutputDir(d)
    if ~isempty(d) && exist(d, 'dir') ~= 7
        mkdir(d);
    end
end
