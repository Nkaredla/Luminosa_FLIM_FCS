function out = flim_bayes_fixed_slb(tcspcPix, irf, pulsePeriodNs, dtNs, tauSlbNs, opts)
%FLIM_BAYES_FIXED_SLB Bayesian FLIM with a fixed short-lifetime SLB state.
%
% out = flim_bayes_fixed_slb(tcspcPix, irf, pulsePeriodNs, dtNs, ...
%     tauSlbNs, opts)
%
% The function compares up to three models independently in every analysed
% pixel:
%   1. fixed SLB only
%   2. fixed SLB + one free membrane lifetime (biexponential)
%   3. fixed SLB + two ordered free membrane lifetimes (triexponential)
%
% opts.maxMembraneStates truncates that list. With the default 2 all three
% models are evaluated. With 1 the triexponential model is never built or
% evaluated: its grid costs nothing, its posterior probability is exactly
% zero, and it cannot contribute to any model-averaged map. Output arrays
% keep their three-model shape so downstream code is unaffected.
%
% tcspcPix may be [nx ny time] or [nx ny time detector]. Detector data are
% summed before inference. The fixed SLB lifetime must be estimated from an
% outside-cell SLB reference decay before calling this function. Its
% amplitude is normally inferred separately per pixel. When
% opts.fixedSlbPhotonCount is supplied, it may be a scalar or an [nx ny]
% spatial calibration. With slbCountPriorNodes == 0 it is a hard expected
% count; with positive prior nodes the count is marginalised over the
% corresponding fixedSlbPhotonCountStd calibration uncertainty.
%
% The likelihood is the TCSPC-shape (multinomial) likelihood conditional on
% the observed photon total. Discrete integration over lifetime, component
% fraction and optional uniform-background grids gives per-pixel marginal
% likelihoods and posterior model probabilities. Equal model priors are
% used by default. These probabilities are conditional on the supplied
% grids and priors; they are not frequentist confidence levels.
%
% Important opts fields (all optional)
%   analysisMask                 [nx ny] logical mask (default all pixels)
%   minPhotons                   minimum photons per fitted pixel (default 1)
%   useGPU                       true/false; [] selects automatically
%   batchSize                    fitted pixels per likelihood batch (128)
%   includeBackground            include a uniform background state (true);
%                                required in fixed-count mode
%   signalGrid                   fluorescence photon fraction grid
%                                (default linspace(0.25,1,4)); in fixed-count
%                                mode, membrane share of residual photons
%   membraneTauGrid              explicit longer-lifetime grid in ns
%   membraneTauBoundsNs          bounds used to generate the default grid
%   membraneTauCount             default log-grid size (12)
%   fractionStep                 fluorescence mixture lattice step (0.125)
%   minimumMembraneFraction      each added state must carry at least this
%                                fluorescence fraction (0.05); in fixed-count
%                                M3, minimum relative membrane-mixture share
%   minimumSlbFraction           minimum fixed-SLB fluorescence fraction (0);
%                                ignored in fixed-count mode
%   minimumTauSeparationNs       separation of the two membrane states
%   minimumTauSeparationFraction relative separation of those states (0.10)
%   modelPrior                   prior for the three models ([1 1 1]/3)
%   maxMembraneStates            highest free membrane-state count that is
%                                evaluated: 2 compares all three models, 1
%                                stops at fixed SLB + one membrane state,
%                                0 leaves only the fixed-SLB model (2)
%   fixedSlbPhotonCount          scalar or [nx ny] outside-calibrated expected
%                                SLB photons; [] leaves amplitude unconstrained
%   fixedSlbPhotonCountStd       scalar or [nx ny] calibration uncertainty
%   slbCountPriorNodes           0 uses a hard count; positive odd values
%                                marginalise over that uncertainty (0)
%   irfShiftBins                 integer shift applied to the IRF (0)
%   convolutionMethod            'auto', 'gui', or 'linear' ('auto')
%
% Principal outputs
%   modelProbability             [nx ny 3], for the models listed above
%   probabilityFixedSlbOnly      one-component probability
%   probabilityBiexponential     fixed SLB + one membrane state
%   probabilityTriexponential    fixed SLB + two membrane states
%   completeExponentialCountMAP  1, 2, or 3 (zero outside validMask)
%   tauMeanArithmetic            model-averaged fluorescence lifetime
%   fixedSlbFraction             model-averaged fluorescence fraction
%   oneMembrane / twoMembrane    model-conditional component maps
%   selected                     component maps selected by model MAP
%   fixedSlbPhotonConstraint      direct-count settings, applied map and clipping
%   fixedSlbCountPrior           backward-compatible alias of that metadata
%
% The helper is intentionally self-contained. If the GUI Convol.m utility
% is on the path, convolutionMethod='auto' uses it; otherwise a causal
% linear convolution is used.

    if nargin < 6 || isempty(opts)
        opts = struct();
    end
    validateattributes(dtNs, {'numeric'}, {'real','finite','scalar','positive'}, ...
        mfilename, 'dtNs');
    validateattributes(tauSlbNs, {'numeric'}, {'real','finite','scalar','positive'}, ...
        mfilename, 'tauSlbNs');

    Ypix = single(tcspcPix);
    if ndims(Ypix) == 4
        Ypix = sum(Ypix, 4);
    end
    if ndims(Ypix) ~= 3
        error('flim_bayes_fixed_slb:TCSPCDimensions', ...
            'tcspcPix must be [nx ny time] or [nx ny time detector].');
    end
    if any(~isfinite(Ypix(:))) || any(Ypix(:) < 0)
        error('flim_bayes_fixed_slb:TCSPCValues', ...
            'tcspcPix must contain finite, nonnegative photon counts.');
    end
    [nx, ny, nt] = size(Ypix);
    modelPeriodNs = double(pulsePeriodNs);
    if ~(isscalar(modelPeriodNs) && isfinite(modelPeriodNs) && ...
            modelPeriodNs >= nt * double(dtNs))
        modelPeriodNs = nt * double(dtNs);
    end

    irf = max(double(irf(:)), 0);
    if numel(irf) ~= nt
        error('flim_bayes_fixed_slb:IRFLength', ...
            'IRF length (%d) must match the TCSPC time dimension (%d).', ...
            numel(irf), nt);
    end
    if sum(irf) <= 0
        irf = zeros(nt, 1);
        irf(1) = 1;
    else
        irf = irf ./ sum(irf);
    end

    opts = fill_options(opts, nx, ny, tauSlbNs, dtNs, modelPeriodNs);
    shiftedIrf = circshift(irf, round(opts.irfShiftBins));
    convolutionMethod = resolve_convolution_method(opts.convolutionMethod);

    tauGrid = resolve_tau_grid(opts, tauSlbNs, dtNs, modelPeriodNs);
    slbPattern = decay_pattern(shiftedIrf, modelPeriodNs, dtNs, ...
        tauSlbNs, convolutionMethod);
    membranePatterns = zeros(nt, numel(tauGrid));
    for tauIndex = 1:numel(tauGrid)
        membranePatterns(:, tauIndex) = decay_pattern(shiftedIrf, ...
            modelPeriodNs, dtNs, tauGrid(tauIndex), convolutionMethod);
    end

    intensity = sum(Ypix, 3);
    analysisMask = logical(opts.analysisMask);
    validMask = analysisMask & double(intensity) >= opts.minPhotons;
    Y = reshape(permute(Ypix, [3 1 2]), nt, []);
    validFlat = validMask(:).';

    fixedCountEnabled = ~isempty(opts.fixedSlbPhotonCount);
    activeModels = false(1, 3);
    activeModels(1:(opts.maxMembraneStates + 1)) = true;
    evaluations = repmat(empty_evaluation(size(Y, 2)), 1, 3);
    modelStateCount = nan(1, 3);
    fractionStepUsed = NaN;
    for modelIndex = 1:3
        if ~activeModels(modelIndex)
            % A skipped model gets zero evidence, so it can never win the
            % comparison and never enters a model average.
            evaluations(modelIndex).logEvidence(:) = -Inf;
            continue;
        end
        membraneCount = modelIndex - 1;
        if fixedCountEnabled
            grid = build_residual_model_grid(membraneCount, ...
                membranePatterns, tauGrid, opts);
            evaluations(modelIndex) = evaluate_fixed_count_grid(Y, ...
                validFlat, grid, slbPattern, tauSlbNs, ...
                opts.fixedSlbPhotonCount, opts.useGPU, opts.batchSize, ...
                opts.slbCountRelTol, opts.fixedSlbPhotonCountStd, ...
                opts.slbCountPriorNodes);
        else
            grid = build_model_grid(membraneCount, slbPattern, ...
                membranePatterns, tauGrid, tauSlbNs, opts);
            evaluations(modelIndex) = evaluate_grid(Y, validFlat, ...
                grid, opts.useGPU, opts.batchSize);
        end
        modelStateCount(modelIndex) = grid.stateCount;
        if modelIndex == 2
            fractionStepUsed = grid.fractionStepUsed;
        end
        % Only one model grid is resident at a time.
        clear grid
    end

    effectivePrior = double(opts.modelPrior(:)).';
    effectivePrior(~activeModels) = 0;
    effectivePrior = effectivePrior ./ sum(effectivePrior);

    [modelProbabilityFlat, modelMapFlat] = combine_models( ...
        evaluations, effectivePrior, validFlat, activeModels);
    conditional = repmat(conditional_maps(evaluations(1), ...
        nx, ny, 0, tauSlbNs), 1, 3);
    for modelIndex = 2:3
        conditional(modelIndex) = conditional_maps(evaluations(modelIndex), ...
            nx, ny, modelIndex - 1, tauSlbNs);
    end

    tauMeanFlat = average_parameter(evaluations, modelProbabilityFlat, 1, ...
        activeModels);
    tauSecondFlat = average_second_moment(evaluations, modelProbabilityFlat, ...
        1, activeModels);
    tauStdFlat = sqrt(max(tauSecondFlat - tauMeanFlat .^ 2, 0));
    slbFractionFlat = average_parameter(evaluations, modelProbabilityFlat, 2, ...
        activeModels);
    backgroundFractionFlat = average_parameter(evaluations, ...
        modelProbabilityFlat, 7, activeModels);
    signalFractionFlat = average_parameter(evaluations, ...
        modelProbabilityFlat, 8, activeModels);
    slbPhotonFractionFlat = average_parameter(evaluations, ...
        modelProbabilityFlat, 9, activeModels);
    membrane1PhotonFractionFlat = average_parameter(evaluations, ...
        modelProbabilityFlat, 10, activeModels);
    membrane2PhotonFractionFlat = average_parameter(evaluations, ...
        modelProbabilityFlat, 11, activeModels);

    selected = selected_maps(evaluations, modelMapFlat, validFlat, ...
        nx, ny, tauSlbNs);
    fixedTauMap = nan(nx, ny, 'single');
    fixedTauMap(validMask) = single(tauSlbNs);

    out = struct();
    out.method = 'fixed-SLB Bayesian TCSPC-shape model comparison';
    out.modelNames = {'fixed SLB only', ...
        'fixed SLB + one membrane lifetime', ...
        'fixed SLB + two membrane lifetimes'};
    out.modelPrior = effectivePrior;
    out.maxMembraneStates = double(opts.maxMembraneStates);
    out.evaluatedModels = activeModels;
    out.modelProbability = reshape(permute(single(modelProbabilityFlat), ...
        [2 3 1]), nx, ny, 3);
    out.probabilityFixedSlbOnly = out.modelProbability(:, :, 1);
    out.probabilityBiexponential = out.modelProbability(:, :, 2);
    out.probabilityTriexponential = out.modelProbability(:, :, 3);
    out.completeExponentialCountMAP = reshape(modelMapFlat, nx, ny);
    out.membraneComponentCountMAP = uint8(max( ...
        double(out.completeExponentialCountMAP) - 1, 0));
    out.intensity = intensity;
    out.analysisMask = analysisMask;
    out.validMask = validMask;
    out.fixedSlbLifetimeNs = double(tauSlbNs);
    out.tauMeanArithmetic = flat_to_map(tauMeanFlat, nx, ny);
    out.tauPosteriorStd = flat_to_map(tauStdFlat, nx, ny);
    out.fixedSlbFraction = flat_to_map(slbFractionFlat, nx, ny);
    out.membraneFraction = 1 - out.fixedSlbFraction;
    out.backgroundFraction = flat_to_map(backgroundFractionFlat, nx, ny);
    out.signalFraction = flat_to_map(signalFractionFlat, nx, ny);
    out.fixedSlbPhotonFraction = flat_to_map(slbPhotonFractionFlat, nx, ny);
    if fixedCountEnabled
        out.fixedSlbExpectedPhotonCount = reshape( ...
            evaluations(1).appliedSlbPhotonCount, nx, ny);
    else
        out.fixedSlbExpectedPhotonCount = single(intensity) .* ...
            out.fixedSlbPhotonFraction;
    end
    out.membrane1PhotonFraction = flat_to_map(membrane1PhotonFractionFlat, nx, ny);
    out.membrane2PhotonFraction = flat_to_map(membrane2PhotonFractionFlat, nx, ny);
    out.oneMembrane = conditional(2);
    out.twoMembrane = conditional(3);
    out.conditionalModel = conditional;
    out.selected = selected;
    out.componentLifetimeNs = cat(3, fixedTauMap, ...
        selected.membraneLifetime1Ns, selected.membraneLifetime2Ns);
    out.componentSpeciesFraction = cat(3, selected.fixedSlbFraction, ...
        selected.membraneFraction1, selected.membraneFraction2);
    out.componentPhotonFraction = cat(3, selected.fixedSlbPhotonFraction, ...
        selected.membrane1PhotonFraction, selected.membrane2PhotonFraction);
    out.logEvidence = cat(3, flat_to_map(evaluations(1).logEvidence, nx, ny), ...
        flat_to_map(evaluations(2).logEvidence, nx, ny), ...
        flat_to_map(evaluations(3).logEvidence, nx, ny));
    out.irf = irf;
    out.irfShifted = shiftedIrf;
    out.irfShiftBins = round(opts.irfShiftBins);
    out.convolutionMethod = convolutionMethod;
    out.pulsePeriodNs = double(pulsePeriodNs);
    out.modelPeriodNs = modelPeriodNs;
    out.dtNs = double(dtNs);
    out.membraneTauGridNs = tauGrid(:).';
    clippedMask = reshape(evaluations(1).fixedSlbCountClipped, nx, ny);
    clippedPixelCount = nnz(clippedMask & validMask);
    clippedPixelFraction = clippedPixelCount / max(nnz(validMask), 1);
    if fixedCountEnabled && opts.slbCountPriorNodes > 0
        constraintMode = 'spatial expected photon-count prior';
        clippingDefinition = ...
            'observed N is below the lower two-sigma calibration support';
    elseif fixedCountEnabled
        constraintMode = 'direct expected photon-count allocation';
        clippingDefinition = ...
            'target count is replaced by observed N when N < target';
    else
        constraintMode = 'unconstrained legacy amplitude';
        clippingDefinition = 'not applicable';
    end
    if fixedCountEnabled
        unsupportedFraction = [0 0 0];
    else
        unsupportedFraction = [NaN NaN NaN];
    end
    out.fixedSlbCountPrior = struct( ...
        'enabled', fixedCountEnabled, ...
        'mode', constraintMode, ...
        'targetPhotonCount', opts.fixedSlbPhotonCount, ...
        'photonCountStd', opts.fixedSlbPhotonCountStd, ...
        'priorNodes', opts.slbCountPriorNodes, ...
        'countMarginalised', fixedCountEnabled && opts.slbCountPriorNodes > 0, ...
        'photonCountStdUsedForInference', ...
            fixedCountEnabled && opts.slbCountPriorNodes > 0, ...
        'appliedExpectedPhotonCount', out.fixedSlbExpectedPhotonCount, ...
        'appliedDetectedPhotonFraction', out.fixedSlbPhotonFraction, ...
        'clippedPixelMask', clippedMask, ...
        'clippedPixelCount', clippedPixelCount, ...
        'clippedPixelFraction', clippedPixelFraction, ...
        'minimumGridDistanceZByModel', [], ...
        'closestGridPhotonCountByModel', [], ...
        'supportThresholdZ', NaN, ...
        'unsupportedPixelFractionByModel', unsupportedFraction, ...
        'statePriorNormalisation', 'uniform over residual states within each model', ...
        'clippingDefinition', clippingDefinition);
    out.fixedSlbPhotonConstraint = out.fixedSlbCountPrior;
    out.gridInfo = struct('modelStateCount', modelStateCount, ...
        'fractionStepRequested', opts.fractionStep, ...
        'fractionStepUsed', fractionStepUsed, ...
        'signalGrid', opts.signalGrid(:).', ...
        'minimumMembraneFraction', opts.minimumMembraneFraction, ...
        'minimumSlbFraction', opts.minimumSlbFraction, ...
        'minimumTauSeparationNs', opts.minimumTauSeparationNs, ...
        'minimumTauSeparationFraction', opts.minimumTauSeparationFraction, ...
        'fixedSlbCountMode', fixedCountEnabled, ...
        'fixedModeSignalGridMeaning', ...
            'membrane share of photons remaining after fixed SLB allocation', ...
        'usedGPU', [evaluations.usedGPU]);
    if fixedCountEnabled && opts.slbCountPriorNodes > 0
        amplitudeAssumption = [ ...
            'The detected SLB count has a spatially calibrated Gaussian prior ' ...
            'that is marginalised during inference.'];
    elseif fixedCountEnabled
        amplitudeAssumption = [ ...
            'The expected detected SLB count is supplied by a scalar or spatial ' ...
            'calibration and applied as a hard constraint.'];
    else
        amplitudeAssumption = 'The SLB amplitude remains pixel-dependent.';
    end
    out.assumptions = { ...
        'tauSlbNs was estimated from an outside-cell SLB reference and is fixed globally.', ...
        amplitudeAssumption, ...
        'Added membrane lifetimes are constrained to be longer than tauSlbNs.', ...
        'Model evidence is conditional on photon total and on the configured discrete priors.', ...
        'The background, when enabled, is uniform across the fitted TCSPC window.', ...
        membrane_state_assumption(opts.maxMembraneStates)};
end

function opts = fill_options(opts, nx, ny, tauSlbNs, dtNs, modelPeriodNs)
    defaults = struct();
    defaults.analysisMask = true(nx, ny);
    defaults.minPhotons = 1;
    defaults.useGPU = [];
    defaults.batchSize = 128;
    defaults.includeBackground = true;
    defaults.signalGrid = linspace(0.25, 1, 4);
    defaults.membraneTauGrid = [];
    tauMinimum = max([1.15 * tauSlbNs, tauSlbNs + 2 * dtNs, 0.05]);
    tauMaximum = min(0.8 * modelPeriodNs, max(5, 4 * tauSlbNs));
    tauMaximum = max(tauMaximum, 1.25 * tauMinimum);
    defaults.membraneTauBoundsNs = [tauMinimum tauMaximum];
    defaults.membraneTauCount = 12;
    defaults.fractionStep = 0.125;
    defaults.minimumMembraneFraction = 0.05;
    defaults.minimumSlbFraction = 0;
    defaults.minimumTauSeparationNs = max(0.10, 2 * dtNs);
    defaults.minimumTauSeparationFraction = 0.10;
    defaults.modelPrior = [1 1 1] / 3;
    defaults.maxMembraneStates = 2;
    defaults.fixedSlbPhotonCount = [];
    defaults.fixedSlbPhotonCountStd = [];
    defaults.slbCountRelTol = 0.02;
    % 0 keeps the historical behaviour: the SLB count is a hard constraint.
    % A positive odd count marginalises over fixedSlbPhotonCountStd instead.
    defaults.slbCountPriorNodes = 0;
    defaults.irfShiftBins = 0;
    defaults.convolutionMethod = 'auto';

    names = fieldnames(defaults);
    for index = 1:numel(names)
        name = names{index};
        if ~isfield(opts, name) || isempty(opts.(name))
            opts.(name) = defaults.(name);
        end
    end
    if ~isequal(size(opts.analysisMask), [nx ny])
        error('flim_bayes_fixed_slb:MaskSize', ...
            'opts.analysisMask must have size [%d %d].', nx, ny);
    end
    opts.analysisMask = logical(opts.analysisMask);
    validateattributes(opts.minPhotons, {'numeric'}, ...
        {'real','finite','scalar','nonnegative'});
    validateattributes(opts.batchSize, {'numeric'}, ...
        {'real','finite','scalar','integer','positive'});
    validateattributes(opts.slbCountRelTol, {'numeric'}, ...
        {'real','finite','scalar','nonnegative','<',1});
    validateattributes(opts.slbCountPriorNodes, {'numeric'}, ...
        {'real','finite','scalar','integer','nonnegative','<=',15});
    if opts.slbCountPriorNodes > 0 && mod(opts.slbCountPriorNodes, 2) == 0
        error('flim_bayes_fixed_slb:SlbPriorNodes', ...
            ['opts.slbCountPriorNodes must be odd so the prior mean is ' ...
             'itself a node (got %d).'], opts.slbCountPriorNodes);
    end
    if opts.slbCountPriorNodes > 0 && isempty(opts.fixedSlbPhotonCountStd)
        error('flim_bayes_fixed_slb:SlbPriorStd', ...
            ['opts.slbCountPriorNodes > 0 requires ' ...
             'opts.fixedSlbPhotonCountStd, which sets the prior width.']);
    end
    validateattributes(opts.includeBackground, {'numeric','logical'}, {'scalar'});
    opts.includeBackground = logical(opts.includeBackground);
    validateattributes(opts.irfShiftBins, {'numeric'}, ...
        {'real','finite','scalar','integer'});
    validateattributes(opts.fractionStep, {'numeric'}, ...
        {'real','finite','scalar','>',0,'<=',0.5});
    validateattributes(opts.minimumMembraneFraction, {'numeric'}, ...
        {'real','finite','scalar','nonnegative','<',0.5});
    validateattributes(opts.minimumSlbFraction, {'numeric'}, ...
        {'real','finite','scalar','nonnegative','<',1});
    validateattributes(opts.minimumTauSeparationNs, {'numeric'}, ...
        {'real','finite','scalar','nonnegative'});
    validateattributes(opts.minimumTauSeparationFraction, {'numeric'}, ...
        {'real','finite','scalar','nonnegative'});
    if opts.includeBackground
        signalGrid = unique(double(opts.signalGrid(:)).', 'sorted');
        if isempty(signalGrid) || any(~isfinite(signalGrid)) || ...
                any(signalGrid < 0 | signalGrid > 1)
            error('flim_bayes_fixed_slb:SignalGrid', ...
                'opts.signalGrid must contain finite values from 0 to 1.');
        end
        opts.signalGrid = signalGrid;
    else
        opts.signalGrid = 1;
    end
    prior = double(opts.modelPrior(:));
    if numel(prior) ~= 3 || any(~isfinite(prior)) || any(prior <= 0)
        error('flim_bayes_fixed_slb:ModelPrior', ...
            'opts.modelPrior must contain three finite positive values.');
    end
    opts.modelPrior = (prior ./ sum(prior)).';
    validateattributes(opts.maxMembraneStates, {'numeric'}, ...
        {'real','finite','scalar','integer','>=',0,'<=',2}, ...
        mfilename, 'opts.maxMembraneStates');
    opts.maxMembraneStates = double(opts.maxMembraneStates);
    fixedCountEnabled = ~isempty(opts.fixedSlbPhotonCount);
    if fixedCountEnabled
        opts.fixedSlbPhotonCount = validate_spatial_calibration( ...
            opts.fixedSlbPhotonCount, nx, ny, ...
            'opts.fixedSlbPhotonCount', false);
        if ~opts.includeBackground
            error('flim_bayes_fixed_slb:FixedSlbCountBackground', ...
                ['opts.fixedSlbPhotonCount requires includeBackground=true ' ...
                 'so model M1 can allocate photons above the fixed SLB count.']);
        end
    end
    if ~isempty(opts.fixedSlbPhotonCountStd)
        opts.fixedSlbPhotonCountStd = validate_spatial_calibration( ...
            opts.fixedSlbPhotonCountStd, nx, ny, ...
            'opts.fixedSlbPhotonCountStd', true);
    end
    if isempty(opts.useGPU)
        opts.useGPU = gpu_available();
    else
        validateattributes(opts.useGPU, {'numeric','logical'}, {'scalar'});
        opts.useGPU = logical(opts.useGPU);
    end
end

function values = validate_spatial_calibration(values, nx, ny, name, allowZero)
    validateattributes(values, {'numeric'}, {'real','finite'});
    if ~(isscalar(values) || numel(values) == nx * ny)
        error('flim_bayes_fixed_slb:SpatialCalibrationSize', ...
            '%s must be scalar or contain exactly %d-by-%d values.', ...
            name, nx, ny);
    end
    values = double(values);
    if allowZero
        invalid = values < 0;
        qualifier = 'nonnegative';
    else
        invalid = values <= 0;
        qualifier = 'positive';
    end
    if any(invalid(:))
        error('flim_bayes_fixed_slb:SpatialCalibrationValue', ...
            '%s must contain only %s values.', name, qualifier);
    end
    if ~isscalar(values)
        values = reshape(values, nx, ny);
    end
end

function tauGrid = resolve_tau_grid(opts, tauSlbNs, dtNs, modelPeriodNs)
    if isempty(opts.membraneTauGrid)
        bounds = double(opts.membraneTauBoundsNs(:)).';
        if numel(bounds) ~= 2 || any(~isfinite(bounds)) || ...
                bounds(1) <= 0 || bounds(2) <= bounds(1)
            error('flim_bayes_fixed_slb:TauBounds', ...
                'opts.membraneTauBoundsNs must be [lower upper] in ns.');
        end
        validateattributes(opts.membraneTauCount, {'numeric'}, ...
            {'real','finite','scalar','integer','>=',3});
        tauGrid = logspace(log10(bounds(1)), log10(bounds(2)), ...
            opts.membraneTauCount);
    else
        tauGrid = double(opts.membraneTauGrid(:)).';
    end
    minimumLongTau = tauSlbNs + max(dtNs, 100 * eps(tauSlbNs));
    tauGrid = unique(tauGrid(isfinite(tauGrid) & tauGrid >= minimumLongTau), ...
        'sorted');
    tauGrid = tauGrid(tauGrid < modelPeriodNs);
    if numel(tauGrid) < 2
        error('flim_bayes_fixed_slb:TauGrid', ...
            ['At least two membrane lifetime grid values must be longer than ' ...
             'tauSlbNs and shorter than the fitted TCSPC window.']);
    end
end

function method = resolve_convolution_method(requested)
    requested = lower(char(requested));
    switch requested
        case 'auto'
            if exist('Convol', 'file') == 2
                method = 'gui';
            else
                method = 'linear';
            end
        case 'gui'
            if exist('Convol', 'file') ~= 2
                error('flim_bayes_fixed_slb:MissingConvol', ...
                    'convolutionMethod=''gui'' requires Convol.m on the MATLAB path.');
            end
            method = 'gui';
        case 'linear'
            method = 'linear';
        otherwise
            error('flim_bayes_fixed_slb:ConvolutionMethod', ...
                'convolutionMethod must be ''auto'', ''gui'', or ''linear''.');
    end
end

function pattern = decay_pattern(irf, modelPeriodNs, dtNs, tauNs, method)
    timeNs = (0:numel(irf)-1).' * dtNs;
    decay = exp(-timeNs ./ tauNs) ./ ...
        max(1 - exp(-modelPeriodNs / tauNs), eps);
    if strcmp(method, 'gui')
        convolved = Convol(irf, decay);
        convolved = double(convolved(:));
        if numel(convolved) < numel(irf)
            error('flim_bayes_fixed_slb:ConvolOutput', ...
                'Convol returned fewer samples than the IRF input.');
        end
        convolved = convolved(1:numel(irf));
    else
        convolved = conv(irf, decay, 'full');
        convolved = convolved(1:numel(irf));
    end
    pattern = max(real(convolved(:)), 0);
    if sum(pattern) <= 0
        error('flim_bayes_fixed_slb:DecayPattern', ...
            'A convolved decay pattern had zero area. Check the IRF.');
    end
    pattern = pattern ./ sum(pattern);
end

function grid = build_model_grid(membraneCount, slbPattern, membranePatterns, ...
        tauGrid, tauSlbNs, opts)
    % Parameter rows: mean tau, SLB species fraction, membrane tau 1/2,
    % membrane species fraction 1/2, background and signal fractions, then
    % SLB and membrane 1/2 fractions of all detected photons.
    parameterCount = 11;
    signalGrid = double(opts.signalGrid(:)).';
    backgroundPattern = ones(size(slbPattern)) / numel(slbPattern);
    latticeCount = max(2, round(1 / opts.fractionStep));
    fractionStepUsed = 1 / latticeCount;

    if membraneCount == 0
        stateCount = numel(signalGrid);
        patterns = zeros(numel(slbPattern), stateCount);
        parameters = zeros(parameterCount, stateCount);
        for signalIndex = 1:numel(signalGrid)
            signal = signalGrid(signalIndex);
            patterns(:, signalIndex) = normalise_pattern( ...
                signal * slbPattern + (1 - signal) * backgroundPattern);
            parameters(:, signalIndex) = [tauSlbNs; 1; 0; 0; 0; 0; ...
                1 - signal; signal; signal; 0; 0];
        end
    elseif membraneCount == 1
        slbFraction = (0:latticeCount) / latticeCount;
        membraneFraction = 1 - slbFraction;
        keep = slbFraction >= opts.minimumSlbFraction & ...
            membraneFraction >= opts.minimumMembraneFraction;
        slbFraction = slbFraction(keep);
        membraneFraction = membraneFraction(keep);
        stateCount = numel(signalGrid) * numel(tauGrid) * numel(slbFraction);
        patterns = zeros(numel(slbPattern), stateCount);
        parameters = zeros(parameterCount, stateCount);
        column = 0;
        for signalIndex = 1:numel(signalGrid)
            signal = signalGrid(signalIndex);
            for tauIndex = 1:numel(tauGrid)
                for fractionIndex = 1:numel(slbFraction)
                    column = column + 1;
                    fSlb = slbFraction(fractionIndex);
                    fMembrane = membraneFraction(fractionIndex);
                    fluorescence = fSlb * slbPattern + ...
                        fMembrane * membranePatterns(:, tauIndex);
                    patterns(:, column) = normalise_pattern(signal * fluorescence + ...
                        (1 - signal) * backgroundPattern);
                    parameters(:, column) = [fSlb * tauSlbNs + ...
                        fMembrane * tauGrid(tauIndex); fSlb; ...
                        tauGrid(tauIndex); 0; fMembrane; 0; ...
                        1 - signal; signal; signal * fSlb; ...
                        signal * fMembrane; 0];
                end
            end
        end
    else
        [tauIndex1, tauIndex2] = lifetime_pairs(tauGrid, ...
            opts.minimumTauSeparationNs, opts.minimumTauSeparationFraction);
        [slbFraction, membraneFraction1, membraneFraction2] = ...
            mixture_lattice(latticeCount, opts.minimumSlbFraction, ...
            opts.minimumMembraneFraction);
        stateCount = numel(signalGrid) * numel(tauIndex1) * numel(slbFraction);
        patterns = zeros(numel(slbPattern), stateCount);
        parameters = zeros(parameterCount, stateCount);
        column = 0;
        for signalIndex = 1:numel(signalGrid)
            signal = signalGrid(signalIndex);
            for pairIndex = 1:numel(tauIndex1)
                t1Index = tauIndex1(pairIndex);
                t2Index = tauIndex2(pairIndex);
                for fractionIndex = 1:numel(slbFraction)
                    column = column + 1;
                    fSlb = slbFraction(fractionIndex);
                    f1 = membraneFraction1(fractionIndex);
                    f2 = membraneFraction2(fractionIndex);
                    fluorescence = fSlb * slbPattern + ...
                        f1 * membranePatterns(:, t1Index) + ...
                        f2 * membranePatterns(:, t2Index);
                    patterns(:, column) = normalise_pattern(signal * fluorescence + ...
                        (1 - signal) * backgroundPattern);
                    parameters(:, column) = [fSlb * tauSlbNs + ...
                        f1 * tauGrid(t1Index) + f2 * tauGrid(t2Index); ...
                        fSlb; tauGrid(t1Index); tauGrid(t2Index); f1; f2; ...
                        1 - signal; signal; signal * fSlb; ...
                        signal * f1; signal * f2];
                end
            end
        end
    end
    if stateCount == 0
        error('flim_bayes_fixed_slb:EmptyGrid', ...
            'Model %d has no states. Relax fraction or lifetime constraints.', ...
            membraneCount + 1);
    end
    grid = struct('patterns', patterns, 'parameters', parameters, ...
        'stateCount', stateCount, 'membraneCount', membraneCount, ...
        'fractionStepUsed', fractionStepUsed);
end

function grid = build_residual_model_grid(membraneCount, membranePatterns, ...
        tauGrid, opts)
    % In fixed-count mode the SLB contribution is injected per photon-total
    % group. These grids describe only the remaining, normalised photons.
    binCount = size(membranePatterns, 1);
    backgroundPattern = ones(binCount, 1) / binCount;
    signalGrid = double(opts.signalGrid(:)).';
    latticeCount = max(2, round(1 / opts.fractionStep));
    fractionStepUsed = 1 / latticeCount;

    if membraneCount == 0
        stateCount = 1;
        residualPatterns = backgroundPattern;
        tau1 = 0;
        tau2 = 0;
        membraneFraction1 = 0;
        membraneFraction2 = 0;
        backgroundFraction = 1;
    elseif membraneCount == 1
        stateCount = numel(signalGrid) * numel(tauGrid);
        residualPatterns = zeros(binCount, stateCount);
        tau1 = zeros(1, stateCount);
        tau2 = zeros(1, stateCount);
        membraneFraction1 = zeros(1, stateCount);
        membraneFraction2 = zeros(1, stateCount);
        backgroundFraction = zeros(1, stateCount);
        column = 0;
        for signalIndex = 1:numel(signalGrid)
            residualSignal = signalGrid(signalIndex);
            for tauIndex = 1:numel(tauGrid)
                column = column + 1;
                residualPatterns(:, column) = normalise_pattern( ...
                    residualSignal * membranePatterns(:, tauIndex) + ...
                    (1 - residualSignal) * backgroundPattern);
                tau1(column) = tauGrid(tauIndex);
                membraneFraction1(column) = residualSignal;
                backgroundFraction(column) = 1 - residualSignal;
            end
        end
    else
        [tauIndex1, tauIndex2] = lifetime_pairs(tauGrid, ...
            opts.minimumTauSeparationNs, opts.minimumTauSeparationFraction);
        relativeFirst = (0:latticeCount) / latticeCount;
        relativeSecond = 1 - relativeFirst;
        keep = relativeFirst >= opts.minimumMembraneFraction & ...
            relativeSecond >= opts.minimumMembraneFraction;
        relativeFirst = relativeFirst(keep);
        relativeSecond = relativeSecond(keep);
        if isempty(relativeFirst)
            error('flim_bayes_fixed_slb:NoResidualMixtures', ...
                ['No two-membrane residual mixtures satisfy ' ...
                 'opts.minimumMembraneFraction.']);
        end
        stateCount = numel(signalGrid) * numel(tauIndex1) * ...
            numel(relativeFirst);
        residualPatterns = zeros(binCount, stateCount);
        tau1 = zeros(1, stateCount);
        tau2 = zeros(1, stateCount);
        membraneFraction1 = zeros(1, stateCount);
        membraneFraction2 = zeros(1, stateCount);
        backgroundFraction = zeros(1, stateCount);
        column = 0;
        for signalIndex = 1:numel(signalGrid)
            residualSignal = signalGrid(signalIndex);
            for pairIndex = 1:numel(tauIndex1)
                firstTauIndex = tauIndex1(pairIndex);
                secondTauIndex = tauIndex2(pairIndex);
                for fractionIndex = 1:numel(relativeFirst)
                    column = column + 1;
                    firstFraction = residualSignal * ...
                        relativeFirst(fractionIndex);
                    secondFraction = residualSignal * ...
                        relativeSecond(fractionIndex);
                    residualPatterns(:, column) = normalise_pattern( ...
                        firstFraction * membranePatterns(:, firstTauIndex) + ...
                        secondFraction * membranePatterns(:, secondTauIndex) + ...
                        (1 - residualSignal) * backgroundPattern);
                    tau1(column) = tauGrid(firstTauIndex);
                    tau2(column) = tauGrid(secondTauIndex);
                    membraneFraction1(column) = firstFraction;
                    membraneFraction2(column) = secondFraction;
                    backgroundFraction(column) = 1 - residualSignal;
                end
            end
        end
    end
    grid = struct('residualPatterns', residualPatterns, ...
        'membraneTau1Ns', tau1, 'membraneTau2Ns', tau2, ...
        'residualMembraneFraction1', membraneFraction1, ...
        'residualMembraneFraction2', membraneFraction2, ...
        'residualBackgroundFraction', backgroundFraction, ...
        'stateCount', stateCount, 'membraneCount', membraneCount, ...
        'fractionStepUsed', fractionStepUsed);
end

function grid = instantiate_fixed_count_grid(residual, slbPattern, ...
        tauSlbNs, fixedSlbPhotonFraction)
    q = min(max(double(fixedSlbPhotonFraction), 0), 1);
    residualFraction = 1 - q;
    patterns = q * slbPattern + residualFraction * residual.residualPatterns;
    stateCount = residual.stateCount;

    fixedSlbPhotonFraction = repmat(q, 1, stateCount);
    membranePhotonFraction1 = residualFraction * ...
        residual.residualMembraneFraction1;
    membranePhotonFraction2 = residualFraction * ...
        residual.residualMembraneFraction2;
    backgroundPhotonFraction = residualFraction * ...
        residual.residualBackgroundFraction;
    signalPhotonFraction = fixedSlbPhotonFraction + ...
        membranePhotonFraction1 + membranePhotonFraction2;
    fixedSlbSpeciesFraction = fixedSlbPhotonFraction ./ signalPhotonFraction;
    membraneSpeciesFraction1 = membranePhotonFraction1 ./ signalPhotonFraction;
    membraneSpeciesFraction2 = membranePhotonFraction2 ./ signalPhotonFraction;
    meanLifetime = (fixedSlbPhotonFraction * tauSlbNs + ...
        membranePhotonFraction1 .* residual.membraneTau1Ns + ...
        membranePhotonFraction2 .* residual.membraneTau2Ns) ./ ...
        signalPhotonFraction;

    parameters = [meanLifetime; fixedSlbSpeciesFraction; ...
        residual.membraneTau1Ns; residual.membraneTau2Ns; ...
        membraneSpeciesFraction1; membraneSpeciesFraction2; ...
        backgroundPhotonFraction; signalPhotonFraction; ...
        fixedSlbPhotonFraction; membranePhotonFraction1; ...
        membranePhotonFraction2];
    grid = struct('patterns', patterns, 'parameters', parameters, ...
        'stateCount', stateCount, 'membraneCount', residual.membraneCount, ...
        'fractionStepUsed', residual.fractionStepUsed);
end

function pattern = normalise_pattern(pattern)
    pattern = max(double(pattern(:)), realmin('single'));
    pattern = pattern ./ sum(pattern);
end

function [firstIndex, secondIndex] = lifetime_pairs(tauGrid, minNs, minFraction)
    firstIndex = [];
    secondIndex = [];
    for first = 1:numel(tauGrid)-1
        for second = first+1:numel(tauGrid)
            difference = tauGrid(second) - tauGrid(first);
            relativeDifference = difference / max(mean(tauGrid([first second])), eps);
            if difference >= minNs && relativeDifference >= minFraction
                firstIndex(end + 1) = first; %#ok<AGROW>
                secondIndex(end + 1) = second; %#ok<AGROW>
            end
        end
    end
    if isempty(firstIndex)
        error('flim_bayes_fixed_slb:NoLifetimePairs', ...
            'No membrane lifetime pairs satisfy the separation constraints.');
    end
end

function [fSlb, f1, f2] = mixture_lattice(latticeCount, minSlb, minMembrane)
    fSlb = [];
    f1 = [];
    f2 = [];
    for slbIndex = 0:latticeCount
        for firstIndex = 0:(latticeCount - slbIndex)
            secondIndex = latticeCount - slbIndex - firstIndex;
            fractions = [slbIndex firstIndex secondIndex] / latticeCount;
            if fractions(1) >= minSlb && ...
                    all(fractions(2:3) >= minMembrane)
                fSlb(end + 1) = fractions(1); %#ok<AGROW>
                f1(end + 1) = fractions(2); %#ok<AGROW>
                f2(end + 1) = fractions(3); %#ok<AGROW>
            end
        end
    end
    if isempty(fSlb)
        error('flim_bayes_fixed_slb:NoMixtures', ...
            'No three-component mixtures satisfy the fraction constraints.');
    end
end

function evaluation = empty_evaluation(pixelCount)
    evaluation = struct('logEvidence', nan(1, pixelCount), ...
        'mean', nan(11, pixelCount, 'single'), ...
        'std', nan(11, pixelCount, 'single'), ...
        'map', nan(11, pixelCount, 'single'), ...
        'appliedSlbPhotonCount', nan(1, pixelCount), ...
        'fixedSlbCountClipped', false(1, pixelCount), ...
        'stateCount', 0, 'usedGPU', false);
end

function evaluation = evaluate_grid(Y, validFlat, grid, useGPU, batchSize)
    pixelCount = size(Y, 2);
    evaluation = empty_evaluation(pixelCount);
    evaluation.stateCount = grid.stateCount;
    logPattern = log(max(double(grid.patterns), realmin('double')));
    parameters = double(grid.parameters);
    parameterSquares = parameters .^ 2;

    gpuUsed = false;
    if useGPU
        try
            logPatternGpu = gpuArray(single(logPattern));
            gpuUsed = true;
        catch
            gpuUsed = false;
        end
    end
    validIndex = find(validFlat);
    for first = 1:batchSize:numel(validIndex)
        index = validIndex(first:min(first + batchSize - 1, numel(validIndex)));
        counts = Y(:, index);
        if gpuUsed
            likelihood = gather(logPatternGpu.' * gpuArray(single(counts)));
            likelihood = double(likelihood);
        else
            likelihood = logPattern.' * double(counts);
        end
        maximum = max(likelihood, [], 1);
        weight = exp(likelihood - maximum);
        weightSum = sum(weight, 1);
        weightSum(weightSum <= 0) = 1;
        weight = weight ./ weightSum;
        meanValue = parameters * weight;
        secondMoment = parameterSquares * weight;
        evaluation.mean(:, index) = single(meanValue);
        evaluation.std(:, index) = single(sqrt(max(secondMoment - meanValue .^ 2, 0)));
        [~, mapIndex] = max(likelihood, [], 1);
        evaluation.map(:, index) = single(parameters(:, mapIndex));
        evaluation.logEvidence(index) = maximum + log(weightSum) - ...
            log(grid.stateCount);
    end
    evaluation.usedGPU = gpuUsed;
end

function evaluation = evaluate_fixed_count_grid(Y, validFlat, residualGrid, ...
        slbPattern, tauSlbNs, fixedSlbPhotonCount, useGPU, batchSize, ...
        slbCountRelTol, slbCountStd, slbCountPriorNodes)
    if nargin < 9 || isempty(slbCountRelTol)
        slbCountRelTol = 0;
    end
    if nargin < 10
        slbCountStd = [];
    end
    if nargin < 11 || isempty(slbCountPriorNodes)
        slbCountPriorNodes = 0;
    end
    pixelCount = size(Y, 2);
    evaluation = empty_evaluation(pixelCount);
    evaluation.stateCount = residualGrid.stateCount;
    validIndex = find(validFlat);
    if isempty(validIndex)
        return;
    end

    photonTotal = sum(double(Y(:, validIndex)), 1);
    targetCountAll = spatial_calibration_values(fixedSlbPhotonCount, ...
        validIndex, pixelCount);
    if isempty(slbCountStd)
        targetStdAll = zeros(size(targetCountAll));
    else
        targetStdAll = spatial_calibration_values(slbCountStd, ...
            validIndex, pixelCount);
    end
    centreFractionAll = min(max(targetCountAll ./ max(photonTotal, eps), 0), 1);
    stdFractionAll = max(targetStdAll ./ max(photonTotal, eps), 0);

    % A spatial prior can give every pixel a different target and uncertainty.
    % Share likelihood grids only when both fractional prior parameters are
    % close. This preserves the calibration map while avoiding one expensive
    % grid construction and upload per pixel.
    groupParameters = [centreFractionAll(:) stdFractionAll(:)];
    if slbCountPriorNodes <= 0
        groupParameters(:, 2) = 0;
    end
    if slbCountRelTol > 0
        denominator = log1p(slbCountRelTol);
        groupKey = floor(log(max(groupParameters, eps)) / denominator);
        groupKey(groupParameters == 0) = -realmax('double');
        [~, ~, totalGroup] = unique(groupKey, 'rows', 'sorted');
    else
        [~, ~, totalGroup] = unique(groupParameters, 'rows', 'sorted');
    end
    totalGroup = reshape(totalGroup, 1, []);
    groupCount = max(totalGroup);

    gpuUsed = false;
    for groupIndex = 1:groupCount
        member = totalGroup == groupIndex;
        groupPixel = validIndex(member);
        % The SLB photon count is a calibration input, not a measurement, so
        % treating it as exact forces every error in it into the free
        % lifetimes or into a spurious extra component. With
        % slbCountPriorNodes > 0 the count is marginalised over a Gaussian
        % prior of width slbCountStd instead, which lets a mis-specified
        % reference be absorbed in the SLB amplitude where it belongs.
        centreFraction = median(centreFractionAll(member));
        stdFraction = median(stdFractionAll(member));
        [candidateFractions, candidateWeights] = slbFractionCandidates( ...
            centreFraction, stdFraction, slbCountPriorNodes);
        candidateCount = numel(candidateFractions);
        pixelCountInGroup = numel(groupPixel);
        logEvidenceAll = nan(candidateCount, pixelCountInGroup);
        meanAll = nan(11, pixelCountInGroup, candidateCount);
        stdByCandidate = nan(11, pixelCountInGroup, candidateCount);
        mapAll = nan(11, pixelCountInGroup, candidateCount);
        for candidateIndex = 1:candidateCount
            groupGrid = instantiate_fixed_count_grid(residualGrid, ...
                slbPattern, tauSlbNs, candidateFractions(candidateIndex));
            [logEvidenceOne, meanOne, stdOne, mapOne, groupUsedGPU] = ...
                evaluate_grid_indices(Y, groupPixel, groupGrid, useGPU, ...
                batchSize);
            logEvidenceAll(candidateIndex, :) = logEvidenceOne;
            meanAll(:, :, candidateIndex) = meanOne;
            stdByCandidate(:, :, candidateIndex) = stdOne;
            mapAll(:, :, candidateIndex) = mapOne;
            gpuUsed = gpuUsed || groupUsedGPU;
        end

        if candidateCount == 1
            logEvidence = logEvidenceAll(1, :);
            meanValue = meanAll(:, :, 1);
            stdValue = stdByCandidate(:, :, 1);
            mapValue = mapAll(:, :, 1);
        else
            [logEvidence, meanValue, stdValue, mapValue] = ...
                mixOverSlbCount(logEvidenceAll, meanAll, stdByCandidate, mapAll, ...
                candidateWeights);
        end
        evaluation.logEvidence(groupPixel) = logEvidence;
        evaluation.mean(:, groupPixel) = meanValue;
        evaluation.std(:, groupPixel) = stdValue;
        evaluation.map(:, groupPixel) = mapValue;
        % Per-pixel bookkeeping stays exact even when the grid is shared.
        memberTotal = photonTotal(member);
        evaluation.appliedSlbPhotonCount(groupPixel) = ...
            double(meanValue(9, :)) .* memberTotal;
        if slbCountPriorNodes > 0
            lowerSupport = max(targetCountAll(member) - ...
                2 * targetStdAll(member), 0);
        else
            lowerSupport = targetCountAll(member);
        end
        evaluation.fixedSlbCountClipped(groupPixel) = ...
            memberTotal < lowerSupport;
        gpuUsed = gpuUsed || groupUsedGPU;
    end
    evaluation.usedGPU = gpuUsed;
end

function values = spatial_calibration_values(calibration, index, pixelCount)
    if isscalar(calibration)
        values = repmat(double(calibration), 1, numel(index));
    else
        if numel(calibration) ~= pixelCount
            error('flim_bayes_fixed_slb:SpatialCalibrationInternalSize', ...
                'Spatial calibration does not match the flattened image size.');
        end
        flattened = double(calibration(:)).';
        values = flattened(index);
    end
end

function [fractions, weights] = slbFractionCandidates(centre, fractionStd, nodes)
%SLBFRACTIONCANDIDATES Gaussian prior nodes for the SLB photon fraction.
% nodes == 0 reproduces the hard constraint exactly: one node, unit weight.
    if nodes <= 0 || isempty(fractionStd) || ...
            ~isfinite(fractionStd) || fractionStd <= 0
        fractions = min(max(centre, 0), 1);
        weights = 1;
        return;
    end
    half = (nodes - 1) / 2;
    offsets = (-half:half);
    % Nodes span +/-2 sigma, which covers the prior without wasting grid
    % instantiations far out in the tails.
    offsets = 2 * offsets / max(half, 1);
    fractions = centre + offsets * fractionStd;
    fractions = min(max(fractions, 0), 1);
    weights = exp(-0.5 * offsets .^ 2);
    weights = weights / sum(weights);
end

function [logEvidence, meanValue, stdValue, mapValue] = ...
        mixOverSlbCount(logEvidenceAll, meanAll, stdAll, mapAll, priorWeights)
%MIXOVERSLBCOUNT Prior-weighted mixture over the candidate SLB counts.
% The marginal likelihood is the prior-weighted average of the per-candidate
% evidences; the reported moments are those of the resulting mixture, so the
% posterior width now includes uncertainty in the SLB calibration itself.
    pixelCount = size(logEvidenceAll, 2);
    candidateCount = size(logEvidenceAll, 1);
    logEvidence = nan(1, pixelCount);
    meanValue = nan(11, pixelCount, 'single');
    stdValue = nan(11, pixelCount, 'single');
    mapValue = nan(11, pixelCount, 'single');

    usable = all(isfinite(logEvidenceAll), 1);
    if ~any(usable)
        return;
    end

    logPrior = reshape(log(priorWeights(:)), candidateCount, 1);
    shifted = logEvidenceAll(:, usable) + logPrior;
    maximum = max(shifted, [], 1);
    weight = exp(shifted - maximum);
    weightSum = sum(weight, 1);
    logEvidence(usable) = maximum + log(weightSum);
    weight = weight ./ weightSum;

    firstMoment = zeros(11, nnz(usable));
    secondMoment = zeros(11, nnz(usable));
    for candidateIndex = 1:candidateCount
        candidateWeight = weight(candidateIndex, :);
        candidateMean = double(meanAll(:, usable, candidateIndex));
        candidateStd = double(stdAll(:, usable, candidateIndex));
        firstMoment = firstMoment + candidateMean .* candidateWeight;
        secondMoment = secondMoment + ...
            (candidateMean .^ 2 + candidateStd .^ 2) .* candidateWeight;
    end
    meanValue(:, usable) = single(firstMoment);
    stdValue(:, usable) = single(sqrt(max(secondMoment - ...
        firstMoment .^ 2, 0)));

    % The MAP is a single parameter vector, so take it from the candidate
    % carrying the most posterior weight rather than averaging.
    [~, bestCandidate] = max(weight, [], 1);
    usableIndex = find(usable);
    for candidateIndex = 1:candidateCount
        selected = usableIndex(bestCandidate == candidateIndex);
        if isempty(selected)
            continue;
        end
        mapValue(:, selected) = mapAll(:, selected, candidateIndex);
    end
end

function [logEvidence, meanOutput, stdOutput, mapOutput, gpuUsed] = ...
        evaluate_grid_indices(Y, pixelIndex, grid, useGPU, batchSize)
    outputCount = numel(pixelIndex);
    logEvidence = nan(1, outputCount);
    meanOutput = nan(11, outputCount, 'single');
    stdOutput = nan(11, outputCount, 'single');
    mapOutput = nan(11, outputCount, 'single');
    logStateCount = log(grid.stateCount);
    logPattern = [];
    parameters = [];
    parameterSquares = [];

    gpuUsed = false;
    if useGPU
        try
            % Take the log on the device. The log of the
            % [bins x stateCount] pattern matrix is evaluated once per
            % group and is the single largest CPU cost in this loop.
            logPatternGpu = log(max(gpuArray(single(grid.patterns)), ...
                realmin('single')));
            parametersGpu = gpuArray(single(grid.parameters));
            parameterSquaresGpu = parametersGpu .^ 2;
            gpuUsed = true;
        catch
            gpuUsed = false;
        end
    end
    if ~gpuUsed
        logPattern = log(max(double(grid.patterns), realmin('double')));
        parameters = double(grid.parameters);
        parameterSquares = parameters .^ 2;
    end

    for first = 1:batchSize:outputCount
        localIndex = first:min(first + batchSize - 1, outputCount);
        counts = Y(:, pixelIndex(localIndex));
        if gpuUsed
            % Reduce on the device and move only the per-pixel summaries
            % back, instead of the full [stateCount x batch] likelihood.
            likelihood = logPatternGpu.' * gpuArray(single(counts));
            maximum = max(likelihood, [], 1);
            weight = exp(likelihood - maximum);
            weightSum = sum(weight, 1);
            weightSum(weightSum <= 0) = 1;
            weight = weight ./ weightSum;
            meanValue = parametersGpu * weight;
            secondMoment = parameterSquaresGpu * weight;
            [~, mapIndex] = max(likelihood, [], 1);
            mapValue = parametersGpu(:, mapIndex);
            meanValue = double(gather(meanValue));
            secondMoment = double(gather(secondMoment));
            mapValue = gather(mapValue);
            maximum = double(gather(maximum));
            weightSum = double(gather(weightSum));
        else
            likelihood = logPattern.' * double(counts);
            maximum = max(likelihood, [], 1);
            weight = exp(likelihood - maximum);
            weightSum = sum(weight, 1);
            weightSum(weightSum <= 0) = 1;
            weight = weight ./ weightSum;
            meanValue = parameters * weight;
            secondMoment = parameterSquares * weight;
            [~, mapIndex] = max(likelihood, [], 1);
            mapValue = parameters(:, mapIndex);
        end
        meanOutput(:, localIndex) = single(meanValue);
        stdOutput(:, localIndex) = single(sqrt(max( ...
            secondMoment - meanValue .^ 2, 0)));
        mapOutput(:, localIndex) = single(mapValue);
        logEvidence(localIndex) = maximum + log(weightSum) - logStateCount;
    end
end

function [probability, modelMap] = combine_models(evaluations, prior, ...
        validFlat, activeModels)
    pixelCount = numel(validFlat);
    probability = nan(3, pixelCount);
    modelMap = zeros(1, pixelCount, 'uint8');
    validIndex = find(validFlat);
    if isempty(validIndex)
        return;
    end
    logPosterior = -inf(3, numel(validIndex));
    for modelIndex = find(activeModels & prior(:).' > 0)
        logPosterior(modelIndex, :) = ...
            evaluations(modelIndex).logEvidence(validIndex) + log(prior(modelIndex));
    end
    maximum = max(logPosterior, [], 1);
    posterior = exp(logPosterior - maximum);
    posterior = posterior ./ sum(posterior, 1);
    probability(:, validIndex) = posterior;
    [~, selected] = max(posterior, [], 1);
    modelMap(validIndex) = uint8(selected);
end

function value = average_parameter(evaluations, probability, row, activeModels)
    value = zeros(1, size(probability, 2));
    for modelIndex = find(activeModels)
        value = value + probability(modelIndex, :) .* ...
            double(evaluations(modelIndex).mean(row, :));
    end
    invalid = any(~isfinite(probability), 1);
    value(invalid) = NaN;
end

function value = average_second_moment(evaluations, probability, row, ...
        activeModels)
    value = zeros(1, size(probability, 2));
    for modelIndex = find(activeModels)
        meanValue = double(evaluations(modelIndex).mean(row, :));
        stdValue = double(evaluations(modelIndex).std(row, :));
        value = value + probability(modelIndex, :) .* ...
            (meanValue .^ 2 + stdValue .^ 2);
    end
    invalid = any(~isfinite(probability), 1);
    value(invalid) = NaN;
end

function maps = conditional_maps(evaluation, nx, ny, membraneCount, tauSlbNs)
    maps = struct();
    maps.membraneComponentCount = membraneCount;
    maps.fixedSlbLifetimeNs = tauSlbNs;
    maps.tauMeanArithmetic = flat_to_map(evaluation.mean(1, :), nx, ny);
    maps.tauMeanStd = flat_to_map(evaluation.std(1, :), nx, ny);
    maps.fixedSlbFraction = flat_to_map(evaluation.mean(2, :), nx, ny);
    maps.membraneLifetime1Ns = flat_to_map(evaluation.mean(3, :), nx, ny);
    maps.membraneLifetime1StdNs = flat_to_map(evaluation.std(3, :), nx, ny);
    maps.membraneLifetime2Ns = flat_to_map(evaluation.mean(4, :), nx, ny);
    maps.membraneLifetime2StdNs = flat_to_map(evaluation.std(4, :), nx, ny);
    maps.membraneFraction1 = flat_to_map(evaluation.mean(5, :), nx, ny);
    maps.membraneFraction2 = flat_to_map(evaluation.mean(6, :), nx, ny);
    maps.backgroundFraction = flat_to_map(evaluation.mean(7, :), nx, ny);
    maps.signalFraction = flat_to_map(evaluation.mean(8, :), nx, ny);
    maps.fixedSlbPhotonFraction = flat_to_map(evaluation.mean(9, :), nx, ny);
    maps.membrane1PhotonFraction = flat_to_map(evaluation.mean(10, :), nx, ny);
    maps.membrane2PhotonFraction = flat_to_map(evaluation.mean(11, :), nx, ny);
    maps.tauMeanMAP = flat_to_map(evaluation.map(1, :), nx, ny);
    maps.membraneLifetime1MAPNs = flat_to_map(evaluation.map(3, :), nx, ny);
    maps.membraneLifetime2MAPNs = flat_to_map(evaluation.map(4, :), nx, ny);
    if membraneCount < 1
        maps.membraneLifetime1Ns(:) = NaN;
        maps.membraneLifetime1StdNs(:) = NaN;
        maps.membraneLifetime1MAPNs(:) = NaN;
    end
    if membraneCount < 2
        maps.membraneLifetime2Ns(:) = NaN;
        maps.membraneLifetime2StdNs(:) = NaN;
        maps.membraneLifetime2MAPNs(:) = NaN;
    end
end

function selected = selected_maps(evaluations, modelMap, validFlat, nx, ny, tauSlbNs)
    selectedFlat = nan(11, numel(validFlat), 'single');
    for modelIndex = 1:3
        use = validFlat & (modelMap == modelIndex);
        selectedFlat(:, use) = evaluations(modelIndex).mean(:, use);
    end
    selected = struct();
    selected.tauMeanArithmetic = flat_to_map(selectedFlat(1, :), nx, ny);
    selected.fixedSlbFraction = flat_to_map(selectedFlat(2, :), nx, ny);
    selected.membraneLifetime1Ns = flat_to_map(selectedFlat(3, :), nx, ny);
    selected.membraneLifetime2Ns = flat_to_map(selectedFlat(4, :), nx, ny);
    selected.membraneFraction1 = flat_to_map(selectedFlat(5, :), nx, ny);
    selected.membraneFraction2 = flat_to_map(selectedFlat(6, :), nx, ny);
    selected.backgroundFraction = flat_to_map(selectedFlat(7, :), nx, ny);
    selected.signalFraction = flat_to_map(selectedFlat(8, :), nx, ny);
    selected.fixedSlbPhotonFraction = flat_to_map(selectedFlat(9, :), nx, ny);
    selected.membrane1PhotonFraction = flat_to_map(selectedFlat(10, :), nx, ny);
    selected.membrane2PhotonFraction = flat_to_map(selectedFlat(11, :), nx, ny);
    noFirstMembrane = reshape(validFlat & (modelMap < 2), nx, ny);
    noSecondMembrane = reshape(validFlat & (modelMap < 3), nx, ny);
    selected.membraneLifetime1Ns(noFirstMembrane) = NaN;
    selected.membraneLifetime2Ns(noSecondMembrane) = NaN;
    selected.membraneFraction1(noFirstMembrane) = 0;
    selected.membraneFraction2(noSecondMembrane) = 0;
    selected.membrane1PhotonFraction(noFirstMembrane) = 0;
    selected.membrane2PhotonFraction(noSecondMembrane) = 0;
    selected.fixedSlbLifetimeNs = tauSlbNs;
end

function map = flat_to_map(values, nx, ny)
    map = reshape(single(values), nx, ny);
end

function text = membrane_state_assumption(maxMembraneStates)
    if maxMembraneStates >= 2
        text = ['All three models were evaluated; no membrane-state ' ...
            'ceiling was applied.'];
    else
        text = sprintf(['The comparison was truncated at %d free membrane ' ...
            'state(s); richer models carry zero posterior probability by ' ...
            'construction, not by evidence.'], maxMembraneStates);
    end
end

function available = gpu_available()
    try
        available = gpuDeviceCount > 0;
    catch
        available = false;
    end
end
