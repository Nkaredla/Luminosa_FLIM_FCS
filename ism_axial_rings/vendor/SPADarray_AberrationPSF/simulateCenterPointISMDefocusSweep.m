function results = simulateCenterPointISMDefocusSweep(varargin)
%SIMULATECENTERPOINTISMDEFOCUSSWEEP Simulate and fit center-point ISM data.
%
%   results = simulateCenterPointISMDefocusSweep()
%
%   Simulates detector micro-images at the bead centre for two planes:
%   focus plus a defocused plane. The defocus distance is swept from
%   100 nm to 1 um in 100 nm steps by default. For each distance, noisy
%   23 x 2 detector data are generated and passed to
%   estimateCenterPointISMWavefront to retrieve the aberration coefficients.
%
%   The output ranks the defocus distances by recovery error and by local
%   center-point information metrics.
%
%   Examples:
%       results = simulateCenterPointISMDefocusSweep();
%
%       results = simulateCenterPointISMDefocusSweep( ...
%           'defocusUm', 0.1:0.1:1.0, ...
%           'truthSet', 'singleModes', ...
%           'nNoiseReplicates', 3, ...
%           'photonsPerPlane', 5e4);
%
%       results = simulateCenterPointISMDefocusSweep( ...
%           'truthCoeffs', struct('coma_x', 0.04, 'spherical', 0.03), ...
%           'defocusUm', [0.3 0.6]);

    opts = parseOptions(varargin{:});
    addRequiredPaths();

    if ~isempty(opts.randomSeed)
        rng(opts.randomSeed);
    end

    baseSim = configureBaseSim(opts);
    cases = buildTruthCases(opts, baseSim);

    trials = repmat(emptyTrial(), 0, 1);
    trialIndex = 0;
    totalTrials = numel(opts.defocusUm) * numel(cases) * opts.nNoiseReplicates;

    for iz = 1:numel(opts.defocusUm)
        dz = opts.defocusUm(iz);
        planeZ = [0 dz];
        sim = prepareSimForPlaneZ(baseSim, planeZ, opts);

        for ic = 1:numel(cases)
            truth = cases(ic);
            [cleanCenter, cleanNorm] = simulateCenterMicroimage(sim, truth.coeffs, planeZ, opts);

            for ir = 1:opts.nNoiseReplicates
                trialIndex = trialIndex + 1;
                observedCounts = addNoiseToCenterData(cleanNorm, opts);

                trial = emptyTrial();
                trial.defocusUm = dz;
                trial.defocusNm = dz * 1000;
                trial.caseName = truth.name;
                trial.replicate = ir;
                trial.planeZ = planeZ;
                trial.trueCoeffVector = truth.coeffVector;
                trial.cleanCenter = cleanCenter;
                trial.cleanCenterNormalized = cleanNorm;
                trial.observedCounts = observedCounts;

                try
                    fitArgs = retrievalArgs(opts, sim, planeZ);
                    fitResult = estimateCenterPointISMWavefront(observedCounts, [], fitArgs{:});
                    trial.fitResult = fitResult;
                    trial.estimatedCoeffVector = fitResult.fit.estCoeffVector(:).';
                    trial.modelCenter = fitResult.fit.modelN;
                    trial.residualNorm = fitResult.fit.residualNorm;
                    trial.rank = fitResult.sufficiency.rank;
                    trial.nParameters = fitResult.sufficiency.nParameters;
                    trial.isFullRank = fitResult.sufficiency.isFullRank;
                    trial.conditionNumber = fitResult.sufficiency.conditionNumber;
                    trial.minSingularValue = fitResult.sufficiency.minSingularValue;
                    trial.errorMessage = '';
                catch err
                    if ~opts.continueOnError
                        rethrow(err);
                    end
                    warning('simulateCenterPointISMDefocusSweep:FitFailed', ...
                        'Fit failed for dz=%.3f um, case %s, replicate %d: %s', ...
                        dz, truth.name, ir, err.message);
                    trial.fitResult = [];
                    trial.estimatedCoeffVector = nan(size(truth.coeffVector));
                    trial.modelCenter = nan(size(cleanNorm));
                    trial.residualNorm = NaN;
                    trial.rank = NaN;
                    trial.nParameters = numel(opts.fitModes) + 2*double(opts.fitXY) + double(opts.fitZ);
                    trial.isFullRank = false;
                    trial.conditionNumber = NaN;
                    trial.minSingularValue = NaN;
                    trial.errorMessage = err.message;
                end

                trial = scoreTrial(trial, baseSim, opts.fitModes);
                trials(end+1, 1) = trial; %#ok<AGROW>

                if opts.verbose
                    fprintf('[simulateCenterPointISMDefocusSweep] %d/%d dz=%.3f um case=%s rep=%d rmse=%.3g waves\n', ...
                        trialIndex, totalTrials, dz, truth.name, ir, trial.coeffRmseWaves);
                end
            end
        end
    end

    trialTable = makeTrialTable(trials);
    summaryTable = makeSummaryTable(trials, opts.defocusUm);
    [bestByError, bestByInformation] = selectBestDefocus(summaryTable);

    results = struct();
    results.options = opts;
    results.sim = baseSim;
    results.truthCases = cases;
    results.trials = trials;
    results.trialTable = trialTable;
    results.summaryTable = summaryTable;
    results.bestByError = bestByError;
    results.bestByInformation = bestByInformation;
    results.outputDir = resolveOutputDir(opts);

    if opts.writeOutputs
        writeSweepOutputs(results);
    end

    if opts.verbose
        printSweepSummary(results);
    end
end

function opts = parseOptions(varargin)
    p = inputParser;
    p.FunctionName = 'simulateCenterPointISMDefocusSweep';

    addParameter(p, 'defocusUm', 0.1:0.1:1.0);
    addParameter(p, 'truthSet', 'mixed');
    addParameter(p, 'truthCoeffs', []);
    addParameter(p, 'singleModeAmplitudeWaves', 0.04);
    addParameter(p, 'mixedCoeffScale', 1.0);
    addParameter(p, 'fitModes', {'defocus','astig_x','astig_y','coma_x','coma_y','spherical'});
    addParameter(p, 'fitXY', true);
    addParameter(p, 'fitZ', false);
    addParameter(p, 'initialCoeffs', struct());
    addParameter(p, 'maxIter', 8);
    addParameter(p, 'fdCoeff', 0.01);
    addParameter(p, 'regCoeff', 1e-5);
    addParameter(p, 'maxCoeffStep', 0.04);
    addParameter(p, 'tolStep', 1e-5);

    addParameter(p, 'photonsPerPlane', 5e4);
    addParameter(p, 'backgroundCountsPerChannel', 0);
    addParameter(p, 'addNoise', true);
    addParameter(p, 'nNoiseReplicates', 3);
    addParameter(p, 'randomSeed', 1);
    addParameter(p, 'centerSampleMode', 'subpixel');
    addParameter(p, 'centerNormalization', 'perPlane');
    addParameter(p, 'modelSampleXY', [0 0]);

    addParameter(p, 'trueXY', [0 0]);
    addParameter(p, 'trueZOffset', 0);
    addParameter(p, 'sim', []);
    addParameter(p, 'objectiveNA', []);
    addParameter(p, 'emissionWavelengthUm', []);
    addParameter(p, 'excitationWavelengthUm', []);
    addParameter(p, 'mediumRefractiveIndex', []);
    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detectorPixelShape', 'hex');
    addParameter(p, 'detectorSubsamples', []);
    addParameter(p, 'modelDz', 0.05);
    addParameter(p, 'modelZPadding', 0.50);

    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'continueOnError', true);
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;

    opts.defocusUm = double(opts.defocusUm(:)).';
    if isempty(opts.defocusUm) || any(~isfinite(opts.defocusUm)) || any(opts.defocusUm <= 0)
        error('simulateCenterPointISMDefocusSweep:BadDefocus', ...
            'defocusUm must contain positive finite distances in micrometers.');
    end
    opts.truthSet = lower(char(opts.truthSet));
    opts.centerSampleMode = lower(char(opts.centerSampleMode));
    opts.centerNormalization = lower(char(opts.centerNormalization));
    opts.modelSampleXY = double(opts.modelSampleXY(:)).';
    if numel(opts.modelSampleXY) ~= 2 || any(~isfinite(opts.modelSampleXY))
        error('simulateCenterPointISMDefocusSweep:BadModelSampleXY', ...
            'modelSampleXY must be finite [x y] coordinates in micrometers.');
    end
    if ischar(opts.fitModes) || isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    end
    opts.fitModes = opts.fitModes(:).';
    opts.nNoiseReplicates = max(1, round(opts.nNoiseReplicates));
    opts.trueXY = double(opts.trueXY(:)).';
    if numel(opts.trueXY) ~= 2
        error('simulateCenterPointISMDefocusSweep:BadTrueXY', ...
            'trueXY must be [x y] in micrometers.');
    end
end

function addRequiredPaths()
    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    parentDir = fileparts(thisDir);
    if exist(parentDir, 'dir') == 7
        addpath(parentDir);
    end
end

function sim = configureBaseSim(opts)
    if ~isempty(opts.sim)
        sim = opts.sim;
    else
        sim = defaultParams();
        sim.detectorLayout = opts.detectorLayout;
        sim.detectorPixelShape = opts.detectorPixelShape;
        [sim.detXY, sim.detectorIndexGrid, sim.detectorLayoutInfo] = ...
            detectorLayout(sim.detectorLayout, sim.detPitch);
        sim.nDet = size(sim.detXY, 1);
        sim.detectorGridSize = size(sim.detectorIndexGrid);
        sim.arrayN = sim.detectorGridSize(1);
        if strcmpi(sim.detectorPixelShape, 'hex')
            sim.detectorHexRadius = sim.detSize / sqrt(3);
        end
    end

    if ~isempty(opts.detectorSubsamples)
        sim.detectorSubsamples = opts.detectorSubsamples;
    end
    sim = applyOpticsOptions(sim, opts);
end

function sim = applyOpticsOptions(sim, opts)
    if isfinitePositiveScalar(opts.objectiveNA)
        sim.NA = double(opts.objectiveNA);
    end
    if isfinitePositiveScalar(opts.excitationWavelengthUm)
        sim.lamExc = double(opts.excitationWavelengthUm);
    end
    if isfinitePositiveScalar(opts.emissionWavelengthUm)
        sim.lamEm = double(opts.emissionWavelengthUm);
        sim.lamRef = double(opts.emissionWavelengthUm);
    end
    if isfinitePositiveScalar(opts.mediumRefractiveIndex)
        sim.nMedium = double(opts.mediumRefractiveIndex);
    end
end

function sim = prepareSimForPlaneZ(sim, planeZ, opts)
    zMin = min([planeZ(:); 0]) - opts.modelZPadding;
    zMax = max([planeZ(:); 0]) + opts.modelZPadding;
    nZ = ceil((zMax - zMin) / opts.modelDz) + 1;
    if mod(nZ, 2) == 0
        nZ = nZ + 1;
    end
    sim.z = linspace(zMin, zMax, nZ);
    sim.nz = numel(sim.z);
    sim.nzRange = zMax - zMin;
    sim.obj = beadObject3D(sim);
end

function cases = buildTruthCases(opts, sim)
    if ~isempty(opts.truthCoeffs)
        cases = truthCasesFromUserInput(opts.truthCoeffs, sim, opts);
        return;
    end

    switch opts.truthSet
        case 'mixed'
            cases = makeCase('mixed', defaultMixedCoeffs(opts), sim);
        case {'singlemodes','single_modes','single'}
            cases = repmat(emptyTruthCase(), 0, 1);
            for k = 1:numel(opts.fitModes)
                coeffs = struct();
                coeffs.(opts.fitModes{k}) = opts.singleModeAmplitudeWaves;
                cases(end+1, 1) = makeCase(opts.fitModes{k}, coeffs, sim); %#ok<AGROW>
            end
        case {'mixedandsinglemodes','mixed_and_single_modes','all'}
            cases = makeCase('mixed', defaultMixedCoeffs(opts), sim);
            for k = 1:numel(opts.fitModes)
                coeffs = struct();
                coeffs.(opts.fitModes{k}) = opts.singleModeAmplitudeWaves;
                cases(end+1, 1) = makeCase(opts.fitModes{k}, coeffs, sim); %#ok<AGROW>
            end
        otherwise
            error('simulateCenterPointISMDefocusSweep:BadTruthSet', ...
                'truthSet must be mixed, singleModes, or mixedAndSingleModes.');
    end
end

function coeffs = defaultMixedCoeffs(opts)
    scale = opts.mixedCoeffScale;
    coeffs = struct();
    coeffs.defocus = 0.030 * scale;
    coeffs.astig_x = 0.025 * scale;
    coeffs.astig_y = -0.020 * scale;
    coeffs.coma_x = 0.020 * scale;
    coeffs.coma_y = -0.015 * scale;
    coeffs.spherical = 0.025 * scale;
end

function cases = truthCasesFromUserInput(inputValue, sim, opts)
    cases = repmat(emptyTruthCase(), 0, 1);
    if isnumeric(inputValue)
        values = double(inputValue);
        if isvector(values)
            values = values(:).';
        end
        if size(values, 2) == numel(opts.fitModes)
            names = opts.fitModes;
        else
            names = sim.modeOrder;
        end
        for k = 1:size(values, 1)
            coeffs = vectorToCoeffStruct(values(k,:), names);
            cases(end+1, 1) = makeCase(sprintf('truth_%d', k), coeffs, sim); %#ok<AGROW>
        end
        return;
    end

    if isstruct(inputValue)
        for k = 1:numel(inputValue)
            coeffs = inputValue(k);
            name = sprintf('truth_%d', k);
            if isfield(coeffs, 'caseName')
                name = char(coeffs.caseName);
                coeffs = rmfield(coeffs, 'caseName');
            elseif isfield(coeffs, 'name')
                name = char(coeffs.name);
                coeffs = rmfield(coeffs, 'name');
            end
            cases(end+1, 1) = makeCase(name, coeffs, sim); %#ok<AGROW>
        end
        return;
    end

    if iscell(inputValue)
        for k = 1:numel(inputValue)
            item = inputValue{k};
            if isnumeric(item)
                item = double(item(:)).';
                if numel(item) == numel(opts.fitModes)
                    names = opts.fitModes;
                else
                    names = sim.modeOrder;
                end
                coeffs = vectorToCoeffStruct(item, names);
            elseif isstruct(item)
                coeffs = item;
            else
                error('simulateCenterPointISMDefocusSweep:BadTruthCoeffs', ...
                    'truthCoeffs cell entries must be numeric vectors or structs.');
            end
            cases(end+1, 1) = makeCase(sprintf('truth_%d', k), coeffs, sim); %#ok<AGROW>
        end
        return;
    end

    error('simulateCenterPointISMDefocusSweep:BadTruthCoeffs', ...
        'truthCoeffs must be a numeric vector/matrix, struct array, or cell array.');
end

function c = emptyTruthCase()
    c = struct('name', '', 'coeffs', struct(), 'coeffVector', []);
end

function c = makeCase(name, coeffs, sim)
    c = emptyTruthCase();
    c.name = char(name);
    c.coeffs = coeffs;
    c.coeffVector = coeffStructToVector(sim, coeffs);
end

function coeffs = vectorToCoeffStruct(values, names)
    values = double(values(:)).';
    coeffs = struct();
    n = min(numel(values), numel(names));
    for k = 1:n
        if abs(values(k)) > 1e-15
            coeffs.(names{k}) = values(k);
        end
    end
end

function [centerValues, centerNorm] = simulateCenterMicroimage(sim, coeffs, planeZ, opts)
    stack = normalizedStackExplicitDetectorZPlanes( ...
        sim, coeffs, planeZ, opts.trueXY(1), opts.trueXY(2), opts.trueZOffset);
    centerValues = sampleModelStackAtXY(stack, sim, opts.modelSampleXY, ...
        opts.centerSampleMode, numel(planeZ));
    centerNorm = normalizeCenterValues(centerValues, opts.centerNormalization);
end

function values = sampleModelStackAtXY(stack, sim, sampleXY, mode, nPlane)
    if nargin < 5
        nPlane = size(stack, 4);
    end

    x = double(sampleXY(1));
    y = double(sampleXY(2));
    x = min(max(x, min(sim.x(:))), max(sim.x(:)));
    y = min(max(y, min(sim.y(:))), max(sim.y(:)));

    switch lower(char(mode))
        case {'nearest','round','pixel'}
            [~, ix] = min(abs(sim.x - x));
            [~, iy] = min(abs(sim.y - y));
            values = squeeze(stack(iy, ix, :, :));
            values = reshape(double(values), sim.nDet, nPlane);
        case {'subpixel','linear','bilinear','interp'}
            values = zeros(sim.nDet, nPlane);
            for ip = 1:nPlane
                for ch = 1:sim.nDet
                    img = double(stack(:,:,ch,ip));
                    values(ch, ip) = interp2(sim.x, sim.y, img, x, y, 'linear', 0);
                end
            end
        otherwise
            error('simulateCenterPointISMDefocusSweep:BadCenterSampleMode', ...
                'centerSampleMode must be subpixel or nearest.');
    end
end

function dataN = normalizeCenterValues(values, mode)
    values = max(double(values), 0);
    switch lower(mode)
        case {'perplane','plane','eachplane'}
            scale = sum(values, 1);
            scale(scale <= 0 | ~isfinite(scale)) = 1;
            dataN = values ./ reshape(scale, 1, []);
        case {'global','all'}
            scale = sum(values(:));
            if scale <= 0 || ~isfinite(scale)
                scale = 1;
            end
            dataN = values / scale;
        case {'none','raw'}
            dataN = values;
        otherwise
            error('simulateCenterPointISMDefocusSweep:BadNormalization', ...
                'centerNormalization must be perPlane, global, or none.');
    end
end

function observed = addNoiseToCenterData(cleanNorm, opts)
    expected = expectedCountsFromNormalizedData(cleanNorm, opts);
    if opts.addNoise
        observed = poissonSample(expected);
    else
        observed = expected;
    end
end

function expected = expectedCountsFromNormalizedData(cleanNorm, opts)
    nCh = size(cleanNorm, 1);
    nPlane = size(cleanNorm, 2);
    photons = double(opts.photonsPerPlane);
    if isscalar(photons)
        photons = repmat(photons, 1, nPlane);
    end
    if numel(photons) ~= nPlane
        error('simulateCenterPointISMDefocusSweep:BadPhotons', ...
            'photonsPerPlane must be scalar or contain one value per plane.');
    end

    bg = double(opts.backgroundCountsPerChannel);
    if isscalar(bg)
        bg = repmat(bg, nCh, nPlane);
    elseif isvector(bg) && numel(bg) == nCh
        bg = repmat(bg(:), 1, nPlane);
    elseif ~isequal(size(bg), [nCh nPlane])
        error('simulateCenterPointISMDefocusSweep:BadBackground', ...
            'backgroundCountsPerChannel must be scalar, 23-vector, or [23 x nPlane].');
    end

    expected = cleanNorm .* reshape(photons, 1, []) + bg;
    expected = max(expected, 0);
end

function args = retrievalArgs(opts, sim, planeZ)
    args = { ...
        'planeZ', planeZ, ...
        'sim', sim, ...
        'fitModes', opts.fitModes, ...
        'fitStrategy', 'joint', ...
        'fitXY', opts.fitXY, ...
        'fitZ', opts.fitZ, ...
        'initialCoeffs', opts.initialCoeffs, ...
        'maxIter', opts.maxIter, ...
        'fdCoeff', opts.fdCoeff, ...
        'regCoeff', opts.regCoeff, ...
        'maxCoeffStep', opts.maxCoeffStep, ...
        'tolStep', opts.tolStep, ...
        'modelDz', opts.modelDz, ...
        'modelZPadding', opts.modelZPadding, ...
        'centerSampleMode', opts.centerSampleMode, ...
        'centerNormalization', opts.centerNormalization, ...
        'modelSampleXY', opts.modelSampleXY, ...
        'darkMode', 'none', ...
        'subtractBoundary', false, ...
        'writeOutputs', false, ...
        'verbose', false};
end

function trial = scoreTrial(trial, sim, fitModes)
    fitMask = ismember(sim.modeOrder, fitModes);
    trueVec = trial.trueCoeffVector(:).';
    estVec = trial.estimatedCoeffVector(:).';
    if numel(estVec) ~= numel(trueVec)
        estVec = nan(size(trueVec));
    end
    err = estVec - trueVec;
    trial.coeffErrorWaves = err;
    trial.coeffRmseWaves = sqrt(mean(err(fitMask).^2, 'omitnan'));
    trial.coeffMaeWaves = mean(abs(err(fitMask)), 'omitnan');
    trial.coeffRmseNm = trial.coeffRmseWaves * sim.lamRef * 1000;
    trial.coeffMaeNm = trial.coeffMaeWaves * sim.lamRef * 1000;
end

function trial = emptyTrial()
    trial = struct();
    trial.defocusUm = NaN;
    trial.defocusNm = NaN;
    trial.caseName = '';
    trial.replicate = NaN;
    trial.planeZ = [];
    trial.trueCoeffVector = [];
    trial.estimatedCoeffVector = [];
    trial.coeffErrorWaves = [];
    trial.coeffRmseWaves = NaN;
    trial.coeffMaeWaves = NaN;
    trial.coeffRmseNm = NaN;
    trial.coeffMaeNm = NaN;
    trial.cleanCenter = [];
    trial.cleanCenterNormalized = [];
    trial.observedCounts = [];
    trial.modelCenter = [];
    trial.residualNorm = NaN;
    trial.rank = NaN;
    trial.nParameters = NaN;
    trial.isFullRank = false;
    trial.conditionNumber = NaN;
    trial.minSingularValue = NaN;
    trial.fitResult = [];
    trial.errorMessage = '';
end

function T = makeTrialTable(trials)
    n = numel(trials);
    defocusUm = zeros(n,1);
    defocusNm = zeros(n,1);
    caseName = cell(n,1);
    replicate = zeros(n,1);
    coeffRmseWaves = zeros(n,1);
    coeffRmseNm = zeros(n,1);
    coeffMaeWaves = zeros(n,1);
    coeffMaeNm = zeros(n,1);
    residualNorm = zeros(n,1);
    rank = zeros(n,1);
    nParameters = zeros(n,1);
    isFullRank = false(n,1);
    conditionNumber = zeros(n,1);
    minSingularValue = zeros(n,1);
    errorMessage = cell(n,1);
    for k = 1:n
        defocusUm(k) = trials(k).defocusUm;
        defocusNm(k) = trials(k).defocusNm;
        caseName{k} = trials(k).caseName;
        replicate(k) = trials(k).replicate;
        coeffRmseWaves(k) = trials(k).coeffRmseWaves;
        coeffRmseNm(k) = trials(k).coeffRmseNm;
        coeffMaeWaves(k) = trials(k).coeffMaeWaves;
        coeffMaeNm(k) = trials(k).coeffMaeNm;
        residualNorm(k) = trials(k).residualNorm;
        rank(k) = trials(k).rank;
        nParameters(k) = trials(k).nParameters;
        isFullRank(k) = trials(k).isFullRank;
        conditionNumber(k) = trials(k).conditionNumber;
        minSingularValue(k) = trials(k).minSingularValue;
        errorMessage{k} = trials(k).errorMessage;
    end
    T = table(defocusUm, defocusNm, caseName, replicate, coeffRmseWaves, ...
        coeffRmseNm, coeffMaeWaves, coeffMaeNm, residualNorm, rank, ...
        nParameters, isFullRank, conditionNumber, minSingularValue, errorMessage);
end

function T = makeSummaryTable(trials, defocusValues)
    n = numel(defocusValues);
    defocusUm = defocusValues(:);
    defocusNm = defocusUm * 1000;
    nTrials = zeros(n,1);
    meanCoeffRmseWaves = nan(n,1);
    medianCoeffRmseWaves = nan(n,1);
    stdCoeffRmseWaves = nan(n,1);
    meanCoeffRmseNm = nan(n,1);
    medianCoeffRmseNm = nan(n,1);
    stdCoeffRmseNm = nan(n,1);
    meanResidualNorm = nan(n,1);
    medianRank = nan(n,1);
    fractionFullRank = nan(n,1);
    medianConditionNumber = nan(n,1);
    medianMinSingularValue = nan(n,1);

    allDz = [trials.defocusUm];
    for k = 1:n
        idx = abs(allDz - defocusValues(k)) < 1e-12;
        subset = trials(idx);
        nTrials(k) = numel(subset);
        if isempty(subset)
            continue;
        end
        rmseW = [subset.coeffRmseWaves].';
        rmseN = [subset.coeffRmseNm].';
        residual = [subset.residualNorm].';
        ranks = [subset.rank].';
        fullRank = [subset.isFullRank].';
        condNum = [subset.conditionNumber].';
        minSing = [subset.minSingularValue].';
        meanCoeffRmseWaves(k) = mean(rmseW, 'omitnan');
        medianCoeffRmseWaves(k) = median(rmseW, 'omitnan');
        stdCoeffRmseWaves(k) = std(rmseW, 'omitnan');
        meanCoeffRmseNm(k) = mean(rmseN, 'omitnan');
        medianCoeffRmseNm(k) = median(rmseN, 'omitnan');
        stdCoeffRmseNm(k) = std(rmseN, 'omitnan');
        meanResidualNorm(k) = mean(residual, 'omitnan');
        medianRank(k) = median(ranks, 'omitnan');
        fractionFullRank(k) = mean(double(fullRank), 'omitnan');
        medianConditionNumber(k) = median(condNum, 'omitnan');
        medianMinSingularValue(k) = median(minSing, 'omitnan');
    end

    T = table(defocusUm, defocusNm, nTrials, meanCoeffRmseWaves, ...
        medianCoeffRmseWaves, stdCoeffRmseWaves, meanCoeffRmseNm, ...
        medianCoeffRmseNm, stdCoeffRmseNm, meanResidualNorm, medianRank, ...
        fractionFullRank, medianConditionNumber, medianMinSingularValue);
end

function [bestByError, bestByInformation] = selectBestDefocus(summaryTable)
    bestByError = summaryTable([], :);
    bestByInformation = summaryTable([], :);
    if isempty(summaryTable)
        return;
    end

    err = summaryTable.meanCoeffRmseNm;
    [~, idx] = min(err);
    if ~isempty(idx) && isfinite(err(idx))
        bestByError = summaryTable(idx, :);
    end

    info = summaryTable.medianMinSingularValue;
    fullRank = summaryTable.fractionFullRank >= 1;
    if any(fullRank & isfinite(info))
        candidates = find(fullRank & isfinite(info));
        [~, ii] = max(info(candidates));
        bestByInformation = summaryTable(candidates(ii), :);
    else
        [~, idx] = max(info);
        if ~isempty(idx) && isfinite(info(idx))
            bestByInformation = summaryTable(idx, :);
        end
    end
end

function outDir = resolveOutputDir(opts)
    if ~isempty(opts.outputDir)
        outDir = char(opts.outputDir);
        return;
    end
    thisDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(thisDir);
    outDir = fullfile(rootDir, 'output_matlab', 'center_point_defocus_sweep');
end

function writeSweepOutputs(results)
    outDir = results.outputDir;
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    writetable(results.trialTable, fullfile(outDir, 'center_point_defocus_sweep_trials.csv'));
    writetable(results.summaryTable, fullfile(outDir, 'center_point_defocus_sweep_summary.csv'));
    writeBestMicroimageCsv(results, fullfile(outDir, 'best_defocus_center_microimages.csv'));
    save(fullfile(outDir, 'center_point_defocus_sweep_results.mat'), 'results', '-v7.3');
    writeSweepFigure(results, fullfile(outDir, 'center_point_defocus_sweep_summary.png'));
end

function writeBestMicroimageCsv(results, outFile)
    if isempty(results.bestByError)
        return;
    end
    bestDz = results.bestByError.defocusUm(1);
    idx = find(abs([results.trials.defocusUm] - bestDz) < 1e-12, 1, 'first');
    if isempty(idx)
        return;
    end
    tr = results.trials(idx);
    nCh = size(tr.observedCounts, 1);
    planeIndex = [ones(nCh,1); 2*ones(nCh,1)];
    zUm = [zeros(nCh,1); bestDz*ones(nCh,1)];
    detectorIndex = repmat((1:nCh).', 2, 1);
    cleanNormalized = [tr.cleanCenterNormalized(:,1); tr.cleanCenterNormalized(:,2)];
    observedCounts = [tr.observedCounts(:,1); tr.observedCounts(:,2)];
    fittedNormalized = [tr.modelCenter(:,1); tr.modelCenter(:,2)];
    T = table(planeIndex, zUm, detectorIndex, cleanNormalized, ...
        observedCounts, fittedNormalized);
    writetable(T, outFile);
end

function writeSweepFigure(results, outFile)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1250 820]);
    tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    S = results.summaryTable;

    ax = nexttile(tl, 1);
    errorbar(ax, S.defocusNm, S.meanCoeffRmseNm, S.stdCoeffRmseNm, ...
        '-o', 'LineWidth', 1.2);
    xlabel(ax, 'defocus from focal plane (nm)');
    ylabel(ax, 'coefficient RMSE (nm RMS)');
    title(ax, 'Retrieval error');
    grid(ax, 'on');

    ax = nexttile(tl, 2);
    plot(ax, S.defocusNm, S.medianMinSingularValue, '-o', 'LineWidth', 1.2);
    xlabel(ax, 'defocus from focal plane (nm)');
    ylabel(ax, 'median min singular value');
    title(ax, 'Local information strength');
    grid(ax, 'on');

    ax = nexttile(tl, 3);
    semilogy(ax, S.defocusNm, S.medianConditionNumber, '-o', 'LineWidth', 1.2);
    xlabel(ax, 'defocus from focal plane (nm)');
    ylabel(ax, 'median condition number');
    title(ax, 'Local conditioning');
    grid(ax, 'on');

    ax = nexttile(tl, 4);
    [caseNames, defocusNm, heat] = caseDefocusErrorMatrix(results.trials);
    if ~isempty(heat)
        imagesc(ax, defocusNm, 1:numel(caseNames), heat);
        set(ax, 'YTick', 1:numel(caseNames), 'YTickLabel', caseNames);
        xlabel(ax, 'defocus from focal plane (nm)');
        ylabel(ax, 'truth case');
        title(ax, 'Mean coefficient RMSE by case (nm)');
        colorbar(ax);
    else
        text(ax, 0.5, 0.5, 'No successful trials', 'HorizontalAlignment', 'center');
        axis(ax, 'off');
    end

    exportFigure(fig, outFile);
end

function [caseNames, defocusNm, heat] = caseDefocusErrorMatrix(trials)
    if isempty(trials)
        caseNames = {};
        defocusNm = [];
        heat = [];
        return;
    end
    caseNames = unique({trials.caseName}, 'stable');
    defocusVals = unique([trials.defocusUm], 'stable');
    defocusNm = defocusVals * 1000;
    heat = nan(numel(caseNames), numel(defocusVals));
    for ic = 1:numel(caseNames)
        for iz = 1:numel(defocusVals)
            idx = strcmp({trials.caseName}, caseNames{ic}) & ...
                abs([trials.defocusUm] - defocusVals(iz)) < 1e-12;
            vals = [trials(idx).coeffRmseNm];
            heat(ic, iz) = mean(vals, 'omitnan');
        end
    end
end

function exportFigure(fig, outFile)
    try
        exportgraphics(fig, outFile, 'Resolution', 180);
    catch
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, outFile, '-dpng', '-r180');
    end
    close(fig);
end

function printSweepSummary(results)
    fprintf('\nCenter-point ISM defocus sweep\n');
    if ~isempty(results.bestByError)
        fprintf('  best by retrieval error: %.0f nm defocus, mean RMSE %.3g nm\n', ...
            results.bestByError.defocusNm(1), results.bestByError.meanCoeffRmseNm(1));
    end
    if ~isempty(results.bestByInformation)
        fprintf('  best by min singular value: %.0f nm defocus, median min singular %.3g\n', ...
            results.bestByInformation.defocusNm(1), ...
            results.bestByInformation.medianMinSingularValue(1));
    end
    fprintf('  outputs: %s\n\n', results.outputDir);
end

function vec = coeffStructToVector(sim, coeffs)
    if isnumeric(coeffs)
        raw = double(coeffs(:)).';
        vec = zeros(1, numel(sim.modeOrder));
        vec(1:min(numel(raw), numel(vec))) = raw(1:min(numel(raw), numel(vec)));
        return;
    end
    vec = zeros(1, numel(sim.modeOrder));
    for k = 1:numel(sim.modeOrder)
        if isfield(coeffs, sim.modeOrder{k})
            vec(k) = coeffs.(sim.modeOrder{k});
        end
    end
end

function tf = isfinitePositiveScalar(v)
    tf = isnumeric(v) && isscalar(v) && isfinite(v) && v > 0;
end
