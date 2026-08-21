function result = simulateRetrieveCenterMicroimageAberrations(varargin)
%SIMULATERETRIEVECENTERMICROIMAGEABERRATIONS Center-microimage aberration test.
%
%   result = simulateRetrieveCenterMicroimageAberrations()
%
%   Simulates the detector microimage at the center of a point emitter placed
%   at a glass-water interface, then attempts to recover coefficients up to
%   primary spherical aberration using only that center detector intensity
%   distribution.
%
%   The default forward model is the repository's layered interface
%   vectorial Richards-Wolf/Bessel model with xyAverage/isotropic collection
%   as the rapidly rotating dipole approximation. This is not an
%   oriented-dipole interface Green-tensor model. Pass
%   'sampleGeometry','homogeneous' to reproduce the homogeneous-medium test.
%   Synthetic observations are SPAD counts by default: the script samples a
%   Poisson model with fixed photonsPerPlane plus darkCountsPerDetector, fits
%   by Poisson deviance, and exports the expected Fisher/CRB precision bounds.
%
%   By default the fit uses supported Zernike modes through primary spherical:
%       tilt_x, tilt_y, defocus, astig_*, coma_*, trefoil_*,
%       spherical.
%
%   Example:
%       result = simulateRetrieveCenterMicroimageAberrations();
%
%       result = simulateRetrieveCenterMicroimageAberrations( ...
%           'planeZ', [0 0.3], 'photonsPerPlane', 2e5);
%
%       result = simulateRetrieveCenterMicroimageAberrations( ...
%           'planeZ', [0 0.3], 'photonsPerPlane', 2e5, ...
%           'monteCarloRepeats', 100);
%
%       result = simulateRetrieveCenterMicroimageAberrations( ...
%           'fitObjective', 'normalizedLS', 'photonsPerPlane', Inf);
%
%       result = simulateRetrieveCenterMicroimageAberrations( ...
%           'sampleGeometry', 'homogeneous');
%
%       result = simulateRetrieveCenterMicroimageAberrations( ...
%           'trueCoeffs', struct('coma_x',0.04,'astig_y',-0.03), ...
%           'fitModes', {'coma_x','astig_y','spherical'});

    addpath(fileparts(mfilename('fullpath')));
    opts = parseInputs(varargin{:});
    sim = configureSimulation(opts);
    fitModes = resolveFitModes(opts);

    trueCoeffVector = coefficientVectorFromInput( ...
        sim, fitModes, opts.trueCoeffs, defaultTrueCoeffs());
    trueCoeffs = coeffsFromFitVector(sim, fitModes, trueCoeffVector);
    initialVector = coefficientVectorFromInput( ...
        sim, fitModes, opts.initialCoeffs, struct());

    truth = centerMicroimageModel(sim, trueCoeffs, opts.planeZ, opts.sampleXY);
    [dataCounts, photonInfo] = makeSyntheticData(truth.values, opts);
    dataN = normalizeMicroimage(dataCounts, opts.normalization);
    trueN = normalizeMicroimage(truth.values, opts.normalization);

    if strcmpi(opts.fitObjective, 'poisson')
        fit = fitCenterMicroimagePoisson( ...
            sim, fitModes, dataCounts, initialVector, opts);
    else
        fit = fitCenterMicroimage(sim, fitModes, dataN, initialVector, opts);
    end
    fitCoeffs = coeffsFromFitVector(sim, fitModes, fit.coeffVector);
    fitRaw = centerMicroimageModel(sim, fitCoeffs, opts.planeZ, opts.sampleXY);
    fitN = normalizeMicroimage(fitRaw.values, opts.normalization);
    if isfinite(opts.photonsPerPlane) && opts.photonsPerPlane > 0
        [fitExpectedCounts, ~] = expectedCountsFromVector( ...
            sim, fitModes, fit.coeffVector, opts);
        poissonFisher = poissonFisherDiagnostics( ...
            sim, fitModes, fit.coeffVector, opts);
    else
        fitExpectedCounts = [];
        poissonFisher = emptyPoissonFisher(numel(fitModes));
    end

    identifiability = centerMicroimageIdentifiability( ...
        sim, fitModes, fit.coeffVector, opts);
    monteCarlo = runMonteCarloCoefficientRecovery( ...
        sim, fitModes, truth.values, trueCoeffVector, initialVector, opts);
    coeffTable = coefficientRecoveryTable(fitModes, trueCoeffVector, ...
        fit.coeffVector, poissonFisher.crbStdWaves, monteCarlo);

    result = struct();
    result.modelDescription = modelDescriptionForSim(sim);
    result.sim = sim;
    result.options = opts;
    result.fitModes = fitModes;
    result.trueCoeffs = trueCoeffs;
    result.trueCoeffVector = trueCoeffVector(:);
    result.trueMicroimage = truth;
    result.trueNormalizedMicroimage = trueN;
    result.syntheticData = dataCounts;
    result.syntheticPhotonInfo = photonInfo;
    result.expectedCountsAtTruth = photonInfo.expectedCounts;
    result.dataNormalized = dataN;
    result.fit = fit;
    result.estimatedCoeffs = fitCoeffs;
    result.estimatedMicroimage = fitRaw;
    result.estimatedNormalizedMicroimage = fitN;
    result.estimatedExpectedCounts = fitExpectedCounts;
    result.coefficientTable = coeffTable;
    result.identifiability = identifiability;
    result.poissonFisher = poissonFisher;
    result.monteCarlo = monteCarlo;
    result.recommendation = recommendationFromIdentifiability(identifiability, opts);
    result.outputDir = resolveOutputDir(opts);

    if opts.writeOutputs
        writeOutputs(result);
    end

    if opts.verbose
        printSummary(result);
    end
end

function opts = parseInputs(varargin)
    p = inputParser;
    p.FunctionName = 'simulateRetrieveCenterMicroimageAberrations';

    addParameter(p, 'maxZernikeOrder', 4);
    addParameter(p, 'fitModes', {});
    addParameter(p, 'includeTilt', true);
    addParameter(p, 'trueCoeffs', []);
    addParameter(p, 'initialCoeffs', []);

    addParameter(p, 'planeZ', 0);
    addParameter(p, 'sampleXY', [0 0]);
    addParameter(p, 'normalization', 'perPlane');
    addParameter(p, 'fitObjective', 'poisson');
    addParameter(p, 'photonsPerPlane', 1e4);
    addParameter(p, 'poissonSampleData', true);
    addParameter(p, 'darkCountsPerDetector', 0);
    addParameter(p, 'minExpectedCount', 1e-9);
    addParameter(p, 'marginalizeFlux', true);
    addParameter(p, 'monteCarloRepeats', 0);
    addParameter(p, 'rngSeed', 1);

    addParameter(p, 'sampleGeometry', 'glassWaterInterface');
    addParameter(p, 'interfaceEmitterHeightUm', 0);
    addParameter(p, 'interfaceRadialWeightMode', 'sampleSolidAngle');
    addParameter(p, 'airInterfaceStageMedium', 'immersion');
    addParameter(p, 'immersionRefractiveIndex', 1.33);
    addParameter(p, 'glassRefractiveIndex', 1.518);
    addParameter(p, 'sampleRefractiveIndex', 1.33);
    addParameter(p, 'designGlassRefractiveIndex', 1.518);
    addParameter(p, 'coverslipThicknessUm', 190);
    addParameter(p, 'designCoverslipThicknessUm', 190);

    addParameter(p, 'fovXY', 1.8);
    addParameter(p, 'nx', 31);
    addParameter(p, 'zRange', 1.2);
    addParameter(p, 'nz', 9);
    addParameter(p, 'Nr', 36);
    addParameter(p, 'Nphi', 72);
    addParameter(p, 'M', 6);

    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detPitchUm', []);
    addParameter(p, 'detFillRatio', 1.0);
    addParameter(p, 'detectorSubsamples', 7);

    addParameter(p, 'objectiveNA', 1.2);
    addParameter(p, 'mediumRefractiveIndex', 1.33);
    addParameter(p, 'excitationWavelengthUm', 0.488);
    addParameter(p, 'emissionWavelengthUm', 0.520);
    addParameter(p, 'referenceWavelengthUm', 0.520);
    addParameter(p, 'excitationPolarizationMode', 'xyAverage');
    addParameter(p, 'collectionPolarizationMode', 'xyAverage');

    addParameter(p, 'maxIter', 8);
    addParameter(p, 'fdCoeff', 0.005);
    addParameter(p, 'regCoeff', 1e-6);
    addParameter(p, 'initialDamping', 1e-3);
    addParameter(p, 'maxCoeffStep', 0.04);
    addParameter(p, 'stepTolerance', 1e-5);
    addParameter(p, 'residualTolerance', 1e-9);
    addParameter(p, 'svdRelativeTolerance', 1e-4);
    addParameter(p, 'maxUsefulCondition', 1e4);

    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;

    opts.maxZernikeOrder = positiveInteger(opts.maxZernikeOrder, 'maxZernikeOrder');
    opts.includeTilt = logical(opts.includeTilt);
    opts.planeZ = double(opts.planeZ(:).');
    if isempty(opts.planeZ) || any(~isfinite(opts.planeZ))
        error('simulateRetrieveCenterMicroimageAberrations:BadPlaneZ', ...
            'planeZ must contain at least one finite axial position in micrometers.');
    end
    opts.sampleXY = double(opts.sampleXY(:).');
    if numel(opts.sampleXY) ~= 2
        error('simulateRetrieveCenterMicroimageAberrations:BadSampleXY', ...
            'sampleXY must be [x y] in micrometers.');
    end
    opts.normalization = lower(char(opts.normalization));
    if ~ismember(opts.normalization, {'perplane','global','none'})
        error('simulateRetrieveCenterMicroimageAberrations:BadNormalization', ...
            'normalization must be perPlane, global, or none.');
    end
    opts.fitObjective = lower(char(opts.fitObjective));
    if ismember(opts.fitObjective, {'ls','leastsquares','normalized','normalizedls'})
        opts.fitObjective = 'normalizedls';
    end
    if ~ismember(opts.fitObjective, {'poisson','normalizedls'})
        error('simulateRetrieveCenterMicroimageAberrations:BadFitObjective', ...
            'fitObjective must be poisson or normalizedLS.');
    end
    opts.photonsPerPlane = double(opts.photonsPerPlane);
    if ~isscalar(opts.photonsPerPlane) || isnan(opts.photonsPerPlane) || ...
            opts.photonsPerPlane < 0
        error('simulateRetrieveCenterMicroimageAberrations:BadPhotons', ...
            'photonsPerPlane must be a nonnegative scalar or Inf.');
    end
    if strcmpi(opts.fitObjective, 'poisson') && ...
            (~isfinite(opts.photonsPerPlane) || opts.photonsPerPlane <= 0)
        error('simulateRetrieveCenterMicroimageAberrations:BadPoissonPhotons', ...
            'Poisson fitting requires finite positive photonsPerPlane.');
    end
    opts.poissonSampleData = logical(opts.poissonSampleData);
    opts.marginalizeFlux = logical(opts.marginalizeFlux);
    opts.monteCarloRepeats = nonnegativeInteger(opts.monteCarloRepeats, ...
        'monteCarloRepeats');
    if opts.monteCarloRepeats > 0 && ...
            (~isfinite(opts.photonsPerPlane) || opts.photonsPerPlane <= 0)
        error('simulateRetrieveCenterMicroimageAberrations:BadMonteCarloPhotons', ...
            'Monte Carlo repeats require finite positive photonsPerPlane.');
    end
    opts.darkCountsPerDetector = double(opts.darkCountsPerDetector(:));
    if isempty(opts.darkCountsPerDetector) || ...
            any(~isfinite(opts.darkCountsPerDetector)) || ...
            any(opts.darkCountsPerDetector < 0)
        error('simulateRetrieveCenterMicroimageAberrations:BadDarkCounts', ...
            'darkCountsPerDetector must be a nonnegative scalar or vector.');
    end
    opts.minExpectedCount = positiveScalar(opts.minExpectedCount, ...
        'minExpectedCount');
    opts.sampleGeometry = normalizeSampleGeometry(opts.sampleGeometry);
    opts.interfaceEmitterHeightUm = nonnegativeScalar( ...
        opts.interfaceEmitterHeightUm, 'interfaceEmitterHeightUm');
    opts.interfaceRadialWeightMode = char(opts.interfaceRadialWeightMode);
    opts.airInterfaceStageMedium = char(opts.airInterfaceStageMedium);
    opts.immersionRefractiveIndex = positiveScalar( ...
        opts.immersionRefractiveIndex, 'immersionRefractiveIndex');
    opts.glassRefractiveIndex = positiveScalar( ...
        opts.glassRefractiveIndex, 'glassRefractiveIndex');
    opts.sampleRefractiveIndex = positiveScalar( ...
        opts.sampleRefractiveIndex, 'sampleRefractiveIndex');
    opts.designGlassRefractiveIndex = positiveScalar( ...
        opts.designGlassRefractiveIndex, 'designGlassRefractiveIndex');
    opts.coverslipThicknessUm = positiveScalar( ...
        opts.coverslipThicknessUm, 'coverslipThicknessUm');
    opts.designCoverslipThicknessUm = positiveScalar( ...
        opts.designCoverslipThicknessUm, 'designCoverslipThicknessUm');
    opts.fovXY = positiveScalar(opts.fovXY, 'fovXY');
    opts.nx = positiveInteger(opts.nx, 'nx');
    opts.zRange = positiveScalar(opts.zRange, 'zRange');
    opts.nz = positiveInteger(opts.nz, 'nz');
    opts.Nr = positiveInteger(opts.Nr, 'Nr');
    opts.Nphi = positiveInteger(opts.Nphi, 'Nphi');
    opts.M = positiveInteger(opts.M, 'M');
    opts.detFillRatio = positiveScalar(opts.detFillRatio, 'detFillRatio');
    opts.detectorSubsamples = positiveInteger(opts.detectorSubsamples, ...
        'detectorSubsamples');
    opts.objectiveNA = positiveScalar(opts.objectiveNA, 'objectiveNA');
    opts.mediumRefractiveIndex = positiveScalar(opts.mediumRefractiveIndex, ...
        'mediumRefractiveIndex');
    opts.excitationWavelengthUm = positiveScalar(opts.excitationWavelengthUm, ...
        'excitationWavelengthUm');
    opts.emissionWavelengthUm = positiveScalar(opts.emissionWavelengthUm, ...
        'emissionWavelengthUm');
    opts.referenceWavelengthUm = positiveScalar(opts.referenceWavelengthUm, ...
        'referenceWavelengthUm');
    opts.maxIter = positiveInteger(opts.maxIter, 'maxIter');
    opts.fdCoeff = positiveScalar(opts.fdCoeff, 'fdCoeff');
    opts.regCoeff = nonnegativeScalar(opts.regCoeff, 'regCoeff');
    opts.initialDamping = positiveScalar(opts.initialDamping, 'initialDamping');
    opts.maxCoeffStep = positiveScalar(opts.maxCoeffStep, 'maxCoeffStep');
    opts.stepTolerance = positiveScalar(opts.stepTolerance, 'stepTolerance');
    opts.residualTolerance = nonnegativeScalar(opts.residualTolerance, ...
        'residualTolerance');
    opts.svdRelativeTolerance = positiveScalar(opts.svdRelativeTolerance, ...
        'svdRelativeTolerance');
    opts.maxUsefulCondition = positiveScalar(opts.maxUsefulCondition, ...
        'maxUsefulCondition');
end

function sim = configureSimulation(opts)
    sim = defaultParams();
    sim.sampleGeometry = opts.sampleGeometry;
    sim.diffractionModel = 'vectorial Richards-Wolf/Bessel';
    sim.includesVectorialPolarization = true;
    sim.vectorialPolarizationMode = 'xyAverage';
    sim.excitationPolarizationMode = char(opts.excitationPolarizationMode);
    sim.collectionPolarizationMode = char(opts.collectionPolarizationMode);

    sim.lamExc = opts.excitationWavelengthUm;
    sim.lamEm = opts.emissionWavelengthUm;
    sim.lamRef = opts.referenceWavelengthUm;
    sim.NA = opts.objectiveNA;
    sim.nMedium = opts.mediumRefractiveIndex;
    if usesInterfaceGeometry(sim)
        sim.interfaceEmitterHeightUm = opts.interfaceEmitterHeightUm;
        sim.interfaceRadialWeightMode = opts.interfaceRadialWeightMode;
        sim.airInterfaceStageMedium = opts.airInterfaceStageMedium;
        sim.nImmersion = opts.immersionRefractiveIndex;
        sim.nGlass = opts.glassRefractiveIndex;
        sim.nSample = opts.sampleRefractiveIndex;
        sim.nDesignGlass = opts.designGlassRefractiveIndex;
        sim.coverslipThicknessUm = opts.coverslipThicknessUm;
        sim.designCoverslipThicknessUm = opts.designCoverslipThicknessUm;
        sim.nMedium = sim.nImmersion;
    end

    sim.fovXY = opts.fovXY;
    sim.nx = opts.nx;
    sim.ny = opts.nx;
    sim.x = linspace(-sim.fovXY/2, sim.fovXY/2, sim.nx);
    sim.y = sim.x;
    sim.dx = abs(sim.x(2) - sim.x(1));

    zMin = min([opts.planeZ(:); -opts.zRange/2]);
    zMax = max([opts.planeZ(:); opts.zRange/2]);
    margin = max(2*sim.dx, 0.05);
    sim.z = linspace(zMin - margin, zMax + margin, opts.nz);
    sim.nz = numel(sim.z);
    sim.nzRange = max(sim.z) - min(sim.z);

    sim.Nr = opts.Nr;
    sim.Nphi = opts.Nphi;
    sim.M = max(opts.M, opts.maxZernikeOrder);
    sim.modeOrder = zernikeModeNames(max(6, opts.maxZernikeOrder));

    if ~isempty(opts.detPitchUm)
        sim.detPitch = positiveScalar(opts.detPitchUm, 'detPitchUm');
    end
    sim.detSize = sim.detPitch * opts.detFillRatio;
    sim.detectorLayout = char(opts.detectorLayout);
    sim.detectorPixelShape = 'hex';
    sim.detectorHexRadius = sim.detSize / sqrt(3);
    [sim.detXY, sim.detectorIndexGrid, sim.detectorLayoutInfo] = ...
        detectorLayout(sim.detectorLayout, sim.detPitch);
    sim.nDet = size(sim.detXY, 1);
    sim.detectorGridSize = size(sim.detectorIndexGrid);
    sim.detectorSubsamples = opts.detectorSubsamples;
    sim.detectorImageInverted = true;
    validateScalarOptics(sim);
end

function geometry = normalizeSampleGeometry(value)
    key = lower(regexprep(char(value), '[\s_\-]', ''));
    switch key
        case 'homogeneous'
            geometry = 'homogeneous';
        case {'aironglass','interface','layeredinterface', ...
                'glasswaterinterface','wateronglass','glasswater', ...
                'glasssampleinterface','sampleonglass'}
            geometry = 'airOnGlass';
        otherwise
            error('simulateRetrieveCenterMicroimageAberrations:BadSampleGeometry', ...
                ['sampleGeometry must be homogeneous or an interface alias ' ...
                'such as glassWaterInterface.']);
    end
end

function fitModes = resolveFitModes(opts)
    if isempty(opts.fitModes)
        if opts.maxZernikeOrder >= 4
            fitModes = primarySphericalModeNames();
        else
            fitModes = zernikeModeNames(opts.maxZernikeOrder);
        end
    else
        fitModes = normalizeModeList(opts.fitModes);
    end
    if ~opts.includeTilt
        fitModes = fitModes(~ismember(fitModes, {'tilt_x','tilt_y'}));
    end
    fitModes = fitModes(:).';
end

function modes = primarySphericalModeNames()
    modes = {'tilt_x','tilt_y','defocus','astig_x','astig_y', ...
        'coma_x','coma_y','trefoil_x','trefoil_y','spherical'};
end

function modes = normalizeModeList(value)
    if ischar(value)
        modes = {value};
    elseif isstring(value)
        modes = cellstr(value(:)).';
    elseif iscell(value)
        modes = value(:).';
    else
        error('simulateRetrieveCenterMicroimageAberrations:BadFitModes', ...
            'fitModes must be a cell array, string array, char, or empty.');
    end
    modes = cellfun(@char, modes, 'UniformOutput', false);
end

function v = coefficientVectorFromInput(sim, fitModes, inputValue, defaultValue)
    if isempty(inputValue)
        inputValue = defaultValue;
    end

    v = zeros(numel(fitModes), 1);
    if isstruct(inputValue)
        for k = 1:numel(fitModes)
            if isfield(inputValue, fitModes{k})
                v(k) = double(inputValue.(fitModes{k}));
            end
        end
        return;
    end

    inputValue = double(inputValue(:));
    if isempty(inputValue)
        return;
    end
    if numel(inputValue) == numel(fitModes)
        v = inputValue(:);
        return;
    end
    if numel(inputValue) <= numel(sim.modeOrder)
        coeffs = coeffStruct(sim, inputValue);
        for k = 1:numel(fitModes)
            if isfield(coeffs, fitModes{k})
                v(k) = coeffs.(fitModes{k});
            end
        end
        return;
    end

    error('simulateRetrieveCenterMicroimageAberrations:BadCoefficientInput', ...
        'Coefficient vector must match fitModes or sim.modeOrder length.');
end

function coeffs = defaultTrueCoeffs()
    coeffs = struct();
    coeffs.tilt_x = 0.010;
    coeffs.tilt_y = -0.008;
    coeffs.defocus = 0.035;
    coeffs.astig_x = -0.030;
    coeffs.astig_y = 0.020;
    coeffs.coma_x = 0.030;
    coeffs.coma_y = -0.022;
    coeffs.trefoil_x = 0.018;
    coeffs.trefoil_y = -0.015;
    coeffs.spherical = -0.028;
end

function coeffs = coeffsFromFitVector(sim, fitModes, p)
    coeffs = struct();
    for k = 1:numel(fitModes)
        idx = find(strcmp(sim.modeOrder, fitModes{k}), 1, 'first');
        if isempty(idx)
            error('simulateRetrieveCenterMicroimageAberrations:UnknownMode', ...
                'Unknown Zernike mode "%s".', fitModes{k});
        end
        if abs(p(k)) > 1e-15
            coeffs.(fitModes{k}) = p(k);
        end
    end
end

function model = centerMicroimageModel(sim, coeffs, planeZ, sampleXY, includeExcitation)
    if nargin < 5 || isempty(includeExcitation)
        includeExcitation = true;
    end
    coeffs = coeffStruct(sim, coeffs);

    if usesInterfaceGeometry(sim)
        emitterHeight = sim.interfaceEmitterHeightUm;
        simDet = setVectorialPolarizationForPath(sim, 'collection');
        simDet.airInterfaceDirection = 'collection';
        hDet = psfBesselAirInterface( ...
            simDet, coeffs, sim.lamEm, planeZ, emitterHeight);
    else
        simDet = setVectorialPolarizationForPath(sim, 'collection');
        hDet = psfBessel(simDet, coeffs, sim.lamEm);
    end
    etaValues = detectorPointCollectionValues(sim, hDet, planeZ, sampleXY);

    nPlane = numel(planeZ);
    values = zeros(sim.nDet, nPlane);
    excitationAtSource = ones(1, nPlane);
    if includeExcitation
        if usesInterfaceGeometry(sim)
            simExc = setVectorialPolarizationForPath(sim, 'excitation');
            simExc.airInterfaceDirection = 'excitation';
            hExc = psfBesselAirInterface( ...
                simExc, coeffs, sim.lamExc, planeZ, emitterHeight);
            for ip = 1:nPlane
                excitationAtSource(ip) = sampleConditionImageAtXY( ...
                    hExc, sim, sampleXY, ip);
            end
        else
            simExc = setVectorialPolarizationForPath(sim, 'excitation');
            hExc = psfBessel(simExc, coeffs, sim.lamExc);
            for ip = 1:nPlane
                excitationAtSource(ip) = sampleVolumeAtXYZ( ...
                    hExc, sim, sampleXY, planeZ(ip));
            end
        end
    end
    for ip = 1:nPlane
        values(:, ip) = excitationAtSource(ip) * etaValues(:, ip);
    end

    model = struct();
    model.values = max(values, 0);
    model.detectorCollectionValues = max(etaValues, 0);
    model.excitationAtSource = excitationAtSource;
    model.planeZ = planeZ(:).';
    model.sampleXY = sampleXY(:).';
end

function etaValues = detectorPointCollectionValues(sim, hDet, planeZ, sampleXY)
    [DU, DV] = detectorPixelOffsets(sim, sim.detectorSubsamples);
    nQ = numel(DU);
    nPlane = numel(planeZ);
    etaValues = zeros(sim.nDet, nPlane);

    if isfield(sim, 'detectorImageInverted') && ~isempty(sim.detectorImageInverted)
        detectorImageInverted = logical(sim.detectorImageInverted);
    else
        detectorImageInverted = true;
    end

    for ip = 1:nPlane
        for k = 1:sim.nDet
            dkx = sim.detXY(k, 1);
            dky = sim.detXY(k, 2);
            acc = 0;
            for q = 1:nQ
                if detectorImageInverted
                    xy = [sampleXY(1) + dkx + DU(q), ...
                        sampleXY(2) + dky + DV(q)];
                else
                    xy = [dkx + DU(q) - sampleXY(1), ...
                        dky + DV(q) - sampleXY(2)];
                end
                if usesInterfaceGeometry(sim)
                    acc = acc + sampleConditionImageAtXY(hDet, sim, xy, ip);
                else
                    acc = acc + sampleVolumeAtXYZ(hDet, sim, xy, planeZ(ip));
                end
            end
            etaValues(k, ip) = acc / nQ;
        end
    end
end

function value = sampleConditionImageAtXY(stack, sim, sampleXY, conditionIndex)
    x = min(max(sampleXY(1), min(sim.x)), max(sim.x));
    y = min(max(sampleXY(2), min(sim.y)), max(sim.y));
    conditionIndex = min(max(round(conditionIndex), 1), size(stack, 3));
    value = interp2(sim.x, sim.y, double(stack(:,:,conditionIndex)), ...
        x, y, 'linear', 0);
end

function value = sampleVolumeAtXYZ(vol, sim, sampleXY, z)
    x = min(max(sampleXY(1), min(sim.x)), max(sim.x));
    y = min(max(sampleXY(2), min(sim.y)), max(sim.y));
    z = min(max(z, min(sim.z)), max(sim.z));

    if numel(sim.z) == 1
        value = interp2(sim.x, sim.y, double(vol(:,:,1)), x, y, 'linear', 0);
        return;
    end

    iz2 = find(sim.z >= z, 1, 'first');
    if isempty(iz2)
        iz2 = numel(sim.z);
    end
    iz1 = max(1, iz2 - 1);
    if iz1 == iz2
        t = 0;
    else
        t = (z - sim.z(iz1)) / max(sim.z(iz2) - sim.z(iz1), eps);
    end

    v1 = interp2(sim.x, sim.y, double(vol(:,:,iz1)), x, y, 'linear', 0);
    v2 = interp2(sim.x, sim.y, double(vol(:,:,iz2)), x, y, 'linear', 0);
    value = (1 - t) * v1 + t * v2;
end

function [dataValues, info] = makeSyntheticData(trueValues, opts)
    if ~strcmpi(opts.fitObjective, 'poisson') && ...
            (~isfinite(opts.photonsPerPlane) || opts.photonsPerPlane <= 0)
        dataValues = trueValues;
        info = struct('photonsPerPlane', opts.photonsPerPlane, ...
            'rngSeed', opts.rngSeed, ...
            'poissonSampled', false, ...
            'expectedCounts', [], ...
            'signalCounts', [], ...
            'darkCounts', []);
        return;
    end

    [mu, signalCounts, darkCounts] = expectedCountsFromMicroimage( ...
        trueValues, opts);
    info = struct('photonsPerPlane', opts.photonsPerPlane, ...
        'rngSeed', opts.rngSeed, ...
        'poissonSampled', false, ...
        'expectedCounts', mu, ...
        'signalCounts', signalCounts, ...
        'darkCounts', darkCounts);

    if strcmpi(opts.fitObjective, 'poisson') || ...
            (isfinite(opts.photonsPerPlane) && opts.photonsPerPlane > 0)
        if opts.poissonSampleData
            rng(opts.rngSeed);
            dataValues = poissonSample(mu);
            info.poissonSampled = true;
        else
            dataValues = mu;
        end
        return;
    end

    dataValues = trueValues;
end

function [mu, signalCounts, darkCounts] = expectedCountsFromMicroimage( ...
        values, opts)
    probs = normalizeMicroimage(values, 'perplane');
    signalCounts = opts.photonsPerPlane * probs;
    darkCounts = darkCountsArray(opts, size(values, 1), size(values, 2));
    mu = max(signalCounts + darkCounts, opts.minExpectedCount);
end

function darkCounts = darkCountsArray(opts, nDet, nPlane)
    dark = double(opts.darkCountsPerDetector(:));
    if isscalar(dark)
        dark = repmat(dark, nDet, 1);
    elseif numel(dark) ~= nDet
        error('simulateRetrieveCenterMicroimageAberrations:BadDarkCounts', ...
            'darkCountsPerDetector must be scalar or contain %d values.', nDet);
    end
    darkCounts = repmat(dark(:), 1, nPlane);
end

function y = normalizeMicroimage(values, mode)
    values = max(double(values), 0);
    switch lower(char(mode))
        case {'none','raw'}
            y = values;
        case {'global','sum'}
            y = values / max(sum(values(:)), eps);
        case {'perplane','plane','eachplane'}
            y = zeros(size(values));
            for ip = 1:size(values, 2)
                y(:,ip) = values(:,ip) / max(sum(values(:,ip)), eps);
            end
        otherwise
            error('simulateRetrieveCenterMicroimageAberrations:BadNormalization', ...
                'Unknown normalization "%s".', mode);
    end
end

function fit = fitCenterMicroimage(sim, fitModes, dataN, initialVector, opts)
    p = initialVector(:);
    damping = opts.initialDamping;
    history = zeros(opts.maxIter, 5);
    converged = false;

    for it = 1:opts.maxIter
        modelN = normalizedModelFromVector(sim, fitModes, p, opts);
        residual = dataN(:) - modelN(:);
        residualNorm = norm(residual);
        J = finiteDifferenceJacobian(sim, fitModes, p, opts);

        H0 = J.' * J;
        g = J.' * residual;
        diagH = max(diag(H0), 1);
        step = (H0 + damping * diag(diagH) + opts.regCoeff * eye(numel(p))) \ g;
        step = limitStep(step, opts.maxCoeffStep);

        accepted = false;
        trialResidualNorm = residualNorm;
        trialStep = step;
        for trial = 1:8
            pTrial = p + trialStep;
            modelTrial = normalizedModelFromVector(sim, fitModes, pTrial, opts);
            trialResidualNorm = norm(dataN(:) - modelTrial(:));
            if trialResidualNorm < residualNorm
                p = pTrial;
                damping = max(damping / 3, 1e-9);
                accepted = true;
                break;
            end
            trialStep = 0.5 * trialStep;
            damping = min(damping * 3, 1e9);
        end

        history(it,:) = [it, residualNorm, trialResidualNorm, ...
            norm(trialStep), damping];

        if opts.verbose
            fprintf(['[simulateRetrieveCenterMicroimageAberrations] iter %d: ' ...
                'residual %.4g -> %.4g, step %.3g, accepted %d\n'], ...
                it, residualNorm, trialResidualNorm, norm(trialStep), accepted);
        end

        if norm(trialStep) < opts.stepTolerance || ...
                trialResidualNorm <= opts.residualTolerance
            converged = true;
            history = history(1:it,:);
            break;
        end
        if ~accepted && norm(trialStep) < opts.stepTolerance
            history = history(1:it,:);
            break;
        end
    end

    modelN = normalizedModelFromVector(sim, fitModes, p, opts);
    fit = struct();
    fit.coeffVector = p;
    fit.coeffs = coeffsFromFitVector(sim, fitModes, p);
    fit.normalizedModel = modelN;
    fit.residual = dataN - modelN;
    fit.residualNorm = norm(fit.residual(:));
    fit.history = history(any(history,2),:);
    fit.converged = converged;
end

function fit = fitCenterMicroimagePoisson( ...
        sim, fitModes, dataCounts, initialVector, opts)
    dataCounts = max(double(dataCounts), 0);
    p = initialVector(:);
    damping = opts.initialDamping;
    history = zeros(opts.maxIter, 6);
    converged = false;

    for it = 1:opts.maxIter
        [mu, modelRaw] = expectedCountsFromVector(sim, fitModes, p, opts);
        dev = poissonDeviance(dataCounts, mu);
        J = finiteDifferenceCountsJacobian(sim, fitModes, p, opts);
        safeMu = max(mu(:), opts.minExpectedCount);
        score = J.' * ((dataCounts(:) - mu(:)) ./ safeMu);
        fisher = J.' * bsxfun(@rdivide, J, safeMu);
        diagF = max(diag(fisher), 1);
        step = (fisher + damping * diag(diagF) + ...
            opts.regCoeff * eye(numel(p))) \ score;
        step = limitStep(step, opts.maxCoeffStep);

        accepted = false;
        trialDev = dev;
        trialStep = step;
        for trial = 1:8
            pTrial = p + trialStep;
            muTrial = expectedCountsFromVector(sim, fitModes, pTrial, opts);
            trialDev = poissonDeviance(dataCounts, muTrial);
            if trialDev < dev
                p = pTrial;
                damping = max(damping / 3, 1e-9);
                accepted = true;
                break;
            end
            trialStep = 0.5 * trialStep;
            damping = min(damping * 3, 1e9);
        end

        history(it,:) = [it, dev, trialDev, norm(trialStep), ...
            damping, double(accepted)];

        if opts.verbose
            fprintf(['[simulateRetrieveCenterMicroimageAberrations] poisson ' ...
                'iter %d: deviance %.4g -> %.4g, step %.3g, accepted %d\n'], ...
                it, dev, trialDev, norm(trialStep), accepted);
        end

        if norm(trialStep) < opts.stepTolerance || ...
                trialDev <= opts.residualTolerance
            converged = true;
            history = history(1:it,:);
            break;
        end
        if ~accepted && norm(trialStep) < opts.stepTolerance
            history = history(1:it,:);
            break;
        end
    end

    [mu, modelRaw] = expectedCountsFromVector(sim, fitModes, p, opts);
    pearson = (dataCounts - mu) ./ sqrt(max(mu, opts.minExpectedCount));
    fit = struct();
    fit.coeffVector = p;
    fit.coeffs = coeffsFromFitVector(sim, fitModes, p);
    fit.expectedCounts = mu;
    fit.rawModel = modelRaw;
    fit.normalizedModel = normalizeMicroimage(mu, opts.normalization);
    fit.residual = dataCounts - mu;
    fit.pearsonResidual = pearson;
    fit.residualNorm = norm(pearson(:));
    fit.deviance = poissonDeviance(dataCounts, mu);
    fit.reducedDeviance = fit.deviance / ...
        max(numel(dataCounts) - numel(p), 1);
    fit.history = history(any(history,2),:);
    fit.converged = converged;
    fit.objective = 'poisson';
end

function mc = runMonteCarloCoefficientRecovery( ...
        sim, fitModes, trueValues, trueCoeffVector, initialVector, opts)
    nRep = opts.monteCarloRepeats;
    nMode = numel(fitModes);
    if nRep <= 0
        mc = emptyMonteCarlo(fitModes);
        return;
    end

    coeffVectors = nan(nMode, nRep);
    deviance = nan(nRep, 1);
    reducedDeviance = nan(nRep, 1);
    converged = false(nRep, 1);

    mcOpts = opts;
    mcOpts.poissonSampleData = true;
    mcOpts.verbose = false;
    for rep = 1:nRep
        mcOpts.rngSeed = opts.rngSeed + rep;
        [dataCounts, ~] = makeSyntheticData(trueValues, mcOpts);
        if strcmpi(opts.fitObjective, 'poisson')
            fit = fitCenterMicroimagePoisson( ...
                sim, fitModes, dataCounts, initialVector, mcOpts);
        else
            dataN = normalizeMicroimage(dataCounts, opts.normalization);
            fit = fitCenterMicroimage( ...
                sim, fitModes, dataN, initialVector, mcOpts);
        end
        coeffVectors(:, rep) = fit.coeffVector(:);
        converged(rep) = fit.converged;
        if isfield(fit, 'deviance')
            deviance(rep) = fit.deviance;
            reducedDeviance(rep) = fit.reducedDeviance;
        end
        if opts.verbose && (rep == 1 || rep == nRep || mod(rep, 10) == 0)
            fprintf('[simulateRetrieveCenterMicroimageAberrations] MC %d/%d\n', ...
                rep, nRep);
        end
    end

    errors = coeffVectors - trueCoeffVector(:);
    mc = struct();
    mc.nRepeats = nRep;
    mc.mode = fitModes(:);
    mc.displayMode = cellfun(@displayModeLabel, fitModes(:), ...
        'UniformOutput', false);
    mc.coeffVectors = coeffVectors;
    mc.errorVectors = errors;
    mc.meanWaves = mean(coeffVectors, 2, 'omitnan');
    mc.stdWaves = std(coeffVectors, 0, 2, 'omitnan');
    mc.biasWaves = mean(errors, 2, 'omitnan');
    mc.rmseWaves = sqrt(mean(errors.^2, 2, 'omitnan'));
    mc.deviance = deviance;
    mc.reducedDeviance = reducedDeviance;
    mc.converged = converged;
    mc.summaryTable = monteCarloSummaryTable(mc, trueCoeffVector);
    mc.samplesTable = monteCarloSamplesTable(mc, trueCoeffVector);
end

function mc = emptyMonteCarlo(fitModes)
    mc = struct();
    mc.nRepeats = 0;
    mc.mode = fitModes(:);
    mc.displayMode = cellfun(@displayModeLabel, fitModes(:), ...
        'UniformOutput', false);
    mc.coeffVectors = [];
    mc.errorVectors = [];
    mc.meanWaves = nan(numel(fitModes), 1);
    mc.stdWaves = nan(numel(fitModes), 1);
    mc.biasWaves = nan(numel(fitModes), 1);
    mc.rmseWaves = nan(numel(fitModes), 1);
    mc.deviance = [];
    mc.reducedDeviance = [];
    mc.converged = [];
    mc.summaryTable = table();
    mc.samplesTable = table();
end

function T = monteCarloSummaryTable(mc, trueCoeffVector)
    mode = mc.mode(:);
    displayMode = mc.displayMode(:);
    trueWaves = trueCoeffVector(:);
    mcMeanWaves = mc.meanWaves(:);
    mcStdWaves = mc.stdWaves(:);
    mcBiasWaves = mc.biasWaves(:);
    mcRmseWaves = mc.rmseWaves(:);
    T = table(mode, displayMode, trueWaves, mcMeanWaves, mcStdWaves, ...
        mcBiasWaves, mcRmseWaves);
end

function T = monteCarloSamplesTable(mc, trueCoeffVector)
    nMode = numel(mc.mode);
    nRep = mc.nRepeats;
    repeatIndex = reshape(repmat(1:nRep, nMode, 1), [], 1);
    mode = repmat(mc.mode(:), nRep, 1);
    displayMode = repmat(mc.displayMode(:), nRep, 1);
    trueWaves = repmat(trueCoeffVector(:), nRep, 1);
    estimatedWaves = mc.coeffVectors(:);
    errorWaves = estimatedWaves - trueWaves;
    T = table(repeatIndex, mode, displayMode, trueWaves, ...
        estimatedWaves, errorWaves);
end

function [mu, model] = expectedCountsFromVector(sim, fitModes, p, opts)
    coeffs = coeffsFromFitVector(sim, fitModes, p);
    model = centerMicroimageModel(sim, coeffs, opts.planeZ, opts.sampleXY, true);
    mu = expectedCountsFromMicroimage(model.values, opts);
end

function modelN = normalizedModelFromVector(sim, fitModes, p, opts)
    coeffs = coeffsFromFitVector(sim, fitModes, p);
    includeExcitation = ~strcmpi(opts.normalization, 'perplane');
    model = centerMicroimageModel(sim, coeffs, opts.planeZ, opts.sampleXY, ...
        includeExcitation);
    modelN = normalizeMicroimage(model.values, opts.normalization);
end

function J = finiteDifferenceJacobian(sim, fitModes, p, opts)
    base = normalizedModelFromVector(sim, fitModes, p, opts);
    J = zeros(numel(base), numel(p));
    for k = 1:numel(p)
        step = opts.fdCoeff;
        pp = p;
        pm = p;
        pp(k) = pp(k) + step;
        pm(k) = pm(k) - step;
        yp = normalizedModelFromVector(sim, fitModes, pp, opts);
        ym = normalizedModelFromVector(sim, fitModes, pm, opts);
        J(:,k) = (yp(:) - ym(:)) / (2 * step);
    end
end

function J = finiteDifferenceCountsJacobian(sim, fitModes, p, opts)
    base = expectedCountsFromVector(sim, fitModes, p, opts);
    J = zeros(numel(base), numel(p));
    for k = 1:numel(p)
        step = opts.fdCoeff;
        pp = p;
        pm = p;
        pp(k) = pp(k) + step;
        pm(k) = pm(k) - step;
        yp = expectedCountsFromVector(sim, fitModes, pp, opts);
        ym = expectedCountsFromVector(sim, fitModes, pm, opts);
        J(:,k) = (yp(:) - ym(:)) / (2 * step);
    end
end

function step = limitStep(step, maxAbsStep)
    scale = max(abs(step)) / maxAbsStep;
    if scale > 1
        step = step / scale;
    end
end

function D = poissonDeviance(dataCounts, mu)
    n = max(double(dataCounts(:)), 0);
    mu = max(double(mu(:)), eps);
    terms = mu - n;
    positive = n > 0;
    terms(positive) = terms(positive) + ...
        n(positive) .* log(n(positive) ./ mu(positive));
    D = 2 * sum(terms);
end

function ident = centerMicroimageIdentifiability(sim, fitModes, p, opts)
    J = finiteDifferenceJacobian(sim, fitModes, p(:), opts);
    s = svd(J, 'econ');
    if isempty(s)
        relative = [];
        rankValue = 0;
        conditionNumber = Inf;
    else
        relative = s / max(s(1), eps);
        rankValue = sum(relative > opts.svdRelativeTolerance);
        if rankValue >= numel(fitModes)
            conditionNumber = s(1) / max(s(numel(fitModes)), eps);
        else
            conditionNumber = Inf;
        end
    end
    ident = struct();
    ident.jacobian = J;
    ident.singularValues = s;
    ident.relativeSingularValues = relative;
    ident.rank = rankValue;
    ident.nModes = numel(fitModes);
    ident.nMeasurements = size(J, 1);
    ident.svdRelativeTolerance = opts.svdRelativeTolerance;
    ident.conditionNumber = conditionNumber;
    ident.isFullRank = rankValue >= numel(fitModes);
    ident.isWellConditioned = ident.isFullRank && ...
        conditionNumber <= opts.maxUsefulCondition;
end

function fisher = poissonFisherDiagnostics(sim, fitModes, p, opts)
    [mu, model] = expectedCountsFromVector(sim, fitModes, p(:), opts);
    J = finiteDifferenceCountsJacobian(sim, fitModes, p(:), opts);
    safeMu = max(mu(:), opts.minExpectedCount);
    weightedJConditional = bsxfun(@rdivide, J, sqrt(safeMu));

    nuisanceRank = 0;
    weightedJ = weightedJConditional;
    if opts.marginalizeFlux
        probs = normalizeMicroimage(model.values, 'perplane');
        fluxSensitivity = perPlaneFluxSensitivity(probs);
        weightedFlux = bsxfun(@rdivide, fluxSensitivity, sqrt(safeMu));
        Q = orthonormalColumns(weightedFlux, opts.svdRelativeTolerance);
        nuisanceRank = size(Q, 2);
        if nuisanceRank > 0
            weightedJ = weightedJConditional - Q * (Q.' * weightedJConditional);
        end
    end

    conditional = fisherStatsFromWeightedJ( ...
        weightedJConditional, numel(fitModes), opts.svdRelativeTolerance);
    marginal = fisherStatsFromWeightedJ( ...
        weightedJ, numel(fitModes), opts.svdRelativeTolerance);

    fisher = struct();
    fisher.information = marginal.information;
    fisher.conditionalInformation = conditional.information;
    fisher.covariance = marginal.covariance;
    fisher.crbStdWaves = marginal.crbStd(:);
    fisher.eigenValues = marginal.eigenValues(:);
    fisher.relativeEigenValues = marginal.relativeEigenValues(:);
    fisher.eigenVectors = marginal.eigenVectors;
    fisher.conditionalCrbStdWaves = conditional.crbStd(:);
    fisher.conditionalEigenValues = conditional.eigenValues(:);
    fisher.rank = marginal.rank;
    fisher.conditionalRank = conditional.rank;
    fisher.nModes = numel(fitModes);
    fisher.nMeasurements = numel(mu);
    fisher.conditionNumber = marginal.conditionNumber;
    fisher.conditionalConditionNumber = conditional.conditionNumber;
    fisher.isFullRank = marginal.isFullRank;
    fisher.conditionalIsFullRank = conditional.isFullRank;
    fisher.marginalizedFlux = logical(opts.marginalizeFlux);
    fisher.nuisanceRank = nuisanceRank;
    fisher.weightedSensitivity = weightedJ;
    fisher.conditionalWeightedSensitivity = weightedJConditional;
    fisher.expectedCounts = mu;
end

function sensitivity = perPlaneFluxSensitivity(probs)
    nDet = size(probs, 1);
    nPlane = size(probs, 2);
    sensitivity = zeros(numel(probs), nPlane);
    for ip = 1:nPlane
        rows = (1:nDet) + (ip - 1) * nDet;
        sensitivity(rows, ip) = probs(:, ip);
    end
end

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

function stats = fisherStatsFromWeightedJ(weightedJ, nModes, relTol)
    I = weightedJ.' * weightedJ;
    I = 0.5 * (I + I.');

    if isempty(I)
        eigVals = [];
        eigVecs = [];
        relative = [];
        rankValue = 0;
        conditionNumber = Inf;
        covar = [];
        crbStd = [];
    else
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
                if ~isempty(nullBasis) && norm(nullBasis(k,:)) > 1e-6
                    crbStd(k) = Inf;
                end
            end
        end
    end

    stats = struct();
    stats.information = I;
    stats.covariance = covar;
    stats.crbStd = crbStd(:);
    stats.eigenValues = eigVals(:);
    stats.relativeEigenValues = relative(:);
    stats.eigenVectors = eigVecs;
    stats.rank = rankValue;
    stats.conditionNumber = conditionNumber;
    stats.isFullRank = rankValue >= nModes;
end

function fisher = emptyPoissonFisher(nModes)
    fisher = struct();
    fisher.information = [];
    fisher.conditionalInformation = [];
    fisher.covariance = [];
    fisher.crbStdWaves = nan(nModes, 1);
    fisher.conditionalCrbStdWaves = nan(nModes, 1);
    fisher.eigenValues = [];
    fisher.relativeEigenValues = [];
    fisher.eigenVectors = [];
    fisher.conditionalEigenValues = [];
    fisher.rank = 0;
    fisher.conditionalRank = 0;
    fisher.nModes = nModes;
    fisher.nMeasurements = 0;
    fisher.conditionNumber = NaN;
    fisher.conditionalConditionNumber = NaN;
    fisher.isFullRank = false;
    fisher.conditionalIsFullRank = false;
    fisher.marginalizedFlux = false;
    fisher.nuisanceRank = 0;
    fisher.weightedSensitivity = [];
    fisher.conditionalWeightedSensitivity = [];
    fisher.expectedCounts = [];
end

function T = coefficientRecoveryTable( ...
        fitModes, trueVector, estimatedVector, crbStdWaves, monteCarlo)
    if nargin < 4 || isempty(crbStdWaves)
        crbStdWaves = nan(numel(fitModes), 1);
    end
    if nargin < 5 || isempty(monteCarlo)
        monteCarlo = emptyMonteCarlo(fitModes);
    end
    mode = fitModes(:);
    displayMode = cellfun(@displayModeLabel, mode, 'UniformOutput', false);
    displayMode = displayMode(:);
    zernikeOrder = zeros(numel(mode), 1);
    for k = 1:numel(mode)
        zernikeOrder(k) = modeOrderNumber(mode{k});
    end
    trueWaves = trueVector(:);
    estimatedWaves = estimatedVector(:);
    errorWaves = estimatedWaves - trueWaves;
    absErrorWaves = abs(errorWaves);
    crbStdWaves = crbStdWaves(:);
    mcMeanWaves = monteCarlo.meanWaves(:);
    mcStdWaves = monteCarlo.stdWaves(:);
    mcBiasWaves = monteCarlo.biasWaves(:);
    mcRmseWaves = monteCarlo.rmseWaves(:);

    fitted = true(numel(mode), 1);
    T = table(mode, displayMode, zernikeOrder, trueWaves, estimatedWaves, ...
        errorWaves, absErrorWaves, crbStdWaves, mcMeanWaves, ...
        mcStdWaves, mcBiasWaves, mcRmseWaves, fitted);
end

function label = displayModeLabel(mode)
    switch char(mode)
        case 'astig_x'
            label = 'astig 0/90';
        case 'astig_y'
            label = 'astig 45/135';
        case 'secondary_astig_x'
            label = 'secondary astig 0/90';
        case 'secondary_astig_y'
            label = 'secondary astig 45/135';
        case 'tertiary_astig_x'
            label = 'tertiary astig 0/90';
        case 'tertiary_astig_y'
            label = 'tertiary astig 45/135';
        otherwise
            label = char(mode);
    end
end

function order = modeOrderNumber(mode)
    all = { ...
        'tilt_x', 1; 'tilt_y', 1; ...
        'defocus', 2; 'astig_x', 2; 'astig_y', 2; ...
        'coma_x', 3; 'coma_y', 3; 'trefoil_x', 3; 'trefoil_y', 3; ...
        'spherical', 4; 'secondary_astig_x', 4; 'secondary_astig_y', 4; ...
        'quadrafoil_x', 4; 'quadrafoil_y', 4; ...
        'secondary_coma_x', 5; 'secondary_coma_y', 5; ...
        'secondary_trefoil_x', 5; 'secondary_trefoil_y', 5; ...
        'pentafoil_x', 5; 'pentafoil_y', 5; ...
        'secondary_spherical', 6; 'tertiary_astig_x', 6; ...
        'tertiary_astig_y', 6; 'secondary_quadrafoil_x', 6; ...
        'secondary_quadrafoil_y', 6; 'hexafoil_x', 6; 'hexafoil_y', 6};
    idx = find(strcmp(all(:,1), mode), 1, 'first');
    if isempty(idx)
        order = NaN;
    else
        order = all{idx,2};
    end
end

function text = recommendationFromIdentifiability(ident, opts)
    if ident.isWellConditioned
        text = sprintf(['Center microimage is locally full rank for these %d ' ...
            'modes at tolerance %.1g. Coefficient estimation is numerically ' ...
            'plausible for this synthetic model.'], ...
            ident.nModes, opts.svdRelativeTolerance);
    elseif ident.isFullRank
        text = sprintf(['Center microimage is formally full rank for these %d ' ...
            'modes, but condition number %.3g exceeds %.3g. Expect unstable ' ...
            'coefficients without stronger priors, more photons, or axial/scan ' ...
            'diversity.'], ...
            ident.nModes, ident.conditionNumber, opts.maxUsefulCondition);
    else
        text = sprintf(['Center microimage is rank deficient: rank %d for %d ' ...
            'fitted modes with %d measurements. Center-only data cannot uniquely ' ...
            'estimate this coefficient set; add a diversity plane or detector-' ...
            'resolved scan images.'], ...
            ident.rank, ident.nModes, ident.nMeasurements);
    end
end

function outDir = resolveOutputDir(opts)
    if ~isempty(opts.outputDir)
        outDir = char(opts.outputDir);
        return;
    end
    thisDir = fileparts(mfilename('fullpath'));
    root = fullfile(fileparts(thisDir), 'output_matlab', ...
        'center_microimage_vectorial_aberration_retrieval');
    outDir = fullfile(root, sprintf('%s_through_spherical_%dplane', ...
        sampleGeometryFolderLabel(opts), numel(opts.planeZ)));
end

function label = sampleGeometryFolderLabel(opts)
    if strcmpi(opts.sampleGeometry, 'homogeneous')
        label = 'homogeneous';
        return;
    end
    if abs(opts.sampleRefractiveIndex - 1.33) < 0.03
        label = 'glass_water_interface';
    elseif opts.sampleRefractiveIndex < 1.1
        label = 'glass_air_interface';
    else
        label = sprintf('glass_sample_n%.3f_interface', ...
            opts.sampleRefractiveIndex);
        label = strrep(label, '.', 'p');
    end
end

function writeOutputs(result)
    outDir = result.outputDir;
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    writetable(result.coefficientTable, ...
        fullfile(outDir, 'center_microimage_true_vs_estimated_coefficients.csv'));

    ident = result.identifiability;
    singularIndex = (1:numel(ident.singularValues)).';
    singularValue = ident.singularValues(:);
    relativeSingularValue = ident.relativeSingularValues(:);
    writetable(table(singularIndex, singularValue, relativeSingularValue), ...
        fullfile(outDir, 'center_microimage_jacobian_singular_values.csv'));

    crbTable = result.coefficientTable(:, ...
        {'mode','displayMode','trueWaves','estimatedWaves','crbStdWaves'});
    crbTable.conditionalCrbStdWaves = ...
        result.poissonFisher.conditionalCrbStdWaves(:);
    writetable(crbTable, fullfile(outDir, 'center_microimage_poisson_crb.csv'));

    writePoissonFisherEigenTable(result, ...
        fullfile(outDir, 'center_microimage_poisson_fisher_eigenmodes.csv'));
    writePoissonCountsTable(result, ...
        fullfile(outDir, 'center_microimage_poisson_counts.csv'));
    if result.monteCarlo.nRepeats > 0
        writetable(result.monteCarlo.summaryTable, ...
            fullfile(outDir, 'center_microimage_monte_carlo_summary.csv'));
        writetable(result.monteCarlo.samplesTable, ...
            fullfile(outDir, 'center_microimage_monte_carlo_samples.csv'));
    end

    recommendation = {result.recommendation};
    poissonDeviance = NaN;
    reducedPoissonDeviance = NaN;
    if isfield(result.fit, 'deviance')
        poissonDeviance = result.fit.deviance;
        reducedPoissonDeviance = result.fit.reducedDeviance;
    end
    fisher = result.poissonFisher;
    summary = table( ...
        ident.nMeasurements, ident.nModes, ident.rank, ident.isFullRank, ...
        ident.conditionNumber, result.fit.residualNorm, poissonDeviance, ...
        reducedPoissonDeviance, fisher.rank, fisher.isFullRank, ...
        fisher.conditionNumber, fisher.marginalizedFlux, ...
        fisher.nuisanceRank, result.monteCarlo.nRepeats, recommendation, ...
        'VariableNames', {'nMeasurements','nModes','rank','isFullRank', ...
        'conditionNumber','residualNorm','poissonDeviance', ...
        'reducedPoissonDeviance','poissonFisherRank', ...
        'poissonFisherFullRank','poissonFisherConditionNumber', ...
        'poissonFisherMarginalizedFlux','poissonFisherNuisanceRank', ...
        'monteCarloRepeats','recommendation'});
    writetable(summary, fullfile(outDir, 'center_microimage_fit_summary.csv'));

    writeMicroimageFigure(result, ...
        fullfile(outDir, 'center_microimage_true_fit_residual_honeycomb.png'));
    writeCoefficientFigure(result, ...
        fullfile(outDir, 'center_microimage_coefficient_recovery.png'));

    save(fullfile(outDir, 'center_microimage_aberration_retrieval.mat'), ...
        'result', '-v7.3');
end

function writePoissonFisherEigenTable(result, outFile)
    fisher = result.poissonFisher;
    eigenIndex = (1:numel(fisher.eigenValues)).';
    eigenValue = fisher.eigenValues(:);
    relativeEigenValue = fisher.relativeEigenValues(:);
    T = table(eigenIndex, eigenValue, relativeEigenValue);
    if ~isempty(fisher.eigenVectors)
        names = matlab.lang.makeValidName(strcat('weight_', result.fitModes));
        W = array2table(fisher.eigenVectors, 'VariableNames', names);
        T = [T W];
    end
    writetable(T, outFile);
end

function writePoissonCountsTable(result, outFile)
    dataCounts = result.syntheticData;
    truthExpectedCounts = result.expectedCountsAtTruth;
    fitExpectedCounts = result.estimatedExpectedCounts;
    signalCounts = result.syntheticPhotonInfo.signalCounts;
    darkCounts = result.syntheticPhotonInfo.darkCounts;
    nDet = size(dataCounts, 1);
    nPlane = size(dataCounts, 2);

    detectorIndex = repmat((1:nDet).', nPlane, 1);
    planeIndex = reshape(repmat(1:nPlane, nDet, 1), [], 1);
    planeZ = reshape(repmat(result.options.planeZ(:).', nDet, 1), [], 1);
    dataCounts = dataCounts(:);
    truthExpectedCounts = vectorOrNaN(truthExpectedCounts, nDet*nPlane);
    fitExpectedCounts = vectorOrNaN(fitExpectedCounts, nDet*nPlane);
    signalCounts = vectorOrNaN(signalCounts, nDet*nPlane);
    darkCounts = vectorOrNaN(darkCounts, nDet*nPlane);

    T = table(detectorIndex, planeIndex, planeZ, dataCounts, ...
        truthExpectedCounts, fitExpectedCounts, signalCounts, darkCounts);
    writetable(T, outFile);
end

function v = vectorOrNaN(value, n)
    if isempty(value)
        v = nan(n, 1);
    else
        v = double(value(:));
    end
end

function writeMicroimageFigure(result, outFile)
    dataN = result.dataNormalized;
    fitN = result.estimatedNormalizedMicroimage;
    residual = dataN - fitN;
    nPlane = size(dataN, 2);
    vmax = max([dataN(:); fitN(:)]);
    if ~isfinite(vmax) || vmax <= 0
        vmax = 1;
    end
    rmax = max(abs(residual(:)));
    if ~isfinite(rmax) || rmax <= 0
        rmax = 1;
    end

    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 1020 max(320, 260*nPlane)]);
    tl = tiledlayout(fig, nPlane, 3, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    for ip = 1:nPlane
        ax = nexttile(tl);
        plotDetectorHexMap(result.sim.detXY, dataN(:,ip), ...
            'Parent', ax, 'CLim', [0 vmax]);
        title(ax, sprintf('data z=%.3g um', result.options.planeZ(ip)));
        colormap(ax, 'parula');

        ax = nexttile(tl);
        plotDetectorHexMap(result.sim.detXY, fitN(:,ip), ...
            'Parent', ax, 'CLim', [0 vmax]);
        title(ax, 'fit');
        colormap(ax, 'parula');

        ax = nexttile(tl);
        plotDetectorHexMap(result.sim.detXY, residual(:,ip), ...
            'Parent', ax, 'CLim', [-rmax rmax]);
        title(ax, 'data - fit');
        colormap(ax, blueWhiteRed(256));
    end
    title(tl, 'Center point detector microimage: data, fitted model, residual');
    exportFigure(fig, outFile);
end

function writeCoefficientFigure(result, outFile)
    T = result.coefficientTable;
    hasMC = isfield(result, 'monteCarlo') && result.monteCarlo.nRepeats > 0 && ...
        ismember('mcMeanWaves', T.Properties.VariableNames);
    if hasMC
        barValues = [T.trueWaves, T.estimatedWaves, T.mcMeanWaves];
        barLabels = {'true','single fit','MC mean'};
        errValues = T.mcStdWaves;
        errBarIndex = 3;
        errLabel = sprintf('MC 1 sigma, n=%d', result.monteCarlo.nRepeats);
    else
        barValues = [T.trueWaves, T.estimatedWaves];
        barLabels = {'true','estimated'};
        errValues = T.crbStdWaves;
        errBarIndex = 2;
        errLabel = 'CRB 1 sigma';
    end

    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 1220 520]);
    b = bar(barValues);
    hold on;
    errPlot = errValues(:);
    errPlot(~isfinite(errPlot)) = NaN;
    xErr = barXEndPoints(b, errBarIndex, size(barValues, 2), height(T));
    e = errorbar(xErr(:), barValues(:,errBarIndex), errPlot, ...
        'k.', 'LineWidth', 1.2, 'CapSize', 8);
    grid on;
    ylabel('coefficient (waves RMS)');
    legend([b(:).' e], [barLabels {errLabel}], 'Location', 'best');
    if ismember('displayMode', T.Properties.VariableNames)
        tickLabels = T.displayMode;
    else
        tickLabels = T.mode;
    end
    set(gca, 'XTick', 1:height(T), 'XTickLabel', tickLabels, ...
        'XTickLabelRotation', 45);
    title(sprintf('Center-microimage coefficient recovery through spherical (%s)', ...
        errLabel));
    exportFigure(fig, outFile);
end

function x = barXEndPoints(b, barIndex, nBars, nGroups)
    try
        x = b(barIndex).XEndPoints;
    catch
        groupWidth = min(0.8, nBars/(nBars + 1.5));
        x = (1:nGroups) - groupWidth/2 + ...
            (2*barIndex - 1) * groupWidth / (2*nBars);
    end
end

function cmap = blueWhiteRed(n)
    if nargin < 1
        n = 256;
    end
    n1 = floor(n/2);
    n2 = n - n1;
    blue = [0.18 0.28 0.72];
    white = [1 1 1];
    red = [0.75 0.16 0.16];
    c1 = [linspace(blue(1), white(1), n1).', ...
        linspace(blue(2), white(2), n1).', ...
        linspace(blue(3), white(3), n1).'];
    c2 = [linspace(white(1), red(1), n2).', ...
        linspace(white(2), red(2), n2).', ...
        linspace(white(3), red(3), n2).'];
    cmap = [c1; c2];
end

function exportFigure(fig, outFile)
    try
        exportgraphics(fig, outFile, 'Resolution', 220);
    catch
        print(fig, outFile, '-dpng', '-r220');
    end
    close(fig);
end

function printSummary(result)
    ident = result.identifiability;
    fisher = result.poissonFisher;
    fprintf('\nCenter microimage vectorial aberration retrieval\n');
    fprintf('  model: %s\n', result.modelDescription);
    fprintf('  objective: %s\n', result.options.fitObjective);
    fprintf('  fit modes: %d through order %d\n', ...
        numel(result.fitModes), result.options.maxZernikeOrder);
    fprintf('  measurements: %d, rank: %d, condition: %.3g\n', ...
        ident.nMeasurements, ident.rank, ident.conditionNumber);
    fprintf('  residual norm: %.4g\n', result.fit.residualNorm);
    if isfield(result.fit, 'deviance')
        fprintf('  poisson deviance: %.4g (reduced %.3g)\n', ...
            result.fit.deviance, result.fit.reducedDeviance);
    end
    fprintf('  poisson Fisher rank: %d/%d, condition: %.3g\n', ...
        fisher.rank, fisher.nModes, fisher.conditionNumber);
    fprintf('  flux nuisance marginalized: %d (rank %d)\n', ...
        fisher.marginalizedFlux, fisher.nuisanceRank);
    if isfield(result, 'monteCarlo') && result.monteCarlo.nRepeats > 0
        fprintf('  Monte Carlo repeats: %d, median sigma %.3g waves\n', ...
            result.monteCarlo.nRepeats, ...
            median(result.monteCarlo.stdWaves, 'omitnan'));
    end
    fprintf('  recommendation: %s\n', result.recommendation);
    if result.options.writeOutputs
        fprintf('  outputs: %s\n', result.outputDir);
    end
end

function tf = usesInterfaceGeometry(sim)
    tf = isfield(sim, 'sampleGeometry') && ...
        strcmpi(sim.sampleGeometry, 'airOnGlass');
end

function text = modelDescriptionForSim(sim)
    if usesInterfaceGeometry(sim)
        if abs(sim.nSample - 1.33) < 0.03
            sampleLabel = 'water';
        elseif sim.nSample < 1.1
            sampleLabel = 'air';
        else
            sampleLabel = sprintf('n=%.4g sample', sim.nSample);
        end
        text = sprintf(['glass-%s interface vectorial Richards-Wolf/Bessel; ' ...
            'emitter height %.4g um above interface; xyAverage/isotropic ' ...
            'collection as rapidly rotating dipole approximation'], ...
            sampleLabel, sim.interfaceEmitterHeightUm);
    else
        text = ['homogeneous vectorial Richards-Wolf/Bessel; ' ...
            'xyAverage/isotropic polarization as rapidly rotating dipole approximation'];
    end
end

function value = positiveInteger(value, name)
    value = double(value);
    if ~isscalar(value) || ~isfinite(value) || value < 1 || round(value) ~= value
        error('simulateRetrieveCenterMicroimageAberrations:BadOption', ...
            '%s must be a positive integer.', name);
    end
end

function value = nonnegativeInteger(value, name)
    value = double(value);
    if ~isscalar(value) || ~isfinite(value) || value < 0 || round(value) ~= value
        error('simulateRetrieveCenterMicroimageAberrations:BadOption', ...
            '%s must be a nonnegative integer.', name);
    end
end

function value = positiveScalar(value, name)
    value = double(value);
    if ~isscalar(value) || ~isfinite(value) || value <= 0
        error('simulateRetrieveCenterMicroimageAberrations:BadOption', ...
            '%s must be a positive scalar.', name);
    end
end

function value = nonnegativeScalar(value, name)
    value = double(value);
    if ~isscalar(value) || ~isfinite(value) || value < 0
        error('simulateRetrieveCenterMicroimageAberrations:BadOption', ...
            '%s must be a nonnegative scalar.', name);
    end
end
