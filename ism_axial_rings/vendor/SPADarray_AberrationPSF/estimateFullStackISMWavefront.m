function result = estimateFullStackISMWavefront(stackInput, varargin)
%ESTIMATEFULLSTACKISMWAVEFRONT Hybrid full-stack Poisson wavefront retrieval.
%
%   result = estimateFullStackISMWavefront(alignmentCsv)
%
%   The likelihood uses every recorded stage position. Selected symmetric
%   planes contribute their full 23-channel XY images; all remaining planes
%   contribute their total photon counts. This partitions the raw data so the
%   axial trace and selected images are not double-counted.
%
%   The optical-axis position z0 is fitted jointly with the requested
%   non-defocus wavefront modes. Zernike defocus is excluded by default because
%   z0 should absorb axial origin, while detector geometry is fixed to a
%   regular PDA23 honeycomb unless estimateDetectorLayout is explicitly enabled.

    if nargin < 1 || isempty(stackInput)
        stackInput = ['D:\Luminosa\Data\PSF_batch_outputs\ISM_Aberation2_73\' ...
            'xz_yz_plots\x0_4uW_0_19collar_80mmlens_20260515_155744_frame_alignment.csv'];
    end
    opts = parseOptions(varargin{:});
    addRequiredPaths(opts);

    data = loadFullStackISMData(stackInput, ...
        'stageZUm', opts.stageZUm, ...
        'xyPixelSizeUm', opts.xyPixelSizeUm, ...
        'channelIDs', opts.channelIDs, ...
        'channelOrder', opts.channelOrder, ...
        'ptuReaderFolder', opts.ptuReaderFolder, ...
        'ptuPhotonsPerChunk', opts.ptuPhotonsPerChunk, ...
        'boundaryWidthPx', opts.boundaryWidthPx, ...
        'darkFile', opts.darkFile, ...
        'darkScale', opts.darkScale, ...
        'darkPerPixel', opts.darkPerPixel, ...
        'darkMeasurementMode', opts.darkMeasurementMode, ...
        'backgroundMode', opts.backgroundMode, ...
        'flatField', opts.flatField, ...
        'flatFieldVariable', opts.flatFieldVariable, ...
        'cacheFile', opts.cacheFile, ...
        'reuseCache', opts.reuseCache, ...
        'verbose', opts.verbose);
    data = applyInterPlaneAlignmentPolicy(data, opts);
    if isempty(opts.xyPixelSizeUm)
        opts.xyPixelSizeUm = data.xyPixelSizeUm;
    end

    [selectedIdx, selection] = selectSymmetricPlanes(data.stageZUm, opts);
    [simBase, geometry] = configureSimulation(data, opts);
    validateValidationPrerequisites(data,simBase,geometry,opts);
    nuisance = nominalNuisance(simBase, opts);
    sim = applyNuisance(simBase, geometry, nuisance);
    obs = buildHybridObservations(data, selectedIdx);
    p0 = initialParameterVector(sim, data, opts);

    if opts.verbose
        fprintf('\n[estimateFullStackISMWavefront] %d total planes, %d full-image planes\n', ...
            numel(data.stageZUm), numel(selectedIdx));
        fprintf('  selected stage z [um]: %s\n', num2str(data.stageZUm(selectedIdx), ' %.3f'));
        fprintf('  raw observations in Poisson objective: %d\n', numel(obs.y));
        fprintf('  flat-field calibrated: %d (%s)\n', ...
            data.flatFieldCalibrated, data.flatFieldSource);
        fprintf('  independent background: %d (%s)\n', ...
            data.backgroundIndependent,data.backgroundSource);
        fprintf('  background exposure calibrated: %d\n', ...
            data.backgroundExposureCalibrated);
        fprintf('  inter-plane alignment: %s\n', ...
            data.interPlaneAlignmentPolicy);
    end

    fit = fitHybridPoisson(obs, data, sim, opts, p0, struct());
    refinedSelectionOpts = opts;
    refinedSelectionOpts.selectionCenterUm = fit.estZ0Um;
    [refinedIdx, refinedSelection] = selectSymmetricPlanes( ...
        data.stageZUm, refinedSelectionOpts);
    if ~isequal(refinedIdx(:),selectedIdx(:))
        if opts.verbose
            fprintf(['[estimateFullStackISMWavefront] Reselecting symmetric planes ' ...
                'around fitted z0 = %.4f um\n'],fit.estZ0Um);
        end
        selectedIdx = refinedIdx;
        selection = refinedSelection;
        obs = buildHybridObservations(data,selectedIdx);
        fit = fitHybridPoisson(obs,data,sim,opts,fit.paramVector,struct());
    end
    residualDiagnostics = diagnoseResiduals(obs, fit, data);
    [baseFitPassed,baseFitReason] = baseFitValidationPassed( ...
        fit,residualDiagnostics,opts);
    validationOpts=opts;
    if opts.skipValidationIfBaseFitRejected && ~baseFitPassed
        validationOpts.runProfiles=false;
        validationOpts.nBootstrap=0;
        validationOpts.runSignProfiles=false;
        if opts.verbose
            fprintf(['[estimateFullStackISMWavefront] Expensive validation ' ...
                'skipped: %s\n'],baseFitReason);
        end
    end
    profiles = runNuisanceProfiles(obs, data, simBase, geometry, ...
        nuisance, fit, validationOpts);
    bootstrap = runBootstrap(obs, data, sim, fit, validationOpts);
    signProfiles = runSignProfiles(obs, data, sim, fit, validationOpts);
    acceptance = buildAcceptance(fit, profiles, bootstrap, ...
        signProfiles, residualDiagnostics, data, validationOpts);

    result = struct();
    result.data = data;
    result.selection = selection;
    result.selectedPlaneIndices = selectedIdx;
    result.sim = sim;
    result.geometry = geometry;
    result.nuisance = nuisance;
    result.fit = fit;
    result.profiles = profiles;
    result.bootstrap = bootstrap;
    result.signProfiles = signProfiles;
    result.residualDiagnostics = residualDiagnostics;
    result.baseFitPassed = baseFitPassed;
    result.baseFitReason = baseFitReason;
    result.acceptance = acceptance;
    result.options = opts;
    result.validationOptions = validationOpts;
    result.outputDir = resolveOutputDir(opts);

    if opts.writeOutputs
        writeOutputs(result);
    end
    if opts.verbose
        printSummary(result);
    end
end

function opts = parseOptions(varargin)
    p = inputParser;
    p.FunctionName = 'estimateFullStackISMWavefront';

    addParameter(p, 'stageZUm', []);
    addParameter(p, 'xyPixelSizeUm', []);
    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'channelOrder', []);
    addParameter(p, 'ptuReaderFolder', defaultPtuReaderFolder());
    addParameter(p, 'ptuPhotonsPerChunk', 1e6);
    addParameter(p, 'cacheFile', '');
    addParameter(p, 'reuseCache', false);
    addParameter(p, 'useInterPlaneAlignment', false);
    addParameter(p, 'interPlaneAlignmentMinPx', 4);
    addParameter(p, 'boundaryWidthPx', 3);
    addParameter(p, 'darkFile', 'D:\Luminosa\Data\ISMdark_counts.ptu');
    addParameter(p, 'darkScale', 1);
    addParameter(p, 'darkPerPixel', []);
    addParameter(p, 'darkMeasurementMode', 'tttr');
    addParameter(p, 'backgroundMode', 'auto');
    addParameter(p, 'flatField', []);
    addParameter(p, 'flatFieldVariable', '');

    addParameter(p, 'selectedOffsetsUm', [-0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8]);
    addParameter(p, 'selectionCenterUm', []);
    addParameter(p, 'sim', []);
    addParameter(p, 'objectiveNA', 1.2);
    addParameter(p, 'objectiveMagnification', 60);
    addParameter(p, 'emissionWavelengthUm', []);
    addParameter(p, 'excitationWavelengthUm', []);
    addParameter(p, 'mediumRefractiveIndex', []);
    addParameter(p, 'sampleGeometry', 'airOnGlass');
    addParameter(p, 'interfaceRadialWeightMode', 'sampleSolidAngle');
    addParameter(p, 'diffractionModel', 'vectorial');
    addParameter(p, 'vectorialPolarizationMode', 'xyAverage');
    addParameter(p, 'excitationPolarizationMode', 'x');
    addParameter(p, 'collectionPolarizationMode', 'xyAverage');
    % Stage/focus defocus reference medium: 'immersion' (piezo moves the
    % coverslip through the immersion gap) or 'glass' (legacy).
    addParameter(p, 'airInterfaceStageMedium', 'immersion');
    addParameter(p, 'immersionRefractiveIndex', 1.33);
    addParameter(p, 'glassRefractiveIndex', 1.518);
    addParameter(p, 'sampleRefractiveIndex', 1.0003);
    addParameter(p, 'designGlassRefractiveIndex', 1.518);
    addParameter(p, 'coverslipThicknessUm', 190);
    addParameter(p, 'designCoverslipThicknessUm', 190);
    addParameter(p, 'beadBottomHeightUm', 0);
    addParameter(p, 'airBeadAxialSamples', 3);
    addParameter(p, 'interfaceGeometryConfirmed', false);
    addParameter(p, 'runAirInterfaceSelfTest', true);
    addParameter(p, 'beadRadiusUm', 0.08);
    addParameter(p, 'beadSubsamples', [3 3 3]);
    addParameter(p, 'modelDzUm', 0.025);
    addParameter(p, 'modelZPaddingUm', 0.30);
    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detectorPixelShape', 'hex');
    addParameter(p, 'estimateDetectorLayout', false);
    addParameter(p, 'detectorXYUm', []);
    addParameter(p, 'detectorSubsamples', 5);
    addParameter(p, 'detectorLayoutScale', 2);
    addParameter(p, 'detectorLayoutPositionSign', -1);
    addParameter(p, 'scanAxisSigns', [1 1]);
    addParameter(p, 'detectorShiftUpsample', 20);
    addParameter(p, 'detectorShiftPlaneIndex', []);
    addParameter(p, 'detectorPitchSampleUm', []);
    addParameter(p, 'detectorHardwarePitchUm', []);
    addParameter(p, 'detectorTotalMagnification', []);

    addParameter(p, 'fitModes', ...
        {'astig_x','astig_y','coma_x','coma_y','spherical'});
    addParameter(p, 'initialCoeffs', struct());
    addParameter(p, 'initialResultFile', '');
    addParameter(p, 'initialZ0Um', []);
    addParameter(p, 'initialXYUm', [0 0]);
    addParameter(p, 'fitXY', false);
    addParameter(p, 'maxIter', 4);
    addParameter(p, 'fdCoeff', 0.01);
    addParameter(p, 'fdZ0Um', 0.03);
    addParameter(p, 'fdXYUm', []);
    addParameter(p, 'maxCoeffStep', 0.04);
    addParameter(p, 'maxZ0StepUm', 0.08);
    addParameter(p, 'maxXYStepUm', 0.04);
    addParameter(p, 'regCoeff', 1e-4);
    addParameter(p, 'regZ0', 1e-5);
    addParameter(p, 'regXY', 1e-4);
    addParameter(p, 'tolStep', 1e-4);
    addParameter(p, 'coefficientBoundsWaves', [-0.30 0.30]);
    addParameter(p, 'z0BoundsUm', []);
    addParameter(p, 'z0SearchHalfRangeUm', 0.60);
    addParameter(p, 'xyBoundsUm', [-0.30 0.30]);

    addParameter(p, 'runProfiles', false);
    addParameter(p, 'profileRefitIterations', 1);
    addParameter(p, 'beadRadiusProfileUm', []);
    addParameter(p, 'detectorScaleProfile', [1.8 2.0 2.2]);
    addParameter(p, 'detectorRotationProfileDeg', [-3 0 3]);
    addParameter(p, 'detectorAnisotropyProfile', [0.95 1 1.05]);
    addParameter(p, 'emissionWavelengthProfileUm', []);
    addParameter(p, 'objectiveNAProfile', []);
    addParameter(p, 'coverslipThicknessProfileUm', []);
    addParameter(p, 'immersionIndexProfile', []);
    addParameter(p, 'glassIndexProfile', []);
    addParameter(p, 'beadBottomHeightProfileUm', []);

    addParameter(p, 'nBootstrap', 0);
    addParameter(p, 'bootstrapMaxIter', 2);
    addParameter(p, 'runSignProfiles', false);
    addParameter(p, 'signProfileMaxIter', 2);
    addParameter(p, 'signRatioThreshold', 1.10);
    addParameter(p, 'bootstrapSignThreshold', 0.95);
    addParameter(p, 'pearsonVarianceRange', [0.5 1.5]);
    addParameter(p, 'maxSpatialResidualCorrelation', 0.10);
    addParameter(p, 'axialPearsonVarianceRange', [0.5 2.0]);
    addParameter(p, 'maxAxialResidualCorrelation', 0.30);
    addParameter(p, 'scalarModelFloorWaves', 0.15);
    addParameter(p, 'knownAberrationCalibration', []);
    addParameter(p, 'requireValidationPrerequisites', false);
    addParameter(p, 'skipValidationIfBaseFitRejected', true);

    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'verbose', true);
    parse(p, varargin{:});
    opts = p.Results;

    opts.channelIDs = double(opts.channelIDs(:)).';
    opts.channelOrder = double(opts.channelOrder(:)).';
    opts.selectedOffsetsUm = unique(double(opts.selectedOffsetsUm(:)).');
    if ischar(opts.fitModes)
        opts.fitModes = {opts.fitModes};
    elseif isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    elseif ~iscell(opts.fitModes)
        error('estimateFullStackISMWavefront:BadFitModes', ...
            'fitModes must be a character vector, string array, or cell array.');
    end
    opts.ptuReaderFolder = char(opts.ptuReaderFolder);
    opts.cacheFile = char(opts.cacheFile);
    opts.flatFieldVariable = char(opts.flatFieldVariable);
    opts.darkFile = char(opts.darkFile);
    opts.darkMeasurementMode = lower(char(opts.darkMeasurementMode));
    opts.backgroundMode = lower(char(opts.backgroundMode));
    opts.detectorLayout = char(opts.detectorLayout);
    opts.detectorPixelShape = char(opts.detectorPixelShape);
    opts.initialResultFile = char(opts.initialResultFile);
    opts.outputDir = char(opts.outputDir);
    opts.sampleGeometry = char(opts.sampleGeometry);
    opts.diffractionModel = char(opts.diffractionModel);
    opts.vectorialPolarizationMode = char(opts.vectorialPolarizationMode);
    opts.excitationPolarizationMode = char(opts.excitationPolarizationMode);
    opts.collectionPolarizationMode = char(opts.collectionPolarizationMode);
    opts.airInterfaceStageMedium = char(opts.airInterfaceStageMedium);
    opts.interPlaneAlignmentMinPx = double(opts.interPlaneAlignmentMinPx);
    if ~isscalar(opts.interPlaneAlignmentMinPx) || ...
            ~isfinite(opts.interPlaneAlignmentMinPx) || opts.interPlaneAlignmentMinPx < 0
        error('estimateFullStackISMWavefront:BadAlignmentThreshold', ...
            'interPlaneAlignmentMinPx must be a finite nonnegative scalar.');
    end
    if ~ismember(lower(opts.sampleGeometry),{'aironglass','homogeneous'})
        error('estimateFullStackISMWavefront:BadSampleGeometry', ...
            'sampleGeometry must be airOnGlass or homogeneous.');
    end
end

function folder = defaultPtuReaderFolder()
    thisDir = fileparts(mfilename('fullpath'));
    folder = fullfile(fileparts(fileparts(thisDir)), 'Luminosa_FLIM_FCS');
end

function addRequiredPaths(opts)
    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    addpath(fileparts(thisDir));
    if ~isempty(opts.ptuReaderFolder) && exist(opts.ptuReaderFolder,'dir') == 7
        addpath(opts.ptuReaderFolder);
    end
end

function data = applyInterPlaneAlignmentPolicy(data, opts)
    nPlane = numel(data.stageZUm);
    data.useInterPlaneAlignment = logical(opts.useInterPlaneAlignment);
    if data.useInterPlaneAlignment
        if ~isfield(data, 'modelShiftToRawXPx') || isempty(data.modelShiftToRawXPx)
            data.modelShiftToRawXPx = zeros(1, nPlane);
        end
        if ~isfield(data, 'modelShiftToRawYPx') || isempty(data.modelShiftToRawYPx)
            data.modelShiftToRawYPx = zeros(1, nPlane);
        end
        shiftMag = hypot(data.modelShiftToRawXPx(:), data.modelShiftToRawYPx(:));
        small = shiftMag < opts.interPlaneAlignmentMinPx;
        data.modelShiftToRawXPx(small) = 0;
        data.modelShiftToRawYPx(small) = 0;
        data.interPlaneAlignmentThresholdPx = opts.interPlaneAlignmentMinPx;
        data.nInterPlaneShiftsApplied = nnz(~small);
        data.nInterPlaneShiftsIgnored = nnz(small);
        data.interPlaneAlignmentPolicy = ...
            sprintf(['using saved model-to-raw shifts from alignment metadata ' ...
            'only when magnitude >= %.3g px'], opts.interPlaneAlignmentMinPx);
        return;
    end

    data.modelShiftToRawXPx = zeros(1, nPlane);
    data.modelShiftToRawYPx = zeros(1, nPlane);
    data.interPlaneAlignmentThresholdPx = opts.interPlaneAlignmentMinPx;
    data.nInterPlaneShiftsApplied = 0;
    data.nInterPlaneShiftsIgnored = nPlane;
    data.interPlaneAlignmentPolicy = ...
        ['ignored saved inter-plane alignment shifts; raw frames are treated ' ...
        'as already aligned in scan coordinates'];
end

function validateValidationPrerequisites(data,sim,geometry,opts)
    if ~opts.requireValidationPrerequisites
        return;
    end
    missing={};
    if ~data.flatFieldCalibrated
        missing{end+1}='independent detector flat field';
    end
    if ~data.backgroundIndependent || ~data.backgroundExposureCalibrated
        missing{end+1}='exposure-calibrated independent dark background';
    end
    if strcmpi(sim.sampleGeometry,'airOnGlass') && ...
            ~opts.interfaceGeometryConfirmed
        missing{end+1}='confirmed glass/interface geometry';
    end
    if ~geometry.physicalPitchCalibrated
        missing{end+1}='detector pitch in sample space';
    end
    oddRequested=any(ismember(opts.fitModes, ...
        {'coma_x','coma_y','tilt_x','tilt_y'}));
    if oddRequested && ~knownCalibrationPassed(opts.knownAberrationCalibration)
        missing{end+1}='successful known-coma sign calibration';
    end
    if ~isempty(missing)
        error('estimateFullStackISMWavefront:ValidationPrerequisites', ...
            ['Validated retrieval was not started because these required ' ...
            'inputs are missing: %s. Run ' ...
            'runDiagnosticAirSurfaceISMWavefront for a non-correction ' ...
            'diagnostic fit, or supply all calibration inputs.'], ...
            strjoin(missing,', '));
    end
end

function [idx, info] = selectSymmetricPlanes(stageZ, opts)
    center = opts.selectionCenterUm;
    if isempty(center)
        center = median(stageZ);
    end
    if isempty(opts.selectedOffsetsUm)
        idx = 1:numel(stageZ);
        info = table(idx(:), stageZ(idx(:)).', ...
            'VariableNames', {'planeIndex','stageZUm'});
        info.requestedCenterUm = repmat(center,height(info),1);
        info.selectionMode = repmat({'allPlanesFullImages'},height(info),1);
        return;
    end
    targets = center + opts.selectedOffsetsUm;
    idx = zeros(size(targets));
    for k = 1:numel(targets)
        [~,idx(k)] = min(abs(stageZ-targets(k)));
    end
    idx = unique(idx,'stable');
    info = table(idx(:), stageZ(idx(:)).', ...
        'VariableNames', {'planeIndex','stageZUm'});
    info.requestedCenterUm = repmat(center,height(info),1);
    info.selectionMode = repmat({'symmetricOffsets'},height(info),1);
end

function [sim, geometry] = configureSimulation(data, opts)
    if isempty(opts.sim)
        sim = defaultParams();
    else
        sim = opts.sim;
    end
    sim = applyOpticalOptions(sim, data.head, opts);
    sim.beadRadius = opts.beadRadiusUm;
    sim.beadSubsamples = opts.beadSubsamples;
    sim.sampleGeometry = opts.sampleGeometry;
    sim.interfaceRadialWeightMode = char(opts.interfaceRadialWeightMode);
    sim.airInterfaceStageMedium = char(opts.airInterfaceStageMedium);
    sim.diffractionModel = char(opts.diffractionModel);
    sim.vectorialPolarizationMode = char(opts.vectorialPolarizationMode);
    sim.excitationPolarizationMode = char(opts.excitationPolarizationMode);
    sim.collectionPolarizationMode = char(opts.collectionPolarizationMode);
    sim.includesVectorialPolarization = usesVectorialPSF(sim);
    sim.nImmersion = opts.immersionRefractiveIndex;
    sim.nGlass = opts.glassRefractiveIndex;
    sim.nSample = opts.sampleRefractiveIndex;
    sim.nDesignGlass = opts.designGlassRefractiveIndex;
    sim.coverslipThicknessUm = opts.coverslipThicknessUm;
    sim.designCoverslipThicknessUm = opts.designCoverslipThicknessUm;
    sim.beadBottomHeightUm = opts.beadBottomHeightUm;
    sim.airBeadAxialSamples = opts.airBeadAxialSamples;
    sim.nMedium = sim.nImmersion;
    if strcmpi(sim.sampleGeometry,'airOnGlass')
        if sim.includesVectorialPolarization
            sim.modelAccuracyNote = ['Vectorial Richards-Wolf/Bessel ' ...
                'water/glass/air model with s/p Fresnel transmission.'];
        else
            sim.modelAccuracyNote = ['Scalar s/p Fresnel Debye-Bessel ' ...
                'water/glass/air model; vectorial polarization disabled.'];
        end
        if opts.runAirInterfaceSelfTest
            if sim.includesVectorialPolarization
                sim.airInterfaceValidation = validateVectorialForwardModel( ...
                    'interfaceRadialWeightMode', sim.interfaceRadialWeightMode, ...
                    'airInterfaceStageMedium', sim.airInterfaceStageMedium, ...
                    'excitationPolarizationMode', sim.excitationPolarizationMode, ...
                    'collectionPolarizationMode', sim.collectionPolarizationMode);
                sim.vectorialValidation = sim.airInterfaceValidation;
                if ~sim.airInterfaceValidation.allPassed
                    error('estimateFullStackISMWavefront:VectorialSelfTest', ...
                        'The vectorial forward-model self-test failed.');
                end
            else
                sim.airInterfaceValidation = validateAirInterfaceForwardModel( ...
                    'interfaceRadialWeightMode', sim.interfaceRadialWeightMode, ...
                    'airInterfaceStageMedium', sim.airInterfaceStageMedium);
                if ~sim.airInterfaceValidation.allPassed
                    error('estimateFullStackISMWavefront:AirInterfaceSelfTest', ...
                        'The air-interface forward-model self-test failed.');
                end
            end
        end
    end
    validateScalarOptics(sim);
    sim.detectorSubsamples = opts.detectorSubsamples;
    sim.detectorLayout = opts.detectorLayout;
    sim.detectorPixelShape = opts.detectorPixelShape;
    [sim.detXY, sim.detectorIndexGrid, sim.detectorLayoutInfo] = ...
        detectorLayout(sim.detectorLayout, sim.detPitch);
    sim.nDet = size(sim.detXY, 1);
    sim.detectorGridSize = size(sim.detectorIndexGrid);
    if strcmpi(sim.detectorPixelShape, 'hex')
        sim.detectorHexRadius = sim.detSize / sqrt(3);
    end
    sim.nx = size(data.rawCounts,2);
    sim.ny = size(data.rawCounts,1);
    sim.xyPixelSizeUm = data.xyPixelSizeUm;
    sim.x = data.xUm;
    sim.y = data.yUm;
    sim.dx = data.xyPixelSizeUm;
    sim.fovX = max(sim.x)-min(sim.x);
    sim.fovY = max(sim.y)-min(sim.y);
    sim.fovXY = max(sim.fovX,sim.fovY);

    z0Bounds = resolveZ0Bounds(data.stageZUm, opts);
    zMin = min(data.stageZUm)-max(z0Bounds)-opts.modelZPaddingUm-sim.beadRadius;
    zMax = max(data.stageZUm)-min(z0Bounds)+opts.modelZPaddingUm+sim.beadRadius;
    nZ = ceil((zMax-zMin)/opts.modelDzUm)+1;
    if mod(nZ,2)==0, nZ=nZ+1; end
    sim.z = linspace(zMin,zMax,nZ);
    sim.nz = nZ;
    sim.nzRange = zMax-zMin;

    shiftPlane = opts.detectorShiftPlaneIndex;
    if isempty(shiftPlane)
        shiftPlane = selectSharpRegistrationPlane(data.rawCounts);
    end
    airyUnitUm = 1.22*sim.lamEm/effectivePropagatingNA(sim);
    if ~isempty(opts.detectorXYUm)
        detXY = validateDetectorXYFullStack(opts.detectorXYUm, size(data.rawCounts,3));
        diagnostics = fixedDetectorDiagnosticsFullStack( ...
            detXY, data.xyPixelSizeUm, airyUnitUm, 'explicit detectorXYUm option');
    elseif opts.estimateDetectorLayout
        layoutOpts = struct();
        layoutOpts.positionSign = opts.detectorLayoutPositionSign;
        layoutOpts.scanAxisSigns = opts.scanAxisSigns;
        layoutOpts.detectorScale = opts.detectorLayoutScale;
        layoutOpts.centerMode = 'reference';
        layoutOpts.upsampleReg = opts.detectorShiftUpsample;
        layoutOpts.airyUnitUm = airyUnitUm;
        [detXY, diagnostics] = estimateDetectorLayoutFromStack( ...
            data.rawCounts(:,:,:,shiftPlane), data.xyPixelSizeUm, layoutOpts);
        diagnostics.source = 'raw ISM phase correlation';
    else
        detXY = sim.detXY;
        diagnostics = fixedDetectorDiagnosticsFullStack( ...
            detXY, data.xyPixelSizeUm, airyUnitUm, ...
            sprintf('fixed regular %s detector layout', sim.detectorLayout));
    end
    diagnostics.detectorPitchSource = detectorPitchSourceFullStack(opts);
    calibratedPitch = resolveDetectorPitchSampleUm(opts);
    recoveredPitch = medianNearestDistance(detXY);
    if isfinite(calibratedPitch) && calibratedPitch > 0 && ...
            isfinite(recoveredPitch) && recoveredPitch > 0
        detXY = detXY*(calibratedPitch/recoveredPitch);
        diagnostics.detXYBeforePitchCalibrationUm = diagnostics.detXYUm;
        diagnostics.detXYUm = detXY;
        diagnostics.detXYNm = 1000*detXY;
        diagnostics.pitchCalibrationScale = calibratedPitch/recoveredPitch;
        diagnostics.calibratedPitchUm = calibratedPitch;
    end

    geometry = struct();
    geometry.shiftPlaneIndex = shiftPlane;
    geometry.shiftPlaneStageZUm = data.stageZUm(shiftPlane);
    geometry.diagnostics = diagnostics;
    geometry.unitScaleXYUm = detXY / opts.detectorLayoutScale;
    geometry.nominalScale = opts.detectorLayoutScale;
    geometry.nominalPitchUm = medianNearestDistance(detXY);
    geometry.physicalPitchCalibrated = isfinite(calibratedPitch);
    geometry.calibratedPitchUm = calibratedPitch;
    geometry.nominalFillRatio = sim.detSize/sim.detPitch;

    sim.detXY = detXY;
    sim.detectorPhysicalPitchCalibrated = geometry.physicalPitchCalibrated;
    sim.detectorCalibratedPitchUm = geometry.calibratedPitchUm;
    sim = updateDetectorSize(sim, geometry);
    if ~strcmpi(sim.sampleGeometry,'airOnGlass')
        sim.obj = beadObject3D(sim);
    end
end

function sim = applyOpticalOptions(sim, head, opts)
    if ~isempty(opts.objectiveNA), sim.NA=opts.objectiveNA; end
    if ~isempty(opts.objectiveMagnification)
        sim.objectiveMagnification=opts.objectiveMagnification;
    end
    if ~isempty(opts.emissionWavelengthUm)
        sim.lamEm=opts.emissionWavelengthUm;
        sim.lamRef=opts.emissionWavelengthUm;
    end
    if ~isempty(opts.excitationWavelengthUm), sim.lamExc=opts.excitationWavelengthUm; end
    if ~isempty(opts.mediumRefractiveIndex)
        sim.nMedium=opts.mediumRefractiveIndex;
        sim.nImmersion=opts.mediumRefractiveIndex;
    end

    if isempty(opts.objectiveNA)
        value = headerNumericByPattern(head, {'Objective.*NA','Numerical.*Aperture','(^|_)NA($|_)'});
        if isfinitePositive(value), sim.NA=value; end
    end
    if isempty(opts.objectiveMagnification)
        value = headerNumericByPattern(head, {'Magnification','Objective.*Mag','Obj.*Mag'});
        if isfinitePositive(value), sim.objectiveMagnification=value; end
    end
end

function idx = selectSharpRegistrationPlane(raw)
    nPlane = size(raw,4);
    score = zeros(1,nPlane);
    for iz = 1:nPlane
        img = sum(raw(:,:,:,iz),3);
        gx = diff(img,1,2);
        gy = diff(img,1,1);
        score(iz) = sum(gx(:).^2)+sum(gy(:).^2);
    end
    [~,idx]=max(score);
end

function nuisance = nominalNuisance(sim, opts)
    nuisance = struct();
    nuisance.beadRadiusUm = sim.beadRadius;
    nuisance.detectorScale = opts.detectorLayoutScale;
    nuisance.detectorRotationDeg = 0;
    nuisance.detectorAnisotropy = 1;
    nuisance.emissionWavelengthUm = sim.lamEm;
    nuisance.objectiveNA = sim.NA;
    nuisance.coverslipThicknessUm = sim.coverslipThicknessUm;
    nuisance.immersionRefractiveIndex = sim.nImmersion;
    nuisance.glassRefractiveIndex = sim.nGlass;
    nuisance.beadBottomHeightUm = sim.beadBottomHeightUm;
end

function sim = applyNuisance(simBase, geometry, nuisance)
    sim = simBase;
    sim.beadRadius = nuisance.beadRadiusUm;
    sim.lamEm = nuisance.emissionWavelengthUm;
    sim.NA = nuisance.objectiveNA;
    sim.coverslipThicknessUm = nuisance.coverslipThicknessUm;
    sim.nImmersion = nuisance.immersionRefractiveIndex;
    sim.nMedium = sim.nImmersion;
    sim.nGlass = nuisance.glassRefractiveIndex;
    sim.beadBottomHeightUm = nuisance.beadBottomHeightUm;

    theta = deg2rad(nuisance.detectorRotationDeg);
    rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    anisotropy = max(nuisance.detectorAnisotropy,eps);
    stretch = diag([sqrt(anisotropy),1/sqrt(anisotropy)]);
    sim.detXY = geometry.unitScaleXYUm * nuisance.detectorScale * ...
        rotation.' * stretch;
    sim = updateDetectorSize(sim, geometry);
    validateScalarOptics(sim);
    if ~strcmpi(sim.sampleGeometry,'airOnGlass')
        sim.obj = beadObject3D(sim);
    end
end

function sim = updateDetectorSize(sim, geometry)
    pitch = medianNearestDistance(sim.detXY);
    if isfinitePositive(pitch)
        sim.detPitch = pitch;
        sim.detSize = geometry.nominalFillRatio*pitch;
        if isfield(sim,'detectorPixelShape') && strcmpi(sim.detectorPixelShape,'hex')
            sim.detectorHexRadius = sim.detSize/sqrt(3);
        end
    end
    sim.nDet = size(sim.detXY,1);
end

function obs = buildHybridObservations(data, selectedIdx)
    nPlane = size(data.rawCounts,4);
    selectedMask = false(1,nPlane);
    selectedMask(selectedIdx)=true;
    unselectedIdx = find(~selectedMask);

    selectedData = data.rawCounts(:,:,:,selectedIdx);
    selectedBg = repmat(data.backgroundPerPixel(:,:,:,selectedIdx), ...
        size(data.rawCounts,1),size(data.rawCounts,2),1,1);
    unselectedData = data.axialRawCounts(unselectedIdx);
    unselectedBg = data.axialBackgroundCounts(unselectedIdx);

    obs = struct();
    obs.selectedIdx = selectedIdx(:).';
    obs.unselectedIdx = unselectedIdx(:).';
    obs.selectedCount = numel(selectedData);
    obs.y = [selectedData(:); unselectedData(:)];
    obs.background = [selectedBg(:); unselectedBg(:)];
    obs.selectedSize = size(selectedData);
    obs.gain = reshape(data.flatFieldGain,1,1,[],1);
end

function p = initialParameterVector(sim, data, opts)
    coeffs = opts.initialCoeffs;
    if ~isempty(opts.initialResultFile)
        coeffs = mergeInitialCoefficients(coeffs,opts.initialResultFile);
    end
    p = zeros(1,numel(opts.fitModes));
    for k=1:numel(opts.fitModes)
        if isfield(coeffs,opts.fitModes{k})
            p(k)=double(coeffs.(opts.fitModes{k}));
        end
    end
    z0 = opts.initialZ0Um;
    if isempty(z0), z0=median(data.stageZUm); end
    p=[p z0];
    if opts.fitXY, p=[p double(opts.initialXYUm(:)).']; end
    p = clampParameters(p,data,opts);
    if ~all(isfinite(p))
        error('estimateFullStackISMWavefront:BadInitialParameters', ...
            'Initial parameter vector contains nonfinite values.');
    end
    if ~isfield(sim,'modeOrder')
        error('estimateFullStackISMWavefront:BadSimulation', ...
            'Simulation structure must contain modeOrder.');
    end
end

function coeffs = mergeInitialCoefficients(coeffs,fileName)
    S=load(fileName);
    value=[];
    names=fieldnames(S);
    for k=1:numel(names)
        candidate=S.(names{k});
        if isstruct(candidate)
            if isfield(candidate,'fit'), candidate=candidate.fit; end
            if isfield(candidate,'estCoeffs')
                value=candidate.estCoeffs;
                break;
            end
        end
    end
    if isempty(value), return; end
    fields=fieldnames(value);
    for k=1:numel(fields), coeffs.(fields{k})=value.(fields{k}); end
end

function fit = fitHybridPoisson(obs,data,sim,opts,p0,constraints)
    p=p0;
    nParam=numel(p);
    steps=finiteDifferenceSteps(opts,nParam);
    reg=regularizationVector(opts,nParam);
    maxStep=maxStepVector(opts,nParam);
    maxIter=opts.maxIter;
    if isfield(constraints,'maxIter'), maxIter=constraints.maxIter; end
    history=nan(maxIter,4);
    converged=false;
    terminationReason='maximum iterations reached';

    [r,deviance,scale,mu,model]=evaluateParameters( ...
        p,obs,data,sim,opts,constraints);
    for it=1:maxIter
        J=zeros(numel(r),nParam);
        for q=1:nParam
            pp=p;
            pp(q)=pp(q)+steps(q);
            pp=applyConstraints(clampParameters(pp,data,opts),constraints);
            rp=evaluateParameters(pp,obs,data,sim,opts,constraints);
            J(:,q)=(rp-r)/max(abs(pp(q)-p(q)),eps);
        end
        lhs=J.'*J+diag(reg);
        rhs=-J.'*r;
        delta=lhs\rhs;
        delta=max(min(delta(:).',maxStep),-maxStep);
        if norm(delta)<opts.tolStep
            converged=true;
            terminationReason='Gauss-Newton step below tolerance';
            history=history(1:max(it-1,0),:);
            break;
        end

        accepted=false;
        lineScale=1;
        for ls=1:8
            trial=applyConstraints(clampParameters(p+lineScale*delta,data,opts),constraints);
            [rt,dt,st,mut,modelt]=evaluateParameters( ...
                trial,obs,data,sim,opts,constraints);
            if dt<deviance
                p=trial; r=rt; deviance=dt; scale=st; mu=mut; model=modelt;
                accepted=true;
                break;
            end
            lineScale=lineScale/2;
        end
        history(it,:)=[deviance,norm(lineScale*delta),scale,accepted];
        if opts.verbose && ~isfield(constraints,'quiet')
            fprintf('[full-stack Poisson] iter %d: deviance %.5g, step %.3g, scale %.4g\n', ...
                it,deviance,norm(lineScale*delta),scale);
        end
        if ~accepted
            terminationReason='line search failed';
            history=history(1:it,:);
            break;
        end
        if norm(lineScale*delta)<opts.tolStep
            converged=true;
            terminationReason='accepted step below tolerance';
            history=history(1:it,:);
            break;
        end
        if it==maxIter, history=history(1:it,:); end
    end

    fit=packFit(p,r,deviance,scale,mu,model,history,sim,opts);
    fit.converged=converged;
    fit.terminationReason=terminationReason;
end

function [residual,deviance,scale,mu,model] = evaluateParameters( ...
        p,obs,data,sim,opts,constraints)
    p=applyConstraints(p,constraints);
    [coeffs,z0,xy]=unpackParameters(p,sim,opts);
    relativeZ=data.stageZUm-z0;
    if strcmpi(sim.sampleGeometry,'airOnGlass')
        model=normalizedStackAirInterfaceZPlanes( ...
            sim,coeffs,relativeZ,xy(1),xy(2),0);
    else
        model=normalizedStackExplicitDetectorZPlanes( ...
            sim,coeffs,relativeZ,xy(1),xy(2),0);
    end

    for iz=1:size(model,4)
        sx=data.modelShiftToRawXPx(iz);
        sy=data.modelShiftToRawYPx(iz);
        if sx~=0 || sy~=0
            for ch=1:size(model,3)
                model(:,:,ch,iz)=max(fourierShift2D( ...
                    model(:,:,ch,iz),sx,sy),0);
            end
        end
    end
    model=model.*obs.gain;
    unit=hybridModelVector(model,obs);
    scale=profileGlobalScale(obs.y,unit,obs.background);
    mu=max(scale*unit+obs.background,realmin);
    residual=poissonDevianceResidual(obs.y,mu);
    deviance=sum(residual.^2);
end

function unit = hybridModelVector(model,obs)
    selected=model(:,:,:,obs.selectedIdx);
    if isempty(obs.unselectedIdx)
        trace=[];
    else
        trace=squeeze(sum(sum(sum(model(:,:,:,obs.unselectedIdx),1),2),3));
    end
    unit=[selected(:);trace(:)];
end

function scale = profileGlobalScale(y,m,b)
    m=max(double(m(:)),0);
    y=double(y(:));
    b=max(double(b(:)),0);
    if sum(m)<=0
        scale=0;
        return;
    end
    signal=max(sum(y-b),0);
    guess=max(signal/sum(m),1);
    derivative=@(a) sum(m.*(1-y./max(a*m+b,realmin)));
    lo=0;
    hi=guess;
    while derivative(hi)<0 && hi<1e18
        hi=2*hi;
    end
    for k=1:60
        mid=(lo+hi)/2;
        if derivative(mid)>0, hi=mid; else, lo=mid; end
    end
    scale=(lo+hi)/2;
end

function r = poissonDevianceResidual(y,mu)
    y=double(y(:));
    mu=max(double(mu(:)),realmin);
    term=mu-y;
    positive=y>0;
    term(positive)=y(positive).*log(y(positive)./mu(positive))- ...
        (y(positive)-mu(positive));
    term=max(term,0);
    r=sign(y-mu).*sqrt(2*term);
end

function fit = packFit(p,r,deviance,scale,mu,model,history,sim,opts)
    [coeffs,z0,xy]=unpackParameters(p,sim,opts);
    vector=zeros(1,numel(sim.modeOrder));
    for k=1:numel(sim.modeOrder)
        if isfield(coeffs,sim.modeOrder{k}), vector(k)=coeffs.(sim.modeOrder{k}); end
    end
    fit=struct();
    fit.paramVector=p;
    fit.estCoeffs=coeffs;
    fit.estCoeffVector=vector;
    fit.estZ0Um=z0;
    fit.estXYUm=xy;
    fit.globalPhotonScale=scale;
    fit.expectedObservation=mu;
    fit.devianceResidual=r;
    fit.deviance=deviance;
    fit.residualNorm=sqrt(deviance);
    fit.model=model;
    fit.history=history;
    fit.sim=sim;
end

function [coeffs,z0,xy] = unpackParameters(p,sim,opts)
    full=zeros(1,numel(sim.modeOrder));
    for k=1:numel(opts.fitModes)
        idx=find(strcmp(sim.modeOrder,opts.fitModes{k}),1);
        if isempty(idx)
            error('estimateFullStackISMWavefront:UnknownMode', ...
                'Unknown mode %s.',opts.fitModes{k});
        end
        full(idx)=p(k);
    end
    coeffs=coeffStruct(sim,full);
    next=numel(opts.fitModes)+1;
    z0=p(next); next=next+1;
    xy=[0 0];
    if opts.fitXY, xy=p(next:next+1); end
end

function steps = finiteDifferenceSteps(opts,nParam)
    steps=opts.fdCoeff*ones(1,numel(opts.fitModes));
    steps=[steps opts.fdZ0Um];
    if opts.fitXY
        fd=opts.fdXYUm;
        if isempty(fd), fd=max(opts.xyPixelSizeUm/4,0.005); end
        if isempty(fd) || ~isfinite(fd), fd=0.01; end
        steps=[steps fd fd];
    end
    steps=steps(1:nParam);
end

function reg = regularizationVector(opts,nParam)
    reg=opts.regCoeff*ones(1,numel(opts.fitModes));
    reg=[reg opts.regZ0];
    if opts.fitXY, reg=[reg opts.regXY opts.regXY]; end
    reg=reg(1:nParam);
end

function values = maxStepVector(opts,nParam)
    values=opts.maxCoeffStep*ones(1,numel(opts.fitModes));
    values=[values opts.maxZ0StepUm];
    if opts.fitXY, values=[values opts.maxXYStepUm opts.maxXYStepUm]; end
    values=values(1:nParam);
end

function p = clampParameters(p,data,opts)
    nMode=numel(opts.fitModes);
    bounds=opts.coefficientBoundsWaves;
    p(1:nMode)=max(min(p(1:nMode),bounds(2)),bounds(1));
    zBounds=resolveZ0Bounds(data.stageZUm,opts);
    p(nMode+1)=max(min(p(nMode+1),zBounds(2)),zBounds(1));
    if opts.fitXY
        xyIdx=nMode+2:nMode+3;
        p(xyIdx)=max(min(p(xyIdx),opts.xyBoundsUm(2)),opts.xyBoundsUm(1));
    end
end

function bounds = resolveZ0Bounds(stageZ,opts)
    bounds=opts.z0BoundsUm;
    if isempty(bounds)
        center=median(stageZ);
        bounds=center+[-1 1]*opts.z0SearchHalfRangeUm;
    end
end

function p = applyConstraints(p,constraints)
    if isfield(constraints,'signIndex') && ~isempty(constraints.signIndex)
        idx=constraints.signIndex;
        signValue=constraints.signValue;
        p(idx)=signValue*max(abs(p(idx)),1e-5);
    end
end

function profiles = runNuisanceProfiles(obs,data,simBase,geometry,nominal,fit,opts)
    profiles=struct('table',table(),'coefficientSamples',[]);
    if ~opts.runProfiles
        return;
    end
    definitions=profileDefinitions(simBase,nominal,opts);
    totalFits=sum(cellfun(@(x) numel(x.values),definitions));
    completed=0;
    totalTimer=tic;
    rows={};
    coeffSamples=[];
    for d=1:numel(definitions)
        def=definitions{d};
        for k=1:numel(def.values)
            completed=completed+1;
            fitTimer=tic;
            if opts.verbose
                fprintf(['[full-stack nuisance] %d/%d: %s = %.6g ' ...
                    '(elapsed %.1f min)\n'],completed,totalFits,def.field, ...
                    def.values(k),toc(totalTimer)/60);
            end
            nuisance=nominal;
            nuisance.(def.field)=def.values(k);
            sim=applyNuisance(simBase,geometry,nuisance);
            localOpts=opts;
            localOpts.maxIter=opts.profileRefitIterations;
            localOpts.verbose=false;
            constraints=struct('maxIter',opts.profileRefitIterations,'quiet',true);
            one=fitHybridPoisson(obs,data,sim,localOpts,fit.paramVector,constraints);
            rows(end+1,:)={def.field,def.values(k),one.deviance}; %#ok<AGROW>
            coeffSamples(end+1,:)=fitModeVector(one,opts.fitModes); %#ok<AGROW>
            if opts.verbose
                fprintf('  completed in %.1f min; deviance %.5g\n', ...
                    toc(fitTimer)/60,one.deviance);
            end
        end
    end
    profiles.table=cell2table(rows, ...
        'VariableNames',{'profileParameter','profileValue','deviance'});
    profiles.coefficientSamples=coeffSamples;
    profiles.coefficientSystematicHalfRange= ...
        0.5*(max(coeffSamples,[],1)-min(coeffSamples,[],1));
    profiles.table.deltaDeviance=profiles.table.deviance-min(profiles.table.deviance);
    coeffTable=array2table(coeffSamples,'VariableNames', ...
        matlab.lang.makeValidName(strcat('coeff_',opts.fitModes)));
    profiles.table=[profiles.table coeffTable];
end

function definitions = profileDefinitions(sim,nominal,opts)
    bead=opts.beadRadiusProfileUm;
    if isempty(bead), bead=max(0.01,nominal.beadRadiusUm+[-0.02 0 0.02]); end
    lam=opts.emissionWavelengthProfileUm;
    if isempty(lam), lam=nominal.emissionWavelengthUm+[-0.01 0 0.01]; end
    na=opts.objectiveNAProfile;
    if isempty(na), na=nominal.objectiveNA+[-0.05 0 0.05]; end
    na=na(na>0 & na<sim.nImmersion);
    cover=opts.coverslipThicknessProfileUm;
    if isempty(cover)
        cover=nominal.coverslipThicknessUm+[-10 0 10];
    end
    definitions={ ...
        struct('field','beadRadiusUm','values',unique(bead)), ...
        struct('field','detectorScale','values',unique(opts.detectorScaleProfile)), ...
        struct('field','detectorRotationDeg','values',unique(opts.detectorRotationProfileDeg)), ...
        struct('field','detectorAnisotropy','values',unique(opts.detectorAnisotropyProfile)), ...
        struct('field','emissionWavelengthUm','values',unique(lam)), ...
        struct('field','objectiveNA','values',unique(na)), ...
        struct('field','coverslipThicknessUm','values',unique(cover))};
    if strcmpi(sim.sampleGeometry,'airOnGlass')
        immersion=opts.immersionIndexProfile;
        if isempty(immersion)
            immersion=nominal.immersionRefractiveIndex+[-0.005 0 0.005];
        end
        glass=opts.glassIndexProfile;
        if isempty(glass)
            glass=nominal.glassRefractiveIndex+[-0.005 0 0.005];
        end
        bottom=opts.beadBottomHeightProfileUm;
        if isempty(bottom)
            bottom=nominal.beadBottomHeightUm+[0 0.01 0.02];
        end
        definitions=[definitions { ...
            struct('field','immersionRefractiveIndex','values',unique(immersion)), ...
            struct('field','glassRefractiveIndex','values',unique(glass)), ...
            struct('field','beadBottomHeightUm','values',unique(bottom))}];
    end
end

function bootstrap = runBootstrap(obs,data,sim,fit,opts)
    n=opts.nBootstrap;
    bootstrap=struct('table',table(),'coefficientStd',nan(1,numel(opts.fitModes)), ...
        'signAgreement',nan(1,numel(opts.fitModes)));
    if n<=0, return; end
    totalTimer=tic;
    coeff=zeros(n,numel(opts.fitModes));
    z0=zeros(n,1);
    for b=1:n
        obsB=obs;
        obsB.y=poissonSample(fit.expectedObservation);
        localOpts=opts;
        localOpts.maxIter=opts.bootstrapMaxIter;
        localOpts.verbose=false;
        constraints=struct('maxIter',opts.bootstrapMaxIter,'quiet',true);
        one=fitHybridPoisson(obsB,data,sim,localOpts,fit.paramVector,constraints);
        coeff(b,:)=fitModeVector(one,opts.fitModes);
        z0(b)=one.estZ0Um;
        if opts.verbose
            fprintf('[full-stack bootstrap] %d/%d (elapsed %.1f min)\n', ...
                b,n,toc(totalTimer)/60);
        end
    end
    bootstrap.coefficientSamples=coeff;
    bootstrap.z0Samples=z0;
    bootstrap.coefficientStd=std(coeff,0,1);
    estimate=fitModeVector(fit,opts.fitModes);
    bootstrap.signAgreement=mean(sign(coeff)==sign(estimate),1);
    names=opts.fitModes(:);
    bootstrap.table=array2table(coeff,'VariableNames',opts.fitModes);
    bootstrap.table.bootstrapIndex=(1:n).';
    bootstrap.table=movevars(bootstrap.table,'bootstrapIndex','Before',1);
    bootstrap.modeSummary=table(names,estimate(:), ...
        bootstrap.coefficientStd(:),bootstrap.signAgreement(:), ...
        'VariableNames',{'mode','estimateWaves','bootstrapStdWaves','signAgreement'});
end

function signProfiles = runSignProfiles(obs,data,sim,fit,opts)
    nMode=numel(opts.fitModes);
    best=fitModeVector(fit,opts.fitModes);
    oppositeDeviance=nan(nMode,1);
    residualRatio=nan(nMode,1);
    if ~opts.runSignProfiles
        signProfiles=table(opts.fitModes(:),best(:),oppositeDeviance,residualRatio, ...
            'VariableNames',{'mode','estimateWaves','oppositeSignDeviance', ...
            'signResidualRatio'});
        return;
    end
    totalTimer=tic;
    for k=1:nMode
        if abs(best(k))<1e-6
            oppositeDeviance(k)=fit.deviance;
            residualRatio(k)=1;
            continue;
        end
        start=fit.paramVector;
        start(k)=-sign(best(k))*max(abs(best(k)),0.01);
        constraints=struct('signIndex',k,'signValue',-sign(best(k)), ...
            'maxIter',opts.signProfileMaxIter,'quiet',true);
        localOpts=opts;
        localOpts.maxIter=opts.signProfileMaxIter;
        localOpts.verbose=false;
        one=fitHybridPoisson(obs,data,sim,localOpts,start,constraints);
        oppositeDeviance(k)=one.deviance;
        residualRatio(k)=sqrt(one.deviance/max(fit.deviance,eps));
        if opts.verbose
            fprintf(['[full-stack sign profile] %d/%d %s ratio %.4f ' ...
                '(elapsed %.1f min)\n'],k,nMode,opts.fitModes{k}, ...
                residualRatio(k),toc(totalTimer)/60);
        end
    end
    signProfiles=table(opts.fitModes(:),best(:),oppositeDeviance,residualRatio, ...
        'VariableNames',{'mode','estimateWaves','oppositeSignDeviance','signResidualRatio'});
end

function diagnostics = diagnoseResiduals(obs,fit,data)
    pearson=(obs.y-fit.expectedObservation)./sqrt(max(fit.expectedObservation,1));
    selected=reshape(pearson(1:obs.selectedCount),obs.selectedSize);
    correlations=[];
    for ip=1:size(selected,4)
        for ch=1:size(selected,3)
            image=selected(:,:,ch,ip);
            correlations(end+1,1)=lagCorrelation(image,1,0); %#ok<AGROW>
            correlations(end+1,1)=lagCorrelation(image,0,1); %#ok<AGROW>
        end
    end
    diagnostics=struct();
    diagnostics.pearsonVariance=mean(pearson.^2);
    diagnostics.meanSpatialLagCorrelation=mean(correlations,'omitnan');
    diagnostics.maxAbsSpatialLagCorrelation=max(abs(correlations),[],'omitnan');
    diagnostics.selectedPearsonResidual=selected;
    diagnostics.allPearsonResidual=pearson;
    modelTrace=squeeze(sum(sum(sum(fit.model,1),2),3)).';
    expectedTrace=fit.globalPhotonScale*modelTrace+data.axialBackgroundCounts;
    axialPearson=(data.axialRawCounts-expectedTrace)./sqrt(max(expectedTrace,1));
    [~,stageOrder]=sort(data.stageZUm);
    diagnostics.axialExpectedCounts=expectedTrace;
    diagnostics.axialPearsonResidual=axialPearson;
    diagnostics.axialPearsonVariance=mean(axialPearson.^2);
    diagnostics.axialLagCorrelation=lagVectorCorrelation(axialPearson(stageOrder));
end

function [passed,reason] = baseFitValidationPassed(fit,residual,opts)
    failures={};
    if ~isfield(fit,'converged') || ~fit.converged
        if isfield(fit,'terminationReason')
            failures{end+1}=sprintf('optimizer did not converge (%s)', ...
                fit.terminationReason);
        else
            failures{end+1}='optimizer did not converge';
        end
    end
    if residual.pearsonVariance<opts.pearsonVarianceRange(1) || ...
            residual.pearsonVariance>opts.pearsonVarianceRange(2)
        failures{end+1}=sprintf('spatial Pearson variance %.3g outside [%.3g %.3g]', ...
            residual.pearsonVariance,opts.pearsonVarianceRange(1), ...
            opts.pearsonVarianceRange(2));
    end
    if residual.maxAbsSpatialLagCorrelation>opts.maxSpatialResidualCorrelation
        failures{end+1}=sprintf('spatial lag correlation %.3g exceeds %.3g', ...
            residual.maxAbsSpatialLagCorrelation, ...
            opts.maxSpatialResidualCorrelation);
    end
    if residual.axialPearsonVariance<opts.axialPearsonVarianceRange(1) || ...
            residual.axialPearsonVariance>opts.axialPearsonVarianceRange(2)
        failures{end+1}=sprintf('axial Pearson variance %.3g outside [%.3g %.3g]', ...
            residual.axialPearsonVariance,opts.axialPearsonVarianceRange(1), ...
            opts.axialPearsonVarianceRange(2));
    end
    if abs(residual.axialLagCorrelation)>opts.maxAxialResidualCorrelation
        failures{end+1}=sprintf('axial lag correlation %.3g exceeds %.3g', ...
            residual.axialLagCorrelation,opts.maxAxialResidualCorrelation);
    end
    passed=isempty(failures);
    if passed
        reason='base-fit convergence and residual gates passed';
    else
        reason=strjoin(failures,'; ');
    end
end

function c = lagCorrelation(image,dy,dx)
    A=image(1:end-dy,1:end-dx);
    B=image(1+dy:end,1+dx:end);
    A=A-mean(A(:),'omitnan');
    B=B-mean(B(:),'omitnan');
    c=sum(A(:).*B(:),'omitnan')/max(sqrt(sum(A(:).^2)*sum(B(:).^2)),eps);
end

function c = lagVectorCorrelation(values)
    A=values(1:end-1);
    B=values(2:end);
    A=A-mean(A,'omitnan');
    B=B-mean(B,'omitnan');
    c=sum(A.*B,'omitnan')/max(sqrt(sum(A.^2)*sum(B.^2)),eps);
end

function acceptance = buildAcceptance(fit,profiles,bootstrap,signProfiles, ...
        residual,data,opts)
    estimate=fitModeVector(fit,opts.fitModes);
    n=numel(estimate);
    bootStd=nan(1,n);
    bootSign=nan(1,n);
    if isfield(bootstrap,'coefficientStd')
        bootStd=bootstrap.coefficientStd;
        bootSign=bootstrap.signAgreement;
    end
    profileStd=nan(1,n);
    if isfield(profiles,'coefficientSystematicHalfRange') && ...
            ~isempty(profiles.coefficientSystematicHalfRange)
        profileStd=profiles.coefficientSystematicHalfRange;
    end
    totalStd=sqrt(bootStd.^2+profileStd.^2);
    aboveUncertainty=abs(estimate)>totalStd;
    aboveModelFloor=abs(estimate)>opts.scalarModelFloorWaves;
    signGate=signProfiles.signResidualRatio.'>opts.signRatioThreshold;
    bootstrapGate=bootSign>opts.bootstrapSignThreshold;
    varianceGate=residual.pearsonVariance>=opts.pearsonVarianceRange(1) && ...
        residual.pearsonVariance<=opts.pearsonVarianceRange(2);
    correlationGate=residual.maxAbsSpatialLagCorrelation<= ...
        opts.maxSpatialResidualCorrelation;
    axialVarianceGate=residual.axialPearsonVariance>= ...
        opts.axialPearsonVarianceRange(1) && ...
        residual.axialPearsonVariance<=opts.axialPearsonVarianceRange(2);
    axialCorrelationGate=abs(residual.axialLagCorrelation)<= ...
        opts.maxAxialResidualCorrelation;
    residualGate=varianceGate && correlationGate && ...
        axialVarianceGate && axialCorrelationGate;
    flatFieldGate=data.flatFieldCalibrated;
    backgroundGate=data.backgroundIndependent && ...
        data.backgroundExposureCalibrated;
    interfaceGeometryGate=~strcmpi(fit.sim.sampleGeometry,'airOnGlass') || ...
        opts.interfaceGeometryConfirmed;
    detectorGeometryGate=isfield(fit.sim,'detectorPhysicalPitchCalibrated') && ...
        fit.sim.detectorPhysicalPitchCalibrated;
    optimizerGate=isfield(fit,'converged') && fit.converged;
    oddCalibration=knownCalibrationPassed(opts.knownAberrationCalibration);
    oddMode=ismember(opts.fitModes,{'coma_x','coma_y','tilt_x','tilt_y'});
    oddGate=true(1,n);
    oddGate(oddMode)=oddCalibration;
    accepted=aboveUncertainty & aboveModelFloor & signGate & bootstrapGate & ...
        residualGate & flatFieldGate & backgroundGate & ...
        interfaceGeometryGate & detectorGeometryGate & optimizerGate & oddGate;

    acceptance=struct();
    acceptance.modeTable=table(opts.fitModes(:),estimate(:),bootStd(:), ...
        profileStd(:),totalStd(:),bootSign(:), ...
        signProfiles.signResidualRatio,aboveUncertainty(:),aboveModelFloor(:),signGate(:), ...
        bootstrapGate(:),repmat(residualGate,n,1), ...
        repmat(flatFieldGate,n,1),repmat(backgroundGate,n,1), ...
        repmat(interfaceGeometryGate,n,1), ...
        repmat(detectorGeometryGate,n,1),repmat(optimizerGate,n,1), ...
        oddGate(:),accepted(:), ...
        'VariableNames',{'mode','estimateWaves','bootstrapStdWaves', ...
        'profileHalfRangeWaves','totalUncertaintyWaves','bootstrapSignAgreement', ...
        'signResidualRatio','aboveUncertainty','aboveScalarModelFloor', ...
        'signGate','bootstrapGate', ...
        'residualGate','flatFieldGate','calibratedBackgroundGate', ...
        'interfaceGeometryGate','detectorGeometryGate', ...
        'optimizerConvergenceGate','oddModeCalibrationGate','accepted'});
    acceptance.residualVarianceGate=varianceGate;
    acceptance.residualCorrelationGate=correlationGate;
    acceptance.axialResidualVarianceGate=axialVarianceGate;
    acceptance.axialResidualCorrelationGate=axialCorrelationGate;
    acceptance.flatFieldGate=flatFieldGate;
    acceptance.calibratedBackgroundGate=backgroundGate;
    acceptance.interfaceGeometryGate=interfaceGeometryGate;
    acceptance.detectorGeometryGate=detectorGeometryGate;
    acceptance.optimizerConvergenceGate=optimizerGate;
    acceptance.knownOddModeCalibrationGate=oddCalibration;
    acceptance.nuisanceProfilesPerformed=opts.runProfiles;
    acceptance.bootstrapPerformed=opts.nBootstrap>0;
    acceptance.signProfilesPerformed=opts.runSignProfiles;
    acceptance.allAccepted=all(accepted);
end

function passed = knownCalibrationPassed(value)
    passed=false;
    if isempty(value), return; end
    if ischar(value) || isstring(value)
        S=load(char(value));
        names=fieldnames(S);
        for k=1:numel(names)
            candidate=S.(names{k});
            if isstruct(candidate)
                passed=comaCalibrationPassed(candidate);
            end
            if passed
                return;
            end
        end
    elseif isstruct(value)
        passed=comaCalibrationPassed(value);
    end
end

function passed = comaCalibrationPassed(report)
    passed=false;
    if ~isfield(report,'table') || ~istable(report.table)
        return;
    end
    T=report.table;
    if ~all(ismember({'mode','signAgreement'},T.Properties.VariableNames))
        return;
    end
    modes=cellstr(string(T.mode));
    isComa=ismember(modes,{'coma_x','coma_y'});
    if ismember('commandWaves',T.Properties.VariableNames)
        isComa=isComa & isfinite(T.commandWaves) & abs(T.commandWaves)>0;
    end
    passed=any(isComa) && all(logical(T.signAgreement(isComa)));
end

function vector = fitModeVector(fit,modes)
    vector=zeros(1,numel(modes));
    for k=1:numel(modes)
        if isfield(fit.estCoeffs,modes{k}), vector(k)=fit.estCoeffs.(modes{k}); end
    end
end

function outDir = resolveOutputDir(opts)
    if ~isempty(opts.outputDir)
        outDir=opts.outputDir;
    else
        root=fileparts(fileparts(mfilename('fullpath')));
        outDir=fullfile(root,'output_matlab','full_stack_ism_wavefront');
    end
end

function writeOutputs(result)
    outDir=result.outputDir;
    if exist(outDir,'dir')~=7, mkdir(outDir); end
    save(fullfile(outDir,'full_stack_ism_wavefront_fit.mat'),'result','-v7.3');
    writetable(result.selection,fullfile(outDir,'selected_symmetric_planes.csv'));
    writetable(coefficientTable(result),fullfile(outDir,'full_stack_coefficients.csv'));
    writetable(result.signProfiles,fullfile(outDir,'sign_profile_likelihood.csv'));
    writetable(result.acceptance.modeTable,fullfile(outDir,'coefficient_acceptance_gates.csv'));
    if ~isempty(result.profiles.table)
        writetable(result.profiles.table,fullfile(outDir,'nuisance_profiles.csv'));
    end
    if ~isempty(result.bootstrap.table)
        writetable(result.bootstrap.table,fullfile(outDir,'bootstrap_coefficients.csv'));
        writetable(result.bootstrap.modeSummary,fullfile(outDir,'bootstrap_summary.csv'));
    end
    writetable(axialTraceTable(result),fullfile(outDir,'axial_trace_fit.csv'));
    darkTable=darkBackgroundTable(result.data);
    if ~isempty(darkTable)
        writetable(darkTable,fullfile(outDir,'dark_background_per_detector.csv'));
    end
    writeSummaryFigure(result,fullfile(outDir,'full_stack_fit_summary.png'));
    try
        plotFullStackISMFitDiagnostics(result,'outputDir',outDir);
    catch diagnosticError
        warning('estimateFullStackISMWavefront:DiagnosticFigures', ...
            'Detector-fit diagnostic figures failed: %s',diagnosticError.message);
    end
end

function T = darkBackgroundTable(data)
    T=table();
    if ~isfield(data,'darkDiagnostics') || ...
            ~isstruct(data.darkDiagnostics) || ...
            ~isfield(data.darkDiagnostics,'expectedCountsPerPixel')
        return;
    end
    D=data.darkDiagnostics;
    n=numel(D.expectedCountsPerPixel);
    channelID=data.channelIDs(:);
    if numel(channelID)~=n, channelID=(1:n).'; end
    totalPhotonCounts=diagnosticColumn(D,'totalPhotonCounts',n);
    countRateHz=diagnosticColumn(D,'countRateHz',n);
    expectedCountsPerPixel=diagnosticColumn(D,'expectedCountsPerPixel',n);
    appliedCountsPerPixel=diagnosticColumn(D,'appliedCountsPerPixel',n);
    T=table(channelID,totalPhotonCounts,countRateHz, ...
        expectedCountsPerPixel,appliedCountsPerPixel);
end

function values = diagnosticColumn(S,name,n)
    values=nan(n,1);
    if isfield(S,name) && numel(S.(name))==n
        values=double(S.(name)(:));
    end
end

function T = coefficientTable(result)
    modes=result.options.fitModes(:);
    estimate=fitModeVector(result.fit,result.options.fitModes).';
    nmRms=estimate*result.fit.sim.lamRef*1000;
    T=table(modes,estimate,nmRms);
end

function T = axialTraceTable(result)
    data=result.data;
    model=result.fit.model;
    unit=squeeze(sum(sum(sum(model,1),2),3)).';
    background=data.axialBackgroundCounts;
    expected=result.fit.globalPhotonScale*unit+background;
    pearsonResidual=(data.axialRawCounts-expected)./sqrt(max(expected,1));
    selected=false(numel(data.stageZUm),1);
    selected(result.selectedPlaneIndices)=true;
    T=table((1:numel(data.stageZUm)).',data.stageZUm(:), ...
        data.axialRawCounts(:),background(:),expected(:), ...
        pearsonResidual(:),selected, ...
        'VariableNames',{'planeIndex','stageZUm','rawCounts', ...
        'measuredBackgroundCounts','expectedCounts','pearsonResidual', ...
        'usedAsFullImage'});
end

function writeSummaryFigure(result,outFile)
    fig=figure('Visible','off','Color','w','Position',[60 60 1300 850]);
    tl=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
    T=axialTraceTable(result);
    ax=nexttile(tl);
    plot(ax,T.stageZUm,T.rawCounts,'ko-','LineWidth',1); hold(ax,'on');
    plot(ax,T.stageZUm,T.expectedCounts,'r-','LineWidth',1.5);
    plot(ax,T.stageZUm(T.usedAsFullImage),T.rawCounts(T.usedAsFullImage), ...
        'bo','MarkerFaceColor','b');
    xline(ax,result.fit.estZ0Um,'--','fitted z0');
    xlabel(ax,'recorded stage z [um]'); ylabel(ax,'raw photons');
    legend(ax,{'measured','fitted','full-image planes'},'Location','best');
    grid(ax,'on'); title(ax,'Unnormalized axial photon trace');

    ax=nexttile(tl);
    A=result.acceptance.modeTable;
    x=1:height(A);
    plot(ax,x,A.estimateWaves,'o','LineWidth',1.2);
    hold(ax,'on');
    finiteUncertainty=isfinite(A.totalUncertaintyWaves);
    if any(finiteUncertainty)
        errorbar(ax,x(finiteUncertainty),A.estimateWaves(finiteUncertainty), ...
            A.totalUncertaintyWaves(finiteUncertainty), ...
            'LineStyle','none','LineWidth',1.2);
    end
    yline(ax,0,'k-');
    set(ax,'XTick',1:height(A),'XTickLabel',A.mode,'XTickLabelRotation',45);
    ylabel(ax,'waves RMS'); grid(ax,'on'); title(ax,'Coefficients and uncertainty');

    ax=nexttile(tl);
    selected=result.residualDiagnostics.selectedPearsonResidual;
    residualMap=squeeze(mean(mean(selected,3),4));
    imagesc(ax,residualMap); axis(ax,'image'); colorbar(ax);
    title(ax,sprintf('Mean Pearson residual; max lag corr %.3f', ...
        result.residualDiagnostics.maxAbsSpatialLagCorrelation));

    ax=nexttile(tl);
    gateNames={'uncertainty','model floor','sign','bootstrap','residual', ...
        'flat field','background','interface','detector','odd calibration', ...
        'optimizer','accepted'};
    gateValues=double([A.aboveUncertainty A.aboveScalarModelFloor ...
        A.signGate A.bootstrapGate A.residualGate A.flatFieldGate ...
        A.calibratedBackgroundGate A.interfaceGeometryGate ...
        A.detectorGeometryGate A.oddModeCalibrationGate ...
        A.optimizerConvergenceGate A.accepted]);
    imagesc(ax,gateValues,[0 1]);
    set(ax,'YTick',1:height(A),'YTickLabel',A.mode, ...
        'XTick',1:numel(gateNames),'XTickLabel',gateNames, ...
        'XTickLabelRotation',45);
    colormap(ax,[0.75 0.08 0.08;0.08 0.55 0.18]);
    colorbar(ax,'Ticks',[0 1],'TickLabels',{'fail','pass'});
    title(ax,'Acceptance gates');
    if result.acceptance.allAccepted
        status='FIT ACCEPTED';
        titleColor=[0 0.45 0];
    else
        status='FIT REJECTED';
        titleColor=[0.75 0 0];
    end
    title(tl,sprintf('%s: Hybrid full-stack Poisson fit; z0 = %.4f um; deviance %.4g', ...
        status,result.fit.estZ0Um,result.fit.deviance), ...
        'FontWeight','bold','Color',titleColor);
    exportgraphics(fig,outFile,'Resolution',180);
    close(fig);
end

function printSummary(result)
    fprintf('\nHybrid full-stack ISM Poisson fit\n');
    fprintf('  sample geometry: %s; forward model: %s\n', ...
        result.fit.sim.sampleGeometry,result.fit.sim.diffractionModel);
    if strcmpi(result.fit.sim.sampleGeometry,'airOnGlass')
        fprintf(['  interface: water %.4g / glass %.4g (%.1f um) / ' ...
            'air %.4g; geometry confirmed %d\n'], ...
            result.fit.sim.nImmersion,result.fit.sim.nGlass, ...
            result.fit.sim.coverslipThicknessUm,result.fit.sim.nSample, ...
            result.acceptance.interfaceGeometryGate);
    end
    fprintf('  calibrated detector pitch: %d', ...
        result.acceptance.detectorGeometryGate);
    if result.acceptance.detectorGeometryGate
        fprintf(' (%.4g um in sample space)', ...
            result.fit.sim.detectorCalibratedPitchUm);
    end
    fprintf('\n');
    if isfield(result.geometry, 'diagnostics') && ...
            isfield(result.geometry.diagnostics, 'source')
        fprintf('  detector geometry source: %s\n', ...
            result.geometry.diagnostics.source);
    end
    fprintf('  fitted optical-axis stage coordinate z0: %.4f um\n',result.fit.estZ0Um);
    fprintf('  global photon scale: %.6g\n',result.fit.globalPhotonScale);
    fprintf('  Poisson deviance: %.6g\n',result.fit.deviance);
    fprintf('  optimizer converged: %d (%s)\n', ...
        result.fit.converged,result.fit.terminationReason);
    fprintf('  spatial Pearson variance: %.3f; max spatial lag correlation: %.3f\n', ...
        result.residualDiagnostics.pearsonVariance, ...
        result.residualDiagnostics.maxAbsSpatialLagCorrelation);
    fprintf('  axial Pearson variance: %.3f; axial lag correlation: %.3f\n', ...
        result.residualDiagnostics.axialPearsonVariance, ...
        result.residualDiagnostics.axialLagCorrelation);
    if ~result.acceptance.nuisanceProfilesPerformed || ...
            ~result.acceptance.bootstrapPerformed || ...
            ~result.acceptance.signProfilesPerformed
        fprintf(['  validation incomplete: runValidatedFullStackISMWavefront ' ...
            'before accepting coefficients\n']);
    end
    disp(result.acceptance.modeTable);
    fprintf('  outputs: %s\n\n',result.outputDir);
end

function pitch = medianNearestDistance(xy)
    n=size(xy,1);
    nearest=nan(n,1);
    for k=1:n
        d=hypot(xy(:,1)-xy(k,1),xy(:,2)-xy(k,2));
        d(k)=inf;
        nearest(k)=min(d);
    end
    pitch=median(nearest(isfinite(nearest)&nearest>0));
end

function detXY = validateDetectorXYFullStack(detXY, nCh)
    detXY = double(detXY);
    if ~ismatrix(detXY) || size(detXY,2) ~= 2 || size(detXY,1) ~= nCh || ...
            any(~isfinite(detXY(:)))
        error('estimateFullStackISMWavefront:BadDetectorXY', ...
            'detectorXYUm must be a finite [%d x 2] array in sample-equivalent um.', nCh);
    end
end

function diagnostics = fixedDetectorDiagnosticsFullStack(detXY, pixelSizeUm, airyUnitUm, source)
    diagnostics = struct();
    diagnostics.source = char(source);
    diagnostics.pixelSizeUm = pixelSizeUm;
    diagnostics.pixelSizeNm = 1000 * pixelSizeUm;
    diagnostics.airyUnitUm = airyUnitUm;
    diagnostics.positionSign = NaN;
    diagnostics.detectorScale = NaN;
    diagnostics.centerMode = 'fixed';
    [~, diagnostics.centerDetectorIndex] = min(sum(detXY.^2, 2));
    diagnostics.signal = nan(size(detXY,1), 1);
    diagnostics.peakValue = nan(size(detXY,1), 1);
    diagnostics.shiftsPx = nan(size(detXY));
    diagnostics.shiftsUm = nan(size(detXY));
    diagnostics.shiftsNm = nan(size(detXY));
    diagnostics.detectorOffsetFromShiftUm = detXY;
    diagnostics.detectorOffsetFromShiftNm = 1000 * detXY;
    diagnostics.detXYUm = detXY;
    diagnostics.detXYNm = 1000 * detXY;
    diagnostics.shiftsAU = nan(size(detXY));
    diagnostics.shiftMagnitudeAU = nan(size(detXY,1), 1);
    diagnostics.detectorOffsetFromShiftAU = detXY / airyUnitUm;
    diagnostics.detXYAU = detXY / airyUnitUm;
    diagnostics.detectorRadiusAU = hypot( ...
        diagnostics.detXYAU(:,1), diagnostics.detXYAU(:,2));
end

function value = headerNumericByPattern(head,patterns)
    value=NaN;
    if ~isstruct(head), return; end
    names=fieldnames(head);
    for p=1:numel(patterns)
        for k=1:numel(names)
            if ~isempty(regexpi(names{k},patterns{p},'once'))
                raw=head.(names{k});
                if isnumeric(raw) && ~isempty(raw)
                    candidate=double(raw(1));
                    if isfinite(candidate)
                        value=candidate;
                        return;
                    end
                end
            end
        end
    end
end

function tf = isfinitePositive(value)
    tf=isnumeric(value)&&isscalar(value)&&isfinite(value)&&value>0;
end

function value = effectivePropagatingNA(sim)
    value = sim.NA;
    if isfield(sim,'sampleGeometry') && strcmpi(sim.sampleGeometry,'airOnGlass')
        value = min(sim.NA,sim.nSample);
    end
end

function pitch = resolveDetectorPitchSampleUm(opts)
    pitch = NaN;
    if isnumeric(opts.detectorPitchSampleUm) && ...
            isscalar(opts.detectorPitchSampleUm) && ...
            isfinite(opts.detectorPitchSampleUm) && ...
            opts.detectorPitchSampleUm > 0
        pitch = double(opts.detectorPitchSampleUm);
        return;
    end
    physical = opts.detectorHardwarePitchUm;
    magnification = opts.detectorTotalMagnification;
    if isnumeric(physical) && isscalar(physical) && isfinite(physical) && ...
            physical > 0 && isnumeric(magnification) && ...
            isscalar(magnification) && isfinite(magnification) && ...
            magnification > 0
        pitch = double(physical)/double(magnification);
    end
end

function source = detectorPitchSourceFullStack(opts)
    if isnumeric(opts.detectorPitchSampleUm) && ...
            isscalar(opts.detectorPitchSampleUm) && ...
            isfinite(opts.detectorPitchSampleUm) && ...
            opts.detectorPitchSampleUm > 0
        source = 'detectorPitchSampleUm';
    elseif isnumeric(opts.detectorHardwarePitchUm) && ...
            isscalar(opts.detectorHardwarePitchUm) && ...
            isfinite(opts.detectorHardwarePitchUm) && ...
            opts.detectorHardwarePitchUm > 0 && ...
            isnumeric(opts.detectorTotalMagnification) && ...
            isscalar(opts.detectorTotalMagnification) && ...
            isfinite(opts.detectorTotalMagnification) && ...
            opts.detectorTotalMagnification > 0
        source = 'detectorHardwarePitchUm / detectorTotalMagnification';
    else
        source = 'defaultParams.detPitch';
    end
end
