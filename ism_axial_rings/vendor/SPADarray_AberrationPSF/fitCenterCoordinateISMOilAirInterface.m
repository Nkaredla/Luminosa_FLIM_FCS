function result = fitCenterCoordinateISMOilAirInterface(stackInput, varargin)
%FITCENTERCOORDINATEISMOILAIRINTERFACE Center-coordinate ISM fit for oil/air data.
%
%   result = fitCenterCoordinateISMOilAirInterface(stackInput)
%
%   Workflow:
%     1) Load the full detector-resolved ISM z stack.
%     2) Find the bead focus plane and center from the detector-summed image
%        using brightness and spatial second moment.
%     3) Fit the detector micro-image sampled at that center using an
%        oil/glass/air-interface scalar model.
%     4) Optionally refine the net wavefront with the multi-plane full-stack
%        Poisson fitter.
%     5) Generate air-interface model volumes and per-detector XY/XZ/YZ
%        measured/model/residual projection figures.
%
%   Defaults deliberately assume the recorded planes are already aligned in
%   scan coordinates. The workflow rereads raw PTU frames by default
%   ('reuseCache', false), ignores saved inter-plane alignment shifts
%   ('useInterPlaneAlignment', false), and does not recenter the measured
%   comparison volume ('recenterMeasuredVolume', false). Set those flags
%   explicitly if you want to reproduce an older alignment-corrected run.
%
%   Single-plane mode uses only the focal center-coordinate detector
%   distribution. This is useful for oil/glass/air data because the interface
%   model is axially asymmetric, so even-mode signs can be partially encoded in
%   one plane. Check result.centerFit.sufficiency before trusting a single-plane
%   estimate.
%
%   By default the bead center is measured from the detector-summed scan image,
%   then held fixed in the model. Zernike defocus is excluded and the PDA23
%   detector coordinates are a regular honeycomb, optionally calibrated by a
%   supplied detector pitch or detectorHardwarePitchUm/detectorTotalMagnification.
%   The old image-derived detector-layout estimate is available only by setting
%   estimateDetectorLayout true.
%
%   Multi-plane mode uses estimateFullStackISMWavefront initialized from the
%   single-plane center-coordinate fit. This is the recommended branch for the
%   final net wavefront estimate.
%
%   Example:
%       res = fitCenterCoordinateISMOilAirInterface( ...
%           'D:\Luminosa\Data\PSF_batch_outputs\ISM_Aberation2_73\xz_yz_plots\x0_4uW_0_19collar_80mmlens_20260515_155744_frame_alignment.csv', ...
%           'useMultiPlane', true, ...
%           'outputDir', 'D:\temp\oil_air_center_fit');

    if nargin < 1 || isempty(stackInput)
        stackInput = ['D:\Luminosa\Data\PSF_batch_outputs\ISM_Aberation2_73\' ...
            'xz_yz_plots\x0_4uW_0_19collar_80mmlens_20260515_155744_frame_alignment.csv'];
    end

    opts = parseOptions(varargin{:});
    addRequiredPaths(opts);

    if isempty(opts.outputDir)
        opts.outputDir = defaultOutputDir();
    end
    if opts.writeOutputs && exist(opts.outputDir, 'dir') ~= 7
        mkdir(opts.outputDir);
    end

    centerOutputDir = '';
    if opts.writeOutputs
        centerOutputDir = fullfile(opts.outputDir, 'center_selection');
    end

    center = findBeadCenterFromFullStack(stackInput, ...
        'channelIDs', opts.channelIDs, ...
        'ptuReaderFolder', opts.ptuReaderFolder, ...
        'darkFile', opts.darkFile, ...
        'darkScale', opts.darkScale, ...
        'backgroundMode', opts.backgroundMode, ...
        'xyPixelSizeUm', opts.xyPixelSizeUm, ...
        'reuseCache', opts.reuseCache, ...
        'cacheFile', opts.cacheFile, ...
        'subtractBackground', true, ...
        'thresholdFraction', opts.centerThresholdFraction, ...
        'makeFigure', opts.writeOutputs, ...
        'figureVisible', opts.figureVisible, ...
        'outputDir', centerOutputDir, ...
        'verbose', opts.verbose);

    data = center.data;
    focusIdx = center.focusIndex;
    focusZUm = center.focusZUm;

    centerFit = runSinglePlaneCenterFit(data, center, opts);
    fitForComparison = centerFit.fit;
    fitKind = 'singlePlaneCenterCoordinate';
    multiPlane = [];

    if opts.useMultiPlane
        multiPlane = runMultiPlaneFullStackFit(stackInput, data, center, centerFit, opts);
        fitForComparison = multiPlane.fit;
        fitKind = 'multiPlaneFullStack';
    end

    comparison = [];
    if opts.generate3DComparison
        comparison = generateAirInterface3DComparison( ...
            data, fitForComparison, focusIdx, focusZUm, fitKind, opts);
    end

    result = struct();
    result.stackInput = stackInput;
    result.options = opts;
    result.centerSelection = center;
    result.centerFit = centerFit;
    result.multiPlaneFit = multiPlane;
    result.fitUsedForComparison = fitKind;
    result.comparison3D = comparison;
    result.outputDir = opts.outputDir;
    result.outputPrefix = opts.outputPrefix;
    result.reportTitle = opts.reportTitle;
    result.centerDetectorIntensity = centerDetectorIntensityTable(result);

    if opts.writeOutputs
        save(fullfile(opts.outputDir, [opts.outputPrefix '_ism_fit.mat']), ...
            'result', '-v7.3');
        writeCoefficientCsv(result, fullfile(opts.outputDir, ...
            [opts.outputPrefix '_coefficients.csv']));
        writeCenterDetectorIntensityOutputs(result, opts.outputDir);
    end

    if opts.verbose
        printSummary(result);
    end
end

function opts = parseOptions(varargin)
    p = inputParser;
    p.FunctionName = 'fitCenterCoordinateISMOilAirInterface';

    addParameter(p, 'useMultiPlane', false);
    addParameter(p, 'generate3DComparison', true);
    addParameter(p, 'selectedOffsetsUm', [-0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8]);
    addParameter(p, 'syntheticDiversityUm', 0.5);
    addParameter(p, 'outputPrefix', 'oil_air_center_coordinate');
    addParameter(p, 'reportTitle', 'Oil/glass/air center-coordinate ISM fit');

    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'ptuReaderFolder', defaultPtuReaderFolder());
    addParameter(p, 'darkFile', 'D:\Luminosa\Data\ISMdark_counts.ptu');
    addParameter(p, 'darkScale', 1);
    addParameter(p, 'backgroundMode', 'auto');
    addParameter(p, 'xyPixelSizeUm', []);
    addParameter(p, 'cacheFile', '');
    addParameter(p, 'reuseCache', false);
    addParameter(p, 'useInterPlaneAlignment', false);

    addParameter(p, 'objectiveNA', 1.45);
    addParameter(p, 'objectiveMagnification', 100);
    addParameter(p, 'emissionWavelengthUm', []);
    addParameter(p, 'excitationWavelengthUm', []);
    addParameter(p, 'immersionRefractiveIndex', 1.515);
    addParameter(p, 'glassRefractiveIndex', 1.518);
    addParameter(p, 'sampleRefractiveIndex', 1.0003);
    addParameter(p, 'interfaceRadialWeightMode', 'sampleSolidAngle');
    addParameter(p, 'diffractionModel', 'vectorial');
    addParameter(p, 'vectorialPolarizationMode', 'xyAverage');
    addParameter(p, 'excitationPolarizationMode', 'x');
    addParameter(p, 'collectionPolarizationMode', 'xyAverage');
    addParameter(p, 'airInterfaceStageMedium', 'immersion');
    addParameter(p, 'designGlassRefractiveIndex', 1.518);
    addParameter(p, 'coverslipThicknessUm', 170);
    addParameter(p, 'designCoverslipThicknessUm', 170);
    addParameter(p, 'beadBottomHeightUm', 0);
    addParameter(p, 'beadDiameterUm', []);
    addParameter(p, 'beadRadiusUm', 0.08);
    addParameter(p, 'airBeadAxialSamples', 5);
    addParameter(p, 'beadSubsamples', [3 3 3]);

    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detectorPixelShape', 'hex');
    addParameter(p, 'detectorSubsamples', 5);
    addParameter(p, 'estimateDetectorLayout', false);
    addParameter(p, 'detectorXYUm', []);
    addParameter(p, 'detectorLayoutScale', 2);
    addParameter(p, 'detectorLayoutPositionSign', -1);
    addParameter(p, 'detectorPitchSampleUm', []);
    addParameter(p, 'detectorHardwarePitchUm', []);
    addParameter(p, 'detectorTotalMagnification', []);
    addParameter(p, 'detectorQE', []);

    addParameter(p, 'fitModes', ...
        {'astig_x','astig_y','coma_x','coma_y','spherical'});
    addParameter(p, 'centerFitMaxIter', 8);
    addParameter(p, 'fullStackMaxIter', 4);
    addParameter(p, 'fitXY', false);
    addParameter(p, 'fitZ', false);
    addParameter(p, 'fitDetectorPitchScale', false);
    addParameter(p, 'initialDetectorPitchScale', 1);
    addParameter(p, 'detectorPitchScaleBounds', [0.7 1.3]);
    addParameter(p, 'centerThresholdFraction', 0.10);
    addParameter(p, 'centerNormalization', 'perPlane');
    addParameter(p, 'centerSampleMode', 'subpixel');
    addParameter(p, 'modelSampleXY', [0 0]);
    addParameter(p, 'initialCoeffs', struct());
    addParameter(p, 'runAirInterfaceSelfTest', true);
    addParameter(p, 'interfaceGeometryConfirmed', false);

    addParameter(p, 'applySavedAlignment', false);
    addParameter(p, 'alignmentShiftMinPx', 4);
    addParameter(p, 'recenterMeasuredVolume', false);
    addParameter(p, 'projectionMode', 'max');
    addParameter(p, 'supportThreshold', 0.01);
    addParameter(p, 'modelBlurSigmaXYUm', 0);
    addParameter(p, 'modelBlurSigmaZUm', 0);
    addParameter(p, 'writePerDetectorFigures', true);
    addParameter(p, 'writeComparisonVolumeMat', true);

    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'figureVisible', 'off');
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;
    opts.channelIDs = double(opts.channelIDs(:)).';
    opts.ptuReaderFolder = char(opts.ptuReaderFolder);
    opts.cacheFile = char(opts.cacheFile);
    opts.darkFile = char(opts.darkFile);
    opts.backgroundMode = lower(char(opts.backgroundMode));
    opts.interfaceRadialWeightMode = char(opts.interfaceRadialWeightMode);
    opts.diffractionModel = char(opts.diffractionModel);
    opts.vectorialPolarizationMode = char(opts.vectorialPolarizationMode);
    opts.excitationPolarizationMode = char(opts.excitationPolarizationMode);
    opts.collectionPolarizationMode = char(opts.collectionPolarizationMode);
    opts.airInterfaceStageMedium = char(opts.airInterfaceStageMedium);
    opts.detectorLayout = char(opts.detectorLayout);
    opts.detectorPixelShape = char(opts.detectorPixelShape);
    opts.centerNormalization = lower(char(opts.centerNormalization));
    opts.centerSampleMode = lower(char(opts.centerSampleMode));
    opts.projectionMode = lower(char(opts.projectionMode));
    opts.figureVisible = char(opts.figureVisible);
    opts.outputDir = char(opts.outputDir);
    opts.outputPrefix = matlab.lang.makeValidName(char(opts.outputPrefix));
    opts.reportTitle = char(opts.reportTitle);
    if ~isempty(opts.beadDiameterUm)
        opts.beadDiameterUm = validatePositiveScalarLocal( ...
            opts.beadDiameterUm, 'beadDiameterUm');
        opts.beadRadiusUm = opts.beadDiameterUm / 2;
    end
    opts.beadRadiusUm = validatePositiveScalarLocal( ...
        opts.beadRadiusUm, 'beadRadiusUm');
    opts.modelBlurSigmaXYUm = validateNonnegativeScalarLocal( ...
        opts.modelBlurSigmaXYUm, 'modelBlurSigmaXYUm');
    opts.modelBlurSigmaZUm = validateNonnegativeScalarLocal( ...
        opts.modelBlurSigmaZUm, 'modelBlurSigmaZUm');
    opts.fitDetectorPitchScale = logical(opts.fitDetectorPitchScale);
    opts.detectorQE = validateDetectorQELocal(opts.detectorQE, numel(opts.channelIDs));
    opts.initialDetectorPitchScale = validatePositiveScalarLocal( ...
        opts.initialDetectorPitchScale, 'initialDetectorPitchScale');
    opts.detectorPitchScaleBounds = validateDetectorPitchScaleBoundsLocal( ...
        opts.detectorPitchScaleBounds);
    opts.alignmentShiftMinPx = double(opts.alignmentShiftMinPx);
    if ~isscalar(opts.alignmentShiftMinPx) || ...
            ~isfinite(opts.alignmentShiftMinPx) || opts.alignmentShiftMinPx < 0
        error('fitCenterCoordinateISMOilAirInterface:BadAlignmentThreshold', ...
            'alignmentShiftMinPx must be a finite nonnegative scalar.');
    end

    if ischar(opts.fitModes) || isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    end
    opts.fitModes = opts.fitModes(:).';
    opts.selectedOffsetsUm = double(opts.selectedOffsetsUm(:)).';
    if ~ismember(opts.projectionMode, {'max','sum'})
        error('fitCenterCoordinateISMOilAirInterface:BadProjectionMode', ...
            'projectionMode must be ''max'' or ''sum''.');
    end
end

function value = validatePositiveScalarLocal(value, name)
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
        error('fitCenterCoordinateISMOilAirInterface:BadPositiveScalar', ...
            '%s must be a finite positive scalar.', name);
    end
    value = double(value);
end

function value = validateNonnegativeScalarLocal(value, name)
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value < 0
        error('fitCenterCoordinateISMOilAirInterface:BadNonnegativeScalar', ...
            '%s must be a finite nonnegative scalar.', name);
    end
    value = double(value);
end

function qe = validateDetectorQELocal(qe, nCh)
    if isempty(qe)
        return;
    end
    qe = double(qe(:));
    if numel(qe) ~= nCh || any(~isfinite(qe)) || any(qe <= 0)
        error('fitCenterCoordinateISMOilAirInterface:BadDetectorQE', ...
            'detectorQE must contain %d finite positive relative gains.', nCh);
    end
    qe = qe / mean(qe);
end

function bounds = validateDetectorPitchScaleBoundsLocal(bounds)
    if ~isnumeric(bounds) || numel(bounds) ~= 2 || ...
            any(~isfinite(bounds(:))) || any(bounds(:) <= 0)
        error('fitCenterCoordinateISMOilAirInterface:BadDetectorPitchScaleBounds', ...
            'detectorPitchScaleBounds must be two finite positive values.');
    end
    bounds = sort(double(bounds(:))).';
    if bounds(1) == bounds(2)
        error('fitCenterCoordinateISMOilAirInterface:BadDetectorPitchScaleBounds', ...
            'detectorPitchScaleBounds must span a nonzero range.');
    end
end

function addRequiredPaths(opts)
    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    parentDir = fileparts(thisDir);
    if exist(parentDir, 'dir') == 7
        addpath(parentDir);
    end
    if isfield(opts, 'ptuReaderFolder') && ~isempty(opts.ptuReaderFolder) && ...
            exist(opts.ptuReaderFolder, 'dir') == 7
        addpath(opts.ptuReaderFolder);
    end
end

function folder = defaultPtuReaderFolder()
    thisDir = fileparts(mfilename('fullpath'));
    luminosaRoot = fileparts(fileparts(thisDir));
    folder = fullfile(luminosaRoot, 'Luminosa_FLIM_FCS');
end

function outDir = defaultOutputDir()
    root = fileparts(fileparts(mfilename('fullpath')));
    stamp = datestr(now, 'yyyymmdd_HHMMSS');
    outDir = fullfile(root, 'output_matlab', ...
        'oil_air_center_coordinate_ism', stamp);
end

function centerFit = runSinglePlaneCenterFit(data, center, opts)
    focusIdx = center.focusIndex;
    focusPlane = correctedPlane(data, focusIdx);
    raw4 = cat(4, focusPlane, focusPlane);

    sim0 = defaultParams();
    sim0.beadRadius = opts.beadRadiusUm;
    sim0.beadSubsamples = opts.beadSubsamples;
    sim0.airBeadAxialSamples = opts.airBeadAxialSamples;
    sim0.detectorLayout = opts.detectorLayout;
    sim0.detectorPixelShape = opts.detectorPixelShape;
    [sim0.detXY, sim0.detectorIndexGrid, sim0.detectorLayoutInfo] = ...
        detectorLayout(sim0.detectorLayout, sim0.detPitch);
    sim0.nDet = size(sim0.detXY, 1);

    centerFit = estimateCenterPointISMWavefront(raw4, [], ...
        'planeZ', [0 opts.syntheticDiversityUm], ...
        'planeWeights', [1 0], ...
        'fitStrategy', 'joint', ...
        'channelIDs', data.channelIDs, ...
        'ptuReaderFolder', opts.ptuReaderFolder, ...
        'darkMode', 'none', ...
        'subtractBoundary', false, ...
        'centerXY', center.centerXY, ...
        'centerNormalization', opts.centerNormalization, ...
        'centerSampleMode', opts.centerSampleMode, ...
        'modelSampleXY', opts.modelSampleXY, ...
        'objectiveNA', opts.objectiveNA, ...
        'objectiveMagnification', opts.objectiveMagnification, ...
        'emissionWavelengthUm', opts.emissionWavelengthUm, ...
        'excitationWavelengthUm', opts.excitationWavelengthUm, ...
        'sampleGeometry', 'airOnGlass', ...
        'interfaceRadialWeightMode', opts.interfaceRadialWeightMode, ...
        'diffractionModel', opts.diffractionModel, ...
        'vectorialPolarizationMode', opts.vectorialPolarizationMode, ...
        'excitationPolarizationMode', opts.excitationPolarizationMode, ...
        'collectionPolarizationMode', opts.collectionPolarizationMode, ...
        'airInterfaceStageMedium', opts.airInterfaceStageMedium, ...
        'immersionRefractiveIndex', opts.immersionRefractiveIndex, ...
        'glassRefractiveIndex', opts.glassRefractiveIndex, ...
        'sampleRefractiveIndex', opts.sampleRefractiveIndex, ...
        'designGlassRefractiveIndex', opts.designGlassRefractiveIndex, ...
        'coverslipThicknessUm', opts.coverslipThicknessUm, ...
        'designCoverslipThicknessUm', opts.designCoverslipThicknessUm, ...
        'beadBottomHeightUm', opts.beadBottomHeightUm, ...
        'airBeadAxialSamples', opts.airBeadAxialSamples, ...
        'beadSubsamples', opts.beadSubsamples, ...
        'sim', sim0, ...
        'xyPixelSizeUm', data.xyPixelSizeUm, ...
        'detectorLayout', opts.detectorLayout, ...
        'detectorPixelShape', opts.detectorPixelShape, ...
        'detectorSubsamples', opts.detectorSubsamples, ...
        'estimateDetectorLayout', opts.estimateDetectorLayout, ...
        'detectorXYUm', opts.detectorXYUm, ...
        'detectorLayoutScale', opts.detectorLayoutScale, ...
        'detectorLayoutPositionSign', opts.detectorLayoutPositionSign, ...
        'detectorPitchSampleUm', opts.detectorPitchSampleUm, ...
        'detectorHardwarePitchUm', opts.detectorHardwarePitchUm, ...
        'detectorTotalMagnification', opts.detectorTotalMagnification, ...
        'detectorQE', opts.detectorQE, ...
        'fitModes', opts.fitModes, ...
        'fitXY', opts.fitXY, ...
        'fitZ', opts.fitZ, ...
        'fitDetectorPitchScale', opts.fitDetectorPitchScale, ...
        'initialDetectorPitchScale', opts.initialDetectorPitchScale, ...
        'detectorPitchScaleBounds', opts.detectorPitchScaleBounds, ...
        'initialCoeffs', opts.initialCoeffs, ...
        'maxIter', opts.centerFitMaxIter, ...
        'writeOutputs', false, ...
        'verbose', opts.verbose);

    centerFit.stackFocusIndex = focusIdx;
    centerFit.stackFocusZUm = center.focusZUm;
    centerFit.singlePlaneOnly = true;
    centerFit.singlePlaneAssumption = ['Only the focal detector micro-image ' ...
        'was weighted. The bead center is taken from the detector-summed scan ' ...
        'image, fitted lateral offset and Zernike defocus are disabled by ' ...
        'default, and detector coordinates are a fixed regular PDA23 honeycomb ' ...
        'unless explicitly overridden. Inspect centerFit.sufficiency.'];
end

function fullResult = runMultiPlaneFullStackFit(stackInput, data, center, centerFit, opts)
    outDir = fullfile(opts.outputDir, 'multi_plane_full_stack');
    if opts.writeOutputs && exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    fullResult = estimateFullStackISMWavefront(stackInput, ...
        'stageZUm', data.stageZUm, ...
        'xyPixelSizeUm', data.xyPixelSizeUm, ...
        'channelIDs', data.channelIDs, ...
        'ptuReaderFolder', opts.ptuReaderFolder, ...
        'darkFile', opts.darkFile, ...
        'darkScale', opts.darkScale, ...
        'backgroundMode', opts.backgroundMode, ...
        'cacheFile', opts.cacheFile, ...
        'reuseCache', opts.reuseCache, ...
        'useInterPlaneAlignment', opts.useInterPlaneAlignment, ...
        'interPlaneAlignmentMinPx', opts.alignmentShiftMinPx, ...
        'selectedOffsetsUm', opts.selectedOffsetsUm, ...
        'selectionCenterUm', center.focusZUm, ...
        'objectiveNA', opts.objectiveNA, ...
        'objectiveMagnification', opts.objectiveMagnification, ...
        'emissionWavelengthUm', opts.emissionWavelengthUm, ...
        'excitationWavelengthUm', opts.excitationWavelengthUm, ...
        'sampleGeometry', 'airOnGlass', ...
        'interfaceRadialWeightMode', opts.interfaceRadialWeightMode, ...
        'diffractionModel', opts.diffractionModel, ...
        'vectorialPolarizationMode', opts.vectorialPolarizationMode, ...
        'excitationPolarizationMode', opts.excitationPolarizationMode, ...
        'collectionPolarizationMode', opts.collectionPolarizationMode, ...
        'airInterfaceStageMedium', opts.airInterfaceStageMedium, ...
        'immersionRefractiveIndex', opts.immersionRefractiveIndex, ...
        'glassRefractiveIndex', opts.glassRefractiveIndex, ...
        'sampleRefractiveIndex', opts.sampleRefractiveIndex, ...
        'designGlassRefractiveIndex', opts.designGlassRefractiveIndex, ...
        'coverslipThicknessUm', opts.coverslipThicknessUm, ...
        'designCoverslipThicknessUm', opts.designCoverslipThicknessUm, ...
        'beadBottomHeightUm', opts.beadBottomHeightUm, ...
        'airBeadAxialSamples', opts.airBeadAxialSamples, ...
        'interfaceGeometryConfirmed', opts.interfaceGeometryConfirmed, ...
        'runAirInterfaceSelfTest', opts.runAirInterfaceSelfTest, ...
        'beadRadiusUm', opts.beadRadiusUm, ...
        'beadSubsamples', opts.beadSubsamples, ...
        'detectorLayout', opts.detectorLayout, ...
        'detectorPixelShape', opts.detectorPixelShape, ...
        'estimateDetectorLayout', opts.estimateDetectorLayout, ...
        'detectorXYUm', opts.detectorXYUm, ...
        'detectorSubsamples', opts.detectorSubsamples, ...
        'detectorLayoutScale', opts.detectorLayoutScale, ...
        'detectorLayoutPositionSign', opts.detectorLayoutPositionSign, ...
        'detectorPitchSampleUm', opts.detectorPitchSampleUm, ...
        'detectorHardwarePitchUm', opts.detectorHardwarePitchUm, ...
        'detectorTotalMagnification', opts.detectorTotalMagnification, ...
        'flatField', opts.detectorQE, ...
        'fitModes', opts.fitModes, ...
        'initialCoeffs', centerFit.fit.estCoeffs, ...
        'initialZ0Um', center.focusZUm, ...
        'fitXY', opts.fitXY, ...
        'maxIter', opts.fullStackMaxIter, ...
        'outputDir', outDir, ...
        'writeOutputs', opts.writeOutputs, ...
        'verbose', opts.verbose);
end

function plane = correctedPlane(data, iz)
    plane = double(data.rawCounts(:,:,:,iz));
    if isfield(data, 'backgroundPerPixel') && ~isempty(data.backgroundPerPixel)
        plane = max(plane - double(data.backgroundPerPixel(:,:,:,iz)), 0);
    end
end

function comparison = generateAirInterface3DComparison( ...
        data, fit, focusIdx, focusZUm, fitKind, opts)
    outDir = fullfile(opts.outputDir, 'air_interface_3d_comparison');
    perDetectorDir = fullfile(outDir, 'per_detector');
    if opts.writeOutputs
        if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
        if opts.writePerDetectorFigures && exist(perDetectorDir, 'dir') ~= 7
            mkdir(perDetectorDir);
        end
    end

    [measured, info] = measuredVolumeForComparison(data, focusIdx, focusZUm, opts);
    [model, modelInfo] = modelAirInterfaceVolume(data, fit, focusZUm, fitKind, opts);

    measuredN = normalizeEachDetector(measured);
    modelN = normalizeEachDetector(model);
    residualN = measuredN - modelN;
    metrics = detectorMetrics(measuredN, modelN, data.channelIDs, opts.supportThreshold);

    comparison = struct();
    comparison.fitKind = fitKind;
    comparison.measuredVolume = measured;
    comparison.modelVolume = model;
    comparison.measuredNormalized = measuredN;
    comparison.modelNormalized = modelN;
    comparison.residualNormalized = residualN;
    comparison.dataInfo = info;
    comparison.modelInfo = modelInfo;
    comparison.metrics = metrics;
    comparison.outputDir = outDir;

    if opts.writeOutputs
        writetable(metrics, fullfile(outDir, 'per_detector_3d_metrics.csv'));
        if opts.writeComparisonVolumeMat
            save(fullfile(outDir, 'air_interface_3d_comparison_volumes.mat'), ...
                'comparison', '-v7.3');
        end
        if opts.writePerDetectorFigures
            nCh = size(measured, 4);
            for ch = 1:nCh
                outFile = fullfile(perDetectorDir, sprintf( ...
                    'detector_%02d_channel_%02d_air_interface_xy_xz_yz.png', ...
                    ch, data.channelIDs(ch)));
                writeDetectorComparisonFigure( ...
                    measuredN(:,:,:,ch), modelN(:,:,:,ch), residualN(:,:,:,ch), ...
                    info, opts, ch, data.channelIDs(ch), outFile);
            end
        end
        writeAllDetectorProjectionFigures(measuredN, modelN, residualN, info, opts, outDir);
        writeMetricSummaryFigure(metrics, fullfile(outDir, ...
            'per_detector_metric_summary.png'));
    end
end

function [volume, info] = measuredVolumeForComparison(data, focusIdx, focusZUm, opts)
    [ny, nx, nCh, nPlane] = size(data.rawCounts);
    volume = zeros(ny, nx, nPlane, nCh);
    for iz = 1:nPlane
        plane = correctedPlane(data, iz);
        if opts.applySavedAlignment && isfield(data, 'alignmentShiftXPx')
            sx = data.alignmentShiftXPx(iz);
            sy = data.alignmentShiftYPx(iz);
            shiftMag = hypot(sx, sy);
            if shiftMag >= opts.alignmentShiftMinPx
                for ch = 1:nCh
                    plane(:,:,ch) = max(fourierShift2D(plane(:,:,ch), sx, sy), 0);
                end
            end
        end
        volume(:,:,iz,:) = permute(plane, [1 2 4 3]);
    end

    recenterShiftXY = [0 0];
    if opts.recenterMeasuredVolume
        focusImage = squeeze(sum(volume(:,:,focusIdx,:), 4));
        centerXY = imageCentroid(focusImage);
        targetXY = [(nx+1)/2, (ny+1)/2];
        recenterShiftXY = targetXY - centerXY;
        for iz = 1:nPlane
            for ch = 1:nCh
                volume(:,:,iz,ch) = max(fourierShift2D( ...
                    volume(:,:,iz,ch), recenterShiftXY(1), recenterShiftXY(2)), 0);
            end
        end
    end

    info = struct();
    info.focusIndex = focusIdx;
    info.focusZUm = focusZUm;
    info.zUm = data.stageZUm(:).';
    info.relativeZUm = data.stageZUm(:).' - focusZUm;
    info.xUm = data.xUm;
    info.yUm = data.yUm;
    info.xyPixelSizeUm = data.xyPixelSizeUm;
    info.channelIDs = data.channelIDs;
    info.recenterShiftXY = recenterShiftXY;
    info.alignmentApplied = opts.applySavedAlignment;
    info.alignmentShiftMinPx = opts.alignmentShiftMinPx;
    if isfield(data, 'alignmentShiftXPx') && isfield(data, 'alignmentShiftYPx')
        shiftMag = hypot(data.alignmentShiftXPx(:), data.alignmentShiftYPx(:));
        info.savedShiftMagnitudePx = shiftMag(:).';
        info.nSavedShiftsAboveThreshold = nnz(shiftMag >= opts.alignmentShiftMinPx);
        info.nSavedShiftsBelowThreshold = nnz(shiftMag < opts.alignmentShiftMinPx);
    end
end

function [model, info] = modelAirInterfaceVolume(data, fit, focusZUm, fitKind, opts)
    sim = fit.sim;
    sim.sampleGeometry = 'airOnGlass';
    sim.interfaceRadialWeightMode = opts.interfaceRadialWeightMode;
    sim.nx = size(data.rawCounts, 2);
    sim.ny = size(data.rawCounts, 1);
    sim.x = data.xUm;
    sim.y = data.yUm;
    sim.dx = data.xyPixelSizeUm;
    sim.fovX = max(sim.x) - min(sim.x);
    sim.fovY = max(sim.y) - min(sim.y);
    sim.fovXY = max(sim.fovX, sim.fovY);

    if isfield(fit, 'estXYUm')
        xy = fit.estXYUm;
    elseif isfield(fit, 'estXY')
        xy = fit.estXY;
    else
        xy = [0 0];
    end

    if isfield(fit, 'estZ0Um')
        planeZForModel = data.stageZUm(:).' - fit.estZ0Um;
        zOffset = 0;
    elseif isfield(fit, 'estZOffset')
        planeZForModel = data.stageZUm(:).' - focusZUm;
        zOffset = fit.estZOffset;
    else
        planeZForModel = data.stageZUm(:).' - focusZUm;
        zOffset = 0;
    end

    model = normalizedStackAirInterfaceZPlanes( ...
        sim, fit.estCoeffs, planeZForModel, xy(1), xy(2), zOffset);
    model = permute(model, [1 2 4 3]); % [y x z detector]
    if ~isempty(opts.detectorQE)
        model = model .* reshape(opts.detectorQE, 1, 1, 1, []);
    end
    blurInfo = struct();
    blurInfo.sigmaXYUm = opts.modelBlurSigmaXYUm;
    blurInfo.sigmaZUm = opts.modelBlurSigmaZUm;
    blurInfo.applied = opts.modelBlurSigmaXYUm > 0 || opts.modelBlurSigmaZUm > 0;
    blurInfo.sigmaXPixel = opts.modelBlurSigmaXYUm / max(data.xyPixelSizeUm, eps);
    blurInfo.sigmaYPixel = blurInfo.sigmaXPixel;
    blurInfo.sigmaZPlane = opts.modelBlurSigmaZUm / max(medianPositiveDiff(data.stageZUm), eps);
    if blurInfo.applied
        model = blurDetectorVolumeModelOnly(model, blurInfo.sigmaXPixel, ...
            blurInfo.sigmaYPixel, blurInfo.sigmaZPlane);
    end

    info = struct();
    info.fitKind = fitKind;
    info.sim = sim;
    info.coeffs = fit.estCoeffs;
    info.xyUm = xy;
    info.planeZForModelUm = planeZForModel;
    info.zOffsetUm = zOffset;
    info.modelBlur = blurInfo;
end

function V = blurDetectorVolumeModelOnly(V, sigmaX, sigmaY, sigmaZ)
    if sigmaY > 0
        V = convolveVolumeAlongDim(V, gaussianKernel1D(sigmaY), 1);
    end
    if sigmaX > 0
        V = convolveVolumeAlongDim(V, gaussianKernel1D(sigmaX), 2);
    end
    if sigmaZ > 0
        V = convolveVolumeAlongDim(V, gaussianKernel1D(sigmaZ), 3);
    end
    V = max(V, 0);
end

function kernel = gaussianKernel1D(sigma)
    sigma = double(sigma);
    if sigma <= 0
        kernel = 1;
        return;
    end
    radius = max(1, ceil(4*sigma));
    x = -radius:radius;
    kernel = exp(-0.5*(x/sigma).^2);
    kernel = kernel / sum(kernel);
end

function Vout = convolveVolumeAlongDim(V, kernel, dim)
    if numel(kernel) == 1
        Vout = V;
        return;
    end
    radius = floor(numel(kernel)/2);
    Vout = zeros(size(V));
    n = size(V, dim);
    for k = 1:numel(kernel)
        offset = k - radius - 1;
        idx = min(max((1:n) + offset, 1), n);
        switch dim
            case 1
                shifted = V(idx,:,:,:);
            case 2
                shifted = V(:,idx,:,:);
            case 3
                shifted = V(:,:,idx,:);
            otherwise
                error('fitCenterCoordinateISMOilAirInterface:BadBlurDim', ...
                    'Blur dimension must be 1, 2, or 3.');
        end
        Vout = Vout + kernel(k) * shifted;
    end
end

function dz = medianPositiveDiff(z)
    z = sort(double(z(:)));
    d = diff(z);
    d = d(isfinite(d) & d > 0);
    if isempty(d)
        dz = 1;
    else
        dz = median(d);
    end
end

function Vn = normalizeEachDetector(V)
    Vn = zeros(size(V));
    for ch = 1:size(V,4)
        s = sum(V(:,:,:,ch), 'all');
        if s > 0 && isfinite(s)
            Vn(:,:,:,ch) = V(:,:,:,ch) / s;
        end
    end
end

function T = detectorMetrics(A, B, channelIDs, threshold)
    nCh = size(A, 4);
    detectorIndex = (1:nCh).';
    channelID = channelIDs(:);
    corr3D = nan(nCh,1);
    cosine3D = nan(nCh,1);
    nrmse3D = nan(nCh,1);
    measuredFraction = nan(nCh,1);
    modelFraction = nan(nCh,1);

    totalA = sum(A(:));
    totalB = sum(B(:));
    for ch = 1:nCh
        a = A(:,:,:,ch);
        b = B(:,:,:,ch);
        corr3D(ch) = maskedCorrelation(a, b, threshold);
        cosine3D(ch) = dot(a(:), b(:)) / max(norm(a(:))*norm(b(:)), eps);
        nrmse3D(ch) = norm(a(:)-b(:)) / max(norm(a(:)), eps);
        measuredFraction(ch) = sum(a(:)) / max(totalA, eps);
        modelFraction(ch) = sum(b(:)) / max(totalB, eps);
    end

    T = table(detectorIndex, channelID, corr3D, cosine3D, nrmse3D, ...
        measuredFraction, modelFraction);
end

function r = maskedCorrelation(a, b, threshold)
    a = double(a);
    b = double(b);
    mask = isfinite(a) & isfinite(b);
    lim = threshold * max([a(:); b(:)]);
    mask = mask & (a >= lim | b >= lim);
    if nnz(mask) < 3
        r = NaN;
        return;
    end
    aa = a(mask) - mean(a(mask));
    bb = b(mask) - mean(b(mask));
    r = dot(aa(:), bb(:)) / max(norm(aa(:))*norm(bb(:)), eps);
end

function centerXY = imageCentroid(img)
    img = double(img);
    img(~isfinite(img)) = 0;
    img = max(img - min(img(:)), 0);
    if sum(img(:)) <= 0
        [~, idx] = max(img(:));
        [iy, ix] = ind2sub(size(img), idx);
        centerXY = [ix iy];
        return;
    end
    [yy, xx] = ndgrid(1:size(img,1), 1:size(img,2));
    centerXY = [sum(xx(:).*img(:)), sum(yy(:).*img(:))] / sum(img(:));
end

function writeDetectorComparisonFigure(A, B, R, info, opts, detIdx, channelID, outFile)
    projections = {'xy','xz','yz'};
    volumes = {A, B, R};
    rowNames = {'measured', 'fitted air-interface model', 'measured - fitted'};
    fig = figure('Visible', opts.figureVisible, 'Color', 'w', ...
        'Position', [50 50 1180 900]);
    tl = tiledlayout(fig, 3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    for row = 1:3
        for col = 1:3
            ax = nexttile(tl);
            img = projectOne(volumes{row}, projections{col}, opts.projectionMode);
            [xAxis, yAxis, xLabel, yLabel] = projectionAxes(info, projections{col});
            imagesc(ax, xAxis, yAxis, img);
            axis(ax, 'image');
            set(ax, 'YDir', 'normal');
            xlabel(ax, xLabel);
            ylabel(ax, yLabel);
            title(ax, sprintf('%s: %s', rowNames{row}, upper(projections{col})));
            if row == 3
                lim = max(abs(img(:)));
                if lim > 0, caxis(ax, [-lim lim]); end
                colormap(ax, redBlueMap(256));
            else
                colormap(ax, hot);
            end
            colorbar(ax);
        end
    end
    title(tl, sprintf('Detector %02d / PTU channel %02d', detIdx, channelID), ...
        'Interpreter', 'none');
    exportgraphics(fig, outFile, 'Resolution', 180);
    close(fig);
end

function writeAllDetectorProjectionFigures(measured, model, residual, info, opts, outDir)
    projections = {'xy','xz','yz'};
    for k = 1:numel(projections)
        planeName = projections{k};
        measuredCanvas = detectorProjectionCanvas(measured, planeName, opts);
        modelCanvas = detectorProjectionCanvas(model, planeName, opts);
        residualCanvas = detectorProjectionCanvas(residual, planeName, opts);
        outFile = fullfile(outDir, sprintf( ...
            'all_detectors_%s_air_interface_projection_comparison.png', upper(planeName)));
        writeAllDetectorProjectionFigure(measuredCanvas, modelCanvas, residualCanvas, ...
            planeName, opts, outFile);
    end
end

function writeAllDetectorProjectionFigure(measuredCanvas, modelCanvas, residualCanvas, ...
        planeName, opts, outFile)
    fig = figure('Visible', opts.figureVisible, 'Color', 'w', ...
        'Position', [60 60 1280 920]);
    tl = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    canvases = {measuredCanvas, modelCanvas, residualCanvas};
    titles = {'measured', 'fitted air-interface model', 'measured - fitted'};
    for row = 1:3
        ax = nexttile(tl);
        imagesc(ax, canvases{row});
        axis(ax, 'image');
        axis(ax, 'off');
        title(ax, titles{row});
        if row == 3
            lim = max(abs(canvases{row}(:)), [], 'omitnan');
            if ~isfinite(lim) || lim <= 0, lim = 1; end
            caxis(ax, [-lim lim]);
            colormap(ax, redBlueMap(256));
        else
            finiteValues = canvases{row}(isfinite(canvases{row}));
            if ~isempty(finiteValues)
                lo = min(finiteValues);
                hi = max(finiteValues);
                if hi > lo
                    caxis(ax, [lo hi]);
                end
            end
            colormap(ax, hot);
        end
        colorbar(ax);
    end
    title(tl, sprintf('All 23 detector pixels: %s projection comparison', upper(planeName)), ...
        'Interpreter', 'none');
    exportgraphics(fig, outFile, 'Resolution', 180);
    close(fig);
end

function canvas = detectorProjectionCanvas(volume, planeName, opts)
    nCh = size(volume, 4);
    try
        [~, idxGrid] = detectorLayout(opts.detectorLayout, 1);
    catch
        idxGrid = detectorIndexGridFallback(nCh);
    end
    if nnz(isfinite(idxGrid)) ~= nCh
        idxGrid = detectorIndexGridFallback(nCh);
    end

    firstImage = projectOne(volume(:,:,:,1), planeName, opts.projectionMode);
    tileSize = size(firstImage);
    gap = max(2, round(min(tileSize) * 0.08));
    [nRow, nCol] = size(idxGrid);
    canvas = nan(nRow * tileSize(1) + (nRow-1) * gap, ...
        nCol * tileSize(2) + (nCol-1) * gap);

    for r = 1:nRow
        for c = 1:nCol
            ch = idxGrid(r,c);
            if ~isfinite(ch) || ch < 1 || ch > nCh
                continue;
            end
            img = projectOne(volume(:,:,:,ch), planeName, opts.projectionMode);
            y0 = (r-1) * (tileSize(1) + gap) + 1;
            x0 = (c-1) * (tileSize(2) + gap) + 1;
            canvas(y0:y0+tileSize(1)-1, x0:x0+tileSize(2)-1) = img;
        end
    end
end

function idxGrid = detectorIndexGridFallback(nCh)
    nCol = ceil(sqrt(nCh));
    nRow = ceil(nCh / nCol);
    idxGrid = nan(nRow, nCol);
    for k = 1:nCh
        r = floor((k-1) / nCol) + 1;
        c = mod(k-1, nCol) + 1;
        idxGrid(r,c) = k;
    end
end

function image = projectOne(vol, planeName, mode)
    switch lower(planeName)
        case 'xy'
            if strcmp(mode, 'sum')
                image = squeeze(sum(vol, 3));
            else
                image = squeeze(max(vol, [], 3));
            end
        case 'xz'
            if strcmp(mode, 'sum')
                image = squeeze(sum(vol, 1)).';
            else
                image = squeeze(max(vol, [], 1)).';
            end
        case 'yz'
            if strcmp(mode, 'sum')
                image = squeeze(sum(vol, 2)).';
            else
                image = squeeze(max(vol, [], 2)).';
            end
        otherwise
            error('fitCenterCoordinateISMOilAirInterface:BadProjection', ...
                'Projection must be XY, XZ, or YZ.');
    end
end

function [xAxis, yAxis, xLabel, yLabel] = projectionAxes(info, planeName)
    switch lower(planeName)
        case 'xy'
            xAxis = info.xUm;
            yAxis = info.yUm;
            xLabel = 'x [um]';
            yLabel = 'y [um]';
        case 'xz'
            xAxis = info.xUm;
            yAxis = info.relativeZUm;
            xLabel = 'x [um]';
            yLabel = 'z from selected focus [um]';
        case 'yz'
            xAxis = info.yUm;
            yAxis = info.relativeZUm;
            xLabel = 'y [um]';
            yLabel = 'z from selected focus [um]';
    end
end

function cmap = redBlueMap(n)
    if nargin < 1, n = 256; end
    x = linspace(-1, 1, n).';
    cmap = [max(-x,0), 1-abs(x), max(x,0)];
    cmap = max(min(cmap, 1), 0);
end

function writeMetricSummaryFigure(metrics, outFile)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 920 320]);
    tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax = nexttile(tl);
    plot(ax, metrics.detectorIndex, metrics.corr3D, '-o', 'LineWidth', 1);
    ylim(ax, [-1 1]);
    grid(ax, 'on');
    xlabel(ax, 'detector index');
    ylabel(ax, '3D correlation');
    ax = nexttile(tl);
    plot(ax, metrics.detectorIndex, metrics.cosine3D, '-o', 'LineWidth', 1);
    ylim(ax, [0 1]);
    grid(ax, 'on');
    xlabel(ax, 'detector index');
    ylabel(ax, 'cosine similarity');
    ax = nexttile(tl);
    plot(ax, metrics.detectorIndex, metrics.nrmse3D, '-o', 'LineWidth', 1);
    grid(ax, 'on');
    xlabel(ax, 'detector index');
    ylabel(ax, 'normalized RMSE');
    exportgraphics(fig, outFile, 'Resolution', 180);
    close(fig);
end

function T = centerDetectorIntensityTable(result)
    data = result.centerSelection.data;
    centerFit = result.centerFit;
    focusIdx = result.centerSelection.focusIndex;
    centerXY = result.centerSelection.centerXY;
    sampleMode = result.options.centerSampleMode;
    nCh = numel(data.channelIDs);

    rawPlane = double(data.rawCounts(:,:,:,focusIdx));
    corrected = correctedPlane(data, focusIdx);
    rawAtCenter = sampleDetectorPlaneAtPixelXY(rawPlane, centerXY, sampleMode);
    correctedAtCenter = sampleDetectorPlaneAtPixelXY(corrected, centerXY, sampleMode);
    backgroundAtCenter = rawAtCenter - correctedAtCenter;

    measuredCorrectedFitCounts = centerFit.centerValues(:,1);
    measuredNormalized = centerFit.normalizedCenterValues(:,1);
    recoveredNormalized = centerFit.fit.modelN(:,1);
    residualNormalized = measuredNormalized - recoveredNormalized;
    scaleCounts = normalizationScaleForPlane(centerFit, 1, measuredCorrectedFitCounts);
    recoveredFitCounts = recoveredNormalized * scaleCounts;
    residualFitCounts = measuredCorrectedFitCounts - recoveredFitCounts;

    detectorIndex = (1:nCh).';
    channelID = data.channelIDs(:);
    fitPositionXPx = repmat(centerXY(1), nCh, 1);
    fitPositionYPx = repmat(centerXY(2), nCh, 1);
    focusPlaneIndex = repmat(focusIdx, nCh, 1);
    focusZUm = repmat(result.centerSelection.focusZUm, nCh, 1);

    detXYUm = nan(nCh,2);
    detXYAU = nan(nCh,2);
    if isfield(centerFit.fit, 'sim') && isfield(centerFit.fit.sim, 'detXY') && ...
            size(centerFit.fit.sim.detXY,1) == nCh
        detXYUm = centerFit.fit.sim.detXY;
        airyUnitUm = centerFit.fit.sim.airyUnitUm;
        if isfinite(airyUnitUm) && airyUnitUm > 0
            detXYAU = detXYUm / airyUnitUm;
        end
    end

    T = table(detectorIndex, channelID, focusPlaneIndex, focusZUm, ...
        fitPositionXPx, fitPositionYPx, ...
        rawAtCenter, backgroundAtCenter, correctedAtCenter, ...
        measuredCorrectedFitCounts, measuredNormalized, ...
        recoveredNormalized, residualNormalized, ...
        recoveredFitCounts, residualFitCounts, ...
        detXYUm(:,1), detXYUm(:,2), detXYAU(:,1), detXYAU(:,2), ...
        'VariableNames', {'detectorIndex','channelID','focusPlaneIndex','focusZUm', ...
        'fitPositionXPx','fitPositionYPx', ...
        'rawCountsAtFitPosition','backgroundCountsAtFitPosition', ...
        'backgroundCorrectedCountsAtFitPosition', ...
        'fitInputCountsAtFitPosition','measuredNormalizedIntensity', ...
        'recoveredNormalizedIntensity','normalizedResidual', ...
        'recoveredCountsAtFitPosition','countsResidual', ...
        'detectorXUm','detectorYUm','detectorXAU','detectorYAU'});
end

function values = sampleDetectorPlaneAtPixelXY(plane, centerXY, mode)
    plane = double(plane);
    nCh = size(plane, 3);
    x = double(centerXY(1));
    y = double(centerXY(2));
    nx = size(plane, 2);
    ny = size(plane, 1);
    values = zeros(nCh, 1);
    switch lower(char(mode))
        case {'nearest','round','pixel'}
            ix = min(max(round(x), 1), nx);
            iy = min(max(round(y), 1), ny);
            values = reshape(plane(iy, ix, :), nCh, 1);
        otherwise
            x = min(max(x, 1), nx);
            y = min(max(y, 1), ny);
            for ch = 1:nCh
                values(ch) = interp2(1:nx, 1:ny, plane(:,:,ch), x, y, 'linear', 0);
            end
    end
end

function scaleCounts = normalizationScaleForPlane(centerFit, planeIndex, fallbackValues)
    scaleCounts = sum(fallbackValues);
    if isfield(centerFit, 'normalization') && isfield(centerFit.normalization, 'scale')
        scale = centerFit.normalization.scale;
        if numel(scale) >= planeIndex && isfinite(scale(planeIndex)) && scale(planeIndex) > 0
            scaleCounts = double(scale(planeIndex));
        elseif isscalar(scale) && isfinite(scale) && scale > 0
            scaleCounts = double(scale);
        end
    end
    if ~isfinite(scaleCounts) || scaleCounts <= 0
        scaleCounts = 1;
    end
end

function writeCenterDetectorIntensityOutputs(result, outDir)
    T = result.centerDetectorIntensity;
    writetable(T, fullfile(outDir, 'center_fit_detector_intensity_at_fit_position.csv'));
    writeCenterDetectorMicroimageFigure(result, fullfile(outDir, ...
        'center_fit_detector_microimages.png'));
end

function writeCenterDetectorMicroimageFigure(result, outFile)
    T = result.centerDetectorIntensity;
    displayXY = detectorDisplayXYForCenterOutput(result, height(T));
    values = { ...
        T.rawCountsAtFitPosition, ...
        T.backgroundCorrectedCountsAtFitPosition, ...
        T.recoveredCountsAtFitPosition, ...
        T.measuredNormalizedIntensity, ...
        T.recoveredNormalizedIntensity, ...
        T.normalizedResidual};
    titles = { ...
        'raw counts at fit position', ...
        'background-corrected counts used for fit', ...
        'recovered fitted counts at fit position', ...
        'measured normalized detector image', ...
        'recovered normalized detector image', ...
        'normalized residual'};

    fig = figure('Visible', result.options.figureVisible, 'Color', 'w', ...
        'Position', [80 80 1180 720]);
    tl = tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    for k = 1:numel(values)
        ax = nexttile(tl);
        signedMap = k == 6;
        if signedMap
            lim = max(abs(values{k}(:)));
            if ~isfinite(lim) || lim <= 0, lim = 1; end
            plotDetectorHexMap(displayXY, values{k}, 'Parent', ax, 'CLim', [-lim lim]);
            colormap(ax, redBlueMap(256));
        else
            plotDetectorHexMap(displayXY, values{k}, 'Parent', ax);
            colormap(ax, parula);
        end
        title(ax, titles{k}, 'Interpreter', 'none');
        colorbar(ax);
        overlayDetectorNumbers(ax, displayXY, T.channelID, values{k});
    end
    title(tl, sprintf('Detector intensity at fitted bead center: x=%.3f px, y=%.3f px; focus z=%.4f um', ...
        result.centerSelection.centerXY(1), result.centerSelection.centerXY(2), ...
        result.centerSelection.focusZUm), 'Interpreter', 'none');
    exportgraphics(fig, outFile, 'Resolution', 180);
    close(fig);
end

function displayXY = detectorDisplayXYForCenterOutput(result, nCh)
    layoutName = 'honeycomb23';
    if isfield(result.options, 'detectorLayout') && ~isempty(result.options.detectorLayout)
        layoutName = result.options.detectorLayout;
    end
    try
        displayXY = detectorLayout(layoutName, 1);
    catch
        displayXY = [];
    end
    if size(displayXY,1) ~= nCh
        displayXY = result.centerFit.fit.sim.detXY;
    end
end

function overlayDetectorNumbers(ax, detXY, channelIDs, values)
    hold(ax, 'on');
    for k = 1:numel(values)
        text(ax, detXY(k,1), detXY(k,2), sprintf('%g\n%.3g', channelIDs(k), values(k)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 7, 'Color', [0 0 0], 'FontWeight', 'bold');
    end
    hold(ax, 'off');
end

function writeCoefficientCsv(result, outFile)
    fit = result.centerFit.fit;
    source = repmat({'singlePlaneCenterCoordinate'}, numel(result.options.fitModes), 1);
    if result.options.useMultiPlane && ~isempty(result.multiPlaneFit)
        fit = result.multiPlaneFit.fit;
        source = repmat({'multiPlaneFullStack'}, numel(result.options.fitModes), 1);
    end
    modes = result.options.fitModes(:);
    waves = zeros(numel(modes), 1);
    for k = 1:numel(modes)
        if isfield(fit.estCoeffs, modes{k})
            waves(k) = fit.estCoeffs.(modes{k});
        end
    end
    nmOPD = waves * fit.sim.lamRef * 1000;
    T = table(source, modes, waves, nmOPD, ...
        'VariableNames', {'fitSource','mode','waves','opdNmAtLambdaRef'});
    writetable(T, outFile);
end

function printSummary(result)
    fprintf('\n%s\n', result.reportTitle);
    fprintf('  focus plane: %d at z = %.4f um\n', ...
        result.centerSelection.focusIndex, result.centerSelection.focusZUm);
    fprintf('  bead center: x = %.2f px, y = %.2f px\n', ...
        result.centerSelection.centerXY(1), result.centerSelection.centerXY(2));
    fprintf('  finite bead model: radius %.1f nm (diameter %.1f nm), axial samples %d\n', ...
        1000*result.options.beadRadiusUm, 2000*result.options.beadRadiusUm, ...
        result.options.airBeadAxialSamples);
    fprintf('  detector positions preserved: estimateDetectorLayout=%d, fitDetectorPitchScale=%d\n', ...
        result.options.estimateDetectorLayout, result.options.fitDetectorPitchScale);
    fprintf('  model-only comparison blur: sigmaXY %.0f nm, sigmaZ %.0f nm\n', ...
        1000*result.options.modelBlurSigmaXYUm, ...
        1000*result.options.modelBlurSigmaZUm);
    fprintf('  single-plane center fit residual: %.4g\n', ...
        result.centerFit.fit.residualNorm);
    if isfield(result, 'centerDetectorIntensity') && ...
            ~isempty(result.centerDetectorIntensity)
        fprintf('  detector intensities at fitted center position:\n');
        T = result.centerDetectorIntensity;
        disp(T(:, {'detectorIndex','channelID', ...
            'rawCountsAtFitPosition', ...
            'backgroundCorrectedCountsAtFitPosition', ...
            'recoveredCountsAtFitPosition', ...
            'measuredNormalizedIntensity', ...
            'recoveredNormalizedIntensity', ...
            'normalizedResidual'}));
    end
    if isfield(result.centerFit, 'sufficiency')
        s = result.centerFit.sufficiency;
        fprintf('  single-plane rank: %d / %d parameters\n', ...
            s.rank, s.nParameters);
    end
    if result.options.useMultiPlane && ~isempty(result.multiPlaneFit)
        fprintf('  multi-plane full-stack deviance: %.4g\n', ...
            result.multiPlaneFit.fit.deviance);
    end
    if ~isempty(result.comparison3D)
        m = result.comparison3D.metrics;
        fprintf('  median per-detector 3D correlation: %.3f\n', ...
            median(m.corr3D, 'omitnan'));
    end
    fprintf('  outputs: %s\n', result.outputDir);
end
