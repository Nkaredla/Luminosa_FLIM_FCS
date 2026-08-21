function result = fitWaterAirCenterCoordinateWithObjectiveCalibration(stackInput, calibrationFile, varargin)
%FITWATERAIRCENTERCOORDINATEWITHOBJECTIVECALIBRATION Fit bead using saved calibration.
%
%   result = fitWaterAirCenterCoordinateWithObjectiveCalibration(stackInput, calibrationFile)
%
%   Loads a calibration produced by calibrateWaterObjectiveISMFromBead3D and
%   applies its fixed nuisance parameters to the water/glass/air center-
%   coordinate fit for another bead:
%       - initial water-objective wavefront coefficients,
%       - fixed regular PDA23 detector pitch / magnification,
%       - relative detector QE / collection efficiency.
%
%   Additional Name/Value pairs override the calibration-derived defaults.

    if nargin < 1
        stackInput = [];
    end
    if nargin < 2 || isempty(calibrationFile)
        calibrationFile = defaultCalibrationFile();
    end

    calibration = loadCalibration(calibrationFile);
    outDir = defaultOutputDir(stackInput, calibration);
    args = calibrationArgs(calibration, outDir);
    result = fitCenterCoordinateISMWaterAirInterface( ...
        stackInput, args{:}, varargin{:});
    result.objectiveCalibration = calibration;
    result.objectiveCalibrationFile = char(calibrationFile);
end

function args = calibrationArgs(calibration, outDir)
    c = calibration;
    detector = c.detector;
    psfSigma = c.psfSigma;
    objective = c.objective;
    defaults = c.centerCoordinateFitDefaults;

    args = { ...
        'objectiveNA', objective.objectiveNA, ...
        'objectiveMagnification', objective.objectiveMagnification, ...
        'immersionRefractiveIndex', objective.immersionRefractiveIndex, ...
        'glassRefractiveIndex', objective.glassRefractiveIndex, ...
        'sampleRefractiveIndex', objective.sampleRefractiveIndex, ...
        'interfaceRadialWeightMode', calibrationInterfaceRadialWeightMode(objective), ...
        'diffractionModel', calibrationDiffractionModel(objective, defaults), ...
        'vectorialPolarizationMode', calibrationVectorialPolarizationMode(objective, defaults), ...
        'excitationPolarizationMode', calibrationPathPolarizationMode( ...
            objective, defaults, 'excitationPolarizationMode', 'x'), ...
        'collectionPolarizationMode', calibrationPathPolarizationMode( ...
            objective, defaults, 'collectionPolarizationMode', 'xyAverage'), ...
        'airInterfaceStageMedium', calibrationAirInterfaceStageMedium(objective, defaults), ...
        'coverslipThicknessUm', objective.coverslipThicknessUm, ...
        'detectorLayout', detector.detectorLayout, ...
        'detectorPixelShape', detector.detectorPixelShape, ...
        'detectorPitchSampleUm', detector.detectorPitchSampleUm, ...
        'detectorHardwarePitchUm', detector.detectorHardwarePitchUm, ...
        'detectorTotalMagnification', detector.detectorTotalMagnification, ...
        'detectorQE', detector.relativeQE, ...
        'estimateDetectorLayout', false, ...
        'fitDetectorPitchScale', false, ...
        'initialCoeffs', defaults.initialCoeffs, ...
        'fitModes', defaults.fitModes, ...
        'modelBlurSigmaXYUm', 0, ...
        'modelBlurSigmaZUm', 0, ...
        'outputPrefix', 'water_air_center_coordinate_calibrated', ...
        'reportTitle', 'Water/glass/air center-coordinate ISM fit with saved objective calibration', ...
        'outputDir', outDir};
end

function mode = calibrationInterfaceRadialWeightMode(objective)
    if isfield(objective, 'interfaceRadialWeightMode') && ...
            ~isempty(objective.interfaceRadialWeightMode)
        mode = objective.interfaceRadialWeightMode;
    else
        mode = 'sampleSolidAngle';
    end
end

function model = calibrationDiffractionModel(objective, defaults)
    if isfield(defaults, 'diffractionModel') && ~isempty(defaults.diffractionModel)
        model = defaults.diffractionModel;
    elseif isfield(objective, 'diffractionModel') && ~isempty(objective.diffractionModel)
        model = objective.diffractionModel;
    else
        model = 'vectorial';
    end
end

function mode = calibrationVectorialPolarizationMode(objective, defaults)
    if isfield(defaults, 'vectorialPolarizationMode') && ...
            ~isempty(defaults.vectorialPolarizationMode)
        mode = defaults.vectorialPolarizationMode;
    elseif isfield(objective, 'vectorialPolarizationMode') && ...
            ~isempty(objective.vectorialPolarizationMode)
        mode = objective.vectorialPolarizationMode;
    else
        mode = 'xyAverage';
    end
end

function mode = calibrationPathPolarizationMode(objective, defaults, fieldName, fallback)
    if isfield(defaults, fieldName) && ~isempty(defaults.(fieldName))
        mode = defaults.(fieldName);
    elseif isfield(objective, fieldName) && ~isempty(objective.(fieldName))
        mode = objective.(fieldName);
    else
        mode = fallback;
    end
end

function medium = calibrationAirInterfaceStageMedium(objective, defaults)
    if isfield(defaults, 'airInterfaceStageMedium') && ...
            ~isempty(defaults.airInterfaceStageMedium)
        medium = defaults.airInterfaceStageMedium;
    elseif isfield(objective, 'airInterfaceStageMedium') && ...
            ~isempty(objective.airInterfaceStageMedium)
        medium = objective.airInterfaceStageMedium;
    else
        medium = 'immersion';
    end
end

function calibration = loadCalibration(fileName)
    fileName = char(fileName);
    if exist(fileName, 'file') ~= 2
        error('fitWaterAirCenterCoordinateWithObjectiveCalibration:MissingCalibration', ...
            'Calibration file was not found: %s', fileName);
    end
    S = load(fileName);
    if isfield(S, 'calibration') && isstruct(S.calibration)
        calibration = S.calibration;
        return;
    end
    names = fieldnames(S);
    for k = 1:numel(names)
        value = S.(names{k});
        if isstruct(value) && isfield(value, 'centerCoordinateFitDefaults') && ...
                isfield(value, 'detector')
            calibration = value;
            return;
        end
    end
    error('fitWaterAirCenterCoordinateWithObjectiveCalibration:BadCalibration', ...
        'No calibration struct was found in %s.', fileName);
end

function outDir = defaultOutputDir(stackInput, calibration)
    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin < 1 || isempty(stackInput)
        stem = 'default_stack';
    else
        [~, stem] = fileparts(char(stackInput));
        stem = regexprep(stem, '_frame_alignment$', '');
    end
    if isfield(calibration, 'createdAt')
        tag = calibration.createdAt;
    else
        tag = 'calibrated';
    end
    outDir = fullfile(root, 'output_matlab', ...
        'water_air_center_coordinate_ism', ...
        [sanitizeName(stem) '_using_calibration_' sanitizeName(tag)]);
end

function fileName = defaultCalibrationFile()
    root = fileparts(fileparts(mfilename('fullpath')));
    fileName = fullfile(root, 'output_matlab', ...
        'water_objective_calibration', 'test_20260515_144001', ...
        'water_objective_ism_calibration.mat');
end

function name = sanitizeName(name)
    name = regexprep(char(name), '[^A-Za-z0-9_\-]+', '_');
    if isempty(name)
        name = 'calibrated';
    end
end
