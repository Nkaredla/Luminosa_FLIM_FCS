function calibration = calibrateWaterObjectiveISMFromBead3D(stackInput, varargin)
%CALIBRATEWATEROBJECTIVEISMFROMBEAD3D Calibrate water-objective ISM nuisance terms.
%
%   calibration = calibrateWaterObjectiveISMFromBead3D(stackInput)
%
%   Uses one bead's full 3-D detector-resolved stack to estimate:
%       - water-objective wavefront coefficients,
%       - in-focus measured/model PSF sigma diagnostic,
%       - detector pitch in sample space and total magnification when the
%         hardware detector pitch is supplied,
%       - relative detector QE/gain, normalized to mean 1.
%
%   The calibration is saved as water_objective_ism_calibration.mat and can
%   be applied to another bead with:
%
%       fitWaterAirCenterCoordinateWithObjectiveCalibration(otherStack, ...
%           calibration.calibrationFile)

    if nargin < 1 || isempty(stackInput)
        stackInput = ['D:\Luminosa\Data\PSF_batch_outputs\ISM_Aberation_72\' ...
            'xz_yz_plots\test_20260515_144001_frame_alignment.csv'];
    end

    opts = parseCalibrationOptions(varargin{:});
    addRequiredPaths();
    opts.outputDir = resolveCalibrationOutputDir(stackInput, opts.outputDir);
    if exist(opts.outputDir, 'dir') ~= 7
        mkdir(opts.outputDir);
    end
    if isempty(opts.cacheFile)
        opts.cacheFile = fullfile(opts.outputDir, 'calibration_raw_detector_cache.mat');
    end
    if isempty(opts.calibrationFile)
        opts.calibrationFile = fullfile(opts.outputDir, ...
            'water_objective_ism_calibration.mat');
    end

    nominalPitchUm = nominalDetectorPitchUm(opts);
    if opts.verbose
        fprintf('\n[water-objective calibration] Initial full-stack fit at pitch %.5g um\n', ...
            nominalPitchUm);
    end
    initialResult = runWaterFullStackFit(stackInput, opts, nominalPitchUm, ...
        struct(), [], fullfile(opts.outputDir, 'initial_no_qe'), ...
        opts.initialMaxIter, false, opts.reuseCache, true);

    [bestPitchUm, pitchProfile] = chooseDetectorPitch( ...
        stackInput, opts, nominalPitchUm, initialResult.fit.estCoeffs);

    if opts.verbose
        fprintf('[water-objective calibration] Final no-QE fit at pitch %.5g um\n', ...
            bestPitchUm);
    end
    noQEResult = runWaterFullStackFit(stackInput, opts, bestPitchUm, ...
        initialResult.fit.estCoeffs, [], fullfile(opts.outputDir, 'final_no_qe'), ...
        opts.maxIter, opts.writeOutputs, true, false);

    [detectorQE, detectorTable] = estimateRelativeDetectorQE(noQEResult);
    finalResult = noQEResult;
    if opts.refitWithDetectorQE
        if opts.verbose
            fprintf('[water-objective calibration] Refitting with estimated relative detector QE\n');
        end
        finalResult = runWaterFullStackFit(stackInput, opts, bestPitchUm, ...
            noQEResult.fit.estCoeffs, detectorQE, ...
            fullfile(opts.outputDir, 'final_with_qe'), ...
            opts.qeRefitMaxIter, opts.writeOutputs, true, false);
    end

    psfSigma = estimateFocusPsfSigma(finalResult);
    totalMagnification = detectorTotalMagnificationFromPitchLocal( ...
        bestPitchUm, opts.detectorHardwarePitchUm, opts.detectorTotalMagnification);

    % Odd modes (coma/tilt) carry an uncertified sign through the
    % rot90/detectorImageInverted/positionSign convention chain. Unless a
    % passing known-aberration sign calibration is supplied, zero them in the
    % saved calibration so one wrong sign cannot propagate to every bead the
    % calibration is later applied to.
    [calibratedCoeffs, oddModeCertification] = ...
        guardOddModeSigns(finalResult.fit.estCoeffs, opts);

    calibration = struct();
    calibration.version = 1;
    calibration.createdAt = datestr(now, 30);
    calibration.calibrationFile = opts.calibrationFile;
    calibration.outputDir = opts.outputDir;
    calibration.sourceStack = char(stackInput);
    calibration.objective = struct( ...
        'description', '60x 1.2 NA water objective, water/glass/air interface', ...
        'objectiveNA', opts.objectiveNA, ...
        'objectiveMagnification', opts.objectiveMagnification, ...
        'immersionRefractiveIndex', opts.immersionRefractiveIndex, ...
        'glassRefractiveIndex', opts.glassRefractiveIndex, ...
        'sampleRefractiveIndex', opts.sampleRefractiveIndex, ...
        'interfaceRadialWeightMode', opts.interfaceRadialWeightMode, ...
        'diffractionModel', opts.diffractionModel, ...
        'vectorialPolarizationMode', opts.vectorialPolarizationMode, ...
        'excitationPolarizationMode', opts.excitationPolarizationMode, ...
        'collectionPolarizationMode', opts.collectionPolarizationMode, ...
        'airInterfaceStageMedium', opts.airInterfaceStageMedium, ...
        'coverslipThicknessUm', opts.coverslipThicknessUm);
    calibration.wavefront = struct( ...
        'coeffs', calibratedCoeffs, ...
        'fittedCoeffsBeforeSignGuard', finalResult.fit.estCoeffs, ...
        'fitModes', {opts.fitModes}, ...
        'lambdaRefUm', finalResult.fit.sim.lamRef, ...
        'z0Um', finalResult.fit.estZ0Um, ...
        'deviance', finalResult.fit.deviance, ...
        'converged', finalResult.fit.converged, ...
        'terminationReason', finalResult.fit.terminationReason);
    calibration.detector = struct( ...
        'channelIDs', finalResult.data.channelIDs(:).', ...
        'relativeQE', detectorQE(:).', ...
        'detectorXYUm', finalResult.fit.sim.detXY, ...
        'detectorPitchSampleUm', bestPitchUm, ...
        'detectorHardwarePitchUm', opts.detectorHardwarePitchUm, ...
        'detectorTotalMagnification', totalMagnification, ...
        'detectorLayout', finalResult.fit.sim.detectorLayout, ...
        'detectorPixelShape', finalResult.fit.sim.detectorPixelShape);
    calibration.psfSigma = psfSigma;
    calibration.centerCoordinateFitDefaults = struct( ...
        'initialCoeffs', calibratedCoeffs, ...
        'fitModes', {opts.fitModes}, ...
        'detectorQE', detectorQE(:).', ...
        'detectorPitchSampleUm', bestPitchUm, ...
        'detectorTotalMagnification', totalMagnification, ...
        'modelBlurSigmaXYUm', 0, ...
        'modelBlurSigmaZUm', 0, ...
        'diffractionModel', opts.diffractionModel, ...
        'interfaceRadialWeightMode', opts.interfaceRadialWeightMode, ...
        'vectorialPolarizationMode', opts.vectorialPolarizationMode, ...
        'excitationPolarizationMode', opts.excitationPolarizationMode, ...
        'collectionPolarizationMode', opts.collectionPolarizationMode, ...
        'airInterfaceStageMedium', opts.airInterfaceStageMedium, ...
        'estimateDetectorLayout', false, ...
        'fitDetectorPitchScale', false);
    calibration.pitchProfile = pitchProfile;
    calibration.oddModeCertification = oddModeCertification;
    calibration.options = opts;

    if opts.writeOutputs
        save(opts.calibrationFile, 'calibration', 'finalResult', ...
            'noQEResult', 'pitchProfile', '-v7.3');
        writetable(detectorTable, fullfile(opts.outputDir, ...
            'water_objective_detector_relative_qe.csv'));
        writetable(wavefrontTable(calibration), fullfile(opts.outputDir, ...
            'water_objective_wavefront_coefficients.csv'));
        writetable(psfSigmaTable(psfSigma), fullfile(opts.outputDir, ...
            'water_objective_psf_sigma.csv'));
        if ~isempty(pitchProfile)
            writetable(pitchProfile, fullfile(opts.outputDir, ...
                'water_objective_detector_pitch_profile.csv'));
        end
    end

    if opts.verbose
        fprintf('\nWater-objective ISM calibration saved:\n  %s\n', ...
            opts.calibrationFile);
        fprintf('  detector pitch %.5g um; total magnification %.5g x\n', ...
            bestPitchUm, totalMagnification);
        fprintf('  relative QE range %.3f to %.3f\n', ...
            min(detectorQE), max(detectorQE));
        fprintf(['  measured-model sigma gap diagnostic: XY %.0f nm, Z %.0f nm ' ...
            '(not applied as blur)\n\n'], ...
            1000*psfSigma.empiricalModelBlurSigmaXYUm, ...
            1000*psfSigma.empiricalModelBlurSigmaZUm);
    end
end

function opts = parseCalibrationOptions(varargin)
    p = inputParser;
    p.FunctionName = 'calibrateWaterObjectiveISMFromBead3D';

    addParameter(p, 'outputDir', '');
    addParameter(p, 'calibrationFile', '');
    addParameter(p, 'cacheFile', '');
    addParameter(p, 'reuseCache', false);
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'verbose', true);

    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'ptuReaderFolder', defaultPtuReaderFolderLocal());
    addParameter(p, 'ptuPhotonsPerChunk', 1e6);
    addParameter(p, 'darkFile', 'D:\Luminosa\Data\ISMdark_counts.ptu');
    addParameter(p, 'darkScale', 1);
    addParameter(p, 'backgroundMode', 'auto');
    addParameter(p, 'xyPixelSizeUm', []);

    addParameter(p, 'objectiveNA', 1.2);
    addParameter(p, 'objectiveMagnification', 60);
    addParameter(p, 'emissionWavelengthUm', []);
    addParameter(p, 'excitationWavelengthUm', []);
    addParameter(p, 'immersionRefractiveIndex', 1.33);
    addParameter(p, 'glassRefractiveIndex', 1.518);
    addParameter(p, 'sampleRefractiveIndex', 1.0003);
    addParameter(p, 'interfaceRadialWeightMode', 'sampleSolidAngle');
    addParameter(p, 'diffractionModel', 'vectorial');
    addParameter(p, 'vectorialPolarizationMode', 'xyAverage');
    addParameter(p, 'excitationPolarizationMode', 'x');
    addParameter(p, 'collectionPolarizationMode', 'xyAverage');
    addParameter(p, 'airInterfaceStageMedium', 'immersion');
    addParameter(p, 'designGlassRefractiveIndex', 1.518);
    addParameter(p, 'coverslipThicknessUm', 190);
    addParameter(p, 'designCoverslipThicknessUm', 190);
    addParameter(p, 'beadDiameterUm', 0.10);
    addParameter(p, 'beadRadiusUm', []);
    addParameter(p, 'airBeadAxialSamples', 7);
    addParameter(p, 'beadSubsamples', [5 5 5]);

    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detectorPixelShape', 'hex');
    addParameter(p, 'detectorSubsamples', 5);
    addParameter(p, 'detectorPitchSampleUm', []);
    addParameter(p, 'detectorHardwarePitchUm', []);
    addParameter(p, 'detectorTotalMagnification', []);
    addParameter(p, 'runPitchProfile', true);
    addParameter(p, 'detectorPitchProfileUm', []);
    addParameter(p, 'pitchProfileScale', 0.55:0.05:1.35);

    addParameter(p, 'fitModes', ...
        {'astig_x','astig_y','coma_x','coma_y', ...
        'trefoil_x','trefoil_y','spherical', ...
        'secondary_astig_x','secondary_astig_y', ...
        'quadrafoil_x','quadrafoil_y'});
    addParameter(p, 'selectedOffsetsUm', []);
    addParameter(p, 'initialMaxIter', 3);
    addParameter(p, 'pitchProfileMaxIter', 2);
    addParameter(p, 'maxIter', 5);
    addParameter(p, 'qeRefitMaxIter', 4);
    addParameter(p, 'refitWithDetectorQE', true);
    addParameter(p, 'interfaceGeometryConfirmed', true);
    addParameter(p, 'runAirInterfaceSelfTest', true);
    addParameter(p, 'knownAberrationCalibration', []);
    addParameter(p, 'requireOddModeSignCertification', true);

    parse(p, varargin{:});
    opts = p.Results;
    opts.outputDir = char(opts.outputDir);
    opts.calibrationFile = char(opts.calibrationFile);
    opts.cacheFile = char(opts.cacheFile);
    opts.ptuReaderFolder = char(opts.ptuReaderFolder);
    opts.darkFile = char(opts.darkFile);
    opts.backgroundMode = lower(char(opts.backgroundMode));
    opts.interfaceRadialWeightMode = char(opts.interfaceRadialWeightMode);
    opts.diffractionModel = char(opts.diffractionModel);
    opts.vectorialPolarizationMode = char(opts.vectorialPolarizationMode);
    opts.excitationPolarizationMode = char(opts.excitationPolarizationMode);
    opts.collectionPolarizationMode = char(opts.collectionPolarizationMode);
    opts.airInterfaceStageMedium = char(opts.airInterfaceStageMedium);
    opts.channelIDs = double(opts.channelIDs(:)).';
    if ischar(opts.fitModes) || isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    end
    opts.fitModes = opts.fitModes(:).';
    if ~isempty(opts.beadDiameterUm)
        opts.beadRadiusUm = double(opts.beadDiameterUm) / 2;
    elseif isempty(opts.beadRadiusUm)
        opts.beadRadiusUm = 0.05;
    end
end

function [coeffs, report] = guardOddModeSigns(coeffs, opts)
    oddModes = {'coma_x','coma_y','tilt_x','tilt_y'};
    presentOdd = intersect(opts.fitModes, oddModes, 'stable');
    certified = oddModeSignCertified(opts.knownAberrationCalibration);
    report = struct();
    report.oddModesInFit = presentOdd;
    report.certified = certified;
    report.certificationSource = describeCertificationSource( ...
        opts.knownAberrationCalibration);
    report.requireCertification = opts.requireOddModeSignCertification;
    report.suppressedModes = {};
    if certified || ~opts.requireOddModeSignCertification
        return;
    end
    suppressed = {};
    for k = 1:numel(presentOdd)
        m = presentOdd{k};
        if isfield(coeffs, m) && coeffs.(m) ~= 0
            coeffs.(m) = 0;
            suppressed{end+1} = m; %#ok<AGROW>
        end
    end
    report.suppressedModes = suppressed;
    if ~isempty(suppressed)
        warning('calibrateWaterObjectiveISMFromBead3D:UncertifiedOddModeSign', ...
            ['Odd modes {%s} were fitted but their SIGN is not certified by a ' ...
             'known-aberration calibration. They have been ZEROED in the saved ' ...
             'calibration so a wrong sign cannot silently propagate to every ' ...
             'downstream bead. Supply a DM-induced known-coma calibration via ' ...
             '''knownAberrationCalibration'' and confirm the sign, or pass ' ...
             '''requireOddModeSignCertification'',false to keep them as fitted.'], ...
            strjoin(suppressed, ', '));
    end
end

function passed = oddModeSignCertified(value)
    passed = false;
    if isempty(value)
        return;
    end
    if ischar(value) || isstring(value)
        fileName = char(value);
        if exist(fileName, 'file') ~= 2
            return;
        end
        S = load(fileName);
        names = fieldnames(S);
        for k = 1:numel(names)
            if isstruct(S.(names{k})) && comaSignReportPasses(S.(names{k}))
                passed = true;
                return;
            end
        end
    elseif isstruct(value)
        passed = comaSignReportPasses(value);
    end
end

function passed = comaSignReportPasses(report)
    passed = false;
    if ~isfield(report, 'table') || ~istable(report.table)
        return;
    end
    T = report.table;
    if ~all(ismember({'mode','signAgreement'}, T.Properties.VariableNames))
        return;
    end
    modes = cellstr(string(T.mode));
    isComa = ismember(modes, {'coma_x','coma_y'});
    if ismember('commandWaves', T.Properties.VariableNames)
        isComa = isComa & isfinite(T.commandWaves) & abs(T.commandWaves) > 0;
    end
    passed = any(isComa) && all(logical(T.signAgreement(isComa)));
end

function s = describeCertificationSource(value)
    if isempty(value)
        s = 'none supplied';
    elseif ischar(value) || isstring(value)
        s = char(value);
    elseif isstruct(value)
        s = 'inline struct report';
    else
        s = class(value);
    end
end

function result = runWaterFullStackFit(stackInput, opts, detectorPitchSampleUm, ...
        initialCoeffs, detectorQE, outputDir, maxIter, writeOutputs, reuseCache, ...
        runSelfTest)
    if nargin < 4 || isempty(initialCoeffs)
        initialCoeffs = struct();
    end
    if nargin < 10 || isempty(runSelfTest)
        runSelfTest = false;
    end
    % The air-interface self-test is invariant across all fits in the
    % calibration, so it is run only on the first (initial) call.
    runSelfTest = runSelfTest && opts.runAirInterfaceSelfTest;
    result = estimateFullStackISMWavefront(stackInput, ...
        'stageZUm', [], ...
        'xyPixelSizeUm', opts.xyPixelSizeUm, ...
        'channelIDs', opts.channelIDs, ...
        'ptuReaderFolder', opts.ptuReaderFolder, ...
        'ptuPhotonsPerChunk', opts.ptuPhotonsPerChunk, ...
        'cacheFile', opts.cacheFile, ...
        'reuseCache', reuseCache, ...
        'useInterPlaneAlignment', false, ...
        'darkFile', opts.darkFile, ...
        'darkScale', opts.darkScale, ...
        'backgroundMode', opts.backgroundMode, ...
        'flatField', detectorQE, ...
        'selectedOffsetsUm', opts.selectedOffsetsUm, ...
        'objectiveNA', opts.objectiveNA, ...
        'objectiveMagnification', opts.objectiveMagnification, ...
        'emissionWavelengthUm', opts.emissionWavelengthUm, ...
        'excitationWavelengthUm', opts.excitationWavelengthUm, ...
        'sampleGeometry', 'airOnGlass', ...
        'immersionRefractiveIndex', opts.immersionRefractiveIndex, ...
        'glassRefractiveIndex', opts.glassRefractiveIndex, ...
        'sampleRefractiveIndex', opts.sampleRefractiveIndex, ...
        'interfaceRadialWeightMode', opts.interfaceRadialWeightMode, ...
        'diffractionModel', opts.diffractionModel, ...
        'vectorialPolarizationMode', opts.vectorialPolarizationMode, ...
        'excitationPolarizationMode', opts.excitationPolarizationMode, ...
        'collectionPolarizationMode', opts.collectionPolarizationMode, ...
        'airInterfaceStageMedium', opts.airInterfaceStageMedium, ...
        'designGlassRefractiveIndex', opts.designGlassRefractiveIndex, ...
        'coverslipThicknessUm', opts.coverslipThicknessUm, ...
        'designCoverslipThicknessUm', opts.designCoverslipThicknessUm, ...
        'beadRadiusUm', opts.beadRadiusUm, ...
        'airBeadAxialSamples', opts.airBeadAxialSamples, ...
        'beadSubsamples', opts.beadSubsamples, ...
        'interfaceGeometryConfirmed', opts.interfaceGeometryConfirmed, ...
        'runAirInterfaceSelfTest', runSelfTest, ...
        'knownAberrationCalibration', opts.knownAberrationCalibration, ...
        'detectorLayout', opts.detectorLayout, ...
        'detectorPixelShape', opts.detectorPixelShape, ...
        'detectorSubsamples', opts.detectorSubsamples, ...
        'estimateDetectorLayout', false, ...
        'detectorPitchSampleUm', detectorPitchSampleUm, ...
        'detectorHardwarePitchUm', opts.detectorHardwarePitchUm, ...
        'detectorTotalMagnification', opts.detectorTotalMagnification, ...
        'fitModes', opts.fitModes, ...
        'initialCoeffs', initialCoeffs, ...
        'fitXY', false, ...
        'maxIter', maxIter, ...
        'runProfiles', false, ...
        'nBootstrap', 0, ...
        'runSignProfiles', false, ...
        'outputDir', outputDir, ...
        'writeOutputs', writeOutputs, ...
        'verbose', opts.verbose);
end

function [bestPitchUm, profile] = chooseDetectorPitch(stackInput, opts, nominalPitchUm, initialCoeffs)
    profile = table();
    bestPitchUm = nominalPitchUm;
    if ~opts.runPitchProfile
        return;
    end
    if ~isempty(opts.detectorPitchProfileUm)
        candidates = double(opts.detectorPitchProfileUm(:)).';
    else
        candidates = nominalPitchUm * double(opts.pitchProfileScale(:)).';
    end
    candidates = unique(candidates(isfinite(candidates) & candidates > 0));
    if isempty(candidates)
        return;
    end

    deviance = nan(numel(candidates),1);
    residualNorm = nan(numel(candidates),1);
    for k = 1:numel(candidates)
        if opts.verbose
            fprintf('[water-objective calibration] pitch profile %d/%d: %.5g um\n', ...
                k, numel(candidates), candidates(k));
        end
        one = runWaterFullStackFit(stackInput, opts, candidates(k), ...
            initialCoeffs, [], fullfile(opts.outputDir, 'pitch_profile'), ...
            opts.pitchProfileMaxIter, false, true, false);
        deviance(k) = one.fit.deviance;
        residualNorm(k) = one.fit.residualNorm;
    end
    profile = table(candidates(:), deviance, residualNorm, ...
        'VariableNames', {'detectorPitchSampleUm','deviance','residualNorm'});
    [~, idx] = min(deviance);
    bestPitchUm = candidates(idx);
end

function [qe, T] = estimateRelativeDetectorQE(result)
    data = result.data;
    signal = backgroundCorrectedCounts(data);
    measured = squeeze(sum(sum(sum(signal,1),2),4));
    model = max(double(result.fit.model), 0);
    predicted = squeeze(sum(sum(sum(model,1),2),4));
    ratio = measured(:) ./ max(predicted(:), realmin);
    valid = isfinite(ratio) & ratio > 0;
    if ~any(valid)
        qe = ones(numel(data.channelIDs),1);
    else
        fill = median(ratio(valid));
        ratio(~valid) = fill;
        qe = ratio / mean(ratio);
    end
    channelID = data.channelIDs(:);
    detectorIndex = (1:numel(channelID)).';
    T = table(detectorIndex, channelID, measured(:), predicted(:), ratio(:), qe(:), ...
        'VariableNames', {'detectorIndex','channelID','measuredSignalCounts', ...
        'modelUnitCounts','measuredOverModel','relativeQE'});
end

function sigma = estimateFocusPsfSigma(result)
    data = result.data;
    signal = backgroundCorrectedCounts(data);
    [~, focusIdx] = min(abs(data.stageZUm - result.fit.estZ0Um));
    measuredImage = sum(signal(:,:,:,focusIdx), 3);
    modelImage = sum(max(double(result.fit.model(:,:,:,focusIdx)),0), 3);
    measuredXY = weightedImageSigma(measuredImage, data.xUm, data.yUm);
    modelXY = weightedImageSigma(modelImage, data.xUm, data.yUm);

    measuredTrace = max(double(data.axialSignalCounts(:)), 0);
    modelTrace = squeeze(sum(sum(sum(max(double(result.fit.model),0),1),2),3));
    measuredZ = weightedVectorSigma(measuredTrace(:), data.stageZUm(:));
    modelZ = weightedVectorSigma(modelTrace(:), data.stageZUm(:));

    sigma = struct();
    sigma.focusPlaneIndex = focusIdx;
    sigma.focusStageZUm = data.stageZUm(focusIdx);
    sigma.fitZ0Um = result.fit.estZ0Um;
    sigma.measuredSigmaXUm = measuredXY.sigmaX;
    sigma.measuredSigmaYUm = measuredXY.sigmaY;
    sigma.measuredSigmaXYUm = measuredXY.sigmaMean;
    sigma.modelSigmaXUm = modelXY.sigmaX;
    sigma.modelSigmaYUm = modelXY.sigmaY;
    sigma.modelSigmaXYUm = modelXY.sigmaMean;
    sigma.measuredSigmaZUm = measuredZ.sigma;
    sigma.modelSigmaZUm = modelZ.sigma;
    sigma.empiricalModelBlurSigmaXYUm = sqrt(max( ...
        sigma.measuredSigmaXYUm.^2 - sigma.modelSigmaXYUm.^2, 0));
    sigma.empiricalModelBlurSigmaZUm = sqrt(max( ...
        sigma.measuredSigmaZUm.^2 - sigma.modelSigmaZUm.^2, 0));
end

function signal = backgroundCorrectedCounts(data)
    raw = double(data.rawCounts);
    bg = double(data.backgroundPerPixel);
    bg = repmat(bg, size(raw,1), size(raw,2), 1, 1);
    signal = max(raw - bg, 0);
end

function s = weightedImageSigma(img, x, y)
    img = max(double(img), 0);
    [X,Y] = meshgrid(double(x(:)).', double(y(:)));
    w = img(:);
    total = sum(w);
    if total <= 0 || ~isfinite(total)
        s = struct('sigmaX', NaN, 'sigmaY', NaN, 'sigmaMean', NaN);
        return;
    end
    x0 = sum(w .* X(:)) / total;
    y0 = sum(w .* Y(:)) / total;
    sx = sqrt(sum(w .* (X(:)-x0).^2) / total);
    sy = sqrt(sum(w .* (Y(:)-y0).^2) / total);
    s = struct('sigmaX', sx, 'sigmaY', sy, ...
        'sigmaMean', sqrt(0.5*(sx.^2 + sy.^2)));
end

function s = weightedVectorSigma(values, axisUm)
    values = max(double(values(:)), 0);
    axisUm = double(axisUm(:));
    total = sum(values);
    if total <= 0 || ~isfinite(total)
        s = struct('centerUm', NaN, 'sigma', NaN);
        return;
    end
    c = sum(values .* axisUm) / total;
    s = struct('centerUm', c, ...
        'sigma', sqrt(sum(values .* (axisUm-c).^2) / total));
end

function T = wavefrontTable(calibration)
    modes = calibration.wavefront.fitModes(:);
    waves = zeros(numel(modes),1);
    for k = 1:numel(modes)
        if isfield(calibration.wavefront.coeffs, modes{k})
            waves(k) = calibration.wavefront.coeffs.(modes{k});
        end
    end
    opdNmAtLambdaRef = waves * calibration.wavefront.lambdaRefUm * 1000;
    T = table(modes, waves, opdNmAtLambdaRef);
end

function T = psfSigmaTable(s)
    name = fieldnames(s);
    value = zeros(numel(name),1);
    for k = 1:numel(name)
        v = s.(name{k});
        if isnumeric(v) && isscalar(v)
            value(k) = double(v);
        else
            value(k) = NaN;
        end
    end
    T = table(name, value);
end

function pitch = nominalDetectorPitchUm(opts)
    pitch = opts.detectorPitchSampleUm;
    if isnumeric(pitch) && isscalar(pitch) && isfinite(pitch) && pitch > 0
        pitch = double(pitch);
        return;
    end
    mag = opts.detectorTotalMagnification;
    physical = opts.detectorHardwarePitchUm;
    if isnumeric(physical) && isscalar(physical) && isfinite(physical) && ...
            physical > 0 && isnumeric(mag) && isscalar(mag) && ...
            isfinite(mag) && mag > 0
        pitch = double(physical) / double(mag);
        return;
    end
    sim = defaultParams();
    pitch = sim.detPitch;
end

function mag = detectorTotalMagnificationFromPitchLocal(pitchUm, hardwarePitchUm, fallback)
    mag = NaN;
    if isnumeric(hardwarePitchUm) && isscalar(hardwarePitchUm) && ...
            isfinite(hardwarePitchUm) && hardwarePitchUm > 0 && ...
            isnumeric(pitchUm) && isscalar(pitchUm) && isfinite(pitchUm) && pitchUm > 0
        mag = double(hardwarePitchUm) / double(pitchUm);
    elseif isnumeric(fallback) && isscalar(fallback) && isfinite(fallback) && fallback > 0
        mag = double(fallback);
    end
end

function outputDir = resolveCalibrationOutputDir(stackInput, outputDir)
    if ~isempty(outputDir)
        return;
    end
    root = fileparts(fileparts(mfilename('fullpath')));
    [~, stem] = fileparts(char(stackInput));
    stem = regexprep(stem, '_frame_alignment$', '');
    outputDir = fullfile(root, 'output_matlab', ...
        'water_objective_calibration', sanitizeName(stem));
end

function name = sanitizeName(name)
    name = regexprep(char(name), '[^A-Za-z0-9_\-]+', '_');
    if isempty(name)
        name = 'calibration';
    end
end

function addRequiredPaths()
    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    addpath(fileparts(thisDir));
end

function folder = defaultPtuReaderFolderLocal()
    thisDir = fileparts(mfilename('fullpath'));
    folder = fullfile(fileparts(fileparts(thisDir)), 'Luminosa_FLIM_FCS');
end
