function biasCorrection = immune_cell_MIET_apply_long_lifetime_bias_correction( ...
        resultMatFile, calibrationMatFile, cfg)
%IMMUNE_CELL_MIET_APPLY_LONG_LIFETIME_BIAS_CORRECTION Correct saved MIET maps.
%
% biasCorrection = immune_cell_MIET_apply_long_lifetime_bias_correction( ...
%     resultMatFile, calibrationMatFile, cfg)
%
% Applies immune_cell_MIET_correct_long_lifetime_bias to the native, 2x2 and
% 4x4 one-membrane results stored by immune_cell_MIET. Raw results are not
% modified. A separate MAT and one diagnostic figure per available binning
% size are written beside the source analysis.
%
% The calibration must match the fixed SLB lifetime, TCSPC bin width,
% lifetime grid and regularized SLB-count estimator. Free-amplitude or hard-
% fixed-amplitude results are rejected because their bias is different.
%
% Optional cfg fields
%   outputDir                 directory beside result MAT by default
%   outputMatFile             immune_cell_MIET_long_lifetime_bias_correction.mat
%   writeFigures              true
%   showFigures               false
%   corrector                 options passed to the inverse corrector
%   maximumIrfL1Difference    0.05 after normalisation
%   maximumPriorRelativeStdDifference 0.15

    if nargin < 3 || isempty(cfg)
        cfg = struct();
    end
    resultMatFile = char(resultMatFile);
    calibrationMatFile = char(calibrationMatFile);
    if ~isfile(resultMatFile)
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:ResultFile', ...
            'Result MAT does not exist: %s', resultMatFile);
    end
    if ~isfile(calibrationMatFile)
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:CalibrationFile', ...
            'Calibration MAT does not exist: %s', calibrationMatFile);
    end
    cfg = fillDefaults(cfg, resultMatFile);
    if ~isfolder(cfg.outputDir)
        mkdir(cfg.outputDir);
    end

    savedResult = load(resultMatFile, 'result');
    if ~isfield(savedResult, 'result') || ~isstruct(savedResult.result)
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:ResultSchema', ...
            'Result MAT must contain a result struct.');
    end
    savedCalibration = load(calibrationMatFile, 'analysis');
    if ~isfield(savedCalibration, 'analysis') || ...
            ~isstruct(savedCalibration.analysis)
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:CalibrationSchema', ...
            'Calibration MAT must contain an analysis struct.');
    end
    result = savedResult.result;
    calibration = savedCalibration.analysis;
    validateCompatibility(result, calibration, cfg);

    correctorCfg = cfg.corrector;
    correctorCfg.tauSlbNs = double(result.slbReference.fixedLifetimeNs);
    stages = struct();
    stages.native = correctStage(result.bayesian, ...
        result.slbReference.amplitudeConstraint, calibration, correctorCfg, ...
        'native pixels');
    stages.sliding2x2 = unavailableStage('2x2 sliding TCSPC was not available');
    if isCompleteStage(result, 'spatialBinning')
        stages.sliding2x2 = correctStage(result.spatialBinning.bayesian, ...
            result.spatialBinning.bayesian.fixedSlbAmplitude, ...
            calibration, correctorCfg, '2x2 sliding TCSPC, one-pixel step');
    end
    stages.sliding4x4 = unavailableStage('4x4 sliding TCSPC was not available');
    if isCompleteStage(result, 'spatialBinning4x4')
        stages.sliding4x4 = correctStage(result.spatialBinning4x4.bayesian, ...
            result.spatialBinning4x4.bayesian.fixedSlbAmplitude, ...
            calibration, correctorCfg, '4x4 sliding TCSPC, one-pixel step');
    end

    biasCorrection = struct();
    biasCorrection.method = ['post-fit conservative inverse calibration; ' ...
        'raw source result left unchanged'];
    biasCorrection.sourceResultMat = resultMatFile;
    biasCorrection.calibrationMat = calibrationMatFile;
    biasCorrection.fixedSlbLifetimeNs = result.slbReference.fixedLifetimeNs;
    biasCorrection.config = cfg;
    biasCorrection.stages = stages;
    biasCorrection.outputFiles = struct('mat', cfg.outputMatFile, ...
        'figures', {{}});

    % Numerical maps are saved before graphics so headless rendering failures
    % cannot lose the correction or reliability masks.
    save(cfg.outputMatFile, 'biasCorrection', '-v7.3');
    figureFiles = {};
    if cfg.writeFigures
        try
            stageNames = fieldnames(stages);
            for index = 1:numel(stageNames)
                name = stageNames{index};
                if ~stages.(name).available
                    continue;
                end
                figureFile = fullfile(cfg.outputDir, sprintf( ...
                    'immune_cell_MIET_%s_long_lifetime_bias_correction.png', name));
                writeCorrectionFigure(stages.(name), figureFile, cfg.showFigures);
                figureFiles{end + 1} = figureFile; %#ok<AGROW>
            end
            biasCorrection.outputFiles.figures = figureFiles;
            save(cfg.outputMatFile, 'biasCorrection', '-v7.3');
        catch figureError
            warning(['immune_cell_MIET_apply_long_lifetime_bias_correction:' ...
                'FigureExport'], ['Numerical correction was saved, but figure ' ...
                'export failed: %s'], figureError.message);
        end
    end

    fprintf('immune_cell_MIET_apply_long_lifetime_bias_correction: saved %s\n', ...
        cfg.outputMatFile);
    printStageSummary(stages);
end

function cfg = fillDefaults(cfg, resultMatFile)
    outputDir = fileparts(resultMatFile);
    defaults = struct('outputDir', outputDir, ...
        'outputMatFile', fullfile(outputDir, ...
            'immune_cell_MIET_long_lifetime_bias_correction.mat'), ...
        'writeFigures', true, 'showFigures', false, 'corrector', struct(), ...
        'maximumIrfL1Difference', 0.05, ...
        'maximumPriorRelativeStdDifference', 0.15);
    cfg = mergeStruct(defaults, cfg);
    cfg.outputDir = char(cfg.outputDir);
    cfg.outputMatFile = char(cfg.outputMatFile);
    validateattributes(cfg.writeFigures, {'numeric','logical'}, {'scalar'});
    validateattributes(cfg.showFigures, {'numeric','logical'}, {'scalar'});
    cfg.writeFigures = logical(cfg.writeFigures);
    cfg.showFigures = logical(cfg.showFigures);
    validateattributes(cfg.maximumIrfL1Difference, {'numeric'}, ...
        {'real','finite','scalar','nonnegative'});
    validateattributes(cfg.maximumPriorRelativeStdDifference, {'numeric'}, ...
        {'real','finite','scalar','nonnegative'});
    if ~(isstruct(cfg.corrector) && isscalar(cfg.corrector))
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:Config', ...
            'cfg.corrector must be a scalar struct.');
    end
end

function validateCompatibility(result, calibration, cfg)
    requiredResult = {'slbReference','channel','irf','bayesian'};
    for index = 1:numel(requiredResult)
        if ~isfield(result, requiredResult{index})
            error('immune_cell_MIET_apply_long_lifetime_bias_correction:ResultSchema', ...
                'Result is missing result.%s.', requiredResult{index});
        end
    end
    if ~isfield(calibration, 'instrument') || ...
            ~isfield(calibration, 'estimator')
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:CalibrationSchema', ...
            'Calibration lacks instrument or estimator metadata.');
    end
    if isfield(result, 'config') && isfield(result.config, 'maxMembraneStates') && ...
            result.config.maxMembraneStates ~= 1
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:ModelOrder', ...
            'This calibration corrects only fixed SLB + one membrane state.');
    end
    compact = result.bayesian.compact;
    if ~isfield(compact, 'fixedSlbCountPrior') || ...
            ~compact.fixedSlbCountPrior.countMarginalised
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:AmplitudeMode', ...
            ['The saved result did not marginalise a regularized spatial SLB-' ...
             'count prior. Generate a matched calibration for its estimator.']);
    end
    calibrationTau = double(calibration.instrument.tauSlbNs);
    resultTau = double(result.slbReference.fixedLifetimeNs);
    if abs(resultTau - calibrationTau) > max(0.01, 0.05 * calibrationTau)
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:TauSlbMismatch', ...
            'Result tau_SLB %.4g ns differs from calibration %.4g ns.', ...
            resultTau, calibrationTau);
    end
    if abs(double(result.channel.dtNs) - double(calibration.instrument.dtNs)) > 1e-9
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:TimingMismatch', ...
            'Result and calibration TCSPC bin widths differ.');
    end
    resultGrid = double(compact.membraneTauGridNs(:));
    calibrationGrid = double(calibration.estimator.lifetimeGridNs(:));
    if numel(resultGrid) ~= numel(calibrationGrid) || ...
            any(abs(resultGrid - calibrationGrid) > 1e-9)
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:GridMismatch', ...
            'Result and calibration lifetime grids differ.');
    end
    resultIrf = normalise(double(result.irf.curve(:)));
    calibrationIrf = normalise(double(calibration.instrument.irf(:)));
    if numel(resultIrf) ~= numel(calibrationIrf) || ...
            sum(abs(resultIrf - calibrationIrf)) > cfg.maximumIrfL1Difference
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:IrfMismatch', ...
            'Result and calibration IRFs differ beyond the configured tolerance.');
    end
    prior = compact.fixedSlbCountPrior;
    target = double(prior.targetPhotonCount(:));
    width = double(prior.photonCountStd(:));
    valid = isfinite(target) & target > 0 & isfinite(width) & width > 0;
    actualRelativeStd = median(width(valid) ./ target(valid));
    requestedRelativeStd = double(calibration.config.slbPriorRelativeStd);
    if isempty(actualRelativeStd) || ~isfinite(actualRelativeStd) || ...
            abs(actualRelativeStd - requestedRelativeStd) > ...
            cfg.maximumPriorRelativeStdDifference
        error(['immune_cell_MIET_apply_long_lifetime_bias_correction:' ...
            'PriorWidthMismatch'], ['Median result SLB-prior relative width ' ...
            'does not match the bias calibration.']);
    end
    if double(prior.priorNodes) ~= double(calibration.config.slbCountPriorNodes)
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:PriorNodeMismatch', ...
            'Result and calibration SLB-prior quadrature node counts differ.');
    end
end

function stage = correctStage(bayesian, amplitude, calibration, cfg, label)
    layers = bayesian.orderedComponents;
    if isfield(bayesian, 'pixelLinearIndices')
        indices = double(bayesian.pixelLinearIndices(:));
    elseif isfield(bayesian, 'membranePixelIndices')
        indices = double(bayesian.membranePixelIndices(:));
    else
        error(['immune_cell_MIET_apply_long_lifetime_bias_correction:' ...
            'PixelIndices'], 'Bayesian result lacks compact pixel indices.');
    end
    imageSize = layers.imageSize;
    slbMap = nan(imageSize);
    if isfield(amplitude, 'priorMeanPhotonCountMapNative') && ...
            isequal(size(amplitude.priorMeanPhotonCountMapNative), imageSize)
        slbMap = double(amplitude.priorMeanPhotonCountMapNative);
    elseif isfield(amplitude, 'priorMeanPhotonCount') && ...
            numel(amplitude.priorMeanPhotonCount) == numel(indices)
        slbMap(indices) = double(amplitude.priorMeanPhotonCount(:));
    else
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:SlbCountMap', ...
            'Could not reconstruct the %s SLB-prior count map.', label);
    end
    rawMap = double(layers.ungated.secondLifetimeNs);
    longMap = double(layers.expectedPhotonCount.secondComponent);
    inverse = immune_cell_MIET_correct_long_lifetime_bias( ...
        rawMap, slbMap, longMap, calibration, cfg);
    inputDisplayMask = logical(layers.masks.secondDisplay) & isfinite(rawMap);
    reliableDisplayMask = inputDisplayMask & inverse.reliableMask;
    displayCorrected = double(inverse.correctedLifetimeNs);
    displayCorrected(~reliableDisplayMask) = NaN;
    inverse.inputDisplayMask = inputDisplayMask;
    inverse.reliableDisplayMask = reliableDisplayMask;
    inverse.displayCorrectedLifetimeNs = single(displayCorrected);
    inverse.summary.inputDisplayPixelCount = nnz(inputDisplayMask);
    inverse.summary.reliableDisplayPixelCount = nnz(reliableDisplayMask);
    inverse.summary.reliableFractionOfDisplayed = ...
        nnz(reliableDisplayMask) / max(nnz(inputDisplayMask), 1);
    inverse.summary.rawDisplayMedianNs = finiteMedian(rawMap(inputDisplayMask));
    inverse.summary.correctedDisplayMedianNs = ...
        finiteMedian(displayCorrected(reliableDisplayMask));
    stage = struct('available', true, 'label', label, ...
        'windowSize', stageWindowSize(amplitude), 'inverse', inverse);
end

function stage = unavailableStage(message)
    stage = struct('available', false, 'message', message);
end

function tf = isCompleteStage(result, fieldName)
    tf = isfield(result, fieldName) && isstruct(result.(fieldName)) && ...
        isfield(result.(fieldName), 'status') && ...
        strcmp(result.(fieldName).status, 'complete') && ...
        isfield(result.(fieldName), 'bayesian');
end

function windowSize = stageWindowSize(amplitude)
    windowSize = [1 1];
    if isfield(amplitude, 'windowSize')
        windowSize = double(amplitude.windowSize(:)).';
    end
end

function writeCorrectionFigure(stage, outputFile, showFigure)
    visibility = 'off';
    if showFigure
        visibility = 'on';
    end
    inverse = stage.inverse;
    displayMask = inverse.inputDisplayMask;
    reliableMask = inverse.reliableDisplayMask;
    raw = double(inverse.rawLifetimeNs);
    corrected = double(inverse.displayCorrectedLifetimeNs);
    applied = double(inverse.appliedCorrectionNs);
    lifetimeValues = [raw(displayMask); corrected(reliableMask)];
    limits = robustLimits(lifetimeValues, [0 1]);
    correctionValues = applied(reliableMask);
    correctionLimit = finitePercentile(abs(correctionValues), 0.99);
    correctionLimit = max(correctionLimit, 0.01);

    h = figure('Color', 'w', 'Visible', visibility, ...
        'Position', [60 60 1500 820]);
    layout = tiledlayout(h, 2, 2, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    showPanel(nexttile(layout), raw, displayMask, turbo(256), limits, ...
        'Raw posterior-mean lifetime', 'Lifetime [ns]');
    showPanel(nexttile(layout), corrected, reliableMask, turbo(256), limits, ...
        'Bias-corrected lifetime (reliable only)', 'Lifetime [ns]');
    showPanel(nexttile(layout), applied, reliableMask, blueWhiteRed(256), ...
        [-correctionLimit correctionLimit], 'Applied correction', 'Correction [ns]');
    reason = double(inverse.reasonCode);
    reason(~displayMask) = NaN;
    showPanel(nexttile(layout), reason, displayMask, lines(11), [0 10], ...
        'Reliability reason code (1 = accepted)', 'Reason code');
    title(layout, sprintf(['Long-lifetime inverse bias calibration | %s | ' ...
        '%d/%d displayed pixels reliable'], stage.label, ...
        inverse.summary.reliableDisplayPixelCount, ...
        inverse.summary.inputDisplayPixelCount), 'FontWeight', 'bold');
    subtitle(layout, ['Raw values are retained in the MAT; correction is NaN ' ...
        'where the simulated inverse is non-identifiable.']);
    exportgraphics(h, outputFile, 'Resolution', 250);
    if ~showFigure
        close(h);
    end
end

function showPanel(ax, data, mask, colourMap, limits, titleText, colourLabel)
    object = imagesc(ax, data.');
    axis(ax, 'image', 'off');
    set(ax, 'YDir', 'normal', 'Color', [0.015 0.015 0.025]);
    set(object, 'AlphaData', mask.');
    colormap(ax, colourMap);
    clim(ax, limits);
    colourBar = colorbar(ax);
    colourBar.Label.String = colourLabel;
    title(ax, titleText);
end

function limits = robustLimits(values, fallback)
    values = double(values(isfinite(values)));
    if isempty(values)
        limits = fallback;
        return;
    end
    limits = [finitePercentile(values, 0.01), ...
        finitePercentile(values, 0.99)];
    if limits(2) <= limits(1)
        pad = max(abs(limits(1)) * 0.1, 0.01);
        limits = limits + [-pad pad];
    end
end

function value = finitePercentile(data, probability)
    data = sort(double(data(isfinite(data))));
    if isempty(data)
        value = NaN;
        return;
    end
    position = 1 + (numel(data) - 1) * probability;
    low = floor(position);
    high = ceil(position);
    fraction = position - low;
    value = data(low) * (1 - fraction) + data(high) * fraction;
end

function value = finiteMedian(data)
    data = double(data(isfinite(data)));
    if isempty(data)
        value = NaN;
    else
        value = median(data);
    end
end

function curve = normalise(curve)
    curve = max(double(curve(:)), 0);
    curve = curve / sum(curve);
end

function map = blueWhiteRed(count)
    x = linspace(-1, 1, count).';
    red = interp1([-1 0 1], [0.10 1.00 0.80], x);
    green = interp1([-1 0 1], [0.30 1.00 0.10], x);
    blue = interp1([-1 0 1], [0.80 1.00 0.10], x);
    map = min(max([red green blue], 0), 1);
end

function printStageSummary(stages)
    names = fieldnames(stages);
    for index = 1:numel(names)
        stage = stages.(names{index});
        if ~stage.available
            continue;
        end
        summary = stage.inverse.summary;
        fprintf('  %s: %d/%d displayed pixels reliable (%.1f%%).\n', ...
            stage.label, summary.reliableDisplayPixelCount, ...
            summary.inputDisplayPixelCount, ...
            100 * summary.reliableFractionOfDisplayed);
    end
end

function out = mergeStruct(base, override)
    if ~(isstruct(override) && isscalar(override))
        error('immune_cell_MIET_apply_long_lifetime_bias_correction:Config', ...
            'cfg must be a scalar struct.');
    end
    out = base;
    names = fieldnames(override);
    for index = 1:numel(names)
        out.(names{index}) = override.(names{index});
    end
end
