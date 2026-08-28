%RUN_BATCH_IMMUNE_CELL_MIET Fixed-SLB-lifetime, regularized-amplitude analysis.
%
% The outside-SLB pixels determine tau_SLB once per acquisition. That
% lifetime is then fixed while each cell-membrane pixel is compared between
% SLB-only and SLB + one longer membrane lifetime. A smooth SLB-amplitude
% surface is fitted outside the cell, extrapolated underneath it, and used as
% a Gaussian prior rather than a fixed amplitude. Native, sliding 2x2 and
% sliding 4x4 results are all analyzed and written as PNG figures.
% dataRoot = 'D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_20260813_1';
dataRoot = 'D:\Luminosa\Data\260813\4deg_Jurkat_CD58_memglow';

% Each acquisition needs its own inverse calibration because tau_SLB, the IRF
% and the lifetime grid can vary between files. The first run generates and
% caches that calibration beside the result; subsequent runs reuse it after
% strict compatibility checks. Raw lifetime maps remain unchanged; separate
% corrected maps and reliability masks are written beside each analysis.
applyLongLifetimeBiasCorrection = true;
biasCalibrationCfg = struct('showFigures', false, 'writeFigures', true);

cfg = struct();
cfg.scanPlanes = {'XY'};          % Cross-sections are catalogued, not fitted.
cfg.minDwellMs = 0.4;             % Excludes fast overview/test scans.
cfg.minCompletionFraction = 0.95; % Excludes interrupted raster scans.
cfg.minRecords = 1e6;
cfg.resume = true;                % Reuse valid per-acquisition MAT results.
cfg.overwrite = false;
cfg.dryRun = false;               % Set true for a fast header-only audit.
cfg.continueOnError = true;
cfg.showFigures = false;          % PNG figures are still saved.
cfg.saveTcspcPix = true;          % Detector-summed, reassigned uint16 cube.
cfg.qc = struct('maxSlbClippedFraction', 0.10);

% cfg.outputDir was silently ignored here - the batch driver derives each
% acquisition's results folder from that acquisition's own location and never
% read this field, which is why repeated runs overwrote one another. Use these
% instead:
cfg.resultsFolderName = ...
    'immune_cell_MIET_slb_regularized_one_membrane_grid48'; % beside each PTU
cfg.versionResults = true;        % Results go in a per-configuration
                                  % subfolder, so changing a setting can no
                                  % longer overwrite an earlier analysis.
cfg.runName = '';                 % '' derives the subfolder from a hash of the
                                  % analysis settings: an UNCHANGED config
                                  % resolves to the same folder and resumes, a
                                  % CHANGED config gets a new folder and leaves
                                  % the old results in place. Set a string here
                                  % to name a run yourself.
cfg.batchOutputDir = '';          % '' puts the batch summary under
                                  % <dataRoot>\immune_cell_MIET_batch\<runName>

pipelineCfg = struct();
pipelineCfg.excitationPulseIndex = 2;
pipelineCfg.excitationNm = 640;
pipelineCfg.ismWavelengthNm = 690;
pipelineCfg.tcspcBinNs = 0.16;    % Exact multiple of both 40-ps and 20-ps data.
pipelineCfg.minPhotonsPerPixel = 3;
pipelineCfg.lifetimeRangeNs = [0.5 1.8];
pipelineCfg.showWaitbar = false;
% Measured on this machine (T1000 Max-Q): the GPU path is ~2x SLOWER than the
% CPU for these matrix shapes, and its single precision flips the selected
% model in ~4.7% of pixels. Keep it off. See
% test_flim_bayes_fixed_slb_gpu_grouping.
pipelineCfg.useGPU = false;
pipelineCfg.bayesMinPhotons = 10;
pipelineCfg.maxMembraneStates = 1; % SLB-only versus SLB + one membrane state.
pipelineCfg.slbAmplitudeMode = 'regularized'; % Smooth outside-SLB spatial prior.
pipelineCfg.slbAmplitudeScale = 1;
pipelineCfg.slbAmplitudePriorRelativeStdFloor = 0.10;
pipelineCfg.slbAmplitudePriorRelativeStdCeiling = 0.50;
pipelineCfg.slbAmplitudeSurfaceClipQuantiles = [0.05 0.95];
pipelineCfg.maxSlbClippedFraction = 0.10;
pipelineCfg.requireValidSegmentation = true;
pipelineCfg.componentMaps = struct( ...
    'enabled', true, ...
    'posteriorThreshold', [0.8 0.95], ...
    'minExpectedPhotons', [10 10], ...
    'probabilityContourLevels', [0.5 0.7 0.9]);
% mergeStruct is shallow, so the whole bayes struct must be restated here.
% The BOUNDS are deliberately left unset. flim_bayes_fixed_slb derives them
% from the fitted SLB lifetime as
%   [max(1.15*tauSlb, tauSlb + 2*dt, 0.05), min(0.8*period, max(5, 4*tauSlb))]
% which for a fitted tauSlb of 0.358 ns at 0.16 ns bins gives a lower bound of
% 0.678 ns. That guard keeps membrane states at least two TCSPC bins away from
% the SLB. Forcing a lower bound (an earlier version used 0.4 ns) puts grid
% states within two bins of the SLB, where they are nearly degenerate with it.
% It is also unphysical for MIET,
% where the SLB is the most quenched layer so membrane lifetimes must be longer.
pipelineCfg.bayes = struct('batchSize', 2048, 'includeBackground', true, ...
    'signalGrid', [0.25 0.5 0.75 1], 'membraneTauCount', 48, ...
    'fractionStep', 0.2, 'minimumMembraneFraction', 0.1, ...
    'slbCountPriorNodes', 5, 'slbCountRelTol', 0.0025);
pipelineCfg.spatialBinning = struct('enabled', true, ...
    'binSize', [2 2], 'step', [1 1]);
pipelineCfg.spatialBinning4x4 = struct('enabled', true, ...
    'binSize', [4 4], 'step', [1 1]);
pipelineCfg.segmentation = struct( ...
    'smoothSigma', 2.0, ...
    'minLifetimeContrastNs', 0.05, ...
    'cellSeparationSigma', 1.5, ...
    'growthMinLifetimeContrastNs', 0.025, ...
    'growthSeparationSigma', 0.75, ...
    'slbLowPopulationQuantile', 0.45, ...
    'minPlausibleCellFraction', 0.005, ...
    'maxCellComponents', 1, ...
    'fallbackUseIntensity', false, ...
    'minCellArea', 500, ...
    'minCellCoreArea', 20, ...
    'fallbackMinContrastNs', 0.025, ...
    'fallbackCloseRadius', 7, ...
    'fallbackMinArea', 250);
cfg.pipeline = pipelineCfg;

[summary, batchInfo] = batch_immune_cell_MIET(dataRoot, cfg);

if applyLongLifetimeBiasCorrection
    correctedCount = 0;
    correctionFailureCount = 0;
    for row = 1:height(summary)
        analysisMat = summary.analysisMat(row);
        if iscell(analysisMat)
            analysisMat = analysisMat{1};
        end
        analysisMat = char(analysisMat);
        if ~isfile(analysisMat)
            continue;
        end
        try
            calibrationDir = fullfile(fileparts(analysisMat), ...
                'long_lifetime_bias');
            calibrationMat = fullfile(calibrationDir, ...
                'long_lifetime_bias_analysis.mat');
            if ~isfile(calibrationMat)
                calibration = immune_cell_MIET_long_lifetime_bias( ...
                    analysisMat, biasCalibrationCfg);
                calibrationMat = calibration.outputFiles.mat;
            end
            try
                immune_cell_MIET_apply_long_lifetime_bias_correction( ...
                    analysisMat, calibrationMat);
            catch staleCalibrationError
                regeneratableIds = { ...
                    'immune_cell_MIET_apply_long_lifetime_bias_correction:CalibrationSchema', ...
                    'immune_cell_MIET_correct_long_lifetime_bias:CalibrationSchema', ...
                    'immune_cell_MIET_apply_long_lifetime_bias_correction:TauSlbMismatch', ...
                    'immune_cell_MIET_apply_long_lifetime_bias_correction:TimingMismatch', ...
                    'immune_cell_MIET_apply_long_lifetime_bias_correction:GridMismatch', ...
                    'immune_cell_MIET_apply_long_lifetime_bias_correction:IrfMismatch', ...
                    'immune_cell_MIET_apply_long_lifetime_bias_correction:PriorWidthMismatch', ...
                    'immune_cell_MIET_apply_long_lifetime_bias_correction:PriorNodeMismatch'};
                if ~ismember(staleCalibrationError.identifier, regeneratableIds)
                    rethrow(staleCalibrationError);
                end
                warning('run_batch_immune_cell_MIET:RegenerateBiasCalibration', ...
                    ['Cached bias calibration is incompatible with %s and ' ...
                     'will be regenerated: %s'], analysisMat, ...
                    staleCalibrationError.message);
                calibration = immune_cell_MIET_long_lifetime_bias( ...
                    analysisMat, biasCalibrationCfg);
                calibrationMat = calibration.outputFiles.mat;
                immune_cell_MIET_apply_long_lifetime_bias_correction( ...
                    analysisMat, calibrationMat);
            end
            correctedCount = correctedCount + 1;
        catch correctionError
            correctionFailureCount = correctionFailureCount + 1;
            warning('run_batch_immune_cell_MIET:BiasCorrectionFailed', ...
                'Bias correction failed for %s: %s', ...
                analysisMat, correctionError.message);
        end
    end
    fprintf(['run_batch_immune_cell_MIET: bias correction complete for %d ' ...
        'results; %d rejected by compatibility/reliability checks.\n'], ...
        correctedCount, correctionFailureCount);
end

fprintf('\nrun_batch_immune_cell_MIET: RESULTS\n');
fprintf('  Batch summary folder: %s\n', batchInfo.outputDir);
fprintf(['  Per-acquisition figures: see the analysisDir column in\n' ...
    '  %s\n'], batchInfo.csvFile);
