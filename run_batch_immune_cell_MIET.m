%RUN_BATCH_IMMUNE_CELL_MIET Analyze all suitable cell acquisitions in a folder.
dataRoot = 'D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_20260813_1';

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
pipelineCfg.fixSlbAmplitude = true;
pipelineCfg.slbAmplitudeScale = 1;
pipelineCfg.maxSlbClippedFraction = 0.10;
pipelineCfg.requireValidSegmentation = true;
pipelineCfg.componentMaps = struct( ...
    'enabled', true, ...
    'posteriorThreshold', [0.8 0.95], ...
    'minExpectedPhotons', [10 10], ...
    'probabilityContourLevels', [0.5 0.7 0.9]);
% mergeStruct is shallow, so the whole bayes struct must be restated here.
% slbCountRelTol groups photon totals so neighbouring totals share one grid.
% 0.0025 is the largest tolerance that leaves model selection identical
% (100% agreement, max |dP| 0.0035). 0.005 shifts it, 0.02 moves dP by 0.07.
% The BOUNDS are deliberately left unset. flim_bayes_fixed_slb derives them
% from the fitted SLB lifetime as
%   [max(1.15*tauSlb, tauSlb + 2*dt, 0.05), min(0.8*period, max(5, 4*tauSlb))]
% which for a fitted tauSlb of 0.358 ns at 0.16 ns bins gives a lower bound of
% 0.678 ns. That guard keeps membrane states at least two TCSPC bins away from
% the SLB. Forcing a lower bound (an earlier version used 0.4 ns) puts grid
% states within two bins of the SLB, where they are nearly degenerate with it -
% and because the SLB amplitude is fixed, such a state acts as a spurious extra
% amplitude degree of freedom. It is also unphysical for MIET, where the SLB is
% the most quenched layer so membrane lifetimes must be longer.
pipelineCfg.bayes = struct('batchSize', 2048, 'includeBackground', true, ...
    'signalGrid', [0.25 0.5 0.75 1], 'membraneTauCount', 48, ...
    'fractionStep', 0.2, 'minimumMembraneFraction', 0.1, ...
    'slbCountRelTol', 0.0025, 'slbCountPriorNodes', 5);
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
