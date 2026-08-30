function [summary, info] = run_batch_immune_cell_MIET_height_maps(dataRoot, calibrationFile, cfg)
%RUN_BATCH_IMMUNE_CELL_MIET_HEIGHT_MAPS Height images for a whole session.
%
% [summary, info] = run_batch_immune_cell_MIET_height_maps()
% [summary, info] = run_batch_immune_cell_MIET_height_maps(dataRoot)
% [summary, info] = run_batch_immune_cell_MIET_height_maps(dataRoot, calibrationFile, cfg)
%
% Walks a session folder, finds every bias-corrected MIET analysis under it,
% and converts each one's native, 2x2 and 4x4 lifetime maps into height maps
% through the MIET calibration curve. One row per acquisition-and-binning is
% returned and written as CSV beside the session.
%
% This deliberately starts from the BIAS-CORRECTED maps rather than re-running
% the fits: the expensive Bayesian analysis is already on disk, the correction
% is already validated against each acquisition's own SLB lifetime and IRF, and
% a height map is a pure post-processing step on top of it. Re-running the fits
% to produce heights would risk a different answer for reasons that have
% nothing to do with the calibration.
%
% ONE CALIBRATION PER ACQUISITION, ANCHORED ON ITS OWN SLB
%
% By default each acquisition gets its OWN calibration curve, with the quantum
% yield solved from that acquisition's measured bare-SLB lifetime. The measured
% tau_SLB varies over the session (0.31-0.39 ns here), and a single shared
% curve would silently convert that variation into a spurious 20-30 nm height
% offset between fields. Anchoring per acquisition puts each height map's zero
% on its own bilayer instead, which is the same reasoning the existing
% per-acquisition inverse bias calibration already follows.
%
% Pass a calibrationFile explicitly to use one shared curve instead.
%
% Defaults: dataRoot  the RT MemGlow session of 2026-08-13
%
% cfg fields
%   resume          true - skip acquisitions whose height MAT is newer than
%                   both its bias correction and the calibration
%   continueOnError true
%   showFigures     false
%   writeFigures    true
%   summaryDir      '' puts the CSV in <dataRoot>\immune_cell_MIET_height_maps
%   heightMaps      struct forwarded to immune_cell_MIET_height_maps
%   anchor          struct controlling the per-acquisition anchoring
%     .enabled              true
%     .slbHeightNm          4
%     .meanCosSquaredTheta  1/3, isotropic wobble
%     .slbLifetimeNsOverride  [] anchors on each acquisition's own fitted SLB
%                           lifetime. Set a scalar to anchor every acquisition
%                           on one session-level value instead - needed when
%                           the SLB lifetime comes from a separate measurement
%                           rather than from the pipeline's fixed component
%                           (see the H2 hypothesis in make_rt_miet_calibration)
%     .tag                  '' derives the output-name tag from
%                           meanCosSquaredTheta; set a string to distinguish
%                           runs that differ in something else, e.g. tau_free
%     .calibration          extra cfg forwarded to make_rt_miet_calibration
%                           (spacer thickness, gold thickness, tauFreeNs, ...)

    projectDir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(dataRoot)
        dataRoot = 'D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_20260813_1';
    end
    if nargin < 3 || isempty(cfg); cfg = struct(); end
    cfg = withDefaults(cfg, struct( ...
        'resume', true, ...
        'continueOnError', true, ...
        'showFigures', false, ...
        'writeFigures', true, ...
        'summaryDir', '', ...
        'heightMaps', struct(), ...
        'anchor', struct()));
    cfg.anchor = withDefaults(cfg.anchor, struct( ...
        'enabled', true, ...
        'slbHeightNm', 4, ...
        'meanCosSquaredTheta', 1/3, ...
        'slbLifetimeNsOverride', [], ...
        'tag', '', ...
        'calibration', struct()));

    if nargin < 2 || isempty(calibrationFile)
        calibrationFile = '';
    end
    calibrationFile = char(calibrationFile);
    if ~isempty(calibrationFile)
        % An explicit file overrides anchoring: the caller asked for one curve.
        cfg.anchor.enabled = false;
    elseif ~cfg.anchor.enabled
        error('run_batch_immune_cell_MIET_height_maps:NoCalibration', ...
            ['With cfg.anchor.enabled = false a calibrationFile must be given ' ...
             'explicitly.']);
    end
    if ~isfolder(dataRoot)
        error('run_batch_immune_cell_MIET_height_maps:DataRoot', ...
            'Session folder does not exist: %s', dataRoot);
    end
    if isempty(cfg.summaryDir)
        cfg.summaryDir = fullfile(dataRoot, ['immune_cell_MIET_height_maps' ...
            variantTag(cfg)]);
    end
    if ~isfolder(cfg.summaryDir); mkdir(cfg.summaryDir); end

    found = dir(fullfile(dataRoot, '**', ...
        'immune_cell_MIET_long_lifetime_bias_correction.mat'));
    found = found(~[found.isdir]);
    if isempty(found)
        error('run_batch_immune_cell_MIET_height_maps:NoInputs', ...
            ['No immune_cell_MIET_long_lifetime_bias_correction.mat under %s. ' ...
             'Run run_batch_immune_cell_MIET with ' ...
             'applyLongLifetimeBiasCorrection = true first.'], dataRoot);
    end
    [~, order] = sort({found.folder});
    found = found(order);

    fprintf('run_batch_immune_cell_MIET_height_maps\n');
    fprintf('  session      %s\n', dataRoot);
    if cfg.anchor.enabled
        fprintf(['  calibration  per acquisition, quantum yield anchored on ' ...
            'its own SLB at %g nm, <cos^2 theta> = %.4f\n'], ...
            cfg.anchor.slbHeightNm, cfg.anchor.meanCosSquaredTheta);
    else
        fprintf('  calibration  %s (shared)\n', calibrationFile);
    end
    fprintf('  inputs       %d bias-corrected analyses\n\n', numel(found));

    rows = {};
    errors = {};
    skipped = 0;
    calibrationInfo = struct([]);
    anchorRows = {};
    for index = 1:numel(found)
        biasMat = fullfile(found(index).folder, found(index).name);
        acquisition = acquisitionLabel(biasMat, dataRoot);
        heightMat = fullfile(found(index).folder, sprintf( ...
            'immune_cell_MIET_height_maps%s.mat', variantTag(cfg)));
        fprintf('[%d/%d] %s\n', index, numel(found), acquisition);

        try
            [acquisitionCalibration, anchorRow] = resolveCalibration( ...
                biasMat, calibrationFile, cfg, acquisition);
        catch calibrationError
            errors{end + 1} = sprintf('%s: %s', acquisition, ...
                calibrationError.message); %#ok<AGROW>
            fprintf(2, '  FAILED (calibration): %s\n', calibrationError.message);
            if ~cfg.continueOnError; rethrow(calibrationError); end
            continue;
        end
        if ~isempty(anchorRow)
            anchorRows{end + 1, 1} = anchorRow; %#ok<AGROW>
        end

        if cfg.resume && isUpToDate(heightMat, biasMat, acquisitionCalibration)
            fprintf('  up to date, reusing %s\n', heightMat);
            skipped = skipped + 1;
            loadedMaps = load(heightMat, 'heightMaps');
            heightMaps = loadedMaps.heightMaps;
        else
            try
                stageCfg = cfg.heightMaps;
                stageCfg.showFigures = cfg.showFigures;
                stageCfg.writeFigures = cfg.writeFigures;
                stageCfg.nameSuffix = variantTag(cfg);
                stageCfg.slbLifetimeNsOverride = cfg.anchor.slbLifetimeNsOverride;
                heightMaps = immune_cell_MIET_height_maps( ...
                    biasMat, acquisitionCalibration, stageCfg);
            catch conversionError
                errors{end + 1} = sprintf('%s: %s', acquisition, ...
                    conversionError.message); %#ok<AGROW>
                fprintf(2, '  FAILED: %s\n', conversionError.message);
                if ~cfg.continueOnError
                    rethrow(conversionError);
                end
                continue;
            end
        end

        if isempty(calibrationInfo)
            calibrationInfo = heightMaps.calibration;
        end
        rows = [rows; stageRows(heightMaps, acquisition, biasMat)]; %#ok<AGROW>
    end

    summary = rowsToTable(rows);
    summaryCsv = fullfile(cfg.summaryDir, 'immune_cell_MIET_height_summary.csv');
    summaryMat = fullfile(cfg.summaryDir, 'immune_cell_MIET_height_summary.mat');
    if ~isempty(summary)
        writetable(summary, summaryCsv);
    end

    info = struct();
    info.dataRoot = dataRoot;
    info.calibrationFile = calibrationFile;
    info.calibration = calibrationInfo;
    info.anchorConfig = cfg.anchor;
    info.anchorTable = rowsToTable(anchorRows);
    info.inputCount = numel(found);
    info.skippedCount = skipped;
    info.errors = {errors};
    info.summaryCsv = summaryCsv;
    info.summaryMat = summaryMat;
    if ~isempty(info.anchorTable)
        info.anchorCsv = fullfile(cfg.summaryDir, ...
            'immune_cell_MIET_height_anchor.csv');
        writetable(info.anchorTable, info.anchorCsv);
    end
    save(summaryMat, 'summary', 'info', '-v7.3');

    printSessionSummary(summary, info, errors);
end

% ------------------------------------------------------------------- helpers

function tag = variantTag(cfg)
% Every output name carries the orientation the curve was built with. Without
% it a second variant would overwrite the first, and two height maps that
% differ only by an assumption would be indistinguishable on disk.
    if ~cfg.anchor.enabled
        tag = '';
        return;
    end
    if ~isempty(cfg.anchor.tag)
        tag = char(cfg.anchor.tag);
        if tag(1) ~= '_'; tag = ['_' tag]; end
        return;
    end
    tag = sprintf('_f%03d', round(cfg.anchor.meanCosSquaredTheta * 100));
end

function [calibrationFile, anchorRow] = resolveCalibration( ...
        biasMat, sharedCalibrationFile, cfg, acquisition)
% Either return the one shared curve, or build this acquisition's own curve
% with the quantum yield anchored on its measured SLB lifetime. The curve is
% written into the acquisition's results folder, so the height map and the
% calibration that produced it stay together.
    anchorRow = [];
    if ~cfg.anchor.enabled
        calibrationFile = sharedCalibrationFile;
        return;
    end

    loaded = load(biasMat, 'biasCorrection');
    pipelineTauSlb = double(loaded.biasCorrection.fixedSlbLifetimeNs);
    if ~isempty(cfg.anchor.slbLifetimeNsOverride)
        % The pipeline's fixed short component is NOT necessarily the bilayer -
        % under the H2 reading it is substrate emission. When the SLB lifetime
        % is measured separately, anchor on that instead, and keep the
        % pipeline's value only for the record.
        tauSlb = double(cfg.anchor.slbLifetimeNsOverride);
    else
        tauSlb = pipelineTauSlb;
    end
    if ~isfinite(tauSlb) || tauSlb <= 0
        error('run_batch_immune_cell_MIET_height_maps:TauSlb', ...
            'No usable SLB lifetime for the anchor.');
    end

    outputDir = fileparts(biasMat);
    calibrationFile = fullfile(outputDir, sprintf( ...
        'miet_calibration_slb_anchored%s.mat', variantTag(cfg)));
    calibrationCfg = cfg.anchor.calibration;
    calibrationCfg.anchorSlbLifetimeNs = tauSlb;
    calibrationCfg.anchorSlbHeightNm = cfg.anchor.slbHeightNm;
    calibrationCfg.meanCosSquaredTheta = cfg.anchor.meanCosSquaredTheta;
    calibrationCfg.outputMat = calibrationFile;
    calibrationCfg.showFigure = false;
    calibrationCfg.writeFigure = cfg.writeFigures;
    calib = make_rt_miet_calibration(calibrationCfg);

    anchorRow = struct();
    anchorRow.acquisition = string(acquisition);
    anchorRow.slbLifetimeNs = tauSlb;
    anchorRow.pipelineFixedShortLifetimeNs = pipelineTauSlb;
    anchorRow.slbLifetimeIsOverride = ~isempty(cfg.anchor.slbLifetimeNsOverride);
    anchorRow.slbHeightNm = cfg.anchor.slbHeightNm;
    anchorRow.meanCosSquaredTheta = cfg.anchor.meanCosSquaredTheta;
    anchorRow.effectiveTiltDeg = acosd(sqrt(cfg.anchor.meanCosSquaredTheta));
    anchorRow.quantumYield = calib.params.quantumYield;
    anchorRow.tauFreeNs = calib.params.tauFreeNs;
    anchorRow.floorLifetimeNs = calib.minLifetimeNs;
    anchorRow.ceilingLifetimeNs = calib.maxLifetimeNs;
    anchorRow.maxHeightNm = calib.maxHeightNm;
    anchorRow.parallelLimitLifetimeNs = calib.anchor.parallelLimitLifetimeNs;
    anchorRow.verticalLimitLifetimeNs = calib.anchor.minAttainableLifetimeNs;
    anchorRow.minimumFeasibleOrder = calib.anchor.minimumFeasibleOrder;
    anchorRow.calibrationMat = string(calibrationFile);
end

function rows = stageRows(heightMaps, acquisition, biasMat)
    rows = {};
    names = fieldnames(heightMaps.stages);
    for index = 1:numel(names)
        stage = heightMaps.stages.(names{index});
        row = struct();
        row.acquisition = string(acquisition);
        row.binning = string(names{index});
        row.available = stage.available;
        row.fixedSlbLifetimeNs = double(heightMaps.fixedSlbLifetimeNs);
        row.slbBelowCalibrationFloor = ...
            logical(heightMaps.slbConsistency.slbLifetimeIsBelowFloor);
        if stage.available
            s = stage.summary;
            row.reliableLifetimePixels = s.reliableLifetimePixelCount;
            row.heightPixels = s.heightPixelCount;
            row.invertedFraction = s.heightFractionOfReliable;
            row.aboveCeilingPixels = s.aboveCalibrationCeilingCount;
            row.belowFloorPixels = s.belowCalibrationFloorCount;
            row.medianHeightNm = s.medianHeightNm;
            row.p05HeightNm = s.p05HeightNm;
            row.p95HeightNm = s.p95HeightNm;
            row.medianHeightStdNm = s.medianHeightStdNm;
            row.medianHeightCorrectionNm = s.medianHeightCorrectionNm;
        else
            row.reliableLifetimePixels = 0;
            row.heightPixels = 0;
            row.invertedFraction = NaN;
            row.aboveCeilingPixels = 0;
            row.belowFloorPixels = 0;
            row.medianHeightNm = NaN;
            row.p05HeightNm = NaN;
            row.p95HeightNm = NaN;
            row.medianHeightStdNm = NaN;
            row.medianHeightCorrectionNm = NaN;
        end
        row.biasCorrectionMat = string(biasMat);
        row.heightMapsMat = string(heightMaps.outputFiles.mat);
        rows{end + 1, 1} = row; %#ok<AGROW>
    end
end

function summary = rowsToTable(rows)
    if isempty(rows); summary = table(); return; end
    summary = struct2table(cell2mat(cellfun(@(r) r, rows, 'UniformOutput', false)));
end

function label = acquisitionLabel(biasMat, dataRoot)
    relative = erase(biasMat, [char(dataRoot) filesep]);
    parts = split(string(relative), filesep);
    if numel(parts) >= 1
        label = char(parts(1));
    else
        label = relative;
    end
    % Keep the configuration subfolder too: one acquisition can hold several
    % analyses that differ only by cfg hash, and collapsing them to the
    % timestamp would silently merge two different runs in the summary.
    if numel(parts) >= 3
        label = sprintf('%s|%s', char(parts(1)), char(parts(end - 1)));
    end
end

function tf = isUpToDate(heightMat, biasMat, calibrationFile)
    tf = false;
    if ~isfile(heightMat); return; end
    target = dir(heightMat);
    sources = [dir(biasMat); dir(calibrationFile)];
    if isempty(sources); return; end
    tf = target.datenum >= max([sources.datenum]);
end

function printSessionSummary(summary, info, errors)
    fprintf('\n--- session summary ---\n');
    fprintf('inputs %d, reused %d, failed %d\n', info.inputCount, ...
        info.skippedCount, numel(errors));
    if ~isempty(info.calibration)
        c = info.calibration;
        fprintf('calibration window %.3f-%.3f ns over %.1f-%.1f nm\n', ...
            c.lifetimeNs(1), c.lifetimeNs(end), c.heightNm(1), c.heightNm(end));
    end
    if isfield(info, 'anchorTable') && ~isempty(info.anchorTable)
        a = info.anchorTable;
        fprintf(['SLB anchor: tau_SLB %.3f-%.3f ns -> quantum yield %.3f-%.3f ' ...
            '(median %.3f) at <cos^2 theta> = %.3f\n'], ...
            min(a.slbLifetimeNs), max(a.slbLifetimeNs), ...
            min(a.quantumYield), max(a.quantumYield), ...
            median(a.quantumYield), a.meanCosSquaredTheta(1));
        fprintf(['            usable height ceiling %.1f-%.1f nm, ' ...
            'lifetime ceiling %.3f-%.3f ns\n'], ...
            min(a.maxHeightNm), max(a.maxHeightNm), ...
            min(a.ceilingLifetimeNs), max(a.ceilingLifetimeNs));
    end
    if ~isempty(summary)
        available = summary(summary.available, :);
        binnings = unique(available.binning);
        for index = 1:numel(binnings)
            rows = available(available.binning == binnings(index), :);
            fprintf(['%-11s median z %.1f nm over %d acquisitions | ' ...
                     'inverted %.1f%% of reliable pixels | above ceiling %.1f%%\n'], ...
                char(binnings(index)), median(rows.medianHeightNm, 'omitnan'), ...
                height(rows), 100 * mean(rows.invertedFraction, 'omitnan'), ...
                100 * sum(rows.aboveCeilingPixels) / ...
                    max(sum(rows.reliableLifetimePixels), 1));
        end
        if any(summary.slbBelowCalibrationFloor)
            fprintf(['\nWARNING: in %d of %d rows the measured SLB lifetime is ' ...
                'shorter than the calibration floor. The absolute height zero ' ...
                'is therefore unverified - see heightMaps.slbConsistency.\n'], ...
                nnz(summary.slbBelowCalibrationFloor), height(summary));
        end
    end
    fprintf('summary CSV %s\n', info.summaryCsv);
    for index = 1:numel(errors)
        fprintf(2, 'error: %s\n', errors{index});
    end
end

function cfg = withDefaults(cfg, defaults)
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(cfg, names{index})
            cfg.(names{index}) = defaults.(names{index});
        end
    end
end
