function [summary, batchInfo] = batch_immune_cell_MIET(dataRoot, cfg)
%BATCH_IMMUNE_CELL_MIET Analyze suitable immune-cell MIET PTUs recursively.
%
% [summary, batchInfo] = batch_immune_cell_MIET(dataRoot, cfg)
%
% Every RawImage.ptu below dataRoot is audited and written to the summary.
% By default, only complete, high-dwell XY scans are sent to
% immune_cell_MIET. The image-derived result is then labelled as a good cell
% experiment only when segmentation, cell/SLB support, lifetime contrast,
% and Bayesian coverage pass configurable QC thresholds.
%
% Important cfg fields (all optional):
%   scanPlanes                 planes to analyze (default {'XY'})
%   minDwellMs                 metadata prefilter (default 0.4 ms/pixel)
%   minCompletionFraction      empirical raster completion (default 0.95)
%   minRecords                 minimum TTTR records (default 1e6)
%   resume                     reuse a valid existing result (default true)
%   overwrite                  recompute existing results (default false)
%   resultsFolderName          name of the per-acquisition results folder
%                              created beside each PTU (default
%                              'immune_cell_MIET')
%   versionResults             put results in a per-configuration subfolder so
%                              a run with different settings cannot overwrite
%                              an earlier one (default true)
%   runName                    name that subfolder explicitly. '' (default)
%                              derives it from a hash of the analysis
%                              configuration, which is the useful behaviour:
%                              an unchanged configuration resolves to the same
%                              folder and RESUMES, while any change to the
%                              pipeline or QC settings resolves to a new folder
%                              and leaves the old results untouched.
%   dryRun                     header audit only (default false)
%   saveTcspcPix               save detector-summed uint16 cube (default true)
%   pipeline                   options passed to immune_cell_MIET
%   qc                         image-derived good-cell QC overrides
%
% Folder-level outputs are checkpointed after every file in
% <dataRoot>/immune_cell_MIET_batch/<runName>. Per-acquisition images and MAT
% files go to <acquisition>/<resultsFolderName>/<runName>, where runName
% defaults to a hash of the analysis configuration - so re-running with changed
% settings writes alongside the previous results instead of over them, while
% re-running with the same settings still resumes. No photon stream is placed
% in a batch output.

    if nargin < 1 || isempty(dataRoot)
        dataRoot = ...
            'E:\Luminosa_data\260813\RT_Jurkat_CD58_memglow_20260813_1';
    end
    if nargin < 2 || isempty(cfg)
        cfg = struct();
    end
    cfg = fillBatchDefaults(cfg);
    dataRoot = char(dataRoot);
    if ~isfolder(dataRoot)
        error('batch_immune_cell_MIET:MissingRoot', ...
            'Data root was not found: %s', dataRoot);
    end
    requireFunction('PTU_Read_Head');
    requireFunction('immune_cell_MIET_scan_geometry');
    if ~cfg.dryRun
        requireFunction('immune_cell_MIET');
    end

    % A batch-level cfg.outputDir has never been read by this function: the
    % per-acquisition folder is derived from each PTU's own location, and the
    % batch summary uses cfg.batchOutputDir. Setting it was therefore silently
    % ineffective, which is exactly how one run overwrote another.
    if isfield(cfg, 'outputDir') && ~isempty(cfg.outputDir)
        warning('batch_immune_cell_MIET:OutputDirIgnored', ...
            ['cfg.outputDir is not used by the batch driver and has been ' ...
             'ignored. Use cfg.resultsFolderName to rename the per-' ...
             'acquisition folder, cfg.runName to name this run, or ' ...
             'cfg.batchOutputDir for the batch summary. Ignored value: %s'], ...
            cfg.outputDir);
    end
    cfg.runName = immune_cell_MIET_run_folder(cfg);
    if isempty(cfg.runName)
        fprintf(['batch_immune_cell_MIET: results go to <acquisition>\\%s ' ...
            '(flat, unversioned - a rerun overwrites)\n'], ...
            cfg.resultsFolderName);
    else
        fprintf('batch_immune_cell_MIET: results go to <acquisition>\\%s\\%s\n', ...
            cfg.resultsFolderName, cfg.runName);
    end
    reportLegacyResults(cfg, dataRoot);

    outputDir = cfg.batchOutputDir;
    if isempty(outputDir)
        outputDir = fullfile(dataRoot, 'immune_cell_MIET_batch', cfg.runName);
    end
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    writeRunManifest(outputDir, cfg, dataRoot);
    csvFile = fullfile(outputDir, 'immune_cell_MIET_batch_summary.csv');
    matFile = fullfile(outputDir, 'immune_cell_MIET_batch_summary.mat');
    errorLogFile = fullfile(outputDir, 'immune_cell_MIET_batch_errors.log');

    batchSourceFile = [mfilename('fullpath') '.m'];
    batchSourceInfo = dir(batchSourceFile);
    if isempty(batchSourceInfo)
        batchSourceStamp = 'unknown';
    else
        batchSourceStamp = batchSourceInfo(1).date;
    end
    requested4x4 = isfield(cfg.pipeline, 'spatialBinning4x4') && ...
        isstruct(cfg.pipeline.spatialBinning4x4) && ...
        isfield(cfg.pipeline.spatialBinning4x4, 'enabled') && ...
        logical(cfg.pipeline.spatialBinning4x4.enabled);
    fprintf(['batch_immune_cell_MIET: 4x4 sliding stage requested = %d ' ...
        '(batch source modified %s).\n'], requested4x4, batchSourceStamp);
    if isfield(cfg.pipeline, 'maxMembraneStates') && ...
            ~isempty(cfg.pipeline.maxMembraneStates)
        maxMembraneStates = double(cfg.pipeline.maxMembraneStates);
    else
        maxMembraneStates = 2;
    end
    fprintf(['batch_immune_cell_MIET: free membrane states per pixel = %d ' ...
        '(results folder %s).\n'], maxMembraneStates, cfg.resultsFolderName);
    if requested4x4 && ~any(strcmp(fieldnames(emptyRow()), ...
            'sliding4x4SecondComponentPixels'))
        error('batch_immune_cell_MIET:Stale4x4Schema', ...
            ['The 4x4 stage was requested but this loaded copy of ' ...
             'batch_immune_cell_MIET has no 4x4 summary columns. Run ' ...
             '"clear functions" or restart MATLAB and try again.']);
    end

    files = dir(fullfile(dataRoot, '**', 'RawImage.ptu'));
    if isempty(files)
        error('batch_immune_cell_MIET:NoFiles', ...
            'No RawImage.ptu files were found below %s.', dataRoot);
    end
    fullNames = arrayfun(@(x) fullfile(x.folder, x.name), files, ...
        'UniformOutput', false);
    [~, order] = sort(fullNames);
    files = files(order);

    rows = repmat(emptyRow(), numel(files), 1);
    fprintf('batch_immune_cell_MIET: auditing %d PTU files below %s\n', ...
        numel(files), dataRoot);
    for fileIndex = 1:numel(files)
        ptuFile = fullfile(files(fileIndex).folder, files(fileIndex).name);
        rows(fileIndex) = initialiseRow(rows(fileIndex), fileIndex, cfg, ...
            files(fileIndex), ptuFile);
        try
            % The repository contains both one- and two-argument copies of
            % PTU_Read_Head. One argument is compatible with both.
            head = PTU_Read_Head(ptuFile);
            geometry = immune_cell_MIET_scan_geometry(head);
            rows(fileIndex) = addHeaderMetadata(rows(fileIndex), head, ...
                geometry, cfg.excitationPulseIndex);
            [rows(fileIndex).preflightEligible, ...
                rows(fileIndex).preflightReason] = preflightDecision( ...
                rows(fileIndex), head, cfg);
            if rows(fileIndex).preflightEligible
                rows(fileIndex).status = 'eligible';
            else
                rows(fileIndex).status = 'skipped_preflight';
            end
        catch headerError
            rows(fileIndex).status = 'header_failed';
            rows(fileIndex).preflightEligible = false;
            rows(fileIndex).preflightReason = cleanText(headerError.message);
            rows(fileIndex).errorIdentifier = headerError.identifier;
            rows(fileIndex).message = cleanText(headerError.message);
            appendErrorLog(errorLogFile, ptuFile, headerError);
        end
    end

    [summary, batchInfo] = checkpoint(rows, cfg, dataRoot, outputDir, ...
        csvFile, matFile, errorLogFile);
    if cfg.dryRun
        fprintf('batch_immune_cell_MIET: dry run complete; no PTU was analyzed.\n');
        printSummary(batchInfo, csvFile);
        return;
    end

    eligibleIndices = find([rows.preflightEligible]);
    if isfinite(cfg.maxFiles)
        eligibleIndices = eligibleIndices(1:min(numel(eligibleIndices), ...
            floor(cfg.maxFiles)));
    end
    for position = 1:numel(eligibleIndices)
        fileIndex = eligibleIndices(position);
        ptuFile = rows(fileIndex).ptuFile;
        analysisMat = rows(fileIndex).analysisMat;
        fprintf('batch_immune_cell_MIET: %d/%d eligible (%d/%d total): %s\n', ...
            position, numel(eligibleIndices), fileIndex, numel(rows), ptuFile);
        started = tic;
        figuresBefore = findall(groot, 'Type', 'figure');
        figureCleanup = onCleanup(@() closeNewPipelineFigures(figuresBefore));
        try
            localCfg = cfg.pipeline;
            localCfg.outputDir = rows(fileIndex).analysisDir;
            localCfg.showFigure = cfg.showFigures;
            localCfg.saveOutputs = true;
            localCfg.saveTcspcPix = cfg.saveTcspcPix;
            reused = false;
            if cfg.resume && ~cfg.overwrite && isfile(analysisMat)
                loaded = load(analysisMat, 'result');
                if isfield(loaded, 'result') && isstruct(loaded.result) && ...
                        existingResultUsable(analysisMat, loaded.result, ...
                        cfg.saveTcspcPix, ptuFile, localCfg)
                    result = loaded.result;
                    reused = true;
                end
            end
            if reused
                rows(fileIndex).action = 'reused_existing_result';
            else
                result = immune_cell_MIET(ptuFile, localCfg);
                rows(fileIndex).action = 'processed';
            end
            rows(fileIndex) = addResultQc(rows(fileIndex), result, cfg.qc);
            if requested4x4
                if isfinite(rows(fileIndex).sliding4x4SecondComponentPixels)
                    if maxMembraneStates == 1
                        fprintf(['batch_immune_cell_MIET: 4x4 stage recorded ' ...
                            '(%d membrane-lifetime pixels; median %.4g ns).\n'], ...
                            rows(fileIndex).sliding4x4SecondComponentPixels, ...
                            rows(fileIndex).sliding4x4MembraneLifetimeMedianNs);
                    else
                        fprintf(['batch_immune_cell_MIET: 4x4 stage recorded ' ...
                            '(second %d, third %d pixels).\n'], ...
                            rows(fileIndex).sliding4x4SecondComponentPixels, ...
                            rows(fileIndex).sliding4x4ThirdComponentPixels);
                    end
                else
                    error('batch_immune_cell_MIET:Missing4x4Result', ...
                        ['The 4x4 stage was requested but no 4x4 summary was ' ...
                         'produced for %s.'], ptuFile);
                end
            end
            if rows(fileIndex).goodCell
                if reused
                    rows(fileIndex).status = 'existing_good_cell';
                else
                    rows(fileIndex).status = 'complete_good_cell';
                end
            else
                if reused
                    rows(fileIndex).status = 'existing_qc_reject';
                else
                    rows(fileIndex).status = 'complete_qc_reject';
                end
            end
            clear result loaded
        catch analysisError
            rows(fileIndex).status = 'analysis_failed';
            rows(fileIndex).action = 'attempted';
            rows(fileIndex).goodCell = false;
            rows(fileIndex).qcReason = cleanText(analysisError.message);
            rows(fileIndex).errorIdentifier = analysisError.identifier;
            rows(fileIndex).message = cleanText(analysisError.message);
            appendErrorLog(errorLogFile, ptuFile, analysisError);
            if ~cfg.continueOnError
                rows(fileIndex).processingSeconds = toc(started);
                clear figureCleanup
                checkpoint(rows, cfg, dataRoot, ...
                    outputDir, csvFile, matFile, errorLogFile);
                rethrow(analysisError);
            end
        end
        rows(fileIndex).processingSeconds = toc(started);
        clear figureCleanup
        checkpoint(rows, cfg, dataRoot, outputDir, ...
            csvFile, matFile, errorLogFile);
    end

    unprocessedEligible = setdiff(find([rows.preflightEligible]), eligibleIndices);
    for k = 1:numel(unprocessedEligible)
        idx = unprocessedEligible(k);
        rows(idx).status = 'eligible_not_run_max_files';
        rows(idx).action = 'not_run';
    end
    [summary, batchInfo] = checkpoint(rows, cfg, dataRoot, outputDir, ...
        csvFile, matFile, errorLogFile);
    printSummary(batchInfo, csvFile);
end

function reportLegacyResults(cfg, dataRoot)
%REPORTLEGACYRESULTS Say plainly when older flat-layout results exist.
% Switching on versioning means results previously written directly into
% resultsFolderName are no longer where the resume check looks. They are not
% deleted or moved - this function makes sure that is stated rather than
% discovered later as an unexplained full re-analysis.
    if isempty(cfg.runName)
        return;
    end
    legacy = dir(fullfile(dataRoot, '**', cfg.resultsFolderName, ...
        'immune_cell_MIET_640nm_red_analysis.mat'));
    if isempty(legacy)
        return;
    end
    fprintf(['batch_immune_cell_MIET: %d result file(s) exist in the older ' ...
        'flat layout\n'], numel(legacy));
    fprintf(['  <acquisition>\\%s\\immune_cell_MIET_640nm_red_analysis.mat' ...
        '\n'], cfg.resultsFolderName);
    fprintf(['  They are left untouched and will NOT be resumed, so this run ' ...
        'recomputes from\n  scratch. To resume them instead, either move ' ...
        'each one into a subfolder and\n  pass its name as cfg.runName, or ' ...
        'set cfg.versionResults = false to keep the\n  old flat layout.\n']);
end

function writeRunManifest(outputDir, cfg, dataRoot)
%WRITERUNMANIFEST Record what produced the results in this folder.
% Without this, two run folders are indistinguishable a month later, which
% defeats the point of keeping both.
    manifest = fullfile(outputDir, 'run_manifest.txt');
    fid = fopen(manifest, 'w');
    if fid < 0
        warning('batch_immune_cell_MIET:Manifest', ...
            'Could not write %s', manifest);
        return;
    end
    closer = onCleanup(@() fclose(fid));
    fprintf(fid, 'run name      : %s\n', cfg.runName);
    fprintf(fid, 'written       : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, 'data root     : %s\n', dataRoot);
    fprintf(fid, 'results folder: %s\n', cfg.resultsFolderName);
    fprintf(fid, 'resume        : %d, overwrite: %d, versioned: %d\n', ...
        cfg.resume, cfg.overwrite, cfg.versionResults);
    fprintf(fid, 'code revision : %s\n', codeRevision());
    fprintf(fid, '\n--- pipeline configuration ---\n');
    printStruct(fid, cfg.pipeline, '');
    fprintf(fid, '\n--- QC thresholds ---\n');
    printStruct(fid, cfg.qc, '');
end

function revision = codeRevision()
    revision = 'unknown (not a git checkout or git unavailable)';
    here = fileparts(mfilename('fullpath'));
    [status, out] = system(sprintf('git -C "%s" rev-parse --short HEAD', here));
    if status == 0 && ~isempty(strtrim(out))
        revision = strtrim(out);
        [dirtyStatus, dirty] = system(sprintf( ...
            'git -C "%s" status --porcelain', here));
        if dirtyStatus == 0 && ~isempty(strtrim(dirty))
            revision = [revision ' (working tree modified)'];
        end
    end
end

function printStruct(fid, value, prefix)
    if ~isstruct(value) || ~isscalar(value)
        fprintf(fid, '%s = %s\n', prefix, compactValue(value));
        return;
    end
    names = fieldnames(value);
    for k = 1:numel(names)
        if isempty(prefix)
            label = names{k};
        else
            label = [prefix '.' names{k}];
        end
        item = value.(names{k});
        if isstruct(item) && isscalar(item)
            printStruct(fid, item, label);
        else
            fprintf(fid, '  %-34s %s\n', label, compactValue(item));
        end
    end
end

function text = compactValue(value)
    if ischar(value)
        text = value;
    elseif isstring(value)
        text = char(strjoin(value, ', '));
    elseif iscell(value)
        parts = cellfun(@compactValue, value, 'UniformOutput', false);
        text = ['{' strjoin(parts, ', ') '}'];
    elseif islogical(value) || isnumeric(value)
        if isempty(value)
            text = '[]';
        else
            text = mat2str(value, 6);
        end
    else
        text = ['<' class(value) '>'];
    end
end

function cfg = fillBatchDefaults(cfg)
    defaults = struct();
    defaults.scanPlanes = {'XY'};
    defaults.minDwellMs = 0.4;
    defaults.minCompletionFraction = 0.95;
    defaults.minRecords = 1e6;
    defaults.requiredPieWindows = 2;
    defaults.excitationPulseIndex = 2;
    defaults.requiredExcitationNm = 640;
    defaults.wavelengthToleranceNm = 10;
    defaults.requireMetadataForPrefilter = true;
    defaults.resume = true;
    defaults.overwrite = false;
    defaults.dryRun = false;
    defaults.continueOnError = true;
    defaults.maxFiles = Inf;
    defaults.showFigures = false;
    defaults.saveTcspcPix = true;
    defaults.batchOutputDir = '';
    defaults.resultsFolderName = 'immune_cell_MIET';
    defaults.versionResults = true;
    defaults.runName = '';
    defaults.pipeline = struct();
    defaults.qc = struct('acceptedSegmentationStatuses', {{'ok'}}, ...
        'minCellFraction', 0.01, 'maxCellFraction', 0.85, ...
        'rejectCellTouchingBorder', true, ...
        'minSlbReferenceFraction', 0.01, ...
        'minSlbReferencePhotons', 500, ...
        'maxSlbRobustSigmaNs', 0.15, ...
        'maxSlbClippedFraction', 0.10, ...
        'minLifetimeContrastNs', 0.05, ...
        'minBayesCoverage', 0.25, ...
        'slbLifetimeRangeNs', [0.03 5]);
    cfg = mergeStruct(defaults, cfg);
    cfg.qc = mergeStruct(defaults.qc, cfg.qc);
    if ischar(cfg.scanPlanes) || (isstring(cfg.scanPlanes) && isscalar(cfg.scanPlanes))
        cfg.scanPlanes = {char(cfg.scanPlanes)};
    elseif isstring(cfg.scanPlanes)
        cfg.scanPlanes = cellstr(cfg.scanPlanes(:));
    end
    if ~iscellstr(cfg.scanPlanes) || isempty(cfg.scanPlanes)
        error('batch_immune_cell_MIET:ScanPlanes', ...
            'cfg.scanPlanes must contain one or more plane names.');
    end
    cfg.scanPlanes = cellfun(@upper, cfg.scanPlanes, 'UniformOutput', false);
    validateattributes(cfg.minDwellMs, {'numeric'}, {'scalar','nonnegative'});
    validateattributes(cfg.minCompletionFraction, {'numeric'}, ...
        {'scalar','nonnegative'});
    validateattributes(cfg.minRecords, {'numeric'}, {'scalar','nonnegative'});
    validateattributes(cfg.maxFiles, {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.qc.maxSlbClippedFraction, {'numeric'}, ...
        {'real','finite','scalar','>=',0,'<=',1});
end

function row = emptyRow()
    row = struct( ...
        'index', NaN, 'acquisition', '', 'ptuFile', '', ...
        'analysisDir', '', 'analysisMat', '', 'preliminaryPng', '', ...
        'segmentationPng', '', 'meanFlimPng', '', 'bayesPng', '', ...
        'secondLifetimePng', '', 'thirdLifetimePng', '', ...
        'layeredLifetimePng', '', ...
        'sliding2x2MeanFlimPng', '', 'sliding2x2BayesPng', '', ...
        'sliding2x2SecondLifetimePng', '', ...
        'sliding2x2ThirdLifetimePng', '', ...
        'sliding2x2LayeredLifetimePng', '', ...
        'sliding4x4MeanFlimPng', '', 'sliding4x4BayesPng', '', ...
        'sliding4x4SecondLifetimePng', '', ...
        'sliding4x4ThirdLifetimePng', '', ...
        'sliding4x4LayeredLifetimePng', '', ...
        'fileBytes', NaN, 'scanDirection', NaN, 'scanPlane', 'unknown', ...
        'scanConfidence', 'low', 'scanReason', '', ...
        'pixelX', NaN, 'pixelY', NaN, 'pixelSizeUm', NaN, ...
        'dwellMs', NaN, 'elapsedMs', NaN, 'expectedRasterMs', NaN, ...
        'completionFraction', NaN, 'bidirectional', NaN, ...
        'headerMaxFrames', NaN, ...
        'numberOfRecords', NaN, 'pieWindows', NaN, ...
        'recordedExcitationNm', NaN, 'tcspcResolutionPs', NaN, ...
        'preflightEligible', false, 'preflightReason', '', ...
        'status', 'not_audited', 'action', 'not_run', ...
        'goodCell', false, 'qcReason', '', 'segmentationStatus', '', ...
        'cellPixels', NaN, 'cellFraction', NaN, ...
        'cellTouchesBorder', NaN, ...
        'slbReferencePixels', NaN, 'slbReferenceFraction', NaN, ...
        'slbReferencePhotons', NaN, 'slbRobustSigmaNs', NaN, ...
        'lifetimeContrastNs', NaN, 'fixedSlbLifetimeNs', NaN, ...
        'slbSignalPhotonsPerPixel', NaN, ...
        'slbAmplitudeStdPhotonsPerPixel', NaN, ...
        'slbClippedPixelFraction', NaN, ...
        'bayesValidPixels', NaN, 'bayesCoverage', NaN, ...
        'secondComponentPixels', NaN, 'thirdComponentPixels', NaN, ...
        'membraneLifetimeMedianNs', NaN, ...
        'sliding2x2SecondComponentPixels', NaN, ...
        'sliding2x2ThirdComponentPixels', NaN, ...
        'sliding2x2MembraneLifetimeMedianNs', NaN, ...
        'sliding4x4SecondComponentPixels', NaN, ...
        'sliding4x4ThirdComponentPixels', NaN, ...
        'sliding4x4MembraneLifetimeMedianNs', NaN, ...
        'processingSeconds', NaN, 'errorIdentifier', '', 'message', '');
end

function row = initialiseRow(row, fileIndex, cfg, fileInfo, ptuFile)
    row.index = fileIndex;
    row.ptuFile = ptuFile;
    row.fileBytes = double(fileInfo.bytes);
    row.acquisition = acquisitionName(fileInfo.folder);
    % The run subfolder is what keeps a new analysis from overwriting an older
    % one. Every artefact path below is derived from analysisDir, so adding one
    % level here is enough - including for the resume check, which looks for
    % analysisMat inside this same folder.
    row.analysisDir = fullfile(fileInfo.folder, cfg.resultsFolderName, ...
        cfg.runName);
    row.analysisMat = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_analysis.mat');
    row.preliminaryPng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_preliminary_mean_FLIM.png');
    row.segmentationPng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_segmentation_check.png');
    row.meanFlimPng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_mean_FLIM.png');
    row.bayesPng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_Bayes_maps.png');
    row.secondLifetimePng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_second_sorted_lifetime.png');
    row.thirdLifetimePng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_third_sorted_lifetime.png');
    row.layeredLifetimePng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_sorted_lifetime_layers.png');
    row.sliding2x2MeanFlimPng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_2x2_step1_mean_FLIM.png');
    row.sliding2x2BayesPng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_2x2_step1_Bayes_maps.png');
    row.sliding2x2SecondLifetimePng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_2x2_step1_second_sorted_lifetime.png');
    row.sliding2x2ThirdLifetimePng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_2x2_step1_third_sorted_lifetime.png');
    row.sliding2x2LayeredLifetimePng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_2x2_step1_sorted_lifetime_layers.png');
    row.sliding4x4MeanFlimPng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_4x4_step1_mean_FLIM.png');
    row.sliding4x4BayesPng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_4x4_step1_Bayes_maps.png');
    row.sliding4x4SecondLifetimePng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_4x4_step1_second_sorted_lifetime.png');
    row.sliding4x4ThirdLifetimePng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_4x4_step1_third_sorted_lifetime.png');
    row.sliding4x4LayeredLifetimePng = fullfile(row.analysisDir, ...
        'immune_cell_MIET_640nm_red_4x4_step1_sorted_lifetime_layers.png');
end

function name = acquisitionName(folder)
    [~, name] = fileparts(folder);
    if isempty(name)
        name = folder;
    end
end

function row = addHeaderMetadata(row, head, geometry, excitationPulseIndex)
    row.scanPlane = geometry.plane;
    row.scanConfidence = geometry.confidence;
    row.scanDirection = numericScalar(head, 'ImgHdr_ScanDirection');
    row.scanReason = cleanText(strjoin(geometry.reasons, ' | '));
    row.pixelX = numericScalar(head, 'ImgHdr_PixX');
    row.pixelY = numericScalar(head, 'ImgHdr_PixY');
    row.pixelSizeUm = numericScalar(head, 'ImgHdr_PixResol');
    row.dwellMs = numericScalar(head, 'ImgHdr_TimePerPixel');
    row.elapsedMs = numericScalar(head, 'TTResult_StopAfter');
    row.bidirectional = numericScalar(head, 'ImgHdr_BiDirect');
    row.headerMaxFrames = numericScalar(head, 'ImgHdr_MaxFrames');
    row.numberOfRecords = numericScalar(head, 'TTResult_NumberOfRecords');
    row.pieWindows = numericScalar(head, 'PIENumPIEWindows');
    resolutionSeconds = numericScalar(head, 'MeasDesc_Resolution');
    if isfinite(resolutionSeconds)
        row.tcspcResolutionPs = resolutionSeconds * 1e12;
    end
    if isfield(head, 'LaserWL') && isnumeric(head.LaserWL)
        wavelengths = double(head.LaserWL(:));
        if numel(wavelengths) >= excitationPulseIndex
            row.recordedExcitationNm = wavelengths(excitationPulseIndex);
        end
    end
    passFactor = 2;
    if isfinite(row.bidirectional) && row.bidirectional ~= 0
        passFactor = 1;
    end
    if all(isfinite([row.pixelX, row.pixelY, row.dwellMs])) && ...
            row.pixelX > 0 && row.pixelY > 0 && row.dwellMs > 0
        % This completion metric is deliberately single-frame. MaxFrames is
        % a configured maximum, not proof of how many frames were acquired.
        row.expectedRasterMs = passFactor * row.pixelX * ...
            row.pixelY * row.dwellMs;
        if isfinite(row.elapsedMs) && row.expectedRasterMs > 0
            row.completionFraction = row.elapsedMs / row.expectedRasterMs;
        end
    end
end

function [eligible, reason] = preflightDecision(row, head, cfg)
    failures = {};
    if ~any(strcmpi(row.scanPlane, cfg.scanPlanes))
        failures{end + 1} = sprintf('scan plane %s is not selected (%s)', ...
            row.scanPlane, strjoin(cfg.scanPlanes, ', '));
    end
    failures = requireThreshold(failures, row.dwellMs, cfg.minDwellMs, ...
        'dwell time', 'ms/pixel', cfg.requireMetadataForPrefilter);
    failures = requireThreshold(failures, row.completionFraction, ...
        cfg.minCompletionFraction, 'empirical raster completion', '', ...
        cfg.requireMetadataForPrefilter);
    failures = requireThreshold(failures, row.numberOfRecords, cfg.minRecords, ...
        'TTTR record count', 'records', cfg.requireMetadataForPrefilter);
    if ~isfinite(row.headerMaxFrames)
        if cfg.requireMetadataForPrefilter
            failures{end + 1} = 'maximum-frame metadata is missing';
        end
    elseif row.headerMaxFrames ~= 1
        failures{end + 1} = sprintf([ ...
            'ImgHdr_MaxFrames=%g is not supported by the single-frame ' ...
            'completion prefilter'], row.headerMaxFrames);
    end
    if ~isfinite(row.pieWindows)
        if cfg.requireMetadataForPrefilter
            failures{end + 1} = 'PIE-window count is missing';
        end
    elseif row.pieWindows < cfg.requiredPieWindows
        failures{end + 1} = sprintf('only %g PIE windows are recorded', ...
            row.pieWindows);
    end
    wavelength = row.recordedExcitationNm;
    if ~isfinite(wavelength)
        if cfg.requireMetadataForPrefilter
            failures{end + 1} = sprintf('pulse-%d wavelength is missing', ...
                cfg.excitationPulseIndex);
        end
    elseif abs(wavelength - cfg.requiredExcitationNm) > cfg.wavelengthToleranceNm
        failures{end + 1} = sprintf('pulse %d is %g nm, not %g nm', ...
            cfg.excitationPulseIndex, wavelength, cfg.requiredExcitationNm);
    end
    if ~(isfield(head, 'ImgHdr_PixX') && isfield(head, 'ImgHdr_PixY'))
        failures{end + 1} = 'image dimensions are missing';
    end
    eligible = isempty(failures);
    if eligible
        reason = sprintf(['eligible: %s, dwell %.4g ms/pixel, completion %.3f, ' ...
            '%.4g million records'], row.scanPlane, row.dwellMs, ...
            row.completionFraction, row.numberOfRecords / 1e6);
    else
        reason = cleanText(strjoin(failures, '; '));
    end
end

function failures = requireThreshold(failures, value, threshold, label, unit, required)
    if ~isfinite(value)
        if required
            failures{end + 1} = sprintf('%s metadata is missing', label);
        end
    elseif value < threshold
        if isempty(unit)
            failures{end + 1} = sprintf('%s %.4g is below %.4g', ...
                label, value, threshold);
        else
            failures{end + 1} = sprintf('%s %.4g %s is below %.4g', ...
                label, value, unit, threshold);
        end
    end
end

function row = addResultQc(row, result, qc)
    reasons = {};
    if isfield(result, 'segmentation') && ...
            isfield(result.segmentation, 'diagnostics')
        diagnostics = result.segmentation.diagnostics;
    else
        error('batch_immune_cell_MIET:ResultSegmentation', ...
            'Existing result has no segmentation diagnostics.');
    end
    row.segmentationStatus = char(diagnostics.status);
    counts = diagnostics.pixelCounts;
    row.cellPixels = double(counts.cellMembrane);
    row.slbReferencePixels = double(counts.slbReference);
    imagePixels = double(counts.image);
    row.cellFraction = row.cellPixels / max(imagePixels, 1);
    row.slbReferenceFraction = row.slbReferencePixels / max(imagePixels, 1);
    row.lifetimeContrastNs = double(diagnostics.cell.contrastOverSlbNs);
    row.cellTouchesBorder = logical(diagnostics.cell.touchesImageBorder);
    row.slbReferencePhotons = double(diagnostics.slb.referencePhotonCount);
    row.slbRobustSigmaNs = double(diagnostics.slb.robustSigmaNs);
    if isfield(result, 'slbReference') && ...
            isfield(result.slbReference, 'fixedLifetimeNs')
        row.fixedSlbLifetimeNs = double(result.slbReference.fixedLifetimeNs);
        if isfield(result.slbReference, 'amplitudeConstraint')
            amplitude = result.slbReference.amplitudeConstraint;
            if isfield(amplitude, 'signalPhotonsPerPixel')
                row.slbSignalPhotonsPerPixel = ...
                    double(amplitude.signalPhotonsPerPixel);
            end
            if isfield(amplitude, 'dispersionPhotonsPerPixel')
                row.slbAmplitudeStdPhotonsPerPixel = ...
                    double(amplitude.dispersionPhotonsPerPixel);
            end
            if isfield(amplitude, 'bayesianClippedPixelFraction')
                row.slbClippedPixelFraction = ...
                    double(amplitude.bayesianClippedPixelFraction);
            end
        end
    end
    if isfield(result, 'bayesian') && isfield(result.bayesian, 'maps') && ...
            isfield(result.bayesian.maps, 'completeExponentialCountMAP')
        orderMap = result.bayesian.maps.completeExponentialCountMAP;
        row.bayesValidPixels = nnz(double(orderMap) > 0);
        row.bayesCoverage = row.bayesValidPixels / max(row.cellPixels, 1);
    end
    if isfield(result, 'bayesian') && ...
            isfield(result.bayesian, 'orderedComponents') && ...
            isfield(result.bayesian.orderedComponents, 'summary')
        orderedComponents = result.bayesian.orderedComponents;
        orderedSummary = orderedComponents.summary;
        row.secondComponentPixels = ...
            double(orderedSummary.secondDisplayPixelCount);
        row.thirdComponentPixels = ...
            double(orderedSummary.thirdDisplayPixelCount);
        row.membraneLifetimeMedianNs = ...
            displayedMembraneLifetimeMedian(orderedComponents);
    end
    if isfield(result, 'spatialBinning') && ...
            isfield(result.spatialBinning, 'bayesian') && ...
            isfield(result.spatialBinning.bayesian, 'orderedComponents') && ...
            isfield(result.spatialBinning.bayesian.orderedComponents, 'summary')
        binnedComponents = ...
            result.spatialBinning.bayesian.orderedComponents;
        binnedSummary = binnedComponents.summary;
        row.sliding2x2SecondComponentPixels = ...
            double(binnedSummary.secondDisplayPixelCount);
        row.sliding2x2ThirdComponentPixels = ...
            double(binnedSummary.thirdDisplayPixelCount);
        row.sliding2x2MembraneLifetimeMedianNs = ...
            displayedMembraneLifetimeMedian(binnedComponents);
    end
    if isfield(result, 'spatialBinning4x4') && ...
            isfield(result.spatialBinning4x4, 'bayesian') && ...
            isfield(result.spatialBinning4x4.bayesian, 'orderedComponents') && ...
            isfield(result.spatialBinning4x4.bayesian.orderedComponents, 'summary')
        binned4x4Components = ...
            result.spatialBinning4x4.bayesian.orderedComponents;
        binned4x4Summary = binned4x4Components.summary;
        row.sliding4x4SecondComponentPixels = ...
            double(binned4x4Summary.secondDisplayPixelCount);
        row.sliding4x4ThirdComponentPixels = ...
            double(binned4x4Summary.thirdDisplayPixelCount);
        row.sliding4x4MembraneLifetimeMedianNs = ...
            displayedMembraneLifetimeMedian(binned4x4Components);
    end

    if ~any(strcmp(row.segmentationStatus, qc.acceptedSegmentationStatuses))
        reasons{end + 1} = sprintf('segmentation status is %s', ...
            row.segmentationStatus);
    end
    if row.cellFraction < qc.minCellFraction
        reasons{end + 1} = sprintf('cell fraction %.3g is below %.3g', ...
            row.cellFraction, qc.minCellFraction);
    elseif row.cellFraction > qc.maxCellFraction
        reasons{end + 1} = sprintf('cell fraction %.3g exceeds %.3g', ...
            row.cellFraction, qc.maxCellFraction);
    end
    if qc.rejectCellTouchingBorder && row.cellTouchesBorder
        reasons{end + 1} = 'segmented cell footprint touches the image border';
    end
    if row.slbReferenceFraction < qc.minSlbReferenceFraction
        reasons{end + 1} = sprintf('SLB-reference fraction %.3g is below %.3g', ...
            row.slbReferenceFraction, qc.minSlbReferenceFraction);
    end
    if ~isfinite(row.slbReferencePhotons) || ...
            row.slbReferencePhotons < qc.minSlbReferencePhotons
        reasons{end + 1} = sprintf('SLB-reference photons %.4g are below %.4g', ...
            row.slbReferencePhotons, qc.minSlbReferencePhotons);
    end
    if ~isfinite(row.slbRobustSigmaNs) || ...
            row.slbRobustSigmaNs > qc.maxSlbRobustSigmaNs
        reasons{end + 1} = sprintf([ ...
            'outside-SLB robust spread %.3g ns exceeds %.3g ns; ' ...
            'the fixed component is not sufficiently homogeneous'], ...
            row.slbRobustSigmaNs, qc.maxSlbRobustSigmaNs);
    end
    if ~isfinite(row.lifetimeContrastNs) || ...
            row.lifetimeContrastNs < qc.minLifetimeContrastNs
        reasons{end + 1} = sprintf('mean-arrival contrast %.3g ns is below %.3g ns', ...
            row.lifetimeContrastNs, qc.minLifetimeContrastNs);
    end
    if ~isfinite(row.bayesCoverage) || row.bayesCoverage < qc.minBayesCoverage
        reasons{end + 1} = sprintf('Bayesian coverage %.3g is below %.3g', ...
            row.bayesCoverage, qc.minBayesCoverage);
    end
    if isfinite(row.slbClippedPixelFraction) && ...
            row.slbClippedPixelFraction > qc.maxSlbClippedFraction
        reasons{end + 1} = sprintf( ...
            'fixed SLB count clipped %.1f%% of cell pixels (maximum %.1f%%)', ...
            100 * row.slbClippedPixelFraction, ...
            100 * qc.maxSlbClippedFraction);
    end
    tauRange = qc.slbLifetimeRangeNs;
    if ~isfinite(row.fixedSlbLifetimeNs) || row.fixedSlbLifetimeNs < tauRange(1) || ...
            row.fixedSlbLifetimeNs > tauRange(2)
        reasons{end + 1} = sprintf('fixed SLB lifetime %.3g ns is outside [%g %g] ns', ...
            row.fixedSlbLifetimeNs, tauRange(1), tauRange(2));
    end
    row.goodCell = isempty(reasons);
    if row.goodCell
        row.qcReason = 'passed image-derived good-cell QC';
    else
        row.qcReason = cleanText(strjoin(reasons, '; '));
    end
end

function value = displayedMembraneLifetimeMedian(orderedComponents)
    value = NaN;
    if ~isfield(orderedComponents, 'display') || ...
            ~isfield(orderedComponents.display, 'secondLifetimeNs')
        return;
    end
    samples = double(orderedComponents.display.secondLifetimeNs(:));
    samples = samples(isfinite(samples));
    if ~isempty(samples)
        value = median(samples);
    end
end

function usable = existingResultUsable(matFile, result, requireTcspcPix, ...
        expectedSourceFile, expectedConfig)
    usable = isfield(result, 'segmentation') && ...
        isfield(result, 'analysisSchemaVersion') && ...
        double(result.analysisSchemaVersion) >= 8 && ...
        isfield(result.segmentation, 'diagnostics') && ...
        isfield(result.segmentation.diagnostics, 'pixelCounts') && ...
        isfield(result, 'slbReference') && ...
        isfield(result.slbReference, 'fixedLifetimeNs') && ...
        isfield(result, 'bayesian') && isfield(result.bayesian, 'maps') && ...
        isfield(result.bayesian.maps, 'completeExponentialCountMAP') && ...
        isfield(result.bayesian, 'orderedComponents');
    if ~usable || ~isfield(result, 'sourceFile') || ...
            ~strcmpi(char(result.sourceFile), char(expectedSourceFile)) || ...
            ~isfield(result, 'config') || ...
            ~structContainsExpectedConfig(result.config, ...
            analysisConfig(expectedConfig))
        usable = false;
        return;
    end
    amplitudeMode = '';
    if isfield(expectedConfig, 'slbAmplitudeMode') && ...
            ~isempty(expectedConfig.slbAmplitudeMode)
        amplitudeMode = lower(char(expectedConfig.slbAmplitudeMode));
    elseif isfield(expectedConfig, 'fixSlbAmplitude') && ...
            logical(expectedConfig.fixSlbAmplitude)
        amplitudeMode = 'fixed';
    end
    fixedCountRequired = ismember(amplitudeMode, {'fixed','regularized'});
    if strcmp(amplitudeMode, 'regularized')
        expectedConstraintMode = 'spatial expected photon-count prior';
    else
        expectedConstraintMode = 'direct expected photon-count allocation';
    end
    if fixedCountRequired && ~(isfield(result.bayesian, 'compact') && ...
            isfield(result.bayesian.compact, 'fixedSlbPhotonConstraint') && ...
            isfield(result.bayesian.compact.fixedSlbPhotonConstraint, 'mode') && ...
            strcmp(result.bayesian.compact.fixedSlbPhotonConstraint.mode, ...
                expectedConstraintMode) && ...
            isfield(result.slbReference, 'amplitudeConstraint') && ...
            isfield(result.slbReference.amplitudeConstraint, ...
                'bayesianClippedPixelFraction'))
        usable = false;
        return;
    end
    componentFiguresRequired = isfield(expectedConfig, 'componentMaps') && ...
        isstruct(expectedConfig.componentMaps) && ...
        isfield(expectedConfig.componentMaps, 'enabled') && ...
        logical(expectedConfig.componentMaps.enabled);
    % A run capped at one free membrane state never writes a third-component
    % figure, so requiring one here would defeat cfg.resume entirely.
    thirdFiguresExpected = ~isfield(expectedConfig, 'maxMembraneStates') || ...
        isempty(expectedConfig.maxMembraneStates) || ...
        double(expectedConfig.maxMembraneStates) >= 2;
    figureFields = {'meanFlimPng', 'bayesPng'};
    if componentFiguresRequired
        figureFields = [figureFields, ...
            {'secondSortedLifetimePng', 'sortedLifetimeLayersPng', ...
             'componentPhotonPng', 'componentLifetimePng'}];
        if thirdFiguresExpected
            figureFields = [figureFields, {'thirdSortedLifetimePng'}];
        end
    end
    if ~isfield(result, 'outputFiles') || ...
            ~all(isfield(result.outputFiles, figureFields))
        usable = false;
        return;
    end
    for figureIndex = 1:numel(figureFields)
        figurePath = result.outputFiles.(figureFields{figureIndex});
        if isempty(figurePath) || ~isfile(figurePath)
            usable = false;
            return;
        end
    end
    spatialBinningRequired = isfield(expectedConfig, 'spatialBinning') && ...
        isstruct(expectedConfig.spatialBinning) && ...
        isfield(expectedConfig.spatialBinning, 'enabled') && ...
        logical(expectedConfig.spatialBinning.enabled);
    if spatialBinningRequired
        usable = isfield(result, 'spatialBinning') && ...
            isfield(result.spatialBinning, 'status') && ...
            strcmp(result.spatialBinning.status, 'complete') && ...
            isfield(result.spatialBinning, 'bayesian') && ...
            isfield(result.spatialBinning.bayesian, 'maps') && ...
            isfield(result.spatialBinning.bayesian.maps, ...
                'completeExponentialCountMAP') && ...
            isfield(result.spatialBinning.bayesian, 'orderedComponents');
        if ~usable
            return;
        end
        binnedFigureFields = {'sliding2x2Step1MeanFlimPng', ...
            'sliding2x2Step1BayesPng'};
        if componentFiguresRequired
            binnedFigureFields = [binnedFigureFields, ...
                 {'sliding2x2Step1SecondSortedLifetimePng', ...
                  'sliding2x2Step1SortedLifetimeLayersPng', ...
                  'sliding2x2Step1ComponentPhotonPng', ...
                  'sliding2x2Step1ComponentLifetimePng'}];
            if thirdFiguresExpected
                binnedFigureFields = [binnedFigureFields, ...
                    {'sliding2x2Step1ThirdSortedLifetimePng'}];
            end
        end
        if ~isfield(result, 'outputFiles') || ...
                ~all(isfield(result.outputFiles, binnedFigureFields))
            usable = false;
            return;
        end
        for figureIndex = 1:numel(binnedFigureFields)
            figurePath = result.outputFiles.(binnedFigureFields{figureIndex});
            if isempty(figurePath) || ~isfile(figurePath)
                usable = false;
                return;
            end
        end
    end
    spatialBinning4x4Required = isfield(expectedConfig, 'spatialBinning4x4') && ...
        isstruct(expectedConfig.spatialBinning4x4) && ...
        isfield(expectedConfig.spatialBinning4x4, 'enabled') && ...
        logical(expectedConfig.spatialBinning4x4.enabled);
    if spatialBinning4x4Required
        usable = isfield(result, 'spatialBinning4x4') && ...
            isfield(result.spatialBinning4x4, 'status') && ...
            strcmp(result.spatialBinning4x4.status, 'complete') && ...
            isfield(result.spatialBinning4x4, 'bayesian') && ...
            isfield(result.spatialBinning4x4.bayesian, 'maps') && ...
            isfield(result.spatialBinning4x4.bayesian.maps, ...
                'completeExponentialCountMAP') && ...
            isfield(result.spatialBinning4x4.bayesian, 'orderedComponents');
        if ~usable
            return;
        end
        binned4x4FigureFields = {'sliding4x4Step1MeanFlimPng', ...
            'sliding4x4Step1BayesPng'};
        if componentFiguresRequired
            binned4x4FigureFields = [binned4x4FigureFields, ...
                 {'sliding4x4Step1SecondSortedLifetimePng', ...
                  'sliding4x4Step1SortedLifetimeLayersPng', ...
                  'sliding4x4Step1ComponentPhotonPng', ...
                  'sliding4x4Step1ComponentLifetimePng'}];
            if thirdFiguresExpected
                binned4x4FigureFields = [binned4x4FigureFields, ...
                    {'sliding4x4Step1ThirdSortedLifetimePng'}];
            end
        end
        if ~isfield(result, 'outputFiles') || ...
                ~all(isfield(result.outputFiles, binned4x4FigureFields))
            usable = false;
            return;
        end
        for figureIndex = 1:numel(binned4x4FigureFields)
            figurePath = result.outputFiles.(binned4x4FigureFields{figureIndex});
            if isempty(figurePath) || ~isfile(figurePath)
                usable = false;
                return;
            end
        end
    end
    if ~requireTcspcPix
        return;
    end
    variables = whos('-file', matFile);
    names = {variables.name};
    cubeIndex = find(strcmp(names, 'tcspc_pix'), 1);
    usable = ~isempty(cubeIndex) && ...
        strcmp(variables(cubeIndex).class, 'uint16');
end

function cfg = analysisConfig(cfg)
    operationalFields = {'outputDir', 'showFigure', 'showWaitbar', ...
        'saveOutputs', 'saveTcspcPix'};
    present = operationalFields(isfield(cfg, operationalFields));
    if ~isempty(present)
        cfg = rmfield(cfg, present);
    end
end

function matches = structContainsExpectedConfig(actual, expected)
    matches = isstruct(actual) && isscalar(actual) && ...
        isstruct(expected) && isscalar(expected);
    if ~matches
        return;
    end
    names = fieldnames(expected);
    for k = 1:numel(names)
        name = names{k};
        if ~isfield(actual, name)
            matches = false;
            return;
        end
        if isstruct(expected.(name)) && isscalar(expected.(name)) && ...
                isstruct(actual.(name)) && isscalar(actual.(name))
            matches = structContainsExpectedConfig(actual.(name), expected.(name));
        else
            matches = isequaln(actual.(name), expected.(name));
        end
        if ~matches
            return;
        end
    end
end

function [summary, batchInfo] = checkpoint(rows, cfg, dataRoot, outputDir, ...
        csvFile, matFile, errorLogFile)
    summary = struct2table(rows);
    batchInfo = buildBatchInfo(rows, dataRoot, outputDir, csvFile, ...
        matFile, errorLogFile);
    % Every per-acquisition MAT and PNG is already on disk by the time a
    % checkpoint runs, so failing to refresh the batch summary must not throw
    % away a completed run. Transient causes are common here: the data lives
    % on an external drive, and the CSV is often open in a viewer. Retry a
    % few times, then warn and carry on - the summary can be rebuilt by
    % rerunning with cfg.resume = true.
    writeWithRetry(@() writetable(summary, csvFile), csvFile, 'summary CSV');
    % save() names its variables as strings, so it must run in a scope that
    % actually holds them. An anonymous function only captures what it
    % lexically references, so wrapping save() directly would look for
    % summary, batchInfo and cfg in an empty workspace and silently write an
    % empty file. saveSummary takes them as arguments instead.
    writeWithRetry(@() saveSummary(matFile, summary, batchInfo, cfg), ...
        matFile, 'summary MAT');
end

function saveSummary(matFile, summary, batchInfo, cfg)
    save(matFile, 'summary', 'batchInfo', 'cfg', '-v7.3');
end

function writeWithRetry(writeFcn, targetFile, label)
    attempts = 4;
    for attempt = 1:attempts
        try
            writeFcn();
            return;
        catch writeError
            if attempt == attempts
                warning('batch_immune_cell_MIET:CheckpointWriteFailed', ...
                    ['Could not update the %s at %s after %d attempts: ' ...
                     '%s. Per-acquisition results are unaffected. Close ' ...
                     'any program holding the file, then rerun with ' ...
                     'cfg.resume = true to rebuild the summary.'], ...
                    label, targetFile, attempts, writeError.message);
                return;
            end
            pause(2^(attempt - 1));
        end
    end
end

function info = buildBatchInfo(rows, dataRoot, outputDir, csvFile, matFile, errorLog)
    statuses = {rows.status};
    info = struct();
    info.pipeline = 'immune_cell_MIET';
    info.generatedAt = datestr(now, 30);
    info.dataRoot = dataRoot;
    info.outputDir = outputDir;
    info.csvFile = csvFile;
    info.matFile = matFile;
    info.errorLogFile = errorLog;
    info.totalFiles = numel(rows);
    info.xyFiles = nnz(strcmp({rows.scanPlane}, 'XY'));
    info.xzFiles = nnz(strcmp({rows.scanPlane}, 'XZ'));
    info.yzFiles = nnz(strcmp({rows.scanPlane}, 'YZ'));
    info.eligibleFiles = nnz([rows.preflightEligible]);
    info.goodCellFiles = nnz([rows.goodCell]);
    info.failedFiles = nnz(contains(statuses, 'failed'));
    info.containsPhotonStream = false;
    info.scanConvention = [ ...
        'Dataset-derived Luminosa mapping: ImgHdr_ScanDirection ' ...
        '0=XY, 1=YZ, 2=XZ.'];
    info.completionDefinition = [ ...
        'TTResult_StopAfter divided by expected raster milliseconds; ' ...
        'unidirectional scans use an empirically validated factor of two.'];
end

function printSummary(info, csvFile)
    fprintf(['batch_immune_cell_MIET: %d total (%d XY, %d XZ, %d YZ); ' ...
        '%d metadata-eligible; %d good cells; %d failures.\n'], ...
        info.totalFiles, info.xyFiles, info.xzFiles, info.yzFiles, ...
        info.eligibleFiles, info.goodCellFiles, info.failedFiles);
    fprintf('batch_immune_cell_MIET: summary saved to %s\n', csvFile);
end

function appendErrorLog(logFile, ptuFile, errorObject)
    fid = fopen(logFile, 'a');
    if fid < 0
        warning('batch_immune_cell_MIET:ErrorLog', ...
            'Could not append the batch error log: %s', logFile);
        return;
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '\n[%s] %s\n%s: %s\n', datestr(now, 31), ptuFile, ...
        errorObject.identifier, errorObject.message);
    for k = 1:numel(errorObject.stack)
        fprintf(fid, '  at %s (line %d)\n', ...
            errorObject.stack(k).name, errorObject.stack(k).line);
    end
    clear cleanup
end

function closeNewPipelineFigures(figuresBefore)
    figuresAfter = findall(groot, 'Type', 'figure');
    for k = 1:numel(figuresAfter)
        if ~any(figuresBefore == figuresAfter(k)) && isgraphics(figuresAfter(k))
            try
                close(figuresAfter(k));
            catch
            end
        end
    end
end

function value = numericScalar(s, fieldName)
    value = NaN;
    if isstruct(s) && isfield(s, fieldName) && ...
            (isnumeric(s.(fieldName)) || islogical(s.(fieldName))) && ...
            ~isempty(s.(fieldName))
        candidate = double(s.(fieldName)(1));
        if isfinite(candidate)
            value = candidate;
        end
    end
end

function text = cleanText(text)
    if isstring(text)
        text = char(text);
    end
    text = strrep(text, sprintf('\r'), ' ');
    text = strrep(text, newline, ' | ');
    text = strtrim(text);
end

function requireFunction(name)
    if exist(name, 'file') ~= 2
        error('batch_immune_cell_MIET:MissingDependency', ...
            'Required function %s.m is not on the MATLAB path.', name);
    end
end

function out = mergeStruct(out, overrides)
    if isempty(overrides)
        return;
    end
    names = fieldnames(overrides);
    for k = 1:numel(names)
        out.(names{k}) = overrides.(names{k});
    end
end
