function batch = fit_tcspc_subfolders_with_fluofit(rootFolder, outputFolder, opts)
%FIT_TCSPC_SUBFOLDERS_WITH_FLUOFIT Batch fit TCSPC files in subfolders.
%
% Usage:
%   batch = fit_tcspc_subfolders_with_fluofit('D:\Luminosa\Data\Natasha')
%
% Expected layout:
%   rootFolder\
%       sample_001\...\one_file.ptu
%       sample_002\...\one_file.pqres
%       ...
%
% The function reads one PTU or PQRES file from inside each main subfolder,
% collects the TCSPC curves, fits every curve with Fluofit or the faster
% tail-fit path, using either a per-file Calc_mIRF IRF, a supplied IRF, or
% an optional global IRF when reconvolution fitting is selected, and
% saves:
%   - one MAT file per fitted curve in the main subfolder result folder
%   - root-level batch_irf.mat
%   - root-level global_irf.mat when opts.irfMode='global'
%   - root-level tcspc_curves_cache.mat for reusing read TCSPC curves
%   - root-level tcspc_batch_fluofit_summary.csv
%   - root-level tcspc_batch_fluofit_all.mat
%
% Main options:
%   opts.tau0                  lifetime starts in ns. In Fluofit mode,
%                              default [] lets DistFluofit choose the
%                              component seeds.
%   opts.limits                Fluofit lifetime limits in ns, default []
%   opts.init                  Fluofit init flag, default 1 when tau0 is
%                              empty, otherwise 0
%   opts.fitMethod             'fluofit' or 'tail'. 'tail' skips IRF
%                              estimation/reconvolution, cuts from the
%                              TCSPC peak plus opts.tailCutAfterPeakNs, and
%                              runs Tailfit, optionally with DistTailfit
%                              when opts.init>0. Default 'fluofit'
%   opts.tailCutAfterPeakNs    tail-fit start after peak, default 0.5 ns
%   opts.tailRejectLastNPoints remove this many final TCSPC bins from the
%                              tail fit after the peak gate, default 0
%   opts.tau0                  in tail mode, Tailfit lifetime initial
%                              guess in ns. opts.tailTau0 is accepted as
%                              a backward-compatible alias when tau0=[]
%   opts.fluofitSolver         'mle', 'ls', or 'pirls', default 'mle'
%   opts.plotFits              draw Fluofit figures, default false
%   opts.saveFitFigures        save each fit figure, default opts.plotFits
%   opts.figureFormats         figure formats, default {'png', 'fig'}
%   opts.copyFitsToBatchFolder copy each per-subfolder fit MAT/figure file
%                              into the common batch results folder,
%                              default true
%   opts.summaryAmplitudeThreshold
%                              omit components below this relative
%                              amplitude from the summary CSV only,
%                              default 0.02. Use 0 to disable.
%   opts.errorMethod           'none', 'bootstrap', or 'chunk'. Bootstrap
%                              Poisson-resamples the fitted/observed TCSPC
%                              histogram and refits each sample. Chunk uses
%                              Poisson-thinned histogram chunks. Default
%                              'none'.
%   opts.nErrorSamples         bootstrap sample count, default 25
%   opts.errorChunks           chunk count when opts.errorMethod='chunk',
%                              default 8
%   opts.errorResampleSource   'fit' or 'counts', default 'fit'
%   opts.irfMode               'per_curve', 'best_per_curve', 'global', or
%                              'supplied', default 'per_curve'.
%                              'best_per_curve' estimates one IRF per file,
%                              selects the candidate with the lowest
%                              bi-exponential screening chi2, then uses
%                              that selected IRF for every final fit.
%   opts.bestIrfTau0           lifetime starts in ns for the best-IRF
%                              screening fit, default [0.4 2.0]
%   opts.bestIrfSolver         Fluofit solver for best-IRF screening,
%                              default opts.fluofitSolver
%   opts.perCurveIrfModel      'spad_exgauss' or 'calc_mirf', default
%                              'spad_exgauss'. The SPAD model is a
%                              Gaussian prompt peak convolved with a
%                              one-sided exponential tail.
%   opts.irfFile               supplied PTU/PQRES IRF file for
%                              opts.irfMode='supplied'
%   opts.globalIrfMethod       'calc_mirf', 'spad_exgauss',
%                              'gamma_shifted_fast', or 'exgauss_fast',
%                              default 'calc_mirf'
%   opts.irfClipFraction       for Calc_mIRF IRFs, set values below this
%                              fraction of max(IRF) to zero before
%                              fitting, default 1e-3
%   opts.searchRecursively     search inside sub-subfolders, default true
%   opts.multipleFileMode      'largest_ptu', 'error', or 'first',
%                              default 'largest_ptu'. With recursive
%                              folder data, 'largest_ptu' selects the
%                              largest PTU candidate under each main
%                              subfolder.
%   opts.resultFolderName      created in each main subfolder, default
%                              'tcspc_fluofit_results'
%   opts.useTcspcCache         load/save read TCSPC curves, default true
%   opts.forceReadTcspc        ignore existing TCSPC cache, default false
%   opts.tcspcCacheFile        default root-level tcspc_curves_cache.mat
%   opts.rejectFirstTimePoint  default true for the Natasha dataset
%   opts.rejectTailAtOrAfterNs default 12.5; set [] to disable
%   opts.rejectTailNPoints     default 8

if nargin < 1 || isempty(rootFolder)
    rootFolder = 'D:\Luminosa\Data\Natasha';
end
if nargin < 2
    outputFolder = '';
end
if nargin < 3 || isempty(opts)
    opts = struct();
end
opts = defaultBatchOptions(opts);
if isempty(outputFolder)
    outputFolder = fullfile(rootFolder, opts.batchResultFolderName);
end
if isempty(opts.tcspcCacheFile)
    opts.tcspcCacheFile = fullfile(outputFolder, 'tcspc_curves_cache.mat');
end

if ~isfolder(rootFolder)
    error('Input root folder does not exist: %s', rootFolder);
end

scriptFolder = fileparts(mfilename('fullpath'));
if ~isempty(scriptFolder)
    addpath(scriptFolder);
end
if strcmp(opts.fitMethod, 'tail')
    if opts.runDistTailfit && exist('DistTailfit', 'file') ~= 2
        error('Could not find DistTailfit.m on the MATLAB path.');
    end
    if exist('Tailfit', 'file') ~= 2
        error('Could not find Tailfit.m on the MATLAB path.');
    end
    if exist('PIRLSnonneg', 'file') ~= 2
        error('Could not find PIRLSnonneg.m on the MATLAB path.');
    end
else
    if exist('Fluofit', 'file') ~= 2
        error('Could not find Fluofit.m on the MATLAB path.');
    end
    if exist('Calc_mIRF', 'file') ~= 2
        error('Could not find Calc_mIRF.m on the MATLAB path.');
    end
end
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

targets = collectSubfolderTargets(rootFolder, opts);
if isempty(targets)
    error('No subfolders containing one .ptu or .pqres file were found in %s', rootFolder);
end
fprintf('Found %d TCSPC target folder(s) under %s\n', numel(targets), rootFolder);

curves = repmat(emptyCurve(), numel(targets), 1);
fits = repmat(emptyFit(), numel(targets), 1);
summary = repmat(emptySummary(), numel(targets), 1);
tcspcCacheLoaded = false;
tcspcCacheMessage = '';

if opts.useTcspcCache
    [tcspcCacheLoaded, curves, summary, tcspcCacheMessage] = ...
        loadTcspcCache(opts.tcspcCacheFile, targets, opts);
    if tcspcCacheLoaded
        fprintf('Loaded TCSPC curves from cache: %s\n', opts.tcspcCacheFile);
    elseif ~isempty(tcspcCacheMessage)
        fprintf('TCSPC cache not used: %s\n', tcspcCacheMessage);
    end
end

if ~tcspcCacheLoaded
    [curves, summary] = readAllTcspcTargets(targets, curves, summary, opts);
    if opts.useTcspcCache
        saveTcspcCache(opts.tcspcCacheFile, rootFolder, outputFolder, targets, curves, summary, opts);
        fprintf('Saved TCSPC curves cache to %s\n', opts.tcspcCacheFile);
    end
end

if strcmp(opts.fitMethod, 'tail')
    batchIrf = emptyNoIrfBatch('Tail fitting does not use an IRF.');
    fprintf('Tail fitting selected; skipping IRF estimation.\n');
else
    batchIrf = loadBatchIrf(curves, opts);
    if strcmp(batchIrf.status, 'ok')
        batchIrfPath = fullfile(outputFolder, 'batch_irf.mat');
        save(batchIrfPath, 'batchIrf');
        fprintf('Saved batch IRF to %s\n', batchIrfPath);
        if strcmp(opts.irfMode, 'global')
            globalIrf = batchIrf; %#ok<NASGU>
            globalIrfPath = fullfile(outputFolder, 'global_irf.mat');
            save(globalIrfPath, 'globalIrf');
            fprintf('Saved global IRF to %s\n', globalIrfPath);
        elseif strcmp(opts.irfMode, 'best_per_curve')
            bestIrf = batchIrf; %#ok<NASGU>
            bestIrfPath = fullfile(outputFolder, 'best_irf.mat');
            save(bestIrfPath, 'bestIrf');
            fprintf('Saved selected best IRF to %s\n', bestIrfPath);
        end
    else
        fprintf('Batch IRF not available: %s\n', batchIrf.message);
    end
end
if strcmp(opts.fitMethod, 'fluofit') && strcmp(opts.irfMode, 'best_per_curve') && ~strcmp(batchIrf.status, 'ok')
    error('Best per-curve IRF selection failed: %s', batchIrf.message);
end

for kk = 1:numel(curves)
    if ~strcmp(curves(kk).status, 'ok')
        continue;
    end

    fprintf('Fitting %d/%d: %s\n', kk, numel(curves), curves(kk).file);
    fits(kk).folder = curves(kk).folder;
    fits(kk).file = curves(kk).file;
    fits(kk).status = 'failed';

    try
        if strcmp(opts.fitMethod, 'tail')
            fitInput = buildBatchTailFitInput(curves(kk), opts);
            [fitResult, fitFigures] = runBatchTailFit(fitInput, opts);
        else
            fitInput = buildBatchFitInput(curves(kk), batchIrf, opts);
            [fitResult, fitFigures] = runBatchFluofit(fitInput, opts);
        end
        fitResult = addFitUncertainty(fitInput, fitResult, opts);

        [~, stem] = fileparts(curves(kk).file);
        outStem = safeFileStem([curves(kk).folderName '_' stem]);
        if ~exist(curves(kk).resultFolder, 'dir')
            mkdir(curves(kk).resultFolder);
        end
        fitResult.figureFiles = {};
        if opts.saveFitFigures
            fitResult.figureFiles = saveFitFigures( ...
                fitFigures, curves(kk).resultFolder, outStem, fitInput, opts);
        end
        if isfield(fitResult, 'fileSuffix') && ~isempty(fitResult.fileSuffix)
            fileSuffix = fitResult.fileSuffix;
        else
            fileSuffix = 'fluofit';
        end
        outPath = fullfile(curves(kk).resultFolder, [outStem '_' fileSuffix '.mat']);
        save(outPath, 'fitInput', 'fitResult', '-v7.3');
        commonFiles = struct('outputFile', '', 'figureFiles', {{}});
        if opts.copyFitsToBatchFolder
            commonFiles = copyFitOutputsToBatchFolder( ...
                outPath, fitResult.figureFiles, outputFolder, outStem);
        end

        fits(kk).fitInput = fitInput;
        fits(kk).fitResult = fitResult;
        fits(kk).outputFile = outPath;
        fits(kk).figureFiles = fitResult.figureFiles;
        fits(kk).commonOutputFile = commonFiles.outputFile;
        fits(kk).commonFigureFiles = commonFiles.figureFiles;
        fits(kk).status = 'ok';
        fits(kk).message = '';

        summary(kk).fitMethod = fitResult.fitMethod;
        summary(kk).irfSource = fitInput.irfSource;
        summary(kk).outputFile = outPath;
        summary(kk).figureFiles = strjoin(fitResult.figureFiles, ';');
        summary(kk).commonOutputFile = commonFiles.outputFile;
        summary(kk).commonFigureFiles = strjoin(commonFiles.figureFiles, ';');
        summary(kk).dtNs = fitInput.dtNs;
        summary(kk).pulsePeriodNs = fitInput.pulsePeriodNs;
        summary(kk).nBins = numel(fitInput.counts);
        summary(kk).nRejected = fitInput.gateInfo.nRejected;
        summaryComponents = filterSummaryComponents( ...
            fitResult.tauNs, fitResult.amplitudes, opts.summaryAmplitudeThreshold);
        summary(kk).tauNs = mat2str(summaryComponents.tauNs(:).', 5);
        summary(kk).amplitudes = mat2str(summaryComponents.amplitudes(:).', 5);
        summary(kk).relativeAmplitudes = mat2str(summaryComponents.relativeAmplitudes(:).', 5);
        summary(kk).nComponents = summaryComponents.nComponents;
        summary(kk).nReportedComponents = summaryComponents.nReportedComponents;
        summary(kk).nRemovedLowAmplitude = summaryComponents.nRemovedLowAmplitude;
        summary(kk).tauErrorNs = mat2str( ...
            filterVectorForSummary(getResampleErrorVector(fitResult, 'tauStdNs'), ...
            summaryComponents.keepMask), 5);
        summary(kk).amplitudeError = mat2str( ...
            filterVectorForSummary(getResampleErrorVector(fitResult, 'amplitudeStd'), ...
            summaryComponents.keepMask), 5);
        summary(kk).relativeAmplitudeError = mat2str( ...
            filterVectorForSummary(getResampleErrorVector(fitResult, 'relativeAmplitudeStd'), ...
            summaryComponents.keepMask), 5);
        [summary(kk).errorMethod, summary(kk).nErrorSamples, summary(kk).nErrorSamplesOk, ...
            summary(kk).errorStatus, summary(kk).errorMessage] = ...
            summarizeResampleErrorStatus(fitResult);
        summary(kk).chi = fitResult.chi;
        summary(kk).status = 'ok';
        summary(kk).message = '';
    catch ME
        fits(kk).message = ME.message;
        summary(kk).status = 'failed';
        summary(kk).message = ME.message;
        warning('Fit failed for %s: %s', curves(kk).file, ME.message);
    end
end

summaryPath = fullfile(outputFolder, 'tcspc_batch_fluofit_summary.csv');
writetable(struct2table(summary, 'AsArray', true), summaryPath);

batch = struct();
batch.rootFolder = rootFolder;
batch.outputFolder = outputFolder;
batch.options = opts;
batch.targets = targets;
batch.curves = curves;
batch.batchIrf = batchIrf;
batch.globalIrf = batchIrf;
batch.bestIrf = batchIrf;
batch.fits = fits;
batch.results = summary;
batch.summaryPath = summaryPath;
batch.tcspcCacheFile = opts.tcspcCacheFile;
batch.tcspcCacheLoaded = tcspcCacheLoaded;
batch.tcspcCacheMessage = tcspcCacheMessage;

allPath = fullfile(outputFolder, 'tcspc_batch_fluofit_all.mat');
save(allPath, 'batch', '-v7.3');
fprintf('Saved summary to %s\n', summaryPath);
fprintf('Saved combined batch data to %s\n', allPath);
end

function opts = defaultBatchOptions(opts)
tau0Provided = isfield(opts, 'tau0') && ~isempty(opts.tau0);
if ~isfield(opts, 'tau0')
    opts.tau0 = [];
end
if ~isfield(opts, 'limits')
    opts.limits = [];
end
if ~isfield(opts, 'init') || isempty(opts.init)
    opts.init = double(~tau0Provided);
end
if ~isfield(opts, 'fitMethod') || isempty(opts.fitMethod)
    opts.fitMethod = 'fluofit';
end
opts.fitMethod = lower(strrep(char(opts.fitMethod), '-', '_'));
if any(strcmp(opts.fitMethod, {'tailfit', 'tail_fit', 'dist_tail', 'disttail'}))
    opts.fitMethod = 'tail';
end
if ~any(strcmp(opts.fitMethod, {'fluofit', 'tail'}))
    error('opts.fitMethod must be ''fluofit'' or ''tail''.');
end
if strcmp(opts.fitMethod, 'fluofit') && isempty(opts.tau0) && opts.init == 0
    opts.init = 1;
end
if ~isfield(opts, 'fluofitSolver') || isempty(opts.fluofitSolver)
    opts.fluofitSolver = 'mle';
end
opts.fluofitSolver = lower(strrep(char(opts.fluofitSolver), '-', '_'));
validSolvers = {'mle', 'ml', 'maximum_likelihood', 'ls', 'least_squares', ...
    'pirls', 'pirlsnonneg', 'poisson_irls'};
if ~any(strcmp(opts.fluofitSolver, validSolvers))
    error('opts.fluofitSolver must be ''mle'', ''ls'', or ''pirls''.');
end
if ~isfield(opts, 'plotFits') || isempty(opts.plotFits)
    opts.plotFits = false;
end
if ~isfield(opts, 'saveFitFigures') || isempty(opts.saveFitFigures)
    opts.saveFitFigures = opts.plotFits;
end
if ~isfield(opts, 'figureFormats') || isempty(opts.figureFormats)
    opts.figureFormats = {'png', 'fig'};
end
opts.figureFormats = normalizeFigureFormats(opts.figureFormats);
if ~isfield(opts, 'copyFitsToBatchFolder') || isempty(opts.copyFitsToBatchFolder)
    opts.copyFitsToBatchFolder = true;
end
if ~isfield(opts, 'summaryAmplitudeThreshold') || isempty(opts.summaryAmplitudeThreshold)
    opts.summaryAmplitudeThreshold = 0.02;
end
opts.summaryAmplitudeThreshold = double(opts.summaryAmplitudeThreshold);
if ~isscalar(opts.summaryAmplitudeThreshold) || ~isfinite(opts.summaryAmplitudeThreshold) || opts.summaryAmplitudeThreshold < 0
    error('opts.summaryAmplitudeThreshold must be a finite non-negative scalar.');
end
if opts.summaryAmplitudeThreshold > 1
    if opts.summaryAmplitudeThreshold <= 100
        opts.summaryAmplitudeThreshold = opts.summaryAmplitudeThreshold / 100;
    else
        error('opts.summaryAmplitudeThreshold must be given as a fraction or percentage.');
    end
end
if ~isfield(opts, 'errorMethod') || isempty(opts.errorMethod)
    opts.errorMethod = 'none';
end
opts.errorMethod = lower(strrep(char(opts.errorMethod), '-', '_'));
if any(strcmp(opts.errorMethod, {'off', 'disable', 'disabled', 'false'}))
    opts.errorMethod = 'none';
elseif any(strcmp(opts.errorMethod, {'poisson', 'poisson_bootstrap', 'boot'}))
    opts.errorMethod = 'bootstrap';
end
if ~any(strcmp(opts.errorMethod, {'none', 'bootstrap', 'chunk'}))
    error('opts.errorMethod must be ''none'', ''bootstrap'', or ''chunk''.');
end
if ~isfield(opts, 'nErrorSamples') || isempty(opts.nErrorSamples)
    opts.nErrorSamples = 25;
end
opts.nErrorSamples = double(opts.nErrorSamples);
if ~isscalar(opts.nErrorSamples) || ~isfinite(opts.nErrorSamples) || opts.nErrorSamples < 0
    error('opts.nErrorSamples must be a finite non-negative scalar.');
end
opts.nErrorSamples = round(opts.nErrorSamples);
if ~isfield(opts, 'errorChunks') || isempty(opts.errorChunks)
    opts.errorChunks = 8;
end
opts.errorChunks = double(opts.errorChunks);
if ~isscalar(opts.errorChunks) || ~isfinite(opts.errorChunks) || opts.errorChunks < 0
    error('opts.errorChunks must be a finite non-negative scalar.');
end
opts.errorChunks = round(opts.errorChunks);
if ~isfield(opts, 'errorResampleSource') || isempty(opts.errorResampleSource)
    opts.errorResampleSource = 'fit';
end
opts.errorResampleSource = lower(strrep(char(opts.errorResampleSource), '-', '_'));
if ~any(strcmp(opts.errorResampleSource, {'fit', 'counts', 'data'}))
    error('opts.errorResampleSource must be ''fit'' or ''counts''.');
end
if strcmp(opts.errorResampleSource, 'data')
    opts.errorResampleSource = 'counts';
end
if ~isfield(opts, 'errorRandomSeed')
    opts.errorRandomSeed = [];
end
if ~isempty(opts.errorRandomSeed)
    opts.errorRandomSeed = double(opts.errorRandomSeed);
    if ~isscalar(opts.errorRandomSeed) || ~isfinite(opts.errorRandomSeed) || opts.errorRandomSeed < 0
        error('opts.errorRandomSeed must be empty or a finite non-negative scalar.');
    end
end
if ~isfield(opts, 'irfMode') || isempty(opts.irfMode)
    opts.irfMode = 'per_curve';
end
opts.irfMode = lower(strrep(char(opts.irfMode), '-', '_'));
if ~any(strcmp(opts.irfMode, {'global', 'supplied', 'per_curve', 'best_per_curve'}))
    error('opts.irfMode must be ''global'', ''supplied'', ''per_curve'', or ''best_per_curve''.');
end
if ~isfield(opts, 'bestIrfTau0') || isempty(opts.bestIrfTau0)
    opts.bestIrfTau0 = [0.4 2.0];
end
opts.bestIrfTau0 = double(opts.bestIrfTau0(:)).';
if numel(opts.bestIrfTau0) ~= 2 || any(~isfinite(opts.bestIrfTau0)) || any(opts.bestIrfTau0 <= 0)
    error('opts.bestIrfTau0 must contain exactly two positive lifetimes in ns.');
end
if ~isfield(opts, 'bestIrfLimits')
    opts.bestIrfLimits = [];
end
if ~isfield(opts, 'bestIrfInit') || isempty(opts.bestIrfInit)
    opts.bestIrfInit = 0;
end
if ~isfield(opts, 'bestIrfSolver') || isempty(opts.bestIrfSolver)
    opts.bestIrfSolver = opts.fluofitSolver;
end
opts.bestIrfSolver = lower(strrep(char(opts.bestIrfSolver), '-', '_'));
if ~any(strcmp(opts.bestIrfSolver, validSolvers))
    error('opts.bestIrfSolver must be ''mle'', ''ls'', or ''pirls''.');
end
if ~isfield(opts, 'tailCutAfterPeakNs') || isempty(opts.tailCutAfterPeakNs)
    opts.tailCutAfterPeakNs = 0.5;
end
opts.tailCutAfterPeakNs = double(opts.tailCutAfterPeakNs);
if ~isscalar(opts.tailCutAfterPeakNs) || ~isfinite(opts.tailCutAfterPeakNs) || opts.tailCutAfterPeakNs < 0
    error('opts.tailCutAfterPeakNs must be a finite non-negative scalar in ns.');
end
if ~isfield(opts, 'tailRejectLastNPoints') || isempty(opts.tailRejectLastNPoints)
    opts.tailRejectLastNPoints = 0;
end
opts.tailRejectLastNPoints = round(double(opts.tailRejectLastNPoints));
if ~isscalar(opts.tailRejectLastNPoints) || ~isfinite(opts.tailRejectLastNPoints) || opts.tailRejectLastNPoints < 0
    error('opts.tailRejectLastNPoints must be a finite non-negative scalar.');
end
if strcmp(opts.fitMethod, 'tail') && ~isempty(opts.tau0)
    opts.tailTau0 = opts.tau0;
elseif ~isfield(opts, 'tailTau0') || isempty(opts.tailTau0)
    opts.tailTau0 = opts.bestIrfTau0;
end
opts.tailTau0 = double(opts.tailTau0(:)).';
if isempty(opts.tailTau0) || any(~isfinite(opts.tailTau0)) || any(opts.tailTau0 <= 0)
    error('opts.tailTau0 must contain positive finite lifetimes in ns.');
end
if ~isfield(opts, 'tailLimits')
    opts.tailLimits = opts.limits;
end
if ~isfield(opts, 'tailSolver') || isempty(opts.tailSolver)
    opts.tailSolver = opts.fluofitSolver;
end
opts.tailSolver = lower(strrep(char(opts.tailSolver), '-', '_'));
if ~any(strcmp(opts.tailSolver, validSolvers))
    error('opts.tailSolver must be ''mle'', ''ls'', or ''pirls''.');
end
if ~isfield(opts, 'tailDistFlag') || isempty(opts.tailDistFlag)
    opts.tailDistFlag = 1;
end
if ~isfield(opts, 'tailDistN') || isempty(opts.tailDistN)
    opts.tailDistN = 200;
end
if ~isfield(opts, 'tailDistTauMaxNs') || isempty(opts.tailDistTauMaxNs)
    opts.tailDistTauMaxNs = [];
end
if ~isfield(opts, 'tailSimplexSteps') || isempty(opts.tailSimplexSteps)
    opts.tailSimplexSteps = [];
end
if ~isfield(opts, 'tailUsePeriodCorrection') || isempty(opts.tailUsePeriodCorrection)
    opts.tailUsePeriodCorrection = true;
end
if ~isfield(opts, 'runDistTailfit') || isempty(opts.runDistTailfit)
    opts.runDistTailfit = opts.init > 0;
end
if ~isfield(opts, 'perCurveIrfModel') || isempty(opts.perCurveIrfModel)
    opts.perCurveIrfModel = 'spad_exgauss';
end
opts.perCurveIrfModel = lower(strrep(char(opts.perCurveIrfModel), '-', '_'));
validPerCurveIrfModels = {'spad_exgauss', 'spad_exgaussian', ...
    'exgauss', 'ex_gauss', 'calc_mirf', 'calcmirf'};
if ~any(strcmp(opts.perCurveIrfModel, validPerCurveIrfModels))
    error('opts.perCurveIrfModel must be ''spad_exgauss'' or ''calc_mirf''.');
end
if ~isfield(opts, 'irfFile') || isempty(opts.irfFile)
    opts.irfFile = '';
end
if isfield(opts, 'suppliedIrfFile') && ~isempty(opts.suppliedIrfFile)
    opts.irfFile = opts.suppliedIrfFile;
end
if strcmp(opts.irfMode, 'supplied') && isempty(opts.irfFile)
    error('opts.irfFile must be set when opts.irfMode is ''supplied''.');
end
if ~isfield(opts, 'globalIrfMethod') || isempty(opts.globalIrfMethod)
    opts.globalIrfMethod = 'calc_mirf';
end
opts.globalIrfMethod = lower(strrep(char(opts.globalIrfMethod), '-', '_'));
validGlobalMethods = {'calc_mirf', 'calcmirf', 'gamma_shifted_fast', ...
    'gamma_fast', 'exgauss_fast', 'ex_gauss_fast', ...
    'spad_exgauss', 'spad_exgaussian'};
if ~any(strcmp(opts.globalIrfMethod, validGlobalMethods))
    error('Unsupported opts.globalIrfMethod: %s', opts.globalIrfMethod);
end
if ~isfield(opts, 'globalIrfOptions') || isempty(opts.globalIrfOptions)
    opts.globalIrfOptions = struct();
end
if ~isfield(opts, 'spadIrfOptions') || isempty(opts.spadIrfOptions)
    opts.spadIrfOptions = struct();
end
if ~isfield(opts, 'globalIrfResample') || isempty(opts.globalIrfResample)
    opts.globalIrfResample = true;
end
if ~isfield(opts, 'globalIrfPeriodToleranceSec') || isempty(opts.globalIrfPeriodToleranceSec)
    opts.globalIrfPeriodToleranceSec = 1e-12;
end
if ~isfield(opts, 'irfClipFraction')
    opts.irfClipFraction = 1e-3;
end
if ~isfield(opts, 'fallbackPerCurveIrf') || isempty(opts.fallbackPerCurveIrf)
    opts.fallbackPerCurveIrf = true;
end
if ~isfield(opts, 'batchResultFolderName') || isempty(opts.batchResultFolderName)
    opts.batchResultFolderName = 'tcspc_batch_fluofit_results';
end
if ~isfield(opts, 'resultFolderName') || isempty(opts.resultFolderName)
    opts.resultFolderName = 'tcspc_fluofit_results';
end
if ~isfield(opts, 'useTcspcCache') || isempty(opts.useTcspcCache)
    opts.useTcspcCache = true;
end
if ~isfield(opts, 'forceReadTcspc') || isempty(opts.forceReadTcspc)
    opts.forceReadTcspc = false;
end
if ~isfield(opts, 'tcspcCacheFile') || isempty(opts.tcspcCacheFile)
    opts.tcspcCacheFile = '';
end
if ~isfield(opts, 'searchRecursively') || isempty(opts.searchRecursively)
    opts.searchRecursively = true;
end
if ~isfield(opts, 'maxSearchDepth') || isempty(opts.maxSearchDepth)
    opts.maxSearchDepth = Inf;
end
if ~isfield(opts, 'includeRootFiles') || isempty(opts.includeRootFiles)
    opts.includeRootFiles = false;
end
if ~isfield(opts, 'multipleFileMode') || isempty(opts.multipleFileMode)
    opts.multipleFileMode = 'largest_ptu';
end
opts.multipleFileMode = lower(strrep(char(opts.multipleFileMode), '-', '_'));
validMultipleFileModes = {'largest_ptu', 'largestptu', 'largest', 'error', 'first'};
if ~any(strcmp(opts.multipleFileMode, validMultipleFileModes))
    error('opts.multipleFileMode must be ''largest_ptu'', ''error'', or ''first''.');
end
if ~isfield(opts, 'ptuResolutionSec')
    opts.ptuResolutionSec = [];
end
if ~isfield(opts, 'ptuChunkRecords') || isempty(opts.ptuChunkRecords)
    opts.ptuChunkRecords = 1e6;
end
if ~isfield(opts, 'detectorChannels')
    opts.detectorChannels = [];
end
if ~isfield(opts, 'rejectFirstTimePoint') || isempty(opts.rejectFirstTimePoint)
    opts.rejectFirstTimePoint = true;
end
if ~isfield(opts, 'rejectTailAtOrAfterNs')
    opts.rejectTailAtOrAfterNs = 12.5;
end
if ~isfield(opts, 'rejectTailNPoints') || isempty(opts.rejectTailNPoints)
    opts.rejectTailNPoints = 8;
end
end

function formats = normalizeFigureFormats(formats)
if ischar(formats) || isstring(formats)
    formats = cellstr(formats);
elseif ~iscell(formats)
    error('opts.figureFormats must be a character vector, string, or cell array.');
end

out = {};
for ii = 1:numel(formats)
    fmt = lower(strrep(strtrim(char(formats{ii})), '.', ''));
    if isempty(fmt)
        continue;
    end
    out{end + 1} = fmt; %#ok<AGROW>
end

if isempty(out)
    error('opts.figureFormats must contain at least one non-empty format.');
end
formats = out;
end

function target = emptyTarget()
target = struct('folder', '', 'folderName', '', 'file', '', 'fileType', '', ...
    'resultFolder', '', 'status', 'pending', 'message', '');
end

function curve = emptyCurve()
curve = struct('folder', '', 'folderName', '', 'file', '', 'fileType', '', ...
    'curveName', '', 'timeSec', [], 'counts', [], 'meta', struct(), ...
    'head', [], 'rawTcspc', [], 'nChannels', NaN, 'status', 'pending', ...
    'resultFolder', '', 'message', '');
end

function batchIrf = emptyNoIrfBatch(message)
batchIrf = struct('status', 'skipped', 'message', message, 'method', 'none', ...
    'timeSec', [], 'counts', [], 'summedCounts', [], 'included', [], ...
    'excluded', [], 'head', [], 'params', []);
end

function fit = emptyFit()
fit = struct('folder', '', 'file', '', 'fitInput', [], 'fitResult', [], ...
    'outputFile', '', 'commonOutputFile', '', 'status', 'pending', 'message', '');
fit.figureFiles = {};
fit.commonFigureFiles = {};
end

function summary = emptySummary()
summary = struct('folder', '', 'file', '', 'fileType', '', 'fitMethod', '', ...
    'irfSource', '', 'resultFolder', '', 'outputFile', '', 'figureFiles', '', ...
    'commonOutputFile', '', 'commonFigureFiles', '', 'dtNs', NaN, ...
    'pulsePeriodNs', NaN, 'nBins', NaN, 'nRejected', NaN, ...
    'totalCounts', NaN, 'peakCounts', NaN, 'tauNs', '', ...
    'amplitudes', '', 'relativeAmplitudes', '', 'nComponents', NaN, ...
    'nReportedComponents', NaN, 'nRemovedLowAmplitude', NaN, ...
    'tauErrorNs', '', 'amplitudeError', '', 'relativeAmplitudeError', '', ...
    'errorMethod', '', 'nErrorSamples', NaN, 'nErrorSamplesOk', NaN, ...
    'errorStatus', '', 'errorMessage', '', ...
    'chi', NaN, 'status', 'pending', 'message', '');
end

function [curves, summary] = readAllTcspcTargets(targets, curves, summary, opts)
for kk = 1:numel(targets)
    target = targets(kk);
    summary(kk).folder = target.folder;
    summary(kk).file = target.file;
    summary(kk).fileType = target.fileType;
    summary(kk).resultFolder = target.resultFolder;

    if ~strcmp(target.status, 'ok')
        curves(kk).status = target.status;
        curves(kk).message = target.message;
        summary(kk).status = target.status;
        summary(kk).message = target.message;
        continue;
    end

    fprintf('Reading %d/%d: %s\n', kk, numel(targets), target.file);
    try
        curves(kk) = readTcspcFile(target, opts, 'decay');
        summary(kk).dtNs = curves(kk).meta.tcspcResolutionSec * 1e9;
        summary(kk).pulsePeriodNs = curves(kk).meta.periodSec * 1e9;
        summary(kk).nBins = numel(curves(kk).counts);
        summary(kk).totalCounts = sum(curves(kk).counts);
        summary(kk).peakCounts = max(curves(kk).counts);
        summary(kk).status = 'read';
        summary(kk).message = '';
    catch ME
        curves(kk).status = 'failed';
        curves(kk).message = ME.message;
        curves(kk).folder = target.folder;
        curves(kk).folderName = target.folderName;
        curves(kk).file = target.file;
        curves(kk).fileType = target.fileType;
        curves(kk).resultFolder = target.resultFolder;
        summary(kk).status = 'failed';
        summary(kk).message = ME.message;
        warning('Read failed for %s: %s', target.file, ME.message);
    end
end
end

function [loaded, curves, summary, message] = loadTcspcCache(cacheFile, targets, opts)
loaded = false;
curves = repmat(emptyCurve(), numel(targets), 1);
summary = repmat(emptySummary(), numel(targets), 1);
message = '';

if opts.forceReadTcspc
    message = 'opts.forceReadTcspc is true.';
    return;
end
if isempty(cacheFile) || exist(cacheFile, 'file') ~= 2
    message = sprintf('cache file does not exist: %s', cacheFile);
    return;
end

try
    cacheData = load(cacheFile, 'tcspcCache');
catch ME
    message = sprintf('could not load cache: %s', ME.message);
    return;
end
if ~isfield(cacheData, 'tcspcCache')
    message = 'cache file does not contain tcspcCache.';
    return;
end

tcspcCache = cacheData.tcspcCache;
requiredFields = {'targetSignatures', 'readOptions', 'curves', 'summary'};
for ii = 1:numel(requiredFields)
    if ~isfield(tcspcCache, requiredFields{ii})
        message = sprintf('cache is missing field %s.', requiredFields{ii});
        return;
    end
end

currentSignatures = tcspcTargetSignatures(targets);
if ~sameTargetSignatures(currentSignatures, tcspcCache.targetSignatures)
    message = 'selected files changed since cache was written.';
    return;
end
if ~isequaln(tcspcReadOptionsForCache(opts), tcspcCache.readOptions)
    message = 'TCSPC read options changed since cache was written.';
    return;
end
if numel(tcspcCache.curves) ~= numel(targets) || numel(tcspcCache.summary) ~= numel(targets)
    message = 'cache target count does not match current target count.';
    return;
end

curves = tcspcCache.curves;
summary = tcspcCache.summary;
loaded = true;
message = sprintf('loaded %d target(s).', numel(targets));
end

function saveTcspcCache(cacheFile, rootFolder, outputFolder, targets, curves, summary, opts)
if isempty(cacheFile)
    return;
end

cacheFolder = fileparts(cacheFile);
if ~isempty(cacheFolder) && exist(cacheFolder, 'dir') ~= 7
    mkdir(cacheFolder);
end

tcspcCache = struct();
tcspcCache.version = 1;
tcspcCache.created = datestr(now, 30);
tcspcCache.rootFolder = rootFolder;
tcspcCache.outputFolder = outputFolder;
tcspcCache.targets = targets;
tcspcCache.targetSignatures = tcspcTargetSignatures(targets);
tcspcCache.readOptions = tcspcReadOptionsForCache(opts);
tcspcCache.curves = curves;
tcspcCache.summary = summary;

save(cacheFile, 'tcspcCache', '-v7.3');
end

function readOptions = tcspcReadOptionsForCache(opts)
readOptions = struct();
readOptions.ptuResolutionSec = opts.ptuResolutionSec;
readOptions.detectorChannels = opts.detectorChannels;
end

function signatures = tcspcTargetSignatures(targets)
signatures = repmat(struct('folder', '', 'resultFolder', '', 'file', '', ...
    'fileType', '', 'status', '', 'message', '', 'bytes', NaN, 'datenum', NaN), numel(targets), 1);

for ii = 1:numel(targets)
    signatures(ii).folder = targets(ii).folder;
    signatures(ii).resultFolder = targets(ii).resultFolder;
    signatures(ii).file = targets(ii).file;
    signatures(ii).fileType = targets(ii).fileType;
    signatures(ii).status = targets(ii).status;
    signatures(ii).message = targets(ii).message;
    if strcmp(targets(ii).status, 'ok') && exist(targets(ii).file, 'file') == 2
        info = dir(targets(ii).file);
        if ~isempty(info)
            signatures(ii).bytes = double(info(1).bytes);
            signatures(ii).datenum = double(info(1).datenum);
        end
    end
end
end

function tf = sameTargetSignatures(a, b)
tf = false;
fields = {'folder', 'resultFolder', 'file', 'fileType', 'status', 'message', 'bytes', 'datenum'};
if ~all(isfield(a, fields)) || ~all(isfield(b, fields))
    return;
end
if numel(a) ~= numel(b)
    return;
end

for ii = 1:numel(a)
    if ~strcmp(a(ii).folder, b(ii).folder) || ...
            ~strcmp(a(ii).resultFolder, b(ii).resultFolder) || ...
            ~strcmp(a(ii).file, b(ii).file) || ...
            ~strcmp(a(ii).fileType, b(ii).fileType) || ...
            ~strcmp(a(ii).status, b(ii).status) || ...
            ~strcmp(a(ii).message, b(ii).message)
        return;
    end
    if ~sameNumberOrNan(a(ii).bytes, b(ii).bytes)
        return;
    end
    if ~sameNumberOrNan(a(ii).datenum, b(ii).datenum)
        return;
    end
end

tf = true;
end

function tf = sameNumberOrNan(a, b)
if isnan(a) && isnan(b)
    tf = true;
else
    tf = abs(double(a) - double(b)) <= 1e-8 * max(1, max(abs(double(a)), abs(double(b))));
end
end

function targets = collectSubfolderTargets(rootFolder, opts)
targets = repmat(emptyTarget(), 0, 1);

if opts.includeRootFiles
    targets = [targets; targetFromFolder(rootFolder, opts)]; %#ok<AGROW>
end

entries = dir(rootFolder);
for kk = 1:numel(entries)
    if ~entries(kk).isdir
        continue;
    end
    name = entries(kk).name;
    if strcmp(name, '.') || strcmp(name, '..') || startsWith(name, '.')
        continue;
    end
    folder = fullfile(entries(kk).folder, name);
    targets = [targets; targetFromFolder(folder, opts)]; %#ok<AGROW>
end
end

function target = targetFromFolder(folder, opts)
target = emptyTarget();
target.folder = folder;
[~, target.folderName] = fileparts(folder);
target.resultFolder = fullfile(folder, opts.resultFolderName);

files = listTcspcFiles(folder, opts);
files = removeExplicitIrfFile(files, opts);
if isempty(files)
    target.status = 'skipped';
    target.message = 'No .ptu or .pqres file found.';
    return;
end
if strcmp(opts.irfMode, 'supplied') && numel(files) > 1
    isIrf = isIrfFileName({files.name});
    decayFiles = files(~isIrf);
    if ~isempty(decayFiles)
        files = decayFiles;
    end
end
if numel(files) > 1
    [files, selectionMessage] = resolveMultipleTcspcFiles(files, opts);
    target.message = selectionMessage;
end
if numel(files) > 1
    target.status = 'failed';
    target.message = sprintf('Expected one .ptu or .pqres file, found %d.', numel(files));
    return;
end

target.file = fullfile(files(1).folder, files(1).name);
[~, ~, ext] = fileparts(target.file);
target.fileType = lower(ext(2:end));
target.status = 'ok';
end

function [files, message] = resolveMultipleTcspcFiles(files, opts)
message = '';
if numel(files) <= 1
    return;
end

switch opts.multipleFileMode
    case 'error'
        return;
    case 'first'
        message = sprintf('Selected first TCSPC file from %d candidate files.', numel(files));
        files = files(1);
    case {'largest_ptu', 'largestptu', 'largest'}
        originalCount = numel(files);
        ptuMask = false(numel(files), 1);
        for ii = 1:numel(files)
            [~, ~, ext] = fileparts(files(ii).name);
            ptuMask(ii) = strcmpi(ext, '.ptu');
        end

        if any(ptuMask)
            candidates = files(ptuMask);
            candidateType = 'PTU';
        else
            candidates = files;
            candidateType = 'TCSPC';
        end

        [~, idx] = max([candidates.bytes]);
        files = candidates(idx);
        message = sprintf('Selected largest %s file from %d candidate files.', ...
            candidateType, originalCount);
    otherwise
        error('Unsupported opts.multipleFileMode: %s', opts.multipleFileMode);
end
end

function files = listTcspcFiles(folder, opts)
if opts.searchRecursively
    files = listTcspcFilesRecursive(folder, opts, 0);
else
    files = listTcspcFilesFlat(folder);
end
files = removeDuplicateDirEntries(files);
end

function files = listTcspcFilesFlat(folder)
files = [dir(fullfile(folder, '*.ptu')); dir(fullfile(folder, '*.PTU')); ...
    dir(fullfile(folder, '*.pqres')); dir(fullfile(folder, '*.PQRES'))];
end

function files = listTcspcFilesRecursive(folder, opts, depth)
files = listTcspcFilesFlat(folder);
if depth >= opts.maxSearchDepth
    return;
end

entries = dir(folder);
for ii = 1:numel(entries)
    if ~entries(ii).isdir
        continue;
    end
    name = entries(ii).name;
    if strcmp(name, '.') || strcmp(name, '..') || startsWith(name, '.')
        continue;
    end
    if strcmpi(name, opts.resultFolderName) || strcmpi(name, opts.batchResultFolderName)
        continue;
    end
    subFolder = fullfile(entries(ii).folder, name);
    files = [files; listTcspcFilesRecursive(subFolder, opts, depth + 1)]; %#ok<AGROW>
end
end

function files = removeExplicitIrfFile(files, opts)
if isempty(files) || ~isfield(opts, 'irfFile') || isempty(opts.irfFile)
    return;
end
irfPath = lower(char(opts.irfFile));
keep = true(numel(files), 1);
for ii = 1:numel(files)
    keep(ii) = ~strcmp(lower(fullfile(files(ii).folder, files(ii).name)), irfPath);
end
files = files(keep);
end

function tf = isIrfFileName(names)
names = lowerCell(names);
tf = containsAny(names, {'irf', 'instrument', 'response', 'prompt'});
tf = tf(:);
end

function curve = readTcspcFile(target, opts, preferredCurve)
if nargin < 3 || isempty(preferredCurve)
    preferredCurve = 'decay';
end
switch target.fileType
    case 'ptu'
        curve = readPtuTcspcFile(target.file, opts);
    case 'pqres'
        curve = readPqresTcspcFile(target.file, preferredCurve);
    otherwise
        error('Unsupported TCSPC file type: %s', target.fileType);
end

curve.folder = target.folder;
curve.folderName = target.folderName;
curve.file = target.file;
curve.fileType = target.fileType;
curve.resultFolder = target.resultFolder;
curve.status = 'ok';
curve.message = '';

curve.timeSec = double(curve.timeSec(:));
curve.counts = double(curve.counts(:));
curve.counts(~isfinite(curve.counts)) = 0;
curve.counts = max(curve.counts, 0);
if numel(curve.timeSec) ~= numel(curve.counts)
    error('Time/count vector length mismatch in %s.', target.file);
end
if ~any(curve.counts > 0)
    error('TCSPC curve has no positive counts in %s.', target.file);
end
end

function curve = readPtuTcspcFile(filePath, opts)
if exist('PTU_Read_Head', 'file') ~= 2 || exist('PTU_Read', 'file') ~= 2
    error('PTU_Read_Head.m and PTU_Read.m must be on the MATLAB path.');
end

try
    head = PTU_Read_Head(filePath, false);
catch
    head = PTU_Read_Head(filePath);
end

rawResolutionSec = getPtuResolutionSec(head);
syncRate = getPtuSyncRate(head);
pulsePeriodSec = 1 / syncRate;

if isempty(opts.ptuResolutionSec)
    chDiv = 1;
else
    chDiv = max(1, ceil(double(opts.ptuResolutionSec) / rawResolutionSec));
end
dtSec = rawResolutionSec * chDiv;
nBins = max(4, ceil(pulsePeriodSec / dtSec));

counts = zeros(nBins, 1);
chunkRecords = max(1, round(opts.ptuChunkRecords));
nRecords = getPtuRecordCount(head);
cnt = 0;

while true
    [~, tcspc, chan, special, num] = PTU_Read(filePath, [cnt + 1, chunkRecords], head);
    if isempty(num) || num <= 0
        break;
    end
    cnt = cnt + double(num);

    keep = ~logical(special);
    if ~isempty(opts.detectorChannels)
        keep = keep & ismember(double(chan), double(opts.detectorChannels(:).'));
    end

    bins = round(double(tcspc(keep)) ./ chDiv);
    valid = isfinite(bins) & bins >= 0 & bins < nBins;
    if any(valid)
        binIdx = bins(valid) + 1;
        counts = counts + accumarray(binIdx(:), ones(numel(binIdx), 1), [nBins, 1], @sum, 0);
    end

    if isfinite(nRecords) && cnt >= nRecords
        break;
    end
    if num < chunkRecords
        break;
    end
end

curve = emptyCurve();
curve.curveName = 'PTU summed TCSPC';
curve.timeSec = (0:nBins-1).' * dtSec;
curve.counts = counts;
curve.meta = struct('periodSec', pulsePeriodSec, 'deltaPulseSec', NaN, ...
    'tcspcResolutionSec', dtSec, 'rawResolutionSec', rawResolutionSec, ...
    'syncRate', syncRate, 'recordsRead', cnt);
curve.head = head;
curve.rawTcspc = counts;
curve.nChannels = NaN;
end

function curve = readPqresTcspcFile(filePath, preferredCurve)
if nargin < 2 || isempty(preferredCurve)
    preferredCurve = 'decay';
end
data = readPqresFile(filePath, preferredCurve);
curve = emptyCurve();
curve.curveName = data.curveName;
curve.timeSec = data.timeSec;
curve.counts = data.counts;
curve.meta = data.meta;
curve.head = [];
curve.rawTcspc = data.counts;
curve.nChannels = 1;
end

function batchIrf = loadBatchIrf(curves, opts)
switch opts.irfMode
    case 'supplied'
        batchIrf = readSuppliedIrf(opts);
    case 'global'
        batchIrf = calculateGlobalIrf(curves, opts);
    case 'best_per_curve'
        batchIrf = calculateBestPerCurveIrf(curves, opts);
    otherwise
        batchIrf = struct('status', 'failed', 'message', 'opts.irfMode is per_curve.', ...
            'method', 'per_curve', 'timeSec', [], 'counts', [], ...
            'summedCounts', [], 'included', [], 'excluded', [], ...
            'head', [], 'params', []);
end
end

function suppliedIrf = readSuppliedIrf(opts)
suppliedIrf = struct('status', 'failed', 'message', '', 'method', 'supplied', ...
    'timeSec', [], 'counts', [], 'summedCounts', [], 'included', [], ...
    'excluded', [], 'head', [], 'params', []);

irfFile = char(opts.irfFile);
if exist(irfFile, 'file') ~= 2
    suppliedIrf.message = sprintf('Supplied IRF file does not exist: %s', irfFile);
    return;
end

[irfFolder, irfStem, ext] = fileparts(irfFile);
target = emptyTarget();
target.folder = irfFolder;
target.folderName = irfStem;
target.file = irfFile;
target.fileType = lower(regexprep(ext, '^\.', ''));
target.resultFolder = irfFolder;
target.status = 'ok';

try
    curve = readTcspcFile(target, opts, 'irf');
catch ME
    suppliedIrf.message = ME.message;
    return;
end

irf = sanitizeIrf(curve.counts);
if ~any(irf > 0)
    suppliedIrf.message = 'Supplied IRF contains no positive counts.';
    return;
end

suppliedIrf.status = 'ok';
suppliedIrf.message = sprintf('Supplied IRF read from %s.', irfFile);
suppliedIrf.method = ['supplied_' curve.fileType];
suppliedIrf.file = irfFile;
suppliedIrf.curveName = curve.curveName;
suppliedIrf.timeSec = curve.timeSec(:);
suppliedIrf.counts = irf(:);
suppliedIrf.summedCounts = curve.counts(:);
suppliedIrf.head = curve.head;
suppliedIrf.params = struct('source', 'supplied', 'file', irfFile, ...
    'fileType', curve.fileType, 'curveName', curve.curveName, ...
    'meta', curve.meta);
end

function bestIrf = calculateBestPerCurveIrf(curves, opts)
method = ['best_per_curve_' opts.perCurveIrfModel];
bestIrf = struct('status', 'failed', 'message', '', 'method', method, ...
    'timeSec', [], 'counts', [], 'summedCounts', [], 'included', [], ...
    'excluded', [], 'head', [], 'params', []);

good = find(strcmp({curves.status}, 'ok'));
if isempty(good)
    bestIrf.message = 'No successfully read TCSPC curves.';
    return;
end

screenOpts = opts;
screenOpts.tau0 = opts.bestIrfTau0;
screenOpts.limits = opts.bestIrfLimits;
screenOpts.init = opts.bestIrfInit;
screenOpts.fluofitSolver = opts.bestIrfSolver;
screenOpts.plotFits = false;
screenOpts.saveFitFigures = false;

candidates = repmat(emptyBestIrfCandidate(), numel(good), 1);
bestLocal = 0;
bestChi = Inf;
selected = struct('timeSec', [], 'irf', [], 'rawCounts', [], ...
    'head', [], 'irfParams', []);

for jj = 1:numel(good)
    idx = good(jj);
    curve = curves(idx);
    candidate = emptyBestIrfCandidate();
    candidate.curveIndex = idx;
    candidate.file = curve.file;
    candidate.folder = curve.folder;
    candidate.curveName = curve.curveName;

    try
        [decayTimeSec, decayCounts, dtSec, pulsePeriodSec] = prepareDecayCurveForIrfSelection(curve);
        [irfOnAxis, irfParams] = estimatePerCurveIrfOnDecayAxis( ...
            decayCounts, dtSec, pulsePeriodSec, screenOpts);
        irfOnAxis = sanitizeIrf(irfOnAxis);
        if ~any(irfOnAxis > 0)
            error('Candidate IRF estimate is all zeros.');
        end

        [~, gatedCounts, gatedIrf, gateInfo] = ...
            applyTimeGate(decayTimeSec, decayCounts, irfOnAxis, screenOpts);
        if numel(gatedCounts) < 4
            error('Time gating left fewer than four samples.');
        end
        if ~any(gatedIrf > 0)
            error('Time gating removed all positive IRF samples.');
        end

        [cshift, ~, amplitudes, tauNs, ~, ~, ~, ~, ~, chi] = ...
            Fluofit(gatedIrf, gatedCounts, pulsePeriodSec * 1e9, ...
            dtSec * 1e9, screenOpts.tau0, screenOpts.limits, ...
            screenOpts.init, screenOpts.fluofitSolver, false);
        chi = realFiniteScalarOrNan(chi);
        if ~isfinite(chi)
            error('Bi-exponential screening fit returned non-finite or complex chi2.');
        end

        candidate.status = 'ok';
        candidate.message = '';
        candidate.chi = chi;
        candidate.colorShiftNs = realFiniteScalarOrNan(cshift);
        candidate.tauNs = real(double(tauNs(:))).';
        candidate.amplitudes = real(double(amplitudes(:))).';
        candidate.nBins = numel(gatedCounts);
        candidate.nRejected = gateInfo.nRejected;
        candidate.totalCounts = sum(gatedCounts);
        candidate.peakCounts = max(gatedCounts);
        candidate.irfSource = perCurveIrfSource(irfParams, opts);
        candidate.irfParams = irfParams;
        candidate.gateInfo = gateInfo;

        if chi < bestChi
            bestChi = chi;
            bestLocal = jj;
            selected.timeSec = decayTimeSec(:);
            selected.irf = irfOnAxis(:);
            selected.rawCounts = decayCounts(:);
            selected.head = struct('Resolution', dtSec * 1e9, ...
                'SyncRate', 1 / pulsePeriodSec);
            selected.irfParams = irfParams;
        end
    catch ME
        candidate.status = 'failed';
        candidate.message = ME.message;
    end

    candidates(jj) = candidate;
end

if bestLocal == 0
    bestIrf.message = 'No candidate IRF passed the bi-exponential screening fit.';
    bestIrf.excluded = candidates;
    return;
end

selectedCandidate = candidates(bestLocal);
okMask = strcmp({candidates.status}, 'ok');
bestIrf.status = 'ok';
bestIrf.message = sprintf(['Selected IRF from %s with lowest ', ...
    'bi-exponential screening chi2 = %.6g.'], selectedCandidate.file, selectedCandidate.chi);
bestIrf.timeSec = selected.timeSec(:);
bestIrf.counts = sanitizeIrf(selected.irf);
bestIrf.summedCounts = selected.rawCounts(:);
bestIrf.included = selectedCandidate.curveIndex;
bestIrf.excluded = candidates(~okMask);
bestIrf.head = selected.head;

params = struct();
params.source = 'best_per_curve';
params.model = opts.perCurveIrfModel;
params.selectionCriterion = 'minimum chi2 from bi-exponential Fluofit screening';
params.screenTau0Ns = screenOpts.tau0;
params.screenLimitsNs = screenOpts.limits;
params.screenInit = screenOpts.init;
params.screenSolver = screenOpts.fluofitSolver;
params.selectedCurveIndex = selectedCandidate.curveIndex;
params.selectedFile = selectedCandidate.file;
params.selectedChi = selectedCandidate.chi;
params.selectedCandidate = selectedCandidate;
params.selectedIrfParams = selected.irfParams;
params.candidates = candidates;
bestIrf.params = params;
end

function candidate = emptyBestIrfCandidate()
candidate = struct('curveIndex', NaN, 'file', '', 'folder', '', ...
    'curveName', '', 'irfSource', '', 'chi', NaN, 'colorShiftNs', NaN, ...
    'tauNs', [], 'amplitudes', [], 'nBins', NaN, 'nRejected', NaN, ...
    'totalCounts', NaN, 'peakCounts', NaN, 'status', 'pending', ...
    'message', '', 'irfParams', [], 'gateInfo', []);
end

function [decayTimeSec, decayCounts, dtSec, pulsePeriodSec] = prepareDecayCurveForIrfSelection(curve)
decayTimeSec = double(curve.timeSec(:));
decayCounts = max(double(curve.counts(:)), 0);

[decayTimeSec, order] = sort(decayTimeSec);
decayCounts = decayCounts(order);
[decayTimeSec, uniqueIdx] = unique(decayTimeSec, 'stable');
decayCounts = decayCounts(uniqueIdx);

dtSec = curve.meta.tcspcResolutionSec;
if ~isfinite(dtSec) || dtSec <= 0
    dtSec = median(diff(decayTimeSec));
end
if ~isfinite(dtSec) || dtSec <= 0
    error('Could not determine TCSPC bin width for %s.', curve.file);
end

pulsePeriodSec = curve.meta.periodSec;
if ~isfinite(pulsePeriodSec) || pulsePeriodSec <= 0
    pulsePeriodSec = (max(decayTimeSec) - min(decayTimeSec)) + dtSec;
end
if ~isfinite(pulsePeriodSec) || pulsePeriodSec <= 0
    error('Could not determine pulse period for %s.', curve.file);
end
end

function src = perCurveIrfSource(params, opts)
if isstruct(params) && isfield(params, 'modelKey') && ~isempty(params.modelKey)
    src = ['per_curve_' char(params.modelKey)];
else
    src = ['per_curve_' opts.perCurveIrfModel];
end
end

function val = realFiniteScalarOrNan(x)
val = NaN;
if isempty(x)
    return;
end
x = x(1);
if ~isfinite(real(x)) || ~isfinite(imag(x))
    return;
end
if abs(imag(x)) > 1e-8 * max(1, abs(real(x)))
    return;
end
val = real(double(x));
end

function globalIrf = calculateGlobalIrf(curves, opts)
globalIrf = struct('status', 'failed', 'message', '', 'method', opts.globalIrfMethod, ...
    'timeSec', [], 'counts', [], 'summedCounts', [], 'included', [], ...
    'excluded', [], 'head', [], 'params', []);

if ~strcmp(opts.irfMode, 'global')
    globalIrf.message = 'opts.irfMode is per_curve.';
    return;
end

good = find(strcmp({curves.status}, 'ok'));
if isempty(good)
    globalIrf.message = 'No successfully read TCSPC curves.';
    return;
end

totals = zeros(numel(good), 1);
for ii = 1:numel(good)
    totals(ii) = sum(curves(good(ii)).counts);
end
[~, refLocal] = max(totals);
refIdx = good(refLocal);
ref = curves(refIdx);
refTime = ref.timeSec(:) - ref.timeSec(1);
refDt = ref.meta.tcspcResolutionSec;
refPeriod = ref.meta.periodSec;

if ~isfinite(refDt) || refDt <= 0 || ~isfinite(refPeriod) || refPeriod <= 0
    globalIrf.message = 'Reference curve has invalid TCSPC resolution or pulse period.';
    return;
end

summed = zeros(size(refTime));
included = [];
excluded = struct('file', {}, 'reason', {});

for ii = 1:numel(good)
    idx = good(ii);
    curve = curves(idx);
    curvePeriod = curve.meta.periodSec;
    curveDt = curve.meta.tcspcResolutionSec;
    if ~isfinite(curvePeriod) || curvePeriod <= 0 || ~isfinite(curveDt) || curveDt <= 0
        excluded(end+1) = struct('file', curve.file, 'reason', 'Invalid dt or period.'); %#ok<AGROW>
        continue;
    end
    if abs(curvePeriod - refPeriod) > opts.globalIrfPeriodToleranceSec
        excluded(end+1) = struct('file', curve.file, 'reason', 'Pulse period differs from reference.'); %#ok<AGROW>
        continue;
    end

    curveTime = curve.timeSec(:) - curve.timeSec(1);
    if numel(curveTime) == numel(refTime) && ...
            max(abs(curveTime - refTime)) <= 0.5 * max(refDt, curveDt)
        y = curve.counts(:);
    elseif opts.globalIrfResample
        [curveTime, uniqueIdx] = unique(curveTime, 'stable');
        y0 = curve.counts(uniqueIdx);
        y = interp1(curveTime, y0, refTime, 'linear', 0);
    else
        excluded(end+1) = struct('file', curve.file, 'reason', 'Time axis differs and resampling is disabled.'); %#ok<AGROW>
        continue;
    end

    y(~isfinite(y)) = 0;
    y = max(double(y(:)), 0);
    summed = summed + y;
    included(end+1) = idx; %#ok<AGROW>
end

if isempty(included) || ~any(summed > 0)
    globalIrf.message = 'No compatible positive curves for global IRF.';
    globalIrf.excluded = excluded;
    return;
end

head = struct();
head.Resolution = refDt * 1e9;
head.SyncRate = 1 / refPeriod;

try
    [irf, params] = estimateGlobalIrfFromSum(head, summed, opts);
catch ME
    globalIrf.message = ME.message;
    globalIrf.included = included;
    globalIrf.excluded = excluded;
    globalIrf.head = head;
    globalIrf.summedCounts = summed;
    return;
end

irf = sanitizeIrf(irf);
if ~any(irf > 0)
    globalIrf.message = 'Global IRF estimate is all zeros.';
    globalIrf.included = included;
    globalIrf.excluded = excluded;
    globalIrf.head = head;
    globalIrf.summedCounts = summed;
    return;
end

globalIrf.status = 'ok';
globalIrf.message = sprintf('Global IRF estimated from %d curve(s).', numel(included));
globalIrf.timeSec = ref.timeSec(:);
globalIrf.counts = irf(:);
globalIrf.summedCounts = summed(:);
globalIrf.included = included;
globalIrf.excluded = excluded;
globalIrf.head = head;
globalIrf.params = params;
end

function [irf, params] = estimateGlobalIrfFromSum(head, summed, opts)
method = lower(opts.globalIrfMethod);
switch method
    case {'calc_mirf', 'calcmirf'}
        irf = Calc_mIRF(head, double(summed(:).'));
        [irf, clipInfo] = clipIrfBelowFraction(irf, opts.irfClipFraction);
        params = struct('source', 'Calc_mIRF', 'model', 'IRF_Fun', ...
            'clipInfo', clipInfo, ...
            'note', 'Estimated by Calc_mIRF from summed collected TCSPC curves.');
    case {'gamma_shifted_fast', 'gamma_fast'}
        if exist('Calc_mIRF_Global_GammaShifted_fast', 'file') ~= 2
            error('Could not find Calc_mIRF_Global_GammaShifted_fast.m.');
        end
        out = Calc_mIRF_Global_GammaShifted_fast(head, double(summed(:)), opts.tau0, opts.globalIrfOptions);
        irf = out.IRF;
        params = out;
    case {'exgauss_fast', 'ex_gauss_fast', 'spad_exgauss', 'spad_exgaussian'}
        if exist('Calc_mIRF_Global_ExGauss_fast', 'file') ~= 2
            error('Could not find Calc_mIRF_Global_ExGauss_fast.m.');
        end
        out = Calc_mIRF_Global_ExGauss_fast(head, double(summed(:)), opts.tau0, opts.globalIrfOptions);
        irf = out.IRF;
        params = out;
    otherwise
        error('Unsupported global IRF method: %s', opts.globalIrfMethod);
end
end

function fitInput = buildBatchFitInput(curve, globalIrf, opts)
decayTimeSec = double(curve.timeSec(:));
decayCounts = max(double(curve.counts(:)), 0);

[decayTimeSec, order] = sort(decayTimeSec);
decayCounts = decayCounts(order);
[decayTimeSec, uniqueIdx] = unique(decayTimeSec, 'stable');
decayCounts = decayCounts(uniqueIdx);

dtSec = curve.meta.tcspcResolutionSec;
if ~isfinite(dtSec) || dtSec <= 0
    dtSec = median(diff(decayTimeSec));
end
if ~isfinite(dtSec) || dtSec <= 0
    error('Could not determine TCSPC bin width for %s.', curve.file);
end

pulsePeriodSec = curve.meta.periodSec;
if ~isfinite(pulsePeriodSec) || pulsePeriodSec <= 0
    pulsePeriodSec = (max(decayTimeSec) - min(decayTimeSec)) + dtSec;
end

irfSource = '';
irfInfo = struct();
if any(strcmp(opts.irfMode, {'global', 'supplied', 'best_per_curve'})) && strcmp(globalIrf.status, 'ok')
    irfOnAxis = mapIrfToDecayAxis(globalIrf.timeSec, globalIrf.counts, decayTimeSec, dtSec);
    if strcmp(opts.irfMode, 'global')
        irfSource = ['global_' globalIrf.method];
    elseif strcmp(opts.irfMode, 'best_per_curve')
        irfSource = globalIrf.method;
    else
        irfSource = globalIrf.method;
    end
    irfInfo = globalIrf;
else
    irfOnAxis = zeros(size(decayCounts));
end

if ~any(irfOnAxis > 0)
    if strcmp(opts.irfMode, 'best_per_curve')
        error('No usable selected best IRF for %s.', curve.file);
    end
    if ~opts.fallbackPerCurveIrf
        error('No usable batch IRF for %s, and fallbackPerCurveIrf is false.', curve.file);
    end
    [irfOnAxis, params] = estimatePerCurveIrfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec, opts);
    irfSource = ['per_curve_' params.modelKey];
    irfInfo = struct('params', params);
end

irfOnAxis = sanitizeIrf(irfOnAxis);
if ~any(irfOnAxis > 0)
    error('IRF is all zeros for %s.', curve.file);
end

[decayTimeSec, decayCounts, irfOnAxis, gateInfo] = ...
    applyTimeGate(decayTimeSec, decayCounts, irfOnAxis, opts);
if numel(decayCounts) < 4
    error('Time gating left fewer than four samples for %s.', curve.file);
end
if ~any(irfOnAxis > 0)
    error('Time gating removed all positive IRF samples for %s.', curve.file);
end

fitInput = struct();
fitInput.file = curve.file;
fitInput.folder = curve.folder;
fitInput.fileType = curve.fileType;
fitInput.curveName = curve.curveName;
fitInput.irfSource = irfSource;
fitInput.irfInfo = irfInfo;
fitInput.gateInfo = gateInfo;
fitInput.timeNs = (decayTimeSec - decayTimeSec(1)) * 1e9;
fitInput.counts = decayCounts(:);
fitInput.irf = irfOnAxis(:);
fitInput.dtNs = dtSec * 1e9;
fitInput.pulsePeriodNs = pulsePeriodSec * 1e9;
fitInput.tau0Ns = double(opts.tau0(:)).';
fitInput.limitsNs = opts.limits;
fitInput.init = opts.init;
fitInput.fluofitSolver = opts.fluofitSolver;
end

function [irfOnAxis, params] = estimatePerCurveIrfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec, opts)
switch opts.perCurveIrfModel
    case {'spad_exgauss', 'spad_exgaussian', 'exgauss', 'ex_gauss'}
        [irfOnAxis, params] = spadExGaussIrfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec, opts);
    case {'calc_mirf', 'calcmirf'}
        [irfOnAxis, params] = calcMirfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec, opts);
    otherwise
        error('Unsupported per-curve IRF model: %s', opts.perCurveIrfModel);
end
end

function [irfOnAxis, params] = spadExGaussIrfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec, opts)
if exist('Calc_mIRF_Global_ExGauss_fast', 'file') ~= 2
    error('Could not find Calc_mIRF_Global_ExGauss_fast.m.');
end

head = struct();
head.Resolution = dtSec * 1e9;
head.SyncRate = 1 / pulsePeriodSec;

fitOpts = opts.spadIrfOptions;
out = Calc_mIRF_Global_ExGauss_fast(head, double(decayCounts(:)), opts.tau0, fitOpts);
[irf, clipInfo] = clipIrfBelowFraction(out.IRF, opts.irfClipFraction);
irfOnAxis = sanitizeIrf(irf);

params = out;
params.source = 'Calc_mIRF_Global_ExGauss_fast';
params.model = 'SPAD ex-Gaussian';
params.modelKey = 'spad_exgauss';
params.head = head;
params.clipInfo = clipInfo;
params.note = ['SPAD IRF model: Gaussian avalanche timing peak convolved ', ...
    'with a one-sided exponential diffusion/timing tail.'];
end

function [irfOnAxis, params] = calcMirfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec, opts)
head = struct();
head.Resolution = dtSec * 1e9;
head.SyncRate = 1 / pulsePeriodSec;
irf = Calc_mIRF(head, double(decayCounts(:).'));
[irf, clipInfo] = clipIrfBelowFraction(irf, opts.irfClipFraction);
irfOnAxis = sanitizeIrf(irf);
params = struct('source', 'Calc_mIRF', 'model', 'IRF_Fun', ...
    'modelKey', 'calc_mirf', 'head', head, ...
    'clipInfo', clipInfo, ...
    'note', 'Estimated per curve by Calc_mIRF(head, tcspc).');
end

function [irf, clipInfo] = clipIrfBelowFraction(irf, fraction)
if nargin < 2 || isempty(fraction)
    fraction = 0;
end
fraction = double(fraction);
if ~isscalar(fraction) || ~isfinite(fraction) || fraction < 0
    error('opts.irfClipFraction must be a finite non-negative scalar.');
end

irf = real(double(squeeze(irf)));
irf = irf(:);
irf(~isfinite(irf)) = 0;
irf = max(irf, 0);

maxVal = max(irf);
threshold = 0;
nClipped = 0;
if fraction > 0 && maxVal > 0
    threshold = fraction * maxVal;
    clipMask = irf < threshold;
    nClipped = sum(clipMask);
    irf(clipMask) = 0;
end

clipInfo = struct('fraction', fraction, 'threshold', threshold, ...
    'maxValue', maxVal, 'nClipped', nClipped);
end

function fitInput = buildBatchTailFitInput(curve, opts)
[decayTimeSec, decayCounts, dtSec, pulsePeriodSec] = prepareDecayCurveForIrfSelection(curve);
dtNs = dtSec * 1e9;
pulsePeriodNs = pulsePeriodSec * 1e9;

[~, peakIdx] = max(decayCounts);
cutBins = round(opts.tailCutAfterPeakNs / dtNs);
tailStartIdx = peakIdx + cutBins;
tailStartIdx = min(max(tailStartIdx, 1), numel(decayCounts));
tailEndIdx = numel(decayCounts) - opts.tailRejectLastNPoints;
if tailEndIdx < tailStartIdx
    error('Tail gate removes all samples after rejecting the final %d point(s) for %s.', ...
        opts.tailRejectLastNPoints, curve.file);
end
if tailEndIdx - tailStartIdx + 1 < 3
    error('Tail gate leaves fewer than three samples for %s.', curve.file);
end

tailIdx = tailStartIdx:tailEndIdx;
tailTimeSec = decayTimeSec(tailIdx);
tailCounts = decayCounts(tailIdx);
if ~any(tailCounts > 0)
    error('Tail gate contains no positive counts for %s.', curve.file);
end

fullTimeNs = (decayTimeSec - decayTimeSec(1)) * 1e9;
tailAbsoluteTimeNs = (tailTimeSec - decayTimeSec(1)) * 1e9;
tailTimeNs = tailAbsoluteTimeNs - tailAbsoluteTimeNs(1);

gateInfo = struct();
gateInfo.method = 'peak_plus_cut';
gateInfo.peakIndex = peakIdx;
gateInfo.peakTimeNs = fullTimeNs(peakIdx);
gateInfo.tailStartIndex = tailStartIdx;
gateInfo.tailStartTimeNs = tailAbsoluteTimeNs(1);
gateInfo.tailEndIndex = tailEndIdx;
gateInfo.tailEndTimeNs = tailAbsoluteTimeNs(end);
gateInfo.tailCutAfterPeakNs = opts.tailCutAfterPeakNs;
gateInfo.tailRejectLastNPoints = opts.tailRejectLastNPoints;
tailRejectedIdx = (tailEndIdx + 1):numel(decayCounts);
earlyRejectedTimeNs = fullTimeNs(1:tailStartIdx-1);
lateRejectedTimeNs = fullTimeNs(tailRejectedIdx);
gateInfo.nRejected = tailStartIdx - 1 + numel(tailRejectedIdx);
gateInfo.rejectedTimeNs = [earlyRejectedTimeNs(:); lateRejectedTimeNs(:)].';
gateInfo.rejectedLastTimeNs = lateRejectedTimeNs(:).';

fitInput = struct();
fitInput.file = curve.file;
fitInput.folder = curve.folder;
fitInput.fileType = curve.fileType;
fitInput.curveName = curve.curveName;
fitInput.fitMethod = 'tail';
fitInput.irfSource = 'none_tail';
fitInput.irfInfo = struct('status', 'skipped', 'message', 'Tail fitting does not use an IRF.');
fitInput.gateInfo = gateInfo;
fitInput.timeNs = tailTimeNs(:);
fitInput.tailAbsoluteTimeNs = tailAbsoluteTimeNs(:);
fitInput.fullTimeNs = fullTimeNs(:);
fitInput.fullCounts = decayCounts(:);
fitInput.counts = tailCounts(:);
fitInput.irf = [];
fitInput.dtNs = dtNs;
fitInput.pulsePeriodNs = pulsePeriodNs;
fitInput.tau0Ns = double(opts.tailTau0(:)).';
fitInput.limitsNs = opts.tailLimits;
fitInput.init = opts.init;
fitInput.fluofitSolver = opts.tailSolver;
end

function [fitResult, fitFigures] = runBatchTailFit(fitInput, opts)
makePlot = opts.plotFits || opts.saveFitFigures;
if makePlot
    figuresBeforeFit = findall(0, 'Type', 'figure');
else
    figuresBeforeFit = [];
end

y = fitInput.counts(:);
dtNs = fitInput.dtNs;
if opts.runDistTailfit
    tailTauMax = opts.tailDistTauMaxNs;
    if isempty(tailTauMax)
        tailTauMax = [];
    end

    if opts.tailUsePeriodCorrection && isfinite(fitInput.pulsePeriodNs) && fitInput.pulsePeriodNs > 0
        [distAmp, distRates, distOffset, distFit, distTimeIndex, distChi] = ...
            DistTailfit(y, dtNs, opts.tailDistFlag, 0, opts.tailDistN, tailTauMax, fitInput.pulsePeriodNs);
    else
        [distAmp, distRates, distOffset, distFit, distTimeIndex, distChi] = ...
            DistTailfit(y, dtNs, opts.tailDistFlag, 0, opts.tailDistN, tailTauMax);
    end
else
    distAmp = [];
    distRates = [];
    distOffset = NaN;
    distFit = NaN(size(y));
    distTimeIndex = [];
    distChi = NaN;
end

[tauNs, amplitudes, offset, z, tAxisNs, chi, coeff] = ...
    Tailfit(y, dtNs, fitInput.tau0Ns, fitInput.limitsNs, ...
    opts.tailSolver, 0, opts.tailSimplexSteps);

fitResult = struct();
if opts.runDistTailfit
    fitResult.signature = 'Tail mode: Tailfit(y,dt,tau0) plus DistTailfit(y,dt)';
else
    fitResult.signature = 'Tail mode: Tailfit(y,dt,tau0); DistTailfit skipped because init=0';
end
fitResult.fitMethod = 'tail';
fitResult.fileSuffix = 'tailfit';
fitResult.fitMode = opts.tailSolver;
fitResult.distTailfitEnabled = opts.runDistTailfit;
fitResult.colorShiftNs = 0;
fitResult.offset = offset;
fitResult.amplitudes = amplitudes(:);
fitResult.tauNs = tauNs(:);
fitResult.colorShiftErrorNs = [];
fitResult.tauErrorNs = [];
fitResult.irfShifted = [];
fitResult.components = tailComponents(tAxisNs, tauNs, coeff);
fitResult.timeAxisNs = tAxisNs(:);
fitResult.tailAbsoluteTimeNs = fitInput.tailAbsoluteTimeNs(:);
fitResult.tailFit = z(:);
fitResult.chi = chi;
fitResult.meanTauNs = weightedMeanTau(tauNs, amplitudes);
fitResult.distAmplitudes = distAmp(:);
fitResult.distRates = distRates(:);
fitResult.distTauNs = 1 ./ distRates(:);
fitResult.distOffset = distOffset;
fitResult.distFit = distFit(:);
fitResult.distTimeIndex = distTimeIndex(:);
fitResult.distChi = distChi;
fitResult.distMeanTauNs = weightedMeanTau(1 ./ distRates(:), distAmp(:));
fitResult.distRateMeanTauNs = distRateMeanTau(distRates(:), distAmp(:));
fitResult.tailGate = fitInput.gateInfo;
fitResult.figureFiles = {};

if makePlot
    plotTailFitFigure(fitInput, fitResult);
    figuresAfterFit = findall(0, 'Type', 'figure');
    fitFigures = newFigureHandles(figuresAfterFit, figuresBeforeFit);
else
    fitFigures = [];
end
end

function components = tailComponents(tNs, tauNs, coeff)
tNs = tNs(:);
tauNs = tauNs(:).';
coeff = coeff(:).';
components = zeros(numel(tNs), numel(coeff));
if isempty(coeff)
    return;
end
components(:, 1) = coeff(1);
for ii = 1:numel(tauNs)
    components(:, ii + 1) = coeff(ii + 1) .* exp(-tNs ./ tauNs(ii));
end
end

function val = weightedMeanTau(tauNs, amplitudes)
tauNs = real(double(tauNs(:)));
amplitudes = real(double(amplitudes(:)));
keep = isfinite(tauNs) & tauNs > 0 & isfinite(amplitudes) & amplitudes > 0;
if ~any(keep)
    val = NaN;
else
    val = sum(tauNs(keep).*amplitudes(keep)) ./ sum(amplitudes(keep));
end
end

function components = filterSummaryComponents(tauNs, amplitudes, relThreshold)
tauNs = real(double(tauNs(:)));
amplitudes = real(double(amplitudes(:)));
n = min(numel(tauNs), numel(amplitudes));
tauNs = tauNs(1:n);
amplitudes = amplitudes(1:n);

valid = isfinite(tauNs) & tauNs > 0 & isfinite(amplitudes);
positiveAmplitudes = amplitudes;
positiveAmplitudes(~isfinite(positiveAmplitudes) | positiveAmplitudes < 0) = 0;
totalPositiveAmplitude = sum(positiveAmplitudes(valid));

relativeAmplitudes = NaN(size(amplitudes));
if totalPositiveAmplitude > 0
    relativeAmplitudes(valid) = positiveAmplitudes(valid) ./ totalPositiveAmplitude;
end

if relThreshold > 0 && totalPositiveAmplitude > 0
    keep = valid & relativeAmplitudes >= relThreshold;
else
    keep = valid;
end

if ~any(keep) && any(valid)
    validIdx = find(valid);
    [~, bestIdx] = max(positiveAmplitudes(validIdx));
    keep(validIdx(bestIdx)) = true;
end

components = struct();
components.tauNs = tauNs(keep);
components.amplitudes = amplitudes(keep);
components.relativeAmplitudes = relativeAmplitudes(keep);
components.nComponents = sum(valid);
components.nReportedComponents = sum(keep);
components.nRemovedLowAmplitude = max(0, components.nComponents - components.nReportedComponents);
components.keepMask = keep(:);
end

function fitResult = addFitUncertainty(fitInput, fitResult, opts)
fitResult.resampleError = emptyResampleError(opts.errorMethod, ...
    'Resampling was not run.');

if strcmp(opts.errorMethod, 'none')
    fitResult.resampleError.message = 'Disabled by opts.errorMethod=''none''.';
    return;
end

if strcmp(opts.errorMethod, 'bootstrap')
    nSamples = opts.nErrorSamples;
    stderrDivisor = 1;
elseif strcmp(opts.errorMethod, 'chunk')
    nSamples = opts.errorChunks;
    stderrDivisor = sqrt(max(nSamples, 1));
else
    error('Unsupported uncertainty method: %s', opts.errorMethod);
end

if nSamples < 2
    fitResult.resampleError.status = 'skipped';
    fitResult.resampleError.message = 'At least two resamples are required.';
    return;
end

if ~isempty(opts.errorRandomSeed)
    rngSeed = mod(round(double(opts.errorRandomSeed)) + stringHash(fitInput.file), 2^31 - 1);
    rng(rngSeed, 'twister');
else
    rngSeed = NaN;
end

lambdaFull = uncertaintyLambda(fitInput, fitResult, opts);
if strcmp(opts.errorMethod, 'chunk')
    lambdaPerSample = lambdaFull ./ nSamples;
    sampleScale = nSamples;
else
    lambdaPerSample = lambdaFull;
    sampleScale = 1;
end

refTauNs = real(double(fitResult.tauNs(:)));
refAmplitudes = real(double(fitResult.amplitudes(:)));
nComp = min(numel(refTauNs), numel(refAmplitudes));
refTauNs = refTauNs(1:nComp);
refAmplitudes = refAmplitudes(1:nComp);

tauSamples = NaN(nSamples, nComp);
amplitudeSamples = NaN(nSamples, nComp);
relativeAmplitudeSamples = NaN(nSamples, nComp);
offsetSamples = NaN(nSamples, 1);
chiSamples = NaN(nSamples, 1);
messages = cell(nSamples, 1);
ok = false(nSamples, 1);

for ii = 1:nSamples
    try
        ySample = poissonRandomCounts(lambdaPerSample) .* sampleScale;
        sample = refitSampleCounts(fitInput, fitResult, opts, ySample);
        [tauAligned, ampAligned] = alignComponentsToReference( ...
            refTauNs, sample.tauNs, sample.amplitudes);
        tauSamples(ii, :) = tauAligned(:).';
        amplitudeSamples(ii, :) = ampAligned(:).';
        relativeAmplitudeSamples(ii, :) = relativeAmplitudesFromAmplitudes(ampAligned).';
        offsetSamples(ii) = sample.offset;
        chiSamples(ii) = sample.chi;
        ok(ii) = true;
    catch ME
        messages{ii} = ME.message;
    end
end

resampleError = emptyResampleError(opts.errorMethod, '');
resampleError.status = 'ok';
if ~any(ok)
    resampleError.status = 'failed';
    resampleError.message = 'All resampled fits failed.';
elseif sum(ok) < max(3, ceil(0.5*nSamples))
    resampleError.status = 'partial';
    resampleError.message = sprintf('%d of %d resampled fits succeeded.', sum(ok), nSamples);
else
    resampleError.message = sprintf('%d of %d resampled fits succeeded.', sum(ok), nSamples);
end

resampleError.method = opts.errorMethod;
resampleError.source = opts.errorResampleSource;
resampleError.randomSeed = rngSeed;
resampleError.nRequested = nSamples;
resampleError.nSucceeded = sum(ok);
resampleError.stderrDivisor = stderrDivisor;
resampleError.referenceTauNs = refTauNs(:);
resampleError.referenceAmplitudes = refAmplitudes(:);
resampleError.referenceRelativeAmplitudes = relativeAmplitudesFromAmplitudes(refAmplitudes);
resampleError.tauSamplesNs = tauSamples;
resampleError.amplitudeSamples = amplitudeSamples;
resampleError.relativeAmplitudeSamples = relativeAmplitudeSamples;
resampleError.offsetSamples = offsetSamples;
resampleError.chiSamples = chiSamples;
resampleError.failedMessages = messages(~ok);
resampleError.tauStdNs = nanStdColumns(tauSamples(ok, :)) ./ stderrDivisor;
resampleError.amplitudeStd = nanStdColumns(amplitudeSamples(ok, :)) ./ stderrDivisor;
resampleError.relativeAmplitudeStd = nanStdColumns(relativeAmplitudeSamples(ok, :)) ./ stderrDivisor;
resampleError.offsetStd = nanStdColumns(offsetSamples(ok)) ./ stderrDivisor;
resampleError.chiStd = nanStdColumns(chiSamples(ok)) ./ stderrDivisor;
resampleError.tauMedianNs = nanPercentileColumns(tauSamples(ok, :), 50);
resampleError.amplitudeMedian = nanPercentileColumns(amplitudeSamples(ok, :), 50);
resampleError.relativeAmplitudeMedian = nanPercentileColumns(relativeAmplitudeSamples(ok, :), 50);
resampleError.tauCI68Ns = [nanPercentileColumns(tauSamples(ok, :), 16); ...
    nanPercentileColumns(tauSamples(ok, :), 84)];
resampleError.amplitudeCI68 = [nanPercentileColumns(amplitudeSamples(ok, :), 16); ...
    nanPercentileColumns(amplitudeSamples(ok, :), 84)];
resampleError.relativeAmplitudeCI68 = [nanPercentileColumns(relativeAmplitudeSamples(ok, :), 16); ...
    nanPercentileColumns(relativeAmplitudeSamples(ok, :), 84)];

fitResult.resampleError = resampleError;
fitResult.tauResampleErrorNs = resampleError.tauStdNs(:);
fitResult.amplitudeResampleError = resampleError.amplitudeStd(:);
fitResult.relativeAmplitudeResampleError = resampleError.relativeAmplitudeStd(:);
end

function resampleError = emptyResampleError(method, message)
resampleError = struct();
resampleError.status = 'skipped';
resampleError.method = method;
resampleError.source = '';
resampleError.message = message;
resampleError.randomSeed = NaN;
resampleError.nRequested = 0;
resampleError.nSucceeded = 0;
resampleError.stderrDivisor = 1;
resampleError.referenceTauNs = [];
resampleError.referenceAmplitudes = [];
resampleError.referenceRelativeAmplitudes = [];
resampleError.tauSamplesNs = [];
resampleError.amplitudeSamples = [];
resampleError.relativeAmplitudeSamples = [];
resampleError.offsetSamples = [];
resampleError.chiSamples = [];
resampleError.failedMessages = {};
resampleError.tauStdNs = [];
resampleError.amplitudeStd = [];
resampleError.relativeAmplitudeStd = [];
resampleError.offsetStd = NaN;
resampleError.chiStd = NaN;
resampleError.tauMedianNs = [];
resampleError.amplitudeMedian = [];
resampleError.relativeAmplitudeMedian = [];
resampleError.tauCI68Ns = [];
resampleError.amplitudeCI68 = [];
resampleError.relativeAmplitudeCI68 = [];
end

function lambda = uncertaintyLambda(fitInput, fitResult, opts)
if strcmp(opts.errorResampleSource, 'counts')
    lambda = fitInput.counts(:);
else
    lambda = [];
    if strcmp(fitResult.fitMethod, 'tail') && isfield(fitResult, 'tailFit')
        lambda = fitResult.tailFit(:);
    elseif isfield(fitResult, 'components') && ~isempty(fitResult.components)
        lambda = sum(real(double(fitResult.components)), 2);
    end
    if numel(lambda) ~= numel(fitInput.counts) || ~any(isfinite(lambda) & lambda > 0)
        lambda = fitInput.counts(:);
    end
end
lambda = real(double(lambda(:)));
lambda(~isfinite(lambda)) = 0;
lambda = max(lambda, 0);
end

function sample = refitSampleCounts(fitInput, fitResult, opts, ySample)
ySample = max(real(double(ySample(:))), 0);
switch fitResult.fitMethod
    case 'tail'
        tauStart = real(double(fitResult.tauNs(:))).';
        limits = resampleLimitsForTau(fitInput.limitsNs, tauStart);
        [tauNs, amplitudes, offset, ~, ~, chi] = Tailfit(ySample, fitInput.dtNs, ...
            tauStart, limits, opts.tailSolver, 0, opts.tailSimplexSteps);
    case 'fluofit'
        tauStart = real(double(fitResult.tauNs(:))).';
        limits = resampleLimitsForTau(fitInput.limitsNs, tauStart);
        [~, offset, amplitudes, tauNs, ~, ~, ~, ~, ~, chi] = ...
            Fluofit(fitInput.irf, ySample, fitInput.pulsePeriodNs, ...
            fitInput.dtNs, tauStart, limits, 0, ...
            opts.fluofitSolver, false);
    otherwise
        error('Unsupported fit method for uncertainty refit: %s', fitResult.fitMethod);
end
sample = struct();
sample.tauNs = real(double(tauNs(:)));
sample.amplitudes = real(double(amplitudes(:)));
sample.offset = real(double(offset));
sample.chi = real(double(chi));
end

function limits = resampleLimitsForTau(limitsIn, tauStart)
limits = [];
if isempty(limitsIn)
    return;
end

tauStart = real(double(tauStart(:)));
nTau = numel(tauStart);
limitsIn = real(double(limitsIn));
if isvector(limitsIn)
    if numel(limitsIn) == 2 || numel(limitsIn) == 2*nTau
        limits = limitsIn(:).';
    end
elseif isequal(size(limitsIn), [nTau 2])
    limits = limitsIn;
end
end

function y = poissonRandomCounts(lambda)
lambda = max(real(double(lambda(:))), 0);
if exist('poissrnd', 'file') == 2
    y = poissrnd(lambda);
    y = double(y(:));
    return;
end

y = zeros(size(lambda));
small = lambda < 35;
idx = find(small & lambda > 0);
for ii = 1:numel(idx)
    lam = lambda(idx(ii));
    limit = exp(-lam);
    count = 0;
    product = 1;
    while product > limit
        count = count + 1;
        product = product * rand;
    end
    y(idx(ii)) = count - 1;
end

large = ~small;
if any(large)
    approx = lambda(large) + sqrt(lambda(large)) .* randn(sum(large), 1);
    y(large) = max(0, round(approx));
end
end

function [tauAligned, ampAligned] = alignComponentsToReference(refTau, tau, amp)
refTau = real(double(refTau(:)));
tau = real(double(tau(:)));
amp = real(double(amp(:)));
nRef = numel(refTau);
tauAligned = NaN(nRef, 1);
ampAligned = NaN(nRef, 1);

validSample = isfinite(tau) & tau > 0 & isfinite(amp);
tau = tau(validSample);
amp = amp(validSample);
unused = true(numel(tau), 1);

for jj = 1:nRef
    if isempty(tau) || ~isfinite(refTau(jj)) || refTau(jj) <= 0 || ~any(unused)
        continue;
    end
    candidates = find(unused);
    distance = abs(log(tau(candidates) ./ refTau(jj)));
    [~, bestLocal] = min(distance);
    bestIdx = candidates(bestLocal);
    tauAligned(jj) = tau(bestIdx);
    ampAligned(jj) = amp(bestIdx);
    unused(bestIdx) = false;
end
end

function relAmp = relativeAmplitudesFromAmplitudes(amplitudes)
amplitudes = real(double(amplitudes(:)));
positiveAmplitudes = amplitudes;
positiveAmplitudes(~isfinite(positiveAmplitudes) | positiveAmplitudes < 0) = 0;
total = sum(positiveAmplitudes);
relAmp = NaN(size(amplitudes));
if total > 0
    relAmp = positiveAmplitudes ./ total;
end
end

function values = getResampleErrorVector(fitResult, fieldName)
values = [];
if isfield(fitResult, 'resampleError') && isfield(fitResult.resampleError, fieldName)
    values = fitResult.resampleError.(fieldName);
end
end

function values = filterVectorForSummary(values, keepMask)
values = real(double(values(:)));
keepMask = logical(keepMask(:));
if isempty(values) || isempty(keepMask)
    values = [];
    return;
end
n = min(numel(values), numel(keepMask));
values = values(1:n);
keepMask = keepMask(1:n);
values = values(keepMask);
end

function [method, nRequested, nSucceeded, status, message] = summarizeResampleErrorStatus(fitResult)
method = '';
nRequested = NaN;
nSucceeded = NaN;
status = '';
message = '';
if isfield(fitResult, 'resampleError') && isstruct(fitResult.resampleError)
    method = fitResult.resampleError.method;
    nRequested = fitResult.resampleError.nRequested;
    nSucceeded = fitResult.resampleError.nSucceeded;
    status = fitResult.resampleError.status;
    message = fitResult.resampleError.message;
end
end

function s = nanStdColumns(x)
x = double(x);
if isempty(x)
    s = [];
    return;
end
if isvector(x)
    x = x(:);
end
s = NaN(1, size(x, 2));
for jj = 1:size(x, 2)
    col = x(:, jj);
    col = col(isfinite(col));
    if numel(col) >= 2
        s(jj) = std(col, 0);
    elseif numel(col) == 1
        s(jj) = 0;
    end
end
end

function p = nanPercentileColumns(x, pct)
x = double(x);
if isempty(x)
    p = [];
    return;
end
if isvector(x)
    x = x(:);
end
p = NaN(1, size(x, 2));
for jj = 1:size(x, 2)
    col = sort(x(isfinite(x(:, jj)), jj));
    if isempty(col)
        continue;
    end
    if numel(col) == 1
        p(jj) = col(1);
    else
        pos = 1 + (numel(col) - 1) * pct / 100;
        lo = floor(pos);
        hi = ceil(pos);
        if lo == hi
            p(jj) = col(lo);
        else
            p(jj) = col(lo) + (pos - lo) * (col(hi) - col(lo));
        end
    end
end
end

function value = stringHash(txt)
txt = char(txt);
value = 0;
for ii = 1:numel(txt)
    value = mod(value * 131 + double(txt(ii)), 2^31 - 1);
end
end

function val = distRateMeanTau(rates, amplitudes)
rates = real(double(rates(:)));
amplitudes = real(double(amplitudes(:)));
keep = isfinite(rates) & rates > 0 & isfinite(amplitudes) & amplitudes > 0;
if ~any(keep)
    val = NaN;
else
    val = sum(amplitudes(keep)) ./ sum(rates(keep).*amplitudes(keep));
end
end

function plotTailFitFigure(fitInput, fitResult)
figure;
hasDist = isfield(fitResult, 'distTailfitEnabled') && fitResult.distTailfitEnabled && ...
    numel(fitResult.distFit) == numel(fitInput.counts) && any(isfinite(fitResult.distFit));

subplot(3,1,1);
hFull = semilogy(fitInput.fullTimeNs, max(fitInput.fullCounts, 0.1), '.', ...
    'Color', [0.65 0.65 0.65]);
hold on
hData = semilogy(fitInput.tailAbsoluteTimeNs, max(fitInput.counts, 0.1), 'ob', ...
    'MarkerSize', 4);
hTail = semilogy(fitInput.tailAbsoluteTimeNs, max(fitResult.tailFit, 0.1), 'r', ...
    'LineWidth', 1.5);
plotHandles = [hFull hData hTail];
plotLabels = {'TCSPC', 'tail data', 'Tailfit'};
if hasDist
    hDist = semilogy(fitInput.tailAbsoluteTimeNs, max(fitResult.distFit, 0.1), 'g', ...
        'LineWidth', 1.2);
    plotHandles = [plotHandles hDist]; %#ok<AGROW>
    plotLabels{end + 1} = 'DistTailfit'; %#ok<AGROW>
end
hold off
axis tight
xlabel('time [ns]');
ylabel('count');
legend(plotHandles, plotLabels, 'Location', 'best');
if hasDist
    titleText = sprintf('Tail start %.3g ns after peak; Tailfit chi2 %.4g, DistTailfit chi2 %.4g', ...
        fitInput.gateInfo.tailCutAfterPeakNs, fitResult.chi, fitResult.distChi);
else
    titleText = sprintf('Tail start %.3g ns after peak; Tailfit chi2 %.4g', ...
        fitInput.gateInfo.tailCutAfterPeakNs, fitResult.chi);
end
if isfield(fitInput.gateInfo, 'tailRejectLastNPoints') && fitInput.gateInfo.tailRejectLastNPoints > 0
    titleText = sprintf('%s; last %d bins rejected', titleText, fitInput.gateInfo.tailRejectLastNPoints);
end
title(titleText, 'Interpreter', 'none');

subplot(3,1,2);
resTail = (fitInput.counts - fitResult.tailFit) ./ sqrt(max(abs(fitResult.tailFit), eps));
if hasDist
    resDist = (fitInput.counts - fitResult.distFit) ./ sqrt(max(abs(fitResult.distFit), eps));
    plot(fitInput.tailAbsoluteTimeNs, resTail, 'r', ...
        fitInput.tailAbsoluteTimeNs, resDist, 'g', ...
        fitInput.tailAbsoluteTimeNs, 0*fitInput.tailAbsoluteTimeNs, 'k:');
    legend({'Tailfit', 'DistTailfit'}, 'Location', 'best');
else
    plot(fitInput.tailAbsoluteTimeNs, resTail, 'r', ...
        fitInput.tailAbsoluteTimeNs, 0*fitInput.tailAbsoluteTimeNs, 'k:');
    legend({'Tailfit'}, 'Location', 'best');
end
axis tight
xlabel('time [ns]');
ylabel('weighted residual');

subplot(3,1,3);
distAmp = fitResult.distAmplitudes(:);
distTau = fitResult.distTauNs(:);
keep = isfinite(distTau) & distTau > 0 & isfinite(distAmp) & distAmp > 0;
distPlotted = false;
if hasDist && any(keep)
    semilogx(distTau(keep), distAmp(keep)./max(distAmp(keep)), 'bo-');
    hold on
    distPlotted = true;
else
    hold on
end
paramAmp = fitResult.amplitudes(:);
paramTau = fitResult.tauNs(:);
keepParam = isfinite(paramTau) & paramTau > 0 & isfinite(paramAmp) & paramAmp > 0;
paramPlotted = false;
if any(keepParam)
    hStem = stem(paramTau(keepParam), paramAmp(keepParam)./max(paramAmp(keepParam)), 'r');
    set(hStem, 'MarkerFaceColor', 'r');
    paramPlotted = true;
end
hold off
axis tight
xlabel('decay time [ns]');
ylabel('relative amplitude');
if distPlotted && paramPlotted
    legend({'DistTailfit', 'Tailfit'}, 'Location', 'best');
elseif paramPlotted
    legend({'Tailfit'}, 'Location', 'best');
elseif distPlotted
    legend({'DistTailfit'}, 'Location', 'best');
end
end

function [fitResult, fitFigures] = runBatchFluofit(fitInput, opts)
makePlot = opts.plotFits || opts.saveFitFigures;
if makePlot
    figuresBeforeFit = findall(0, 'Type', 'figure');
else
    figuresBeforeFit = [];
end

[cshift, offset, amplitudes, tauNs, dc, dtau, irfShifted, components, tAxisNs, chi] = ...
    Fluofit(fitInput.irf, fitInput.counts, fitInput.pulsePeriodNs, ...
    fitInput.dtNs, fitInput.tau0Ns, fitInput.limitsNs, fitInput.init, ...
    opts.fluofitSolver, makePlot);

if makePlot
    figuresAfterFit = findall(0, 'Type', 'figure');
    fitFigures = newFigureHandles(figuresAfterFit, figuresBeforeFit);
    if isempty(fitFigures) && ~isempty(figuresAfterFit)
        currentFigure = gcf;
        if ~isempty(currentFigure) && ishandle(currentFigure)
            fitFigures = currentFigure;
        end
    end
else
    fitFigures = [];
end

fitResult = struct();
fitResult.signature = 'Fluofit(irf, y, p, dt, tau, lim, init, fitMode, plotFlag)';
fitResult.fitMethod = 'fluofit';
fitResult.fileSuffix = 'fluofit';
fitResult.fitMode = opts.fluofitSolver;
fitResult.colorShiftNs = cshift;
fitResult.offset = offset;
fitResult.amplitudes = amplitudes;
fitResult.tauNs = tauNs;
fitResult.colorShiftErrorNs = dc;
fitResult.tauErrorNs = dtau;
fitResult.irfShifted = irfShifted;
fitResult.components = components;
fitResult.timeAxisNs = tAxisNs;
fitResult.chi = chi;
fitResult.figureFiles = {};
end

function newFigures = newFigureHandles(figuresAfterFit, figuresBeforeFit)
newFigures = figuresAfterFit([]);
for ii = 1:numel(figuresAfterFit)
    isOld = false;
    for jj = 1:numel(figuresBeforeFit)
        if figuresAfterFit(ii) == figuresBeforeFit(jj)
            isOld = true;
            break;
        end
    end
    if ~isOld
        newFigures(end + 1) = figuresAfterFit(ii); %#ok<AGROW>
    end
end
end
function figureFiles = saveFitFigures(figures, resultFolder, outStem, fitInput, opts)
figureFiles = {};
figures = figures(:).';
figures = figures(ishandle(figures));
hasAxes = false(size(figures));
for ii = 1:numel(figures)
    hasAxes(ii) = ~isempty(findall(figures(ii), 'Type', 'axes'));
end
figures = figures(hasAxes);

if isempty(figures)
    error('No Fluofit figure is available to save for %s.', fitInput.file);
end

for jj = 1:numel(figures)
    fig = figures(jj);
    set(fig, 'Name', sprintf('Fluofit: %s', fitInput.curveName), 'NumberTitle', 'off');
    drawnow;

    if numel(figures) == 1
        basePath = fullfile(resultFolder, [outStem '_fit']);
    else
        basePath = fullfile(resultFolder, sprintf('%s_fit_%02d', outStem, jj));
    end

    for ii = 1:numel(opts.figureFormats)
        fmt = opts.figureFormats{ii};
        figPath = [basePath '.' fmt];
        switch fmt
            case 'fig'
                savefig(fig, figPath);
            case 'png'
                if exist('exportgraphics', 'file') == 2
                    exportgraphics(fig, figPath, 'Resolution', 200);
                else
                    saveas(fig, figPath);
                end
            otherwise
                saveas(fig, figPath);
        end
        figureFiles{end + 1} = figPath; %#ok<AGROW>
    end
end
end

function commonFiles = copyFitOutputsToBatchFolder(outputFile, figureFiles, outputFolder, outStem)
commonFiles = struct('outputFile', '', 'figureFiles', {{}});
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

if ~isempty(outputFile) && exist(outputFile, 'file') == 2
    [~, outputName, outputExt] = fileparts(outputFile);
    if isempty(outputName)
        outputName = [outStem '_fit'];
        outputExt = '.mat';
    end
    dst = fullfile(outputFolder, [outputName outputExt]);
    copyFileIfDifferent(outputFile, dst);
    commonFiles.outputFile = dst;
end

for ii = 1:numel(figureFiles)
    src = figureFiles{ii};
    if isempty(src) || exist(src, 'file') ~= 2
        continue;
    end
    [~, name, ext] = fileparts(src);
    dst = fullfile(outputFolder, [name ext]);
    copyFileIfDifferent(src, dst);
    commonFiles.figureFiles{end + 1} = dst; %#ok<AGROW>
end
end

function copyFileIfDifferent(src, dst)
src = char(src);
dst = char(dst);
if strcmpi(src, dst)
    return;
end
[ok, msg] = copyfile(src, dst, 'f');
if ~ok
    error('Could not copy %s to %s: %s', src, dst, msg);
end
end

function irf = sanitizeIrf(irf)
irf = real(double(squeeze(irf)));
irf = irf(:);
irf(~isfinite(irf)) = 0;
irf = max(irf, 0);
s = sum(irf);
if s > 0
    irf = irf ./ s;
end
end

function irfOnAxis = mapIrfToDecayAxis(irfTimeSec, irfCounts, decayTimeSec, dtSec)
irfTimeSec = double(irfTimeSec(:));
irfCounts = double(irfCounts(:));
decayTimeSec = double(decayTimeSec(:));

if numel(irfCounts) == numel(decayTimeSec) && ...
        max(abs((irfTimeSec - irfTimeSec(1)) - (decayTimeSec - decayTimeSec(1)))) <= 0.5 * dtSec
    irfOnAxis = irfCounts;
    return;
end

irfRel = irfTimeSec - irfTimeSec(1);
decayRel = decayTimeSec - decayTimeSec(1);
[irfRel, uniqueIdx] = unique(irfRel, 'stable');
irfCounts = irfCounts(uniqueIdx);
irfOnAxis = interp1(irfRel, irfCounts, decayRel, 'linear', 0);
irfOnAxis = max(double(irfOnAxis(:)), 0);
end

function [timeSec, counts, irfOnAxis, gateInfo] = applyTimeGate(timeSec, counts, irfOnAxis, opts)
timeSec = double(timeSec(:));
counts = double(counts(:));
irfOnAxis = double(irfOnAxis(:));
timeNs = (timeSec - timeSec(1)) * 1e9;

keep = true(size(timeSec));
firstRejected = false;
tailRejected = [];

if opts.rejectFirstTimePoint && ~isempty(keep)
    keep(1) = false;
    firstRejected = true;
end

if ~isempty(opts.rejectTailAtOrAfterNs) && opts.rejectTailNPoints > 0
    tailCandidates = find(timeNs >= opts.rejectTailAtOrAfterNs);
    if ~isempty(tailCandidates)
        nTail = min(opts.rejectTailNPoints, numel(tailCandidates));
        tailRejected = tailCandidates(end - nTail + 1:end);
        keep(tailRejected) = false;
    end
end

gateInfo = struct();
gateInfo.rejectFirstTimePoint = opts.rejectFirstTimePoint;
gateInfo.rejectTailAtOrAfterNs = opts.rejectTailAtOrAfterNs;
gateInfo.rejectTailNPoints = opts.rejectTailNPoints;
gateInfo.firstRejected = firstRejected;
gateInfo.tailRejected = tailRejected(:).';
gateInfo.rejectedTimeNs = timeNs(~keep).';
gateInfo.nRejected = sum(~keep);

timeSec = timeSec(keep);
counts = counts(keep);
irfOnAxis = irfOnAxis(keep);
end

function rawResolutionSec = getPtuResolutionSec(head)
if isfield(head, 'MeasDesc_Resolution') && ~isempty(head.MeasDesc_Resolution)
    rawResolutionSec = double(head.MeasDesc_Resolution(1));
elseif isfield(head, 'Resolution') && ~isempty(head.Resolution)
    rawResolutionSec = double(head.Resolution(1)) * 1e-9;
else
    error('PTU header does not contain MeasDesc_Resolution.');
end
if ~isfinite(rawResolutionSec) || rawResolutionSec <= 0
    error('Invalid PTU TCSPC resolution.');
end
end

function syncRate = getPtuSyncRate(head)
if isfield(head, 'TTResult_SyncRate') && ~isempty(head.TTResult_SyncRate)
    syncRate = double(head.TTResult_SyncRate(1));
elseif isfield(head, 'SyncRate') && ~isempty(head.SyncRate)
    syncRate = double(head.SyncRate(1));
else
    error('PTU header does not contain TTResult_SyncRate.');
end
if ~isfinite(syncRate) || syncRate <= 0
    error('Invalid PTU sync rate.');
end
end

function nRecords = getPtuRecordCount(head)
nRecords = Inf;
if isfield(head, 'TTResult_NumberOfRecords') && ~isempty(head.TTResult_NumberOfRecords)
    nRecords = double(head.TTResult_NumberOfRecords(1));
end
end

function data = readPqresFile(filePath, preferredCurve)
if nargin < 2 || isempty(preferredCurve)
    preferredCurve = 'decay';
end

if isPqresTaggedFile(filePath)
    [tags, version] = readPqresTags(filePath);
    [timeSec, counts, curveName] = selectPqresCurve(tags, preferredCurve);
    meta = parsePqresMeta(tags, timeSec);
else
    [timeSec, counts] = readNumericTextCurve(filePath);
    version = '';
    curveName = 'text';
    meta = parsePqresMeta(struct([]), timeSec);
end

counts = double(counts(:));
counts(~isfinite(counts)) = 0;
counts = max(counts, 0);

data = struct();
data.file = filePath;
data.version = version;
data.curveName = curveName;
data.timeSec = double(timeSec(:));
data.counts = counts;
data.meta = meta;
end

function tf = isPqresTaggedFile(filePath)
fid = fopen(filePath, 'r');
if fid < 0
    error('Could not open %s', filePath);
end
cleanup = onCleanup(@() fclose(fid));
magic = fread(fid, 8, '*char')';
tf = strcmp(stripNulls(magic), 'PQRESLT');
end

function [tags, version] = readPqresTags(filePath)
fid = fopen(filePath, 'r', 'ieee-le');
if fid < 0
    error('Could not open %s', filePath);
end
cleanup = onCleanup(@() fclose(fid));

magic = fread(fid, 8, '*char')';
if ~strcmp(stripNulls(magic), 'PQRESLT')
    error('%s is not a PQRESLT file.', filePath);
end
version = stripNulls(fread(fid, 8, '*char')');

TY_EMPTY8 = uint32(hex2dec('FFFF0008'));
TY_BOOL8 = uint32(hex2dec('00000008'));
TY_INT8 = uint32(hex2dec('10000008'));
TY_BITSET64 = uint32(hex2dec('11000008'));
TY_COLOR8 = uint32(hex2dec('12000008'));
TY_FLOAT8 = uint32(hex2dec('20000008'));
TY_DATETIME = uint32(hex2dec('21000008'));
TY_FLOAT8_ARRAY = uint32(hex2dec('2001FFFF'));
TY_ANSI_STRING = uint32(hex2dec('4001FFFF'));
TY_WIDE_STRING = uint32(hex2dec('4002FFFF'));
TY_BINARY_BLOB = uint32(hex2dec('FFFFFFFF'));

tags = repmat(struct('name', '', 'idx', -1, 'type', uint32(0), 'value', []), 0, 1);

while ~feof(fid)
    identBytes = fread(fid, 32, '*uint8')';
    if numel(identBytes) < 32
        break;
    end

    tag.name = stripNulls(char(identBytes));
    tag.idx = double(fread(fid, 1, 'int32=>int32'));
    tag.type = fread(fid, 1, 'uint32=>uint32');
    payloadBytes = fread(fid, 8, '*uint8')';
    if numel(payloadBytes) < 8
        break;
    end

    rawInt = typecast(uint8(payloadBytes), 'int64');
    rawDouble = typecast(uint8(payloadBytes), 'double');

    switch tag.type
        case TY_EMPTY8
            tag.value = [];
        case {TY_INT8, TY_BITSET64, TY_COLOR8}
            tag.value = double(rawInt);
        case TY_BOOL8
            tag.value = rawInt ~= 0;
        case {TY_FLOAT8, TY_DATETIME}
            tag.value = double(rawDouble);
        case TY_ANSI_STRING
            nBytes = double(rawInt);
            bytes = fread(fid, nBytes, '*uint8')';
            tag.value = stripNulls(char(bytes));
        case TY_FLOAT8_ARRAY
            nBytes = double(rawInt);
            if rem(nBytes, 8) ~= 0
                error('Invalid float array byte count in tag %s of %s.', tag.name, filePath);
            end
            tag.value = fread(fid, nBytes / 8, 'double=>double');
        case {TY_BINARY_BLOB, TY_WIDE_STRING}
            nBytes = double(rawInt);
            tag.value = [];
            fseek(fid, nBytes, 'cof');
        otherwise
            tag.value = double(rawInt);
    end

    tags(end + 1, 1) = tag; %#ok<AGROW>
    if strcmp(tag.name, 'Header_End')
        break;
    end
end
end

function txt = stripNulls(txt)
if isempty(txt)
    txt = '';
    return;
end
zeroPos = find(txt == char(0), 1, 'first');
if ~isempty(zeroPos)
    txt = txt(1:zeroPos - 1);
end
txt = strtrim(txt);
end

function [timeSec, counts, curveName] = selectPqresCurve(tags, preferredCurve)
switch lower(preferredCurve)
    case 'irf'
        prefixes = {'VarIRF', 'VarDecay', 'VarOverallDecay', 'VarSmoothedDecay'};
    otherwise
        prefixes = {'VarDecay', 'VarOverallDecay', 'VarSmoothedDecay', 'VarIRF'};
end

for ii = 1:numel(prefixes)
    prefix = prefixes{ii};
    x = tagValue(tags, [prefix, 'X'], []);
    y = tagValue(tags, [prefix, 'Y'], []);
    if ~isempty(x) && ~isempty(y)
        n = min(numel(x), numel(y));
        if n < 4
            continue;
        end
        timeSec = double(x(1:n));
        counts = double(y(1:n));
        curveName = prefix;
        return;
    end
end

error('No usable PQRES curve found for preferred curve type: %s', preferredCurve);
end

function meta = parsePqresMeta(tags, timeSec)
meta = struct();
meta.periodSec = tagValue(tags, 'VarPeriod', NaN);
meta.deltaPulseSec = tagValue(tags, 'VarDeltaPulse', NaN);
meta.tcspcResolutionSec = tcspcResolutionFromTags(tags);

if isnan(meta.tcspcResolutionSec) && numel(timeSec) > 1
    meta.tcspcResolutionSec = median(diff(double(timeSec(:))));
end
if (isnan(meta.periodSec) || meta.periodSec <= 0) && numel(timeSec) > 1
    meta.periodSec = (max(double(timeSec(:))) - min(double(timeSec(:)))) + meta.tcspcResolutionSec;
end
end

function value = tcspcResolutionFromTags(tags)
value = NaN;
if isempty(tags)
    return;
end

for ii = 1:numel(tags)
    if ~isempty(regexp(tags(ii).name, '^VarTGR_N\d+$', 'once')) && ischar(tags(ii).value)
        if strcmp(tags(ii).value, 'TCSPCResol')
            suffix = regexprep(tags(ii).name, '^VarTGR_N', '');
            value = tagValue(tags, ['VarTGR_V', suffix], NaN);
            return;
        end
    end
end
end

function value = tagValue(tags, name, defaultValue)
value = defaultValue;
if isempty(tags)
    return;
end
idx = find(strcmp({tags.name}, name), 1, 'last');
if ~isempty(idx)
    value = tags(idx).value;
end
end

function [timeSec, counts] = readNumericTextCurve(filePath)
try
    numericData = readmatrix(filePath, 'FileType', 'text');
catch ME
    error(['Could not read numeric data from %s. The file is neither a ', ...
        'PQRESLT tagged file nor readable numeric text. MATLAB error: %s'], filePath, ME.message);
end

numericData = double(numericData);
numericData = numericData(all(isfinite(numericData), 2), :);
numericData = numericData(:, any(isfinite(numericData), 1));

if isempty(numericData)
    error('No numeric text data found in %s.', filePath);
end

if size(numericData, 2) == 1
    counts = numericData(:, 1);
    timeSec = (0:numel(counts) - 1).';
else
    timeSec = numericData(:, 1);
    counts = numericData(:, end);
end
end

function files = removeDuplicateDirEntries(files)
if isempty(files)
    return;
end

paths = cell(numel(files), 1);
for ii = 1:numel(files)
    paths{ii} = lower(fullfile(files(ii).folder, files(ii).name));
end

[~, keep] = unique(paths, 'stable');
files = files(keep);
end

function mask = containsAny(values, needles)
mask = false(size(values));
for ii = 1:numel(values)
    for jj = 1:numel(needles)
        if ~isempty(strfind(values{ii}, needles{jj}))
            mask(ii) = true;
            break;
        end
    end
end
end

function values = lowerCell(values)
for ii = 1:numel(values)
    values{ii} = lower(values{ii});
end
end

function stem = safeFileStem(stem)
stem = regexprep(stem, '[^\w.-]', '_');
if isempty(stem)
    stem = 'tcspc';
end
end
