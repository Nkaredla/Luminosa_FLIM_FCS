function batch = fit_tcspc_subfolders_with_fluofit(rootFolder, outputFolder, opts)
%FIT_TCSPC_SUBFOLDERS_WITH_FLUOFIT Batch fit TCSPC files in subfolders.
%
% Usage:
%   batch = fit_tcspc_subfolders_with_fluofit('D:\Luminosa\Data\Natasha')
%
% Expected layout:
%   rootFolder\
%       sample_001\one_file.ptu
%       sample_002\one_file.pqres
%       ...
%
% The function reads one PTU or PQRES file from each immediate subfolder,
% collects the TCSPC curves, estimates one global IRF from the collected
% curves when possible, fits every curve with Fluofit, and saves:
%   - one MAT file per fitted curve
%   - batch_irf.mat
%   - global_irf.mat when opts.irfMode='global'
%   - tcspc_batch_fluofit_summary.csv
%   - tcspc_batch_fluofit_all.mat
%
% Main options:
%   opts.tau0                  lifetime starts in ns, default [0.3 1.0 2.0 4.0]
%   opts.limits                Fluofit lifetime limits in ns, default []
%   opts.init                  Fluofit init flag, default 0
%   opts.fluofitSolver         'mle', 'ls', or 'pirls', default 'mle'
%   opts.plotFits              draw Fluofit figures, default false
%   opts.irfMode               'global', 'supplied', or 'per_curve',
%                              default 'global'
%   opts.irfFile               supplied PTU/PQRES IRF file for
%                              opts.irfMode='supplied'
%   opts.globalIrfMethod       'calc_mirf', 'gamma_shifted_fast', or
%                              'exgauss_fast', default 'calc_mirf'
%   opts.rejectFirstTimePoint  default true for the Natasha dataset
%   opts.rejectTailAtOrAfterNs default 12.5; set [] to disable
%   opts.rejectTailNPoints     default 8

if nargin < 1 || isempty(rootFolder)
    rootFolder = 'D:\Luminosa\Data\Natasha';
end
if nargin < 2 || isempty(outputFolder)
    outputFolder = fullfile(rootFolder, 'tcspc_batch_fluofit_results');
end
if nargin < 3 || isempty(opts)
    opts = struct();
end
opts = defaultBatchOptions(opts);

if ~isfolder(rootFolder)
    error('Input root folder does not exist: %s', rootFolder);
end

scriptFolder = fileparts(mfilename('fullpath'));
if ~isempty(scriptFolder)
    addpath(scriptFolder);
end
if exist('Fluofit', 'file') ~= 2
    error('Could not find Fluofit.m on the MATLAB path.');
end
if exist('Calc_mIRF', 'file') ~= 2
    error('Could not find Calc_mIRF.m on the MATLAB path.');
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

for kk = 1:numel(targets)
    target = targets(kk);
    summary(kk).folder = target.folder;
    summary(kk).file = target.file;
    summary(kk).fileType = target.fileType;

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
        summary(kk).status = 'failed';
        summary(kk).message = ME.message;
        warning('Read failed for %s: %s', target.file, ME.message);
    end
end

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
    end
else
    fprintf('Batch IRF not available: %s\n', batchIrf.message);
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
        fitInput = buildBatchFitInput(curves(kk), batchIrf, opts);
        fitResult = runBatchFluofit(fitInput, opts);

        [~, stem] = fileparts(curves(kk).file);
        outStem = safeFileStem([curves(kk).folderName '_' stem]);
        outPath = fullfile(outputFolder, [outStem '_fluofit.mat']);
        save(outPath, 'fitInput', 'fitResult', '-v7.3');

        fits(kk).fitInput = fitInput;
        fits(kk).fitResult = fitResult;
        fits(kk).outputFile = outPath;
        fits(kk).status = 'ok';
        fits(kk).message = '';

        summary(kk).irfSource = fitInput.irfSource;
        summary(kk).outputFile = outPath;
        summary(kk).dtNs = fitInput.dtNs;
        summary(kk).pulsePeriodNs = fitInput.pulsePeriodNs;
        summary(kk).nBins = numel(fitInput.counts);
        summary(kk).nRejected = fitInput.gateInfo.nRejected;
        summary(kk).tauNs = mat2str(fitResult.tauNs(:).', 5);
        summary(kk).amplitudes = mat2str(fitResult.amplitudes(:).', 5);
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
batch.fits = fits;
batch.results = summary;
batch.summaryPath = summaryPath;

allPath = fullfile(outputFolder, 'tcspc_batch_fluofit_all.mat');
save(allPath, 'batch', '-v7.3');
fprintf('Saved summary to %s\n', summaryPath);
fprintf('Saved combined batch data to %s\n', allPath);
end

function opts = defaultBatchOptions(opts)
if ~isfield(opts, 'tau0') || isempty(opts.tau0)
    opts.tau0 = [0.3 1.0 2.0 4.0];
end
if ~isfield(opts, 'limits')
    opts.limits = [];
end
if ~isfield(opts, 'init') || isempty(opts.init)
    opts.init = 0;
end
if isempty(opts.tau0) && opts.init == 0
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
if ~isfield(opts, 'irfMode') || isempty(opts.irfMode)
    opts.irfMode = 'global';
end
opts.irfMode = lower(strrep(char(opts.irfMode), '-', '_'));
if ~any(strcmp(opts.irfMode, {'global', 'supplied', 'per_curve'}))
    error('opts.irfMode must be ''global'', ''supplied'', or ''per_curve''.');
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
    'gamma_fast', 'exgauss_fast', 'ex_gauss_fast'};
if ~any(strcmp(opts.globalIrfMethod, validGlobalMethods))
    error('Unsupported opts.globalIrfMethod: %s', opts.globalIrfMethod);
end
if ~isfield(opts, 'globalIrfOptions') || isempty(opts.globalIrfOptions)
    opts.globalIrfOptions = struct();
end
if ~isfield(opts, 'globalIrfResample') || isempty(opts.globalIrfResample)
    opts.globalIrfResample = true;
end
if ~isfield(opts, 'globalIrfPeriodToleranceSec') || isempty(opts.globalIrfPeriodToleranceSec)
    opts.globalIrfPeriodToleranceSec = 1e-12;
end
if ~isfield(opts, 'fallbackPerCurveIrf') || isempty(opts.fallbackPerCurveIrf)
    opts.fallbackPerCurveIrf = true;
end
if ~isfield(opts, 'includeRootFiles') || isempty(opts.includeRootFiles)
    opts.includeRootFiles = false;
end
if ~isfield(opts, 'multipleFileMode') || isempty(opts.multipleFileMode)
    opts.multipleFileMode = 'error';
end
opts.multipleFileMode = lower(char(opts.multipleFileMode));
if ~any(strcmp(opts.multipleFileMode, {'error', 'first'}))
    error('opts.multipleFileMode must be ''error'' or ''first''.');
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

function target = emptyTarget()
target = struct('folder', '', 'folderName', '', 'file', '', 'fileType', '', ...
    'status', 'pending', 'message', '');
end

function curve = emptyCurve()
curve = struct('folder', '', 'folderName', '', 'file', '', 'fileType', '', ...
    'curveName', '', 'timeSec', [], 'counts', [], 'meta', struct(), ...
    'head', [], 'rawTcspc', [], 'nChannels', NaN, 'status', 'pending', ...
    'message', '');
end

function fit = emptyFit()
fit = struct('folder', '', 'file', '', 'fitInput', [], 'fitResult', [], ...
    'outputFile', '', 'status', 'pending', 'message', '');
end

function summary = emptySummary()
summary = struct('folder', '', 'file', '', 'fileType', '', 'irfSource', '', ...
    'outputFile', '', 'dtNs', NaN, 'pulsePeriodNs', NaN, 'nBins', NaN, ...
    'nRejected', NaN, 'totalCounts', NaN, 'peakCounts', NaN, 'tauNs', '', ...
    'amplitudes', '', 'chi', NaN, 'status', 'pending', 'message', '');
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

files = listTcspcFiles(folder);
files = removeExplicitIrfFile(files, opts);
if isempty(files)
    target.status = 'skipped';
    target.message = 'No .ptu or .pqres file found.';
    return;
end
if strcmp(opts.irfMode, 'supplied') && numel(files) > 1
    isIrf = isIrfFileName({files.name});
    decayFiles = files(~isIrf);
    if numel(decayFiles) == 1
        files = decayFiles;
    end
end
if numel(files) > 1 && strcmp(opts.multipleFileMode, 'error')
    target.status = 'failed';
    target.message = sprintf('Expected one .ptu or .pqres file, found %d.', numel(files));
    return;
end

target.file = fullfile(files(1).folder, files(1).name);
[~, ~, ext] = fileparts(target.file);
target.fileType = lower(ext(2:end));
target.status = 'ok';
target.message = '';
end

function files = listTcspcFiles(folder)
files = [dir(fullfile(folder, '*.ptu')); dir(fullfile(folder, '*.PTU')); ...
    dir(fullfile(folder, '*.pqres')); dir(fullfile(folder, '*.PQRES'))];
files = removeDuplicateDirEntries(files);
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
        params = struct('source', 'Calc_mIRF', 'model', 'IRF_Fun', ...
            'note', 'Estimated by Calc_mIRF from summed collected TCSPC curves.');
    case {'gamma_shifted_fast', 'gamma_fast'}
        if exist('Calc_mIRF_Global_GammaShifted_fast', 'file') ~= 2
            error('Could not find Calc_mIRF_Global_GammaShifted_fast.m.');
        end
        out = Calc_mIRF_Global_GammaShifted_fast(head, double(summed(:)), opts.tau0, opts.globalIrfOptions);
        irf = out.IRF;
        params = out;
    case {'exgauss_fast', 'ex_gauss_fast'}
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
if any(strcmp(opts.irfMode, {'global', 'supplied'})) && strcmp(globalIrf.status, 'ok')
    irfOnAxis = mapIrfToDecayAxis(globalIrf.timeSec, globalIrf.counts, decayTimeSec, dtSec);
    if strcmp(opts.irfMode, 'global')
        irfSource = ['global_' globalIrf.method];
    else
        irfSource = globalIrf.method;
    end
    irfInfo = globalIrf;
else
    irfOnAxis = zeros(size(decayCounts));
end

if ~any(irfOnAxis > 0)
    if ~opts.fallbackPerCurveIrf
        error('No usable batch IRF for %s, and fallbackPerCurveIrf is false.', curve.file);
    end
    [irfOnAxis, params] = calcMirfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec);
    irfSource = 'per_curve_calc_mirf';
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

function [irfOnAxis, params] = calcMirfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec)
head = struct();
head.Resolution = dtSec * 1e9;
head.SyncRate = 1 / pulsePeriodSec;
irf = Calc_mIRF(head, double(decayCounts(:).'));
irfOnAxis = sanitizeIrf(irf);
params = struct('source', 'Calc_mIRF', 'model', 'IRF_Fun', 'head', head, ...
    'note', 'Estimated per curve by Calc_mIRF(head, tcspc).');
end

function fitResult = runBatchFluofit(fitInput, opts)
[cshift, offset, amplitudes, tauNs, dc, dtau, irfShifted, components, tAxisNs, chi] = ...
    Fluofit(fitInput.irf, fitInput.counts, fitInput.pulsePeriodNs, ...
    fitInput.dtNs, fitInput.tau0Ns, fitInput.limitsNs, fitInput.init, ...
    opts.fluofitSolver, opts.plotFits);

fitResult = struct();
fitResult.signature = 'Fluofit(irf, y, p, dt, tau, lim, init, fitMode, plotFlag)';
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
