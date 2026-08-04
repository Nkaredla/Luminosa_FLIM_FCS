function results = fit_natasha_pqres_with_fluofit(dataFolder, outputFolder, opts)
%FIT_NATASHA_PQRES_WITH_FLUOFIT Read PQRES IRF/TCSPC files and fit decays.
%
% Usage:
%   results = fit_natasha_pqres_with_fluofit;
%   results = fit_natasha_pqres_with_fluofit('D:\Luminosa\Data\Natasha');
%
% IRF files are identified by filename. Include one of these tokens in the
% IRF filename: IRF, instrument, response, or prompt. Other .pqres files are
% treated as TCSPC decay files.
%
% Optional fit settings:
%   opts.tau0        - initial lifetime guesses in ns. Default [] lets
%                      DistFluofit choose the component seeds.
%   opts.limits      - Fluofit lifetime limits in ns, default []
%   opts.init        - Fluofit init flag, default 1 when tau0 is empty,
%                      otherwise 0
%   opts.fluofitSolver - 'mle', 'ls', or 'pirls', default 'mle'
%   opts.plotFits    - draw Fluofit figures, default false
%   opts.irfMode     - 'supplied' or 'parametric', default 'supplied'
%   opts.irfModel    - parametric IRF model, default 'calc_mirf'.
%                      Use 'spad_exgauss' for a SPAD Gaussian peak plus
%                      one-sided exponential tail model.
%   opts.irfClipFraction - for Calc_mIRF IRFs, set values below this
%                      fraction of max(IRF) to zero before fitting,
%                      default 1e-3
%   opts.rejectFirstTimePoint - default true for this Natasha dataset
%   opts.rejectTailAtOrAfterNs - default 12.5; set [] to disable
%   opts.rejectTailNPoints - default 8; used with rejectTailAtOrAfterNs
%   opts.summaryAmplitudeThreshold - omit components below this relative
%                      amplitude from the summary CSV only, default 0.02.
%                      Use 0 to disable.

if nargin < 1 || isempty(dataFolder)
    dataFolder = 'D:\Luminosa\Data\Natasha';
end

if nargin < 2 || isempty(outputFolder)
    outputFolder = fullfile(dataFolder, 'fluofit_results');
end

if nargin < 3 || isempty(opts)
    opts = struct();
end
opts = defaultFitOptions(opts);

if ~isfolder(dataFolder)
    error('Input folder does not exist: %s', dataFolder);
end

scriptFolder = fileparts(mfilename('fullpath'));
if ~isempty(scriptFolder)
    addpath(scriptFolder);
end

if exist('Fluofit', 'file') ~= 2
    error(['Could not find Fluofit.m on the MATLAB path. ', ...
        'Add the folder containing Fluofit.m, then rerun this function.']);
end

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

files = dir(fullfile(dataFolder, '*.pqres'));
files = [files; dir(fullfile(dataFolder, '*.PQRES'))];
files = removeDuplicateDirEntries(files);

if isempty(files)
    error('No .pqres files found in %s', dataFolder);
end

names = {files.name};
isIrf = containsAny(lowerCell(names), {'irf', 'instrument', 'response', 'prompt'});
isIrf = isIrf(:);

irfFiles = files(isIrf);
decayFiles = files(~isIrf);
needsIrfFile = irfFileRequired(opts);

if isempty(irfFiles) && needsIrfFile
    error(['No IRF .pqres file was identified. Rename the IRF file so its ', ...
        'filename contains IRF, instrument, response, or prompt.']);
end

if isempty(decayFiles)
    error('No TCSPC decay .pqres files were found after excluding IRF files.');
end

irfData = cell(numel(irfFiles), 1);
if needsIrfFile
    fprintf('Reading %d IRF file(s) from %s\n', numel(irfFiles), dataFolder);
    for ii = 1:numel(irfFiles)
        irfPath = fullfile(irfFiles(ii).folder, irfFiles(ii).name);
        irfData{ii} = readPqresFile(irfPath, 'irf');
    end
else
    fprintf('Using Calc_mIRF parametric IRF; supplied IRF files are not required.\n');
end

results = repmat(struct( ...
    'decayFile', '', ...
    'irfFile', '', ...
    'decayCurve', '', ...
    'irfCurve', '', ...
    'irfMode', '', ...
    'irfModel', '', ...
    'outputFile', '', ...
    'dtNs', NaN, ...
    'pulsePeriodNs', NaN, ...
    'nBins', NaN, ...
    'nRejected', NaN, ...
    'peakCounts', NaN, ...
    'tauNs', '', ...
    'amplitudes', '', ...
    'relativeAmplitudes', '', ...
    'nComponents', NaN, ...
    'nReportedComponents', NaN, ...
    'nRemovedLowAmplitude', NaN, ...
    'chi', NaN, ...
    'status', '', ...
    'message', ''), numel(decayFiles), 1);

for kk = 1:numel(decayFiles)
    decayPath = fullfile(decayFiles(kk).folder, decayFiles(kk).name);
    if needsIrfFile
        irfIdx = chooseIrfForDecay(decayFiles(kk).name, irfFiles);
        irf = irfData{irfIdx};
        irfName = irfFiles(irfIdx).name;
    else
        irf = emptyParametricIrf();
        irfName = 'Calc_mIRF';
    end

    fprintf('Fitting %s with IRF %s (%s mode)\n', ...
        decayFiles(kk).name, irfName, opts.irfMode);

    results(kk).decayFile = decayPath;
    results(kk).irfFile = irf.file;
    results(kk).irfCurve = irf.curveName;

    try
        decay = readPqresFile(decayPath, 'decay');
        fitInput = buildFitInput(decay, irf, opts);
        fitResult = runLocalFluofit(fitInput, opts);

        [~, stem] = fileparts(decayFiles(kk).name);
        outStem = safeFileStem(stem);
        if strcmpi(fitInput.irfMode, 'parametric')
            outStem = [outStem, '_parametric_irf'];
        end
        outPath = fullfile(outputFolder, [outStem, '_fluofit.mat']);
        save(outPath, 'fitInput', 'fitResult');

        results(kk).decayCurve = decay.curveName;
        results(kk).irfMode = fitInput.irfMode;
        results(kk).irfModel = fitInput.irfModel;
        results(kk).outputFile = outPath;
        results(kk).dtNs = fitInput.dtNs;
        results(kk).pulsePeriodNs = fitInput.pulsePeriodNs;
        results(kk).nBins = numel(fitInput.counts);
        results(kk).nRejected = fitInput.gateInfo.nRejected;
        results(kk).peakCounts = max(fitInput.counts);
        summaryComponents = filterSummaryComponents( ...
            fitResult.tauNs, fitResult.amplitudes, opts.summaryAmplitudeThreshold);
        results(kk).tauNs = mat2str(summaryComponents.tauNs(:).', 5);
        results(kk).amplitudes = mat2str(summaryComponents.amplitudes(:).', 5);
        results(kk).relativeAmplitudes = mat2str(summaryComponents.relativeAmplitudes(:).', 5);
        results(kk).nComponents = summaryComponents.nComponents;
        results(kk).nReportedComponents = summaryComponents.nReportedComponents;
        results(kk).nRemovedLowAmplitude = summaryComponents.nRemovedLowAmplitude;
        results(kk).chi = fitResult.chi;
        results(kk).status = 'ok';
        results(kk).message = '';
    catch ME
        results(kk).status = 'failed';
        results(kk).message = ME.message;
        warning('Fit failed for %s: %s', decayFiles(kk).name, ME.message);
    end
end

summaryPath = fullfile(outputFolder, 'fluofit_summary.csv');
writetable(struct2table(results, 'AsArray', true), summaryPath);
fprintf('Saved summary to %s\n', summaryPath);
end

function opts = defaultFitOptions(opts)
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
if ~isfield(opts, 'irfMode') || isempty(opts.irfMode)
    opts.irfMode = 'supplied';
%     opts.irfMode = 'parametric';
end
opts.irfMode = lower(char(opts.irfMode));
if ~any(strcmp(opts.irfMode, {'supplied', 'parametric'}))
    error('opts.irfMode must be ''supplied'' or ''parametric''.');
end
if ~isfield(opts, 'irfModel') || isempty(opts.irfModel)
    opts.irfModel = 'calc_mirf';
end
opts.irfModel = lower(char(opts.irfModel));
opts.irfModel = strrep(opts.irfModel, '-', '_');
validIrfModels = {'calc_mirf', 'calcirf', 'walther', 'irf_fun', ...
    'gaussian', 'spad_exgauss', 'spad_exgaussian', 'exgauss', 'ex_gauss'};
if ~any(strcmp(opts.irfModel, validIrfModels))
    error('opts.irfModel must be ''calc_mirf'', ''gaussian'', or ''spad_exgauss''.');
end
if ~isfield(opts, 'irfParams')
    opts.irfParams = [];
end
if ~isfield(opts, 'spadIrfOptions') || isempty(opts.spadIrfOptions)
    opts.spadIrfOptions = struct();
end
if ~isfield(opts, 'irfClipFraction')
    opts.irfClipFraction = 1e-3;
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

function tf = irfFileRequired(opts)
tf = strcmpi(opts.irfMode, 'supplied') || ...
    (strcmpi(opts.irfMode, 'parametric') && strcmpi(opts.irfModel, 'gaussian'));
end

function irf = emptyParametricIrf()
irf = struct();
irf.file = '';
irf.version = '';
irf.curveName = 'Calc_mIRF';
irf.timeSec = [];
irf.counts = [];
irf.meta = struct('periodSec', NaN, 'deltaPulseSec', NaN, 'tcspcResolutionSec', NaN);
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

function fitInput = buildFitInput(decay, irf, opts)
decayTimeSec = double(decay.timeSec(:));
decayCounts = double(decay.counts(:));
irfTimeSec = double(irf.timeSec(:));
irfCounts = double(irf.counts(:));

if numel(decayTimeSec) ~= numel(decayCounts)
    error('Decay time/count vectors have different lengths.');
end
if numel(irfTimeSec) ~= numel(irfCounts)
    error('IRF time/count vectors have different lengths.');
end

[decayTimeSec, order] = sort(decayTimeSec);
decayCounts = decayCounts(order);
[decayTimeSec, uniqueIdx] = unique(decayTimeSec, 'stable');
decayCounts = decayCounts(uniqueIdx);

dtSec = decay.meta.tcspcResolutionSec;
if isnan(dtSec) || dtSec <= 0
    dtSec = median(diff(decayTimeSec));
end
if isnan(dtSec) || dtSec <= 0
    error('Could not determine TCSPC bin width.');
end

pulsePeriodSec = decay.meta.periodSec;
if isnan(pulsePeriodSec) || pulsePeriodSec <= 0
    pulsePeriodSec = irf.meta.periodSec;
end
if isnan(pulsePeriodSec) || pulsePeriodSec <= 0
    pulsePeriodSec = (max(decayTimeSec) - min(decayTimeSec)) + dtSec;
end

irfInfo = struct();
[irfOnDecayAxis, irfInfo] = makeIrfOnDecayAxis( ...
    irfTimeSec, irfCounts, decayTimeSec, decayCounts, dtSec, pulsePeriodSec, opts);
if ~any(irfOnDecayAxis > 0)
    error('Mapped IRF is all zeros; IRF and decay time axes do not overlap.');
end

[decayTimeSec, decayCounts, irfOnDecayAxis, gateInfo] = ...
    applyTimeGate(decayTimeSec, decayCounts, irfOnDecayAxis, opts);
if numel(decayCounts) < 4
    error('Time gating left fewer than four decay samples.');
end
if ~any(irfOnDecayAxis > 0)
    error('Time gating removed all positive IRF samples.');
end

fitInput = struct();
fitInput.decayFile = decay.file;
fitInput.irfFile = irf.file;
fitInput.decayCurve = decay.curveName;
fitInput.irfCurve = irf.curveName;
fitInput.irfMode = opts.irfMode;
fitInput.irfModel = '';
if strcmpi(opts.irfMode, 'parametric')
    fitInput.irfModel = opts.irfModel;
end
fitInput.irfInfo = irfInfo;
fitInput.gateInfo = gateInfo;
fitInput.timeNs = (decayTimeSec - decayTimeSec(1)) * 1e9;
fitInput.counts = decayCounts(:);
fitInput.irf = irfOnDecayAxis(:);
fitInput.dtNs = dtSec * 1e9;
fitInput.pulsePeriodNs = pulsePeriodSec * 1e9;
fitInput.tau0Ns = double(opts.tau0(:)).';
fitInput.limitsNs = opts.limits;
fitInput.init = opts.init;
end

function [irfOnAxis, info] = makeIrfOnDecayAxis(irfTimeSec, irfCounts, ...
    decayTimeSec, decayCounts, dtSec, pulsePeriodSec, opts)
switch lower(opts.irfMode)
    case 'supplied'
        irfOnAxis = mapIrfToDecayAxis(irfTimeSec, irfCounts, decayTimeSec, dtSec);
        info = struct('mode', 'supplied', 'model', '', 'params', []);
    case 'parametric'
        [irfOnAxis, params] = parametricIrfOnDecayAxis( ...
            irfTimeSec, irfCounts, decayTimeSec, decayCounts, ...
            dtSec, pulsePeriodSec, opts);
        info = struct('mode', 'parametric', 'model', opts.irfModel, 'params', params);
    otherwise
        error('Unsupported IRF mode: %s', opts.irfMode);
end
end

function irfOnAxis = mapIrfToDecayAxis(irfTimeSec, irfCounts, decayTimeSec, dtSec)
irfTimeSec = double(irfTimeSec(:));
irfCounts = double(irfCounts(:));
decayTimeSec = double(decayTimeSec(:));

if numel(irfCounts) == numel(decayTimeSec) && max(abs(irfTimeSec - decayTimeSec)) <= 0.5 * dtSec
    irfOnAxis = irfCounts;
    return;
end

bin = round((irfTimeSec - decayTimeSec(1)) ./ dtSec) + 1;
valid = isfinite(bin) & bin >= 1 & bin <= numel(decayTimeSec);

if any(valid)
    irfOnAxis = accumarray(bin(valid), irfCounts(valid), [numel(decayTimeSec), 1], @sum, 0);
else
    [irfTimeSec, uniqueIdx] = unique(irfTimeSec, 'stable');
    irfCounts = irfCounts(uniqueIdx);
    irfOnAxis = interp1(irfTimeSec, irfCounts, decayTimeSec, 'linear', 0);
end

irfOnAxis = max(double(irfOnAxis(:)), 0);
end

function [irfOnAxis, params] = parametricIrfOnDecayAxis(irfTimeSec, irfCounts, ...
    decayTimeSec, decayCounts, dtSec, pulsePeriodSec, opts)
switch lower(opts.irfModel)
    case {'calc_mirf', 'calcirf', 'walther', 'irf_fun'}
        [irfOnAxis, params] = calcMirfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec, opts);
    case {'spad_exgauss', 'spad_exgaussian', 'exgauss', 'ex_gauss'}
        [irfOnAxis, params] = spadExGaussIrfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec, opts);
    case 'gaussian'
        if ~isempty(opts.irfParams)
            params = parseIrfParams(opts.irfParams);
        else
            params = fitGaussianIrfParams(irfTimeSec, irfCounts, decayTimeSec);
        end

        axisNs = (double(decayTimeSec(:)) - double(decayTimeSec(1))) * 1e9;
        irfOnAxis = gaussianIrf(axisNs, params);
    otherwise
        error('Unsupported parametric IRF model: %s', opts.irfModel);
end

irfOnAxis = max(double(irfOnAxis(:)), 0);
end

function [irfOnAxis, params] = calcMirfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec, opts)
if exist('Calc_mIRF', 'file') ~= 2
    error('Could not find Calc_mIRF.m on the MATLAB path.');
end

head = struct();
head.Resolution = dtSec * 1e9;
head.SyncRate = 1 / pulsePeriodSec;

tcspc = double(decayCounts(:).');
irf = Calc_mIRF(head, tcspc);
[irf, clipInfo] = clipIrfBelowFraction(irf, opts.irfClipFraction);
irfOnAxis = squeeze(irf);
irfOnAxis = double(irfOnAxis(:));

params = struct();
params.source = 'Calc_mIRF';
params.model = 'IRF_Fun';
params.head = head;
params.clipInfo = clipInfo;
params.note = 'Estimated by existing Calc_mIRF(head, tcspc) code.';
end

function [irfOnAxis, params] = spadExGaussIrfOnDecayAxis(decayCounts, dtSec, pulsePeriodSec, opts)
if exist('Calc_mIRF_Global_ExGauss_fast', 'file') ~= 2
    error('Could not find Calc_mIRF_Global_ExGauss_fast.m.');
end

head = struct();
head.Resolution = dtSec * 1e9;
head.SyncRate = 1 / pulsePeriodSec;

out = Calc_mIRF_Global_ExGauss_fast(head, double(decayCounts(:)), opts.tau0, opts.spadIrfOptions);
[irf, clipInfo] = clipIrfBelowFraction(out.IRF, opts.irfClipFraction);
irfOnAxis = double(irf(:));

params = out;
params.source = 'Calc_mIRF_Global_ExGauss_fast';
params.model = 'SPAD ex-Gaussian';
params.modelKey = 'spad_exgauss';
params.head = head;
params.clipInfo = clipInfo;
params.note = ['SPAD IRF model: Gaussian avalanche timing peak convolved ', ...
    'with a one-sided exponential diffusion/timing tail.'];
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

function params = fitGaussianIrfParams(irfTimeSec, irfCounts, decayTimeSec)
xNs = (double(irfTimeSec(:)) - double(decayTimeSec(1))) * 1e9;
y = max(double(irfCounts(:)), 0);

finiteMask = isfinite(xNs) & isfinite(y);
xNs = xNs(finiteMask);
y = y(finiteMask);
if numel(y) < 4 || ~any(y > 0)
    error('Cannot fit parametric IRF: not enough positive IRF samples.');
end

ys = sort(y);
nBaseline = max(1, ceil(0.2 * numel(ys)));
baseline0 = median(ys(1:nBaseline));
y0 = max(y - baseline0, 0);
if ~any(y0 > 0)
    y0 = y;
    baseline0 = min(y);
end

amp0 = max(y0);
mu0 = sum(xNs .* y0) / sum(y0);
sigma0 = sqrt(max(sum(((xNs - mu0) .^ 2) .* y0) / sum(y0), eps));
theta0 = [log(max(amp0, eps)), mu0, log(max(sigma0, eps)), log(max(baseline0, eps))];

objective = @(theta) sum((sqrt(gaussianIrf(xNs, thetaToGaussianParams(theta)) + 1) - sqrt(y + 1)) .^ 2);
fitOpts = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 4000);
theta = fminsearch(objective, theta0, fitOpts);
params = thetaToGaussianParams(theta);
end

function params = thetaToGaussianParams(theta)
params = struct();
params.amplitude = exp(theta(1));
params.centerNs = theta(2);
params.sigmaNs = exp(theta(3));
params.baseline = exp(theta(4));
end

function params = parseIrfParams(rawParams)
if isstruct(rawParams)
    params = rawParams;
else
    vals = double(rawParams(:));
    if numel(vals) < 3
        error('opts.irfParams must contain at least [amplitude centerNs sigmaNs].');
    end
    params = struct();
    params.amplitude = vals(1);
    params.centerNs = vals(2);
    params.sigmaNs = vals(3);
    if numel(vals) >= 4
        params.baseline = vals(4);
    else
        params.baseline = 0;
    end
end

required = {'amplitude', 'centerNs', 'sigmaNs'};
for ii = 1:numel(required)
    if ~isfield(params, required{ii}) || isempty(params.(required{ii}))
        error('opts.irfParams is missing field %s.', required{ii});
    end
end
if ~isfield(params, 'baseline') || isempty(params.baseline)
    params.baseline = 0;
end
end

function y = gaussianIrf(xNs, params)
y = params.baseline + params.amplitude .* ...
    exp(-0.5 .* ((xNs - params.centerNs) ./ max(params.sigmaNs, eps)) .^ 2);
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

function fitResult = runLocalFluofit(fitInput, opts)
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
end

function idx = chooseIrfForDecay(decayName, irfFiles)
if numel(irfFiles) == 1
    idx = 1;
    return;
end

decayTokens = fileTokens(decayName);
bestScore = -Inf;
idx = 1;

for ii = 1:numel(irfFiles)
    irfTokens = fileTokens(irfFiles(ii).name);
    score = sum(ismember(decayTokens, irfTokens));
    if score > bestScore
        bestScore = score;
        idx = ii;
    end
end
end

function tokens = fileTokens(fileName)
[~, stem] = fileparts(lower(fileName));
tokens = regexp(stem, '[a-z0-9]+', 'match');
drop = {'irf', 'instrument', 'response', 'prompt', 'decay', 'tcspc'};
tokens = setdiff(tokens, drop, 'stable');
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

function stem = safeFileStem(stem)
stem = regexprep(stem, '[^\w.-]', '_');
if isempty(stem)
    stem = 'decay';
end
end
