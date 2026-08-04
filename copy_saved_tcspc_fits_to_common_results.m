function copied = copy_saved_tcspc_fits_to_common_results(rootFolder, outputFolder, opts)
%COPY_SAVED_TCSPC_FITS_TO_COMMON_RESULTS Collect saved TCSPC fit outputs.
%
% Usage:
%   copied = copy_saved_tcspc_fits_to_common_results(rootFolder)
%
% The function scans each main subfolder for the per-sample result folder
% created by fit_tcspc_subfolders_with_fluofit, then copies saved fit MAT
% files and fit figures into the common batch results folder.

if nargin < 1 || isempty(rootFolder)
    rootFolder = 'D:\Luminosa\Data\Natasha';
end
if nargin < 2
    outputFolder = '';
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

if ~isfield(opts, 'resultFolderName') || isempty(opts.resultFolderName)
    opts.resultFolderName = 'tcspc_fluofit_results';
end
if ~isfield(opts, 'batchResultFolderName') || isempty(opts.batchResultFolderName)
    opts.batchResultFolderName = 'tcspc_batch_fluofit_results';
end
if ~isfield(opts, 'filePatterns') || isempty(opts.filePatterns)
    opts.filePatterns = {'*_fluofit.mat', '*_tailfit.mat', '*_fit*.png', '*_fit*.fig'};
end
if ischar(opts.filePatterns) || isstring(opts.filePatterns)
    opts.filePatterns = cellstr(opts.filePatterns);
end

if isempty(outputFolder)
    outputFolder = fullfile(rootFolder, opts.batchResultFolderName);
end
if ~isfolder(rootFolder)
    error('Input root folder does not exist: %s', rootFolder);
end
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

resultFolders = findPerSampleResultFolders(rootFolder, opts);
copied = repmat(emptyCopyRecord(), 0, 1);

for kk = 1:numel(resultFolders)
    resultFolder = resultFolders{kk};
    for pp = 1:numel(opts.filePatterns)
        files = dir(fullfile(resultFolder, opts.filePatterns{pp}));
        for ii = 1:numel(files)
            src = fullfile(files(ii).folder, files(ii).name);
            dst = fullfile(outputFolder, files(ii).name);
            if strcmpi(src, dst)
                continue;
            end
            [ok, msg] = copyfile(src, dst, 'f');
            if ~ok
                error('Could not copy %s to %s: %s', src, dst, msg);
            end

            rec = emptyCopyRecord();
            rec.sourceFile = src;
            rec.commonFile = dst;
            rec.bytes = files(ii).bytes;
            rec.status = 'ok';
            copied(end + 1, 1) = rec; %#ok<AGROW>
        end
    end
end

if ~isempty(copied)
    copiedTable = struct2table(copied);
    summaryPath = fullfile(outputFolder, 'tcspc_common_fit_files.csv');
    writetable(copiedTable, summaryPath);
    fprintf('Copied %d TCSPC fit file(s) to %s\n', numel(copied), outputFolder);
    fprintf('Saved copy summary to %s\n', summaryPath);
else
    fprintf('No saved TCSPC fit files found under %s\n', rootFolder);
end
end

function resultFolders = findPerSampleResultFolders(rootFolder, opts)
resultFolders = {};
entries = dir(rootFolder);
for ii = 1:numel(entries)
    if ~entries(ii).isdir
        continue;
    end
    name = entries(ii).name;
    if strcmp(name, '.') || strcmp(name, '..') || startsWith(name, '.')
        continue;
    end
    if strcmpi(name, opts.batchResultFolderName)
        continue;
    end

    resultFolder = fullfile(entries(ii).folder, name, opts.resultFolderName);
    if isfolder(resultFolder)
        resultFolders{end + 1, 1} = resultFolder; %#ok<AGROW>
    end
end
end

function rec = emptyCopyRecord()
rec = struct('sourceFile', '', 'commonFile', '', 'bytes', NaN, 'status', '');
end
