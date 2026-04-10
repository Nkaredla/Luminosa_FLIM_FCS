folderName = 'D:\Luminosa\Data\270326\nk_b7h6_miet_64\';
resultFolderName = 'CellSeg2Color_results';

cfg = struct();
cfg.segDetectorIndex = 1;   % d1
cfg.fitDetectorIndex = 2;   % d2
cfg.segPulseIndex = 1;      % channel 1
cfg.fitPulseIndex = 2;      % channel 2
cfg.minCellDiameterUm = 3.0;
cfg.maxPixelSizeUmToProcess = 0.1;
cfg.segmentationMode = 'variance';
cfg.segVarianceWindow = 3;
cfg.segVarianceThresholdScale = 2;
cfg.segVarianceBackgroundQuantile = 0.6;
cfg.segRemoveBorderObjects = false;
cfg.segBorderMarginPixels = 1;
cfg.segThresholdScale = 0.9;
cfg.segAdaptiveSensitivity = 0.52;
cfg.segTophatRadiusUm = 1.0;
cfg.segCloseRadiusUm = 0.8;
cfg.segCloseRadiusPixelsMax = 5;
cfg.segFlatFilterSize = 9;
cfg.segNormalizedClusterPeak = 250;
cfg.segPrefilterMinArea = 30;
cfg.tau0 = [0.4 2.5];
cfg.useGuiStyleMultiExpFit = true;
cfg.usePixelwiseIntensityLifetimeMap = true;
cfg.includeBackgroundInFit = true;
cfg.pixelwiseIncludeBackground = true;
cfg.pixelwiseMinPhotons = 10;
cfg.showFigures = false;
cfg.saveOutputs = true;
cfg.saveMat = true;
cfg.saveCsv = true;
cfg.saveFigures = true;
cfg.closeFiguresAfterSave = true;
cfg.verbose = true;
cfg.plotProgress = true;
cfg.closeProgressFigure = true;

targets = collectBatchTargets(folderName);
if isempty(targets)
    error('No PTU files were found inside subfolders of %s', folderName);
end

batchSummary = repmat(struct( ...
    'ptuPath', '', ...
    'outputDir', '', ...
    'pixelSizeUm', NaN, ...
    'nCells', NaN, ...
    'backgroundTauNs', NaN, ...
    'status', '', ...
    'message', ''), numel(targets), 1);

for k = 1:numel(targets)
    target = targets(k);
    fprintf('Processing %d/%d: %s\n', k, numel(targets), target.ptuPath);

    localCfg = cfg;
    localCfg.outputDir = fullfile(target.folderPath, resultFolderName);
    localCfg.outputStem = target.baseName;

    batchSummary(k).ptuPath = target.ptuPath;
    batchSummary(k).outputDir = localCfg.outputDir;

    try
        head = PTU_Read_Head(target.ptuPath);
        pixelSizeUm = localResolvePixelSizeUmFromHead(head);
        batchSummary(k).pixelSizeUm = pixelSizeUm;
        if isfinite(pixelSizeUm) && isfield(localCfg, 'maxPixelSizeUmToProcess') && ...
                isfinite(localCfg.maxPixelSizeUmToProcess) && localCfg.maxPixelSizeUmToProcess > 0 && ...
                pixelSizeUm > localCfg.maxPixelSizeUmToProcess
            batchSummary(k).status = 'skipped';
            batchSummary(k).message = sprintf('Skipped: pixel size %.4f um is above max threshold %.4f um.', ...
                pixelSizeUm, localCfg.maxPixelSizeUmToProcess);
            fprintf('%s\n', batchSummary(k).message);
            continue;
        end
        res = PTU_CellSeg2ColorLifetime(target.ptuPath, localCfg);
        batchSummary(k).pixelSizeUm = res.pixelSizeUm;
        batchSummary(k).nCells = numel(res.rois);
        batchSummary(k).backgroundTauNs = res.background.fitCh2D2.tauIntNs;
        batchSummary(k).status = 'ok';
        batchSummary(k).message = '';
    catch ME
        batchSummary(k).status = 'failed';
        batchSummary(k).message = ME.message;
        fprintf('Failed: %s\n', ME.message);
        for s = 1:numel(ME.stack)
            fprintf('  at %s (line %d)\n', ME.stack(s).name, ME.stack(s).line);
        end
    end
end

summaryTable = struct2table(batchSummary);
summaryPath = fullfile(folderName, 'CellSeg2Color_batch_summary.csv');
if exist('writetable', 'file') == 2
    writetable(summaryTable, summaryPath);
else
    writecell([summaryTable.Properties.VariableNames; table2cell(summaryTable)], summaryPath);
end


function targets = collectBatchTargets(folderName)
targets = struct('ptuPath', {}, 'folderPath', {}, 'baseName', {});

entries = dir(folderName);
for k = 1:numel(entries)
    if ~entries(k).isdir
        continue;
    end
    entryName = entries(k).name;
    if strcmp(entryName, '.') || strcmp(entryName, '..')
        continue;
    end

    subfolder = fullfile(folderName, entryName);
    ptuFiles = dir(fullfile(subfolder, '*.ptu'));
    for j = 1:numel(ptuFiles)
        ptuPath = fullfile(subfolder, ptuFiles(j).name);
        [~, stem] = fileparts(ptuFiles(j).name);
        targets(end+1).ptuPath = ptuPath; %#ok<AGROW>
        targets(end).folderPath = subfolder;
        targets(end).baseName = stem;
    end
end
end

function pixelSizeUm = localResolvePixelSizeUmFromHead(head)
pixelSizeUm = NaN;
if ~isstruct(head) || ~isfield(head, 'ImgHdr_PixResol') || isempty(head.ImgHdr_PixResol)
    return;
end
pixelSizeUm = double(head.ImgHdr_PixResol);
if ~isfinite(pixelSizeUm) || pixelSizeUm <= 0
    pixelSizeUm = NaN;
    return;
end
if pixelSizeUm < 1e-3
    pixelSizeUm = pixelSizeUm * 1e6;
end
end
