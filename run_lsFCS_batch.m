% Batch lsFCS correlation over subfolders.
%
% Each subfolder under rootDir is expected to contain exactly one .ptu file.
% The output for each file is saved by lsFCS as *_lsFCS.mat next to the .ptu.

% === User settings ===
rootDir   = 'D:\Luminosa\Data\270226\POPC_DiO_5pc_NTA_B7H6_Atto488_20\'; % <-- change this
cnum      = 2;     % number of pulses in one sync (set to your experiment)
maxtime   = 10;    % max correlation time (seconds)
timegates = [];    % [] = auto-detect; or Nx2 array [start len] per pulse
overwrite = false; % set true to re-run even if output exists

% === Batch processing ===
dirInfo = dir(rootDir);
for k = 1:numel(dirInfo)
    if ~dirInfo(k).isdir
        continue;
    end
    name = dirInfo(k).name;
    if strcmp(name,'.') || strcmp(name,'..')
        continue;
    end

    subDir = fullfile(rootDir, name);
    ptuFiles = dir(fullfile(subDir, '*.ptu'));

    if isempty(ptuFiles)
        fprintf('No PTU found in %s\n', subDir);
        continue;
    end
    if numel(ptuFiles) > 1
        fprintf('Multiple PTUs found in %s (skipping)\n', subDir);
        continue;
    end

    ptuPath = fullfile(subDir, ptuFiles(1).name);
    outMat = [ptuPath(1:end-4), '_lsFCS.mat'];
    if ~overwrite && exist(outMat, 'file')
        fprintf('Output exists, skipping: %s\n', outMat);
        continue;
    end

    fprintf('Processing %s\n', ptuPath);
    try
        lsFCS(ptuPath, cnum, maxtime, timegates);
    catch ME
        fprintf('Failed %s: %s\n', ptuPath, ME.message);
    end
end
