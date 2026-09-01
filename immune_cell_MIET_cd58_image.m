function out = immune_cell_MIET_cd58_image(ptuFile, cfg)
%IMMUNE_CELL_MIET_CD58_IMAGE Mean intensity image of the 485-excited channel.
%
% out = immune_cell_MIET_cd58_image(ptuFile)
% out = immune_cell_MIET_cd58_image(ptuFile, cfg)
%
% Builds the CD58-Atto488 intensity image from the FIRST PIE window. The
% immune-cell pipeline reads only the 640 window (piePulseIndex = 2) and throws
% the 485 photons away at read time, so nothing downstream holds this image and
% the raw PTU has to be re-read.
%
% WHY THE READER IS GATED RATHER THAN THE PHOTON LIST FILTERED
%
% The 485 window carries about 2.5% of the photons in these files. Reading the
% whole photon stream and discarding 97.5% of it afterwards would need several
% gigabytes of photon lists for a 3e8-record file. PTU_FLIM_GPU can gate in
% NATIVE TCSPC bins while reading, so only the wanted window is ever
% materialised - the same mechanism immune_cell_MIET uses to select the 640
% window.
%
% The window is chosen from the header (PIETimeGate paired with LaserWL), not
% by peak-finding: with 2.5% of the photons the 485 peak is not reliably the
% brightest feature, and a mis-identified window would silently return the
% MemGlow image instead of the CD58 one.
%
% cfg fields
%   excitationNm    485
%   imageSize       [] validates the result against the analysis grid
%   photonsPerChunk 5e6
%   maxNgate        4096, high enough for a 100 ns period at 40 ps
%   cacheFile       '' caches beside the PTU as immune_cell_MIET_cd58_image.mat
%   overwrite       false - reuse the cache when it is newer than the PTU
%
% out fields
%   intensity       nRow-by-nCol photon counts in the 485 window
%   photons, gate, imageSize, ptuFile, elapsedSeconds
%
% See also IMMUNE_CELL_MIET_TOPOGRAPHY_FIGURE, EXTRACT_FULL_PERIOD_TCSPC.

    if nargin < 2 || isempty(cfg); cfg = struct(); end
    cfg = withDefaults(cfg, struct( ...
        'excitationNm', 485, 'imageSize', [], 'photonsPerChunk', 5e6, ...
        'maxNgate', 4096, 'cacheFile', '', 'overwrite', false, ...
        'photonBlockSize', 2e6));
    ptuFile = char(ptuFile);
    if ~isfile(ptuFile)
        error('immune_cell_MIET_cd58_image:MissingFile', ...
            'PTU not found: %s', ptuFile);
    end
    if isempty(cfg.cacheFile)
        folder = fileparts(ptuFile);
        cfg.cacheFile = fullfile(folder, sprintf( ...
            'immune_cell_MIET_pie%dnm_image.mat', round(cfg.excitationNm)));
        % The 485 window was cached under an older name before this function
        % handled both PIE windows. Reuse it rather than spend 80 s per
        % acquisition re-reading the PTU for a file that already exists.
        legacy = fullfile(folder, 'immune_cell_MIET_cd58_image.mat');
        if round(cfg.excitationNm) == 485 && isfile(legacy) && ~isfile(cfg.cacheFile)
            cfg.cacheFile = legacy;
        end
    end
    if ~cfg.overwrite && isfile(cfg.cacheFile)
        cached = dir(cfg.cacheFile);
        source = dir(ptuFile);
        if cached.datenum >= source.datenum
            loaded = load(cfg.cacheFile, 'out');
            % Version 1 caches were written with rows and columns transposed.
            % Recompute rather than transpose in place, so the stored image can
            % never disagree with the code that produced it.
            if isfield(loaded, 'out') && isfield(loaded.out, 'formatVersion') ...
                    && loaded.out.formatVersion >= 2
                out = loaded.out;
                fprintf('immune_cell_MIET_cd58_image: reusing %s\n', cfg.cacheFile);
                return;
            end
        end
    end

    started = tic;
    head = PTU_Read_Head(ptuFile);
    gate = pieGateFromHeader(head, cfg.excitationNm);

    ptu = PTU_FLIM_GPU(ptuFile, struct( ...
        'photonsPerChunk', cfg.photonsPerChunk, 'computePerFrame', false, ...
        'storeTcspcPix', false, 'storePhotonLists', true, ...
        'storeTimeCell', false, 'showWaitbar', false, ...
        'maxNgate', cfg.maxNgate, 'computeGlobalMaps', false, ...
        'storePhotonSync', false, 'storePhotonFrame', false, ...
        'lowMemoryPhotonLists', true, 'photonBlockSize', cfg.photonBlockSize, ...
        'photonGateStartNative', gate.startNative, ...
        'photonGateLengthNative', gate.lengthNative));

    line = double(ptu.im_line(:));
    col = double(ptu.im_col(:));
    keep = isfinite(line) & isfinite(col) & line >= 1 & col >= 1;
    line = round(line(keep));
    col = round(col(keep));
    if isempty(line)
        error('immune_cell_MIET_cd58_image:NoPhotons', ...
            ['The %g nm window holds no photons with pixel coordinates in ' ...
             '%s.'], cfg.excitationNm, ptuFile);
    end

    % im_col indexes the analysis grid's ROWS and im_line its COLUMNS, not the
    % other way round. Verified against the pipeline's own intensity map on a
    % square acquisition: this assignment correlates 0.982, the transpose only
    % 0.814. A square scan hides the error - a roughly round, centred cell still
    % correlates well when transposed - so it only surfaced on the 558x560
    % acquisitions, where the photons fell outside the grid.
    imageSize = [max(col) max(line)];
    if ~isempty(cfg.imageSize)
        target = double(cfg.imageSize(:)).';
        if any(imageSize > target)
            error('immune_cell_MIET_cd58_image:GridMismatch', ...
                ['Photons land outside the analysis grid (%dx%d against ' ...
                 '%dx%d), so the CD58 image cannot be overlaid on the ' ...
                 'height map.'], imageSize(1), imageSize(2), target(1), target(2));
        end
        imageSize = target;
    end
    intensity = accumarray([col line], 1, imageSize);

    out = struct();
    out.formatVersion = 2;   % 1 had rows and columns transposed
    out.intensity = intensity;
    out.photons = numel(line);
    out.gate = gate;
    out.imageSize = imageSize;
    out.ptuFile = ptuFile;
    out.excitationNm = cfg.excitationNm;
    out.config = cfg;
    out.elapsedSeconds = toc(started);
    save(cfg.cacheFile, 'out', '-v7.3');
    fprintf(['immune_cell_MIET_cd58_image: %g nm window, %.3e photons over ' ...
        '%dx%d (%.1f per pixel) in %.0f s -> %s\n'], cfg.excitationNm, ...
        out.photons, imageSize(1), imageSize(2), ...
        out.photons / prod(imageSize), out.elapsedSeconds, cfg.cacheFile);
end

function gate = pieGateFromHeader(head, wantedNm)
% PIETimeGate holds each window's start in native TCSPC bins and LaserWL the
% wavelength per laser, index for index. Windows run from one start to the
% next, and the last one to the end of the period.
    if ~isfield(head, 'PIETimeGate') || isempty(head.PIETimeGate)
        error('immune_cell_MIET_cd58_image:NoPieHeader', ...
            'The PTU header has no PIETimeGate, so the 485 window is unknown.');
    end
    starts = double(head.PIETimeGate(:)).';
    windowCount = numel(starts);
    if isfield(head, 'PIENumPIEWindows') && ~isempty(head.PIENumPIEWindows)
        windowCount = min(windowCount, round(double(head.PIENumPIEWindows(1))));
    end
    starts = starts(1:windowCount);

    nativeDtNs = 1e9 * double(head.MeasDesc_Resolution);
    periodNs = 1e9 / double(head.TTResult_SyncRate);
    periodNative = round(periodNs / nativeDtNs);
    stops = [starts(2:end), periodNative];

    laserNm = nan(1, windowCount);
    if isfield(head, 'LaserWL') && ~isempty(head.LaserWL)
        available = double(head.LaserWL(:)).';
        laserNm(1:min(end, numel(available))) = ...
            available(1:min(windowCount, numel(available)));
    end
    if all(isnan(laserNm))
        index = 1;   % first window by convention
    else
        [~, index] = min(abs(laserNm - wantedNm));
    end

    gate = struct();
    gate.index = index;
    gate.windowCount = windowCount;
    gate.laserNm = laserNm;
    gate.selectedNm = laserNm(index);
    gate.startNative = starts(index);
    gate.lengthNative = stops(index) - starts(index);
    gate.startNs = starts(index) * nativeDtNs;
    gate.stopNs = stops(index) * nativeDtNs;
    gate.periodNs = periodNs;
    gate.nativeDtNs = nativeDtNs;
end

function cfg = withDefaults(cfg, defaults)
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(cfg, names{index}) || isempty(cfg.(names{index}))
            cfg.(names{index}) = defaults.(names{index});
        end
    end
end
