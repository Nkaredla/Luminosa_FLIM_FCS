function out = run_batch_immune_cell_MIET_topography(dataRoot, cfg)
%RUN_BATCH_IMMUNE_CELL_MIET_TOPOGRAPHY Topography figures for a whole session.
%
% out = run_batch_immune_cell_MIET_topography()
% out = run_batch_immune_cell_MIET_topography(dataRoot, cfg)
%
% Renders one 3D topography figure per cell, per binning, per colormap, each
% with the MIET height projection and the CD58-Atto488 bilayer image as flat
% layers underneath. The PIE-window images are extracted from the raw PTU on
% first use and cached beside the acquisition, so a second pass over the same
% session only re-renders.
%
% WHY ALL THREE BINNINGS GET THE SAME TREATMENT
%
% Native, 2x2 and 4x4 sliding results are rendered with identical smoothing,
% downsampling and colour limits. That makes the three directly comparable, and
% the comparison is the point: 4x4 buys reliable pixels at the cost of lateral
% detail, native keeps the detail but inverts far fewer pixels. Tuning the
% rendering per binning would hide exactly that trade-off. Expect the native
% figures to look noisier and to carry a larger interpolated fraction - both
% are real, and both are stated on each figure.
%
% cfg fields
%   variantTag   'h2_rnd_t220' - which calibration variant's height maps to use
%   binnings     {'native','sliding2x2','sliding4x4'}
%   colormaps    {'turbo','parula','copper'}
%   figure       struct forwarded to immune_cell_MIET_topography_figure
%
% See also IMMUNE_CELL_MIET_TOPOGRAPHY_FIGURE, IMMUNE_CELL_MIET_CD58_IMAGE.

    if nargin < 1 || isempty(dataRoot)
        dataRoot = 'D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_20260813_1';
    end
    if nargin < 2 || isempty(cfg); cfg = struct(); end
    if ~isfield(cfg, 'variantTag') || isempty(cfg.variantTag)
        cfg.variantTag = 'h2_rnd_t220';
    end
    % cfg.binning stays accepted as a single-binning alias for older callers.
    if ~isfield(cfg, 'binnings') || isempty(cfg.binnings)
        if isfield(cfg, 'binning') && ~isempty(cfg.binning)
            cfg.binnings = cellstr(cfg.binning);
        else
            cfg.binnings = {'native', 'sliding2x2', 'sliding4x4'};
        end
    end
    cfg.binnings = cellstr(cfg.binnings);
    if ~isfield(cfg, 'colormaps') || isempty(cfg.colormaps)
        cfg.colormaps = {'turbo', 'parula', 'copper'};
    end
    if ~isfield(cfg, 'figure'); cfg.figure = struct(); end

    pattern = sprintf('immune_cell_MIET_height_maps_%s.mat', cfg.variantTag);
    found = dir(fullfile(dataRoot, '**', pattern));
    found = found(~[found.isdir]);
    if isempty(found)
        error('run_batch_immune_cell_MIET_topography:NoInputs', ...
            'No %s under %s', pattern, dataRoot);
    end
    [~, order] = sort({found.folder});
    found = found(order);

    fprintf('run_batch_immune_cell_MIET_topography\n');
    fprintf('  session   %s\n', dataRoot);
    fprintf('  variant   %s\n', cfg.variantTag);
    fprintf('  binnings  %s\n', strjoin(cfg.binnings, ', '));
    fprintf('  %d cell(s) x %d binning(s) x %d colormap(s)\n\n', ...
        numel(found), numel(cfg.binnings), numel(cfg.colormaps));

    out = struct('results', {{}}, 'failures', {{}}, 'skipped', {{}});
    for index = 1:numel(found)
        matFile = fullfile(found(index).folder, found(index).name);
        fprintf('[%d/%d] %s\n', index, numel(found), found(index).folder);
        for b = 1:numel(cfg.binnings)
            binning = cfg.binnings{b};
            figureCfg = cfg.figure;
            figureCfg.binning = binning;
            figureCfg.colormaps = cfg.colormaps;
            try
                out.results{end + 1} = ...
                    immune_cell_MIET_topography_figure(matFile, figureCfg);
            catch err
                % An unavailable or too-sparse binning is expected for some
                % acquisitions, especially native at 4 C, and must not stop the
                % other binnings of the same cell.
                switch err.identifier
                    case {'immune_cell_MIET_topography_figure:StageUnavailable', ...
                          'immune_cell_MIET_topography_figure:TooFewPixels', ...
                          'immune_cell_MIET_topography_figure:NoBinning'}
                        fprintf('  %-11s skipped: %s\n', binning, err.message);
                        out.skipped{end + 1} = struct('file', matFile, ...
                            'binning', binning, 'message', err.message);
                    otherwise
                        fprintf(2, '  %-11s FAILED: %s\n', binning, err.message);
                        out.failures{end + 1} = struct('file', matFile, ...
                            'binning', binning, 'message', err.message);
                end
            end
        end
    end
    fprintf('\nrendered %d figure set(s), %d skipped, %d failure(s)\n', ...
        numel(out.results), numel(out.skipped), numel(out.failures));
end
