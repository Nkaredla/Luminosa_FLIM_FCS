function out = run_batch_immune_cell_MIET_profiles(dataRoot, cfg)
%RUN_BATCH_IMMUNE_CELL_MIET_PROFILES Cross-section profile figures for a session.
%
% out = run_batch_immune_cell_MIET_profiles()
% out = run_batch_immune_cell_MIET_profiles(dataRoot, cfg)
%
% Renders one twin-axis cross-section figure per cell per binning: MIET height
% and CD58-Atto488 intensity along a cut in x and a cut in y, both in the same
% axes. The CD58 image is taken from the cached 485 nm PIE extraction.
%
% cfg fields
%   variantTag   'h2_rnd_t220' - which calibration variant's height maps to use
%   binnings     {'native','sliding2x2','sliding4x4'}
%   figure       struct forwarded to immune_cell_MIET_profile_figure
%
% See also IMMUNE_CELL_MIET_PROFILE_FIGURE, RUN_BATCH_IMMUNE_CELL_MIET_TOPOGRAPHY.

    if nargin < 1 || isempty(dataRoot)
        dataRoot = 'D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_20260813_1';
    end
    if nargin < 2 || isempty(cfg); cfg = struct(); end
    if ~isfield(cfg, 'variantTag') || isempty(cfg.variantTag)
        cfg.variantTag = 'h2_rnd_t220';
    end
    if ~isfield(cfg, 'binnings') || isempty(cfg.binnings)
        cfg.binnings = {'native', 'sliding2x2', 'sliding4x4'};
    end
    cfg.binnings = cellstr(cfg.binnings);
    if ~isfield(cfg, 'figure'); cfg.figure = struct(); end

    pattern = sprintf('immune_cell_MIET_height_maps_%s.mat', cfg.variantTag);
    found = dir(fullfile(dataRoot, '**', pattern));
    found = found(~[found.isdir]);
    if isempty(found)
        error('run_batch_immune_cell_MIET_profiles:NoInputs', ...
            'No %s under %s', pattern, dataRoot);
    end
    [~, order] = sort({found.folder});
    found = found(order);

    fprintf('run_batch_immune_cell_MIET_profiles\n');
    fprintf('  session  %s\n  variant  %s\n  binnings %s\n', dataRoot, ...
        cfg.variantTag, strjoin(cfg.binnings, ', '));
    fprintf('  %d cell(s)\n\n', numel(found));

    out = struct('results', {{}}, 'skipped', {{}}, 'failures', {{}});
    for index = 1:numel(found)
        matFile = fullfile(found(index).folder, found(index).name);
        fprintf('[%d/%d] %s\n', index, numel(found), found(index).folder);
        for b = 1:numel(cfg.binnings)
            figureCfg = cfg.figure;
            figureCfg.binning = cfg.binnings{b};
            try
                out.results{end + 1} = ...
                    immune_cell_MIET_profile_figure(matFile, figureCfg);
            catch err
                % A binning with too few inverted pixels is expected on the
                % sparser 4 C acquisitions and must not stop the others.
                switch err.identifier
                    case {'immune_cell_MIET_profile_figure:StageUnavailable', ...
                          'immune_cell_MIET_profile_figure:TooFewPixels', ...
                          'immune_cell_MIET_profile_figure:NoBinning'}
                        fprintf('  %-11s skipped: %s\n', cfg.binnings{b}, err.message);
                        out.skipped{end + 1} = struct('file', matFile, ...
                            'binning', cfg.binnings{b}, 'message', err.message);
                    otherwise
                        fprintf(2, '  %-11s FAILED: %s\n', cfg.binnings{b}, err.message);
                        out.failures{end + 1} = struct('file', matFile, ...
                            'binning', cfg.binnings{b}, 'message', err.message);
                end
            end
        end
    end
    fprintf('\nrendered %d, skipped %d, failed %d\n', numel(out.results), ...
        numel(out.skipped), numel(out.failures));
end
