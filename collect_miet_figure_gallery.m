function out = collect_miet_figure_gallery(kind, cfg)
%COLLECT_MIET_FIGURE_GALLERY Gather per-cell MIET figures into one flat folder.
%
% out = collect_miet_figure_gallery('contacts')
% out = collect_miet_figure_gallery('profiles')
% out = collect_miet_figure_gallery('topography')
% out = collect_miet_figure_gallery(kind, cfg)
%
% The per-cell figures are written deep inside each acquisition, beside the
% height maps that produced them, which is right for provenance and useless for
% browsing: comparing seventeen cells means opening seventeen folders. This
% copies them into one directory per figure kind, renamed so the session, the
% acquisition and the binning are all visible in the file name and the whole set
% sorts sensibly.
%
% WHY COPY RATHER THAN MOVE OR LINK
%
% The originals stay where the pipeline put them, next to the .mat they came
% from, so nothing downstream breaks and a stale gallery can always be deleted
% and rebuilt. The gallery is derived data - treat it as disposable.
%
% NAMING
%
%   <sessionTag>_<acquisition>_<binning>[_<colormap>].png
%
% e.g. 4degC_20260813-133830_sliding4x4.png. This is the convention the
% topography gallery already used; contacts and cross-sections now match it, so
% the three folders sort the same way and a cell's three views line up by name.
%
% cfg fields
%   sessions    n-by-2 cellstr, {sessionFolder, sessionTag}. Defaults to the
%               260813 pair.
%   outputDir   '' places the gallery beside the sessions, named after the kind
%   overwrite   true - copy even when a same-named file is already there
%   dryRun      false - list what would be copied without touching the disk
%
% See also RUN_BATCH_IMMUNE_CELL_MIET_CONTACTS, RUN_BATCH_IMMUNE_CELL_MIET_PROFILES.

    if nargin < 1 || isempty(kind); kind = 'contacts'; end
    if nargin < 2 || isempty(cfg); cfg = struct(); end

    switch lower(char(kind))
        case 'contacts'
            prefix = 'immune_cell_MIET_contacts_';
            folderName = 'MIET_contact_figures';
        case {'profiles', 'crosssections', 'cross_sections'}
            prefix = 'immune_cell_MIET_profiles_';
            folderName = 'MIET_crosssection_figures';
        case 'topography'
            prefix = 'immune_cell_MIET_topography_';
            folderName = 'MIET_topography_figures';
        otherwise
            error('collect_miet_figure_gallery:UnknownKind', ...
                'kind must be contacts, profiles or topography, not %s', kind);
    end

    if ~isfield(cfg, 'sessions') || isempty(cfg.sessions)
        cfg.sessions = { ...
            'D:\Luminosa\Data\260813\4deg_Jurkat_CD58_memglow',            '4degC'; ...
            'D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_20260813_1',   'RT'};
    end
    if ~isfield(cfg, 'overwrite') || isempty(cfg.overwrite); cfg.overwrite = true; end
    if ~isfield(cfg, 'dryRun') || isempty(cfg.dryRun); cfg.dryRun = false; end
    if ~isfield(cfg, 'outputDir') || isempty(cfg.outputDir)
        parent = fileparts(cfg.sessions{1, 1});
        cfg.outputDir = fullfile(parent, folderName);
    end
    if ~cfg.dryRun && ~isfolder(cfg.outputDir); mkdir(cfg.outputDir); end

    fprintf('collect_miet_figure_gallery  (%s)\n', kind);
    fprintf('  into  %s\n', cfg.outputDir);

    out = struct('copied', {{}}, 'skipped', {{}}, 'outputDir', cfg.outputDir, ...
        'kind', char(kind));
    for s = 1:size(cfg.sessions, 1)
        root = cfg.sessions{s, 1};
        tag = cfg.sessions{s, 2};
        if ~isfolder(root)
            warning('collect_miet_figure_gallery:NoSession', ...
                'Session folder not found: %s', root);
            continue;
        end
        found = dir(fullfile(root, '**', [prefix '*.png']));
        found = found(~[found.isdir]);
        fprintf('  %-6s %3d file(s) under %s\n', tag, numel(found), root);
        for k = 1:numel(found)
            source = fullfile(found(k).folder, found(k).name);
            stem = regexprep(found(k).name, ['^' regexptranslate('escape', prefix)], '');
            target = fullfile(cfg.outputDir, sprintf('%s_%s', tag, stem));
            if ~cfg.overwrite && isfile(target)
                out.skipped{end + 1} = target; %#ok<AGROW>
                continue;
            end
            if cfg.dryRun
                out.copied{end + 1} = target; %#ok<AGROW>
                continue;
            end
            [ok, message] = copyfile(source, target, 'f');
            if ok
                out.copied{end + 1} = target; %#ok<AGROW>
            else
                warning('collect_miet_figure_gallery:CopyFailed', ...
                    'Could not copy %s (%s)', source, message);
                out.skipped{end + 1} = source; %#ok<AGROW>
            end
        end
    end
    fprintf('  copied %d, skipped %d\n', numel(out.copied), numel(out.skipped));
end
