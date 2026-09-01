function out = run_batch_immune_cell_MIET_contacts(dataRoot, cfg)
%RUN_BATCH_IMMUNE_CELL_MIET_CONTACTS Contact-region analysis for a session.
%
% out = run_batch_immune_cell_MIET_contacts()
% out = run_batch_immune_cell_MIET_contacts(dataRoot, cfg)
%
% Runs immune_cell_MIET_contact_regions over every cell in a session and writes
% one summary CSV per session alongside the per-cell figures.
%
% WHY THE SUMMARY TABLE MATTERS MORE THAN THE INDIVIDUAL FIGURES
%
% Each cell's randomisation gives a p-value for that cell alone, and with a
% handful of contact regions per cell none of them will be decisive. The
% question - whether membrane contacts and CD58 clusters coincide - is answered
% across cells, by asking whether the enrichment ratios sit above 1 more often
% than the nulls do. The CSV is what supports that; the figures are for looking
% at individual cells afterwards.
%
% The same fixed height threshold is used for every cell, so the "0 regions"
% rows are meaningful data and are written out rather than skipped.
%
% cfg fields
%   variantTag  'h2_rnd_t220' - which calibration variant's height maps to use
%   binnings    {'sliding4x4'} - contacts need the reliable pixels; native and
%               2x2 invert too few to define an extended region on the sparser
%               acquisitions. Pass more if you want the comparison.
%   summaryFile '' writes <dataRoot>/immune_cell_MIET_contacts_summary.csv
%   contacts    struct forwarded to immune_cell_MIET_contact_regions
%
% See also IMMUNE_CELL_MIET_CONTACT_REGIONS, RUN_BATCH_IMMUNE_CELL_MIET_PROFILES.

    if nargin < 1 || isempty(dataRoot)
        dataRoot = 'D:\Luminosa\Data\260813\4deg_Jurkat_CD58_memglow';
    end
    if nargin < 2 || isempty(cfg); cfg = struct(); end
    if ~isfield(cfg, 'variantTag') || isempty(cfg.variantTag)
        cfg.variantTag = 'h2_rnd_t220';
    end
    if ~isfield(cfg, 'binnings') || isempty(cfg.binnings)
        cfg.binnings = {'sliding4x4'};
    end
    cfg.binnings = cellstr(cfg.binnings);
    if ~isfield(cfg, 'contacts'); cfg.contacts = struct(); end
    if ~isfield(cfg, 'summaryFile') || isempty(cfg.summaryFile)
        cfg.summaryFile = fullfile(dataRoot, ...
            sprintf('immune_cell_MIET_contacts_summary_%s.csv', cfg.variantTag));
    end

    pattern = sprintf('immune_cell_MIET_height_maps_%s.mat', cfg.variantTag);
    found = dir(fullfile(dataRoot, '**', pattern));
    found = found(~[found.isdir]);
    if isempty(found)
        error('run_batch_immune_cell_MIET_contacts:NoInputs', ...
            'No %s under %s', pattern, dataRoot);
    end
    [~, order] = sort({found.folder});
    found = found(order);

    fprintf('run_batch_immune_cell_MIET_contacts\n');
    fprintf('  session  %s\n  variant  %s\n  binnings %s\n', dataRoot, ...
        cfg.variantTag, strjoin(cfg.binnings, ', '));
    fprintf('  %d cell(s)\n\n', numel(found));

    out = struct('results', {{}}, 'skipped', {{}}, 'failures', {{}});
    rows = cell(0, 1);
    for index = 1:numel(found)
        matFile = fullfile(found(index).folder, found(index).name);
        fprintf('[%d/%d] %s\n', index, numel(found), found(index).folder);
        for b = 1:numel(cfg.binnings)
            contactCfg = cfg.contacts;
            contactCfg.binning = cfg.binnings{b};
            try
                result = immune_cell_MIET_contact_regions(matFile, contactCfg);
                out.results{end + 1} = result;
                rows{end + 1, 1} = summaryRow(result); %#ok<AGROW>
            catch err
                switch err.identifier
                    case {'immune_cell_MIET_contact_regions:StageUnavailable', ...
                          'immune_cell_MIET_contact_regions:TooFewPixels', ...
                          'immune_cell_MIET_contact_regions:NoBinning'}
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

    out.summaryFile = '';
    if ~isempty(rows)
        writeSummary(cfg.summaryFile, rows);
        out.summaryFile = cfg.summaryFile;
        fprintf('\nsummary -> %s\n', cfg.summaryFile);
    end
    fprintf('analysed %d, skipped %d, failed %d\n', numel(out.results), ...
        numel(out.skipped), numel(out.failures));
end

function row = summaryRow(result)
    row = struct();
    row.acquisition = regexprep(result.acquisition, '^_', '');
    row.binning = result.binning;
    row.thresholdNm = result.contactThresholdNm;
    row.footprintUm2 = result.footprintAreaUm2;
    row.measuredFraction = result.measuredFraction;
    row.nRegions = result.nRegions;
    row.contactUm2 = result.contactAreaUm2;
    row.contactFraction = result.contactAreaFraction;
    row.nClusters = result.nClusters;
    row.clustersInContact = result.stats.clustersInContact;
    row.cd58MeanInside = result.stats.meanInside;
    row.cd58MeanOutside = result.stats.meanOutside;
    row.cd58FootprintMean = result.stats.footprintMean;
    row.enrichment = result.stats.enrichment;
    row.pEnriched = result.p.enrichedGreater;
    row.pDepleted = result.p.enrichedLess;
    row.nullMedianEnrichment = medianOf(result.null.enrichment);
    row.meanRatioInOut = result.stats.meanRatio;
    row.nullMaskAreaRatio = fieldOr(result.null, 'maskPixelRatio', NaN);
    row.nullTrials = result.null.trials;
    row.nullInsideFraction = fieldOr(result.null, 'meanInsideFraction', NaN);
end

function writeSummary(file, rows)
% The table is assembled by hand rather than through struct2table so the column
% order is fixed and a rerun produces a byte-identical header.
    names = fieldnames(rows{1});
    fid = fopen(file, 'w');
    if fid < 0
        warning('run_batch_immune_cell_MIET_contacts:NoWrite', ...
            'Could not open %s for writing.', file);
        return;
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', strjoin(names(:).', ','));
    for index = 1:numel(rows)
        values = cell(1, numel(names));
        for k = 1:numel(names)
            value = rows{index}.(names{k});
            if ischar(value)
                values{k} = value;
            elseif ~isfinite(value)
                values{k} = 'NaN';
            else
                values{k} = sprintf('%.6g', value);
            end
        end
        fprintf(fid, '%s\n', strjoin(values, ','));
    end
end

function value = medianOf(values)
    values = values(isfinite(values));
    if isempty(values); value = NaN; else; value = median(values); end
end

function value = fieldOr(structure, name, fallback)
    value = fallback;
    if isstruct(structure) && isfield(structure, name)
        value = structure.(name);
    end
end
