function runName = immune_cell_MIET_run_folder(cfg)
%IMMUNE_CELL_MIET_RUN_FOLDER Name the results subfolder for one analysis run.
%
% runName = immune_cell_MIET_run_folder(cfg)
%
% Returns the subfolder name that batch_immune_cell_MIET creates inside each
% acquisition's results folder, so that a run with different settings cannot
% overwrite an earlier one.
%
% cfg fields used
%   versionResults   false -> '' (results sit directly in the results folder,
%                    the historical layout, and a rerun overwrites)
%   runName          non-empty -> used as given
%   pipeline, qc and the prefilter thresholds -> hashed when runName is empty
%
% Kept as its own function, rather than as a subfunction of the batch driver,
% so test_immune_cell_MIET_run_folder can exercise it without running an
% analysis. The properties it has to satisfy are easy to state and easy to get
% wrong, which is exactly what wants a test:
%
%   unchanged settings                  -> identical name (so results RESUME)
%   changed pipeline or QC settings     -> different name (so nothing is lost)
%   changed output path or figure flags -> identical name (cosmetic only)
%   permuted struct field order         -> identical name

    if ~cfg.versionResults
        runName = '';
        return;
    end
    runName = strtrim(char(cfg.runName));
    if ~isempty(runName)
        if ~isempty(regexp(runName, '[\\/:*?"<>|]', 'once'))
            error('batch_immune_cell_MIET:RunName', ...
                'cfg.runName must be a single folder name, got: %s', runName);
        end
        return;
    end
    runName = ['cfg_' configurationHash(cfg)];
end

function digest = configurationHash(cfg)
%CONFIGURATIONHASH Stable 8-hex-digit digest of the settings that affect results.
% FNV-1a over a canonical text rendering. Implemented locally rather than via
% Java or a toolbox so it works in -batch runs with no display and no extra
% licences, and so the same input always gives the same digest across MATLAB
% versions.
    scoring = struct();
    scoring.pipeline = stripNonNumericFields(cfg.pipeline);
    scoring.qc = cfg.qc;
    scoring.saveTcspcPix = cfg.saveTcspcPix;
    scoring.scanPlanes = cfg.scanPlanes;
    scoring.minDwellMs = cfg.minDwellMs;
    scoring.minCompletionFraction = cfg.minCompletionFraction;
    scoring.minRecords = cfg.minRecords;
    scoring.requiredPieWindows = cfg.requiredPieWindows;
    scoring.excitationPulseIndex = cfg.excitationPulseIndex;
    scoring.requiredExcitationNm = cfg.requiredExcitationNm;
    scoring.wavelengthToleranceNm = cfg.wavelengthToleranceNm;

    text = canonicalText(scoring, '');
    hash = 2166136261;
    bytes = double(unicode2native(text, 'UTF-8'));
    for k = 1:numel(bytes)
        hash = double(bitxor(uint32(hash), uint32(bytes(k))));
        hash = fnvMultiply(hash);
    end
    digest = lower(dec2hex(hash, 8));
end

function product = fnvMultiply(hash)
%FNVMULTIPLY hash * 16777619 mod 2^32, computed exactly in double precision.
% The direct product overflows the exact-integer range of a double: the hash
% reaches 2^32 and the FNV prime is 2^24-ish, so their product approaches
% 7.2e16 while doubles are only exact to 2^53 = 9.0e15. Multiplying the two
% 16-bit halves separately keeps every intermediate under 1.1e12, which is
% exact, and the high half's carry beyond 2^32 is discarded by the modulus
% anyway.
    prime = 16777619;
    modulus = 4294967296;              % 2^32
    low = mod(hash, 65536);
    high = floor(hash / 65536);
    product = mod(mod(low * prime, modulus) + ...
        mod(mod(high * prime, 65536) * 65536, modulus), modulus);
end

function value = stripNonNumericFields(value)
% Remove settings that change presentation or destinations but not results.
    cosmetic = {'outputDir', 'showFigure', 'showWaitbar', 'saveOutputs', ...
        'figureVisible', 'verbose', 'saveTcspcPix'};
    if ~isstruct(value) || ~isscalar(value)
        return;
    end
    for k = 1:numel(cosmetic)
        if isfield(value, cosmetic{k})
            value = rmfield(value, cosmetic{k});
        end
    end
end

function text = canonicalText(value, prefix)
% Field names sorted so that struct field ORDER cannot change the digest.
    if ~isstruct(value) || ~isscalar(value)
        text = sprintf('%s=%s;', prefix, compactValue(value));
        return;
    end
    names = sort(fieldnames(value));
    parts = cell(numel(names), 1);
    for k = 1:numel(names)
        if isempty(prefix)
            label = names{k};
        else
            label = [prefix '.' names{k}];
        end
        parts{k} = canonicalText(value.(names{k}), label);
    end
    text = strjoin(parts, '');
end

function text = compactValue(value)
% Local copy: the batch driver has its own for the run manifest. Duplicating a
% ten-line formatter is preferable to either file depending on the other's
% private helpers.
    if ischar(value)
        text = value;
    elseif isstring(value)
        text = char(strjoin(value, ', '));
    elseif iscell(value)
        parts = cellfun(@compactValue, value, 'UniformOutput', false);
        text = ['{' strjoin(parts, ', ') '}'];
    elseif islogical(value) || isnumeric(value)
        if isempty(value)
            text = '[]';
        else
            text = mat2str(value, 6);
        end
    else
        text = ['<' class(value) '>'];
    end
end
