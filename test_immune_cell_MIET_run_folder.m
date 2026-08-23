function test_immune_cell_MIET_run_folder()
%TEST_IMMUNE_CELL_MIET_RUN_FOLDER Contract test for the results-run naming.
%
% test_immune_cell_MIET_run_folder()
%
% The run folder exists so that re-analysing with different settings cannot
% overwrite an earlier analysis. That only works if the name satisfies four
% properties at once, and getting any one of them wrong fails silently - either
% by orphaning hours of finished work, or by overwriting it. Hence a test.
%
%   RESUMES     identical settings must give an identical name
%   PROTECTS    any change that changes the numbers must give a different name
%   IGNORES     changes to output paths and figure flags must NOT change it
%   STABLE      struct field order must not change it

    total = 0; failed = 0;
    fprintf('test_immune_cell_MIET_run_folder\n');

    base = struct( ...
        'versionResults', true, 'runName', '', ...
        'saveTcspcPix', true, 'scanPlanes', {{'XY'}}, ...
        'minDwellMs', 0.4, 'minCompletionFraction', 0.95, ...
        'minRecords', 1e6, 'requiredPieWindows', 2, ...
        'excitationPulseIndex', 2, 'requiredExcitationNm', 640, ...
        'wavelengthToleranceNm', 10, ...
        'pipeline', struct('outputDir', 'D:\\somewhere', 'showFigure', true, ...
            'saveOutputs', true, ...
            'bayes', struct('membraneTauCount', 48, 'fractionStep', 0.2, ...
                'slbCountRelTol', 0.0025, 'slbCountPriorNodes', 5)), ...
        'qc', struct('maxSlbClippedFraction', 0.10, 'minBayesCoverage', 0.25));

    reference = immune_cell_MIET_run_folder(base);
    [total, failed] = check(total, failed, 'name has the cfg_ prefix', ...
        startsWith(reference, 'cfg_') && numel(reference) == 12, reference);

    % ---- RESUMES ----------------------------------------------------------
    [total, failed] = check(total, failed, ...
        'identical settings give an identical name', ...
        strcmp(immune_cell_MIET_run_folder(base), reference), reference);

    % ---- IGNORES ----------------------------------------------------------
    cosmetic = base;
    cosmetic.pipeline.outputDir = 'E:\\a\\completely\\different\\path';
    cosmetic.pipeline.showFigure = false;
    cosmetic.pipeline.saveOutputs = false;
    [total, failed] = check(total, failed, ...
        'output path and figure flags do not change the name', ...
        strcmp(immune_cell_MIET_run_folder(cosmetic), reference), ...
        immune_cell_MIET_run_folder(cosmetic));

    % ---- STABLE -----------------------------------------------------------
    permuted = base;
    permuted.pipeline.bayes = orderfields(base.pipeline.bayes, ...
        {'slbCountPriorNodes', 'fractionStep', 'slbCountRelTol', ...
         'membraneTauCount'});
    [total, failed] = check(total, failed, ...
        'permuting struct field order does not change the name', ...
        strcmp(immune_cell_MIET_run_folder(permuted), reference), ...
        immune_cell_MIET_run_folder(permuted));

    % ---- PROTECTS ---------------------------------------------------------
    variants = { ...
        'membraneTauCount', @(c) setBayes(c, 'membraneTauCount', 10); ...
        'fractionStep',     @(c) setBayes(c, 'fractionStep', 0.1); ...
        'slbCountRelTol',   @(c) setBayes(c, 'slbCountRelTol', 0.02); ...
        'qc threshold',     @(c) setQc(c, 'maxSlbClippedFraction', 0.2); ...
        'minDwellMs',       @(c) setField(c, 'minDwellMs', 0.8); ...
        'saveTcspcPix',     @(c) setField(c, 'saveTcspcPix', false); ...
        'scanPlanes',       @(c) setField(c, 'scanPlanes', {{'XZ'}})};
    names = containers.Map('KeyType', 'char', 'ValueType', 'char');
    names(reference) = 'baseline';
    for k = 1:size(variants, 1)
        candidate = immune_cell_MIET_run_folder(variants{k, 2}(base));
        distinct = ~strcmp(candidate, reference) && ~isKey(names, candidate);
        [total, failed] = check(total, failed, ...
            sprintf('changing %s gives a new name', variants{k, 1}), ...
            distinct, candidate);
        if ~isKey(names, candidate)
            names(candidate) = variants{k, 1};
        end
    end

    % ---- explicit and disabled modes --------------------------------------
    explicit = base; explicit.runName = 'grid48_rerun';
    [total, failed] = check(total, failed, 'an explicit runName is used verbatim', ...
        strcmp(immune_cell_MIET_run_folder(explicit), 'grid48_rerun'), ...
        immune_cell_MIET_run_folder(explicit));

    flat = base; flat.versionResults = false;
    [total, failed] = check(total, failed, ...
        'versionResults=false gives the flat legacy layout', ...
        isempty(immune_cell_MIET_run_folder(flat)), '<empty>');

    bad = base; bad.runName = 'has\\a\\separator';
    threw = false;
    try
        immune_cell_MIET_run_folder(bad);
    catch err
        threw = strcmp(err.identifier, 'batch_immune_cell_MIET:RunName');
    end
    [total, failed] = check(total, failed, ...
        'a runName containing a path separator is rejected', threw, '');

    fprintf('\n  %d/%d checks passed\n', total - failed, total);
    if failed > 0
        error('test_immune_cell_MIET_run_folder:Failed', ...
            '%d check(s) failed.', failed);
    end
    fprintf('test_immune_cell_MIET_run_folder: PASS\n');
end

function [total, failed] = check(total, failed, label, condition, detail)
    total = total + 1;
    if condition
        fprintf('  ok    %-52s %s\n', label, detail);
    else
        failed = failed + 1;
        fprintf('  FAIL  %-52s %s\n', label, detail);
    end
end

function cfg = setBayes(cfg, field, value)
    cfg.pipeline.bayes.(field) = value;
end

function cfg = setQc(cfg, field, value)
    cfg.qc.(field) = value;
end

function cfg = setField(cfg, field, value)
    if iscell(value) && isscalar(value) && iscell(value{1})
        cfg.(field) = value{1};
    else
        cfg.(field) = value;
    end
end
