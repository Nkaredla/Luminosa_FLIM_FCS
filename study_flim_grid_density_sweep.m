function sweep = study_flim_grid_density_sweep(opts)
%STUDY_FLIM_GRID_DENSITY_SWEEP How fine must the lifetime grid be?
%
% sweep = study_flim_grid_density_sweep()
% sweep = study_flim_grid_density_sweep(opts)
%
% Finds the smallest membraneTauCount that keeps the third-component
% false-positive rate under control across a range of photon counts, and
% prints a recommended value.
%
% WHY
%
% Calibration showed that with the production grid (10 points over 0.62-5 ns)
% the false three-component rate at a fixed 0.95 threshold climbs from 0.3% to
% 63% as photons rise from 500 to 64000. Snapping the simulated truth onto the
% fit grid removes that entirely - false positives stay at 0.000-0.004 across
% the whole range, while detection power is unchanged. So the drift is caused
% by lifetime discretisation: an off-grid single exponential is approximated
% better by two on-grid exponentials than by one, and enough photons make that
% improvement decisive.
%
% Real lifetimes are continuous, so they are always off-grid by up to half a
% step. The question is therefore not "can we sit on the grid" but "how fine
% must the grid be for the residual off-grid error to stay below photon
% noise". Shape error scales with the grid step and the KL gap with its
% square, while the evidence for an extra component scales with N, so the
% required step should fall as N^(-1/2): roughly 3x finer grid per 10x
% photons. This measures that directly instead of trusting the argument.
%
% WHAT IT DOES NOT FIX
%
% A finer grid removes discretisation only. The other mis-specification
% sources are untouched: the hard SLB amplitude (see slbCountPriorNodes),
% continuously distributed MIET lifetimes, and spatial mixing of different
% lifetimes inside a sliding window. The last two mean the true decay is not a
% sum of three exponentials at all, which no grid can repair.
%
% opts fields (all optional)
%   gridCounts          membraneTauCount values, default [10 16 24 32 48]
%   photonTotals        default [2000 8000 32000]
%   pixelsPerCondition  default 500
%   membraneTauBoundsNs default [0.4 5.5]. Note the pipeline default lower
%                       bound is max(1.15*tauSlb, tauSlb+2*dt, 0.05) = 0.62 ns
%                       for tauSlb = 0.3 and dt = 0.16, which EXCLUDES any
%                       membrane component below 0.62 ns.
%   maxFalsePositive    acceptance criterion, default 0.05
%   referenceThreshold  fixed threshold under test, default 0.95
%   outputDir           default pwd
%   seed                default 11

    if nargin < 1 || isempty(opts)
        opts = struct();
    end
    defaults = struct( ...
        'gridCounts', [10 16 24 32 48], ...
        'photonTotals', [2000 8000 32000], ...
        'pixelsPerCondition', 500, ...
        'membraneTauBoundsNs', [0.4 5.5], ...
        'maxFalsePositive', 0.05, ...
        'referenceThreshold', 0.95, ...
        'outputDir', pwd, ...
        'seed', 11);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if ~isfolder(opts.outputDir)
        mkdir(opts.outputDir);
    end

    fprintf('\nstudy_flim_grid_density_sweep\n');
    fprintf('  grid counts  : %s\n', mat2str(opts.gridCounts));
    fprintf('  photons      : %s\n', mat2str(opts.photonTotals));
    fprintf('  criterion    : false positives <= %.3g at fixed threshold %.2f\n\n', ...
        opts.maxFalsePositive, opts.referenceThreshold);

    rows = struct([]);
    for k = 1:numel(opts.gridCounts)
        gridCount = opts.gridCounts(k);
        inner = struct( ...
            'membraneTauCount', gridCount, ...
            'membraneTauBoundsNs', opts.membraneTauBoundsNs, ...
            'photonTotals', opts.photonTotals, ...
            'pixelsPerCondition', opts.pixelsPerCondition, ...
            'referenceThreshold', opts.referenceThreshold, ...
            'placeTruthOnGrid', false, ...
            'seed', opts.seed, ...
            'outputDir', fullfile(opts.outputDir, sprintf('grid%02d', gridCount)));
        started = tic;
        calibration = study_flim_component_thresholds(inner);
        elapsed = toc(started);

        summary = calibration.summary;
        for r = 1:height(summary)
            entry = struct( ...
                'membraneTauCount', gridCount, ...
                'photonTotal', summary.photonTotal(r), ...
                'falsePositiveAtFixed', summary.falsePositiveAtReference(r), ...
                'powerAtFixed', summary.powerAtReference(r), ...
                'calibratedThreshold', summary.calibratedThreshold(r), ...
                'thresholdSaturated', summary.thresholdSaturated(r), ...
                'nullMedianQ3', summary.nullMedianQ3(r), ...
                'secondsForGrid', elapsed);
            if isempty(rows)
                rows = entry;
            else
                rows(end+1) = entry; %#ok<AGROW>
            end
        end

        worst = max(summary.falsePositiveAtReference);
        fprintf(['  grid %2d : worst false-positive %.4f across photon ' ...
            'range (%.0f s)\n'], gridCount, worst, elapsed);
    end

    sweepTable = struct2table(rows);
    sweep = struct('summary', sweepTable, 'opts', opts);
    csvFile = fullfile(opts.outputDir, 'flim_grid_density_sweep.csv');
    writetable(sweepTable, csvFile);

    fprintf('\n  grid   %s\n', sprintf('  N=%-7d', opts.photonTotals));
    for k = 1:numel(opts.gridCounts)
        gridCount = opts.gridCounts(k);
        fprintf('  %4d ', gridCount);
        for t = 1:numel(opts.photonTotals)
            selected = sweepTable.membraneTauCount == gridCount & ...
                sweepTable.photonTotal == opts.photonTotals(t);
            fprintf('  FP=%.4f ', sweepTable.falsePositiveAtFixed(selected));
        end
        fprintf('\n');
    end

    % Smallest grid whose worst false-positive rate meets the criterion.
    recommended = NaN;
    for k = 1:numel(opts.gridCounts)
        gridCount = opts.gridCounts(k);
        selected = sweepTable.membraneTauCount == gridCount;
        if max(sweepTable.falsePositiveAtFixed(selected)) <= opts.maxFalsePositive
            recommended = gridCount;
            break;
        end
    end

    fprintf('\n');
    if isnan(recommended)
        fprintf(['  NO tested grid met the criterion. Either extend ' ...
            'gridCounts upward, or accept\n  that discretisation is not the ' ...
            'only mis-specification source here and treat the\n  remaining ' ...
            'false-positive rate as irreducible for this model family.\n']);
    else
        fprintf('  RECOMMENDED membraneTauCount = %d\n\n', recommended);
        fprintf('  pipelineCfg.bayes = struct( ...\n');
        fprintf('      ''batchSize'', 2048, ''includeBackground'', true, ...\n');
        fprintf('      ''signalGrid'', [0.25 0.5 0.75 1], ...\n');
        fprintf('      ''membraneTauCount'', %d, ...\n', recommended);
        fprintf('      ''membraneTauBoundsNs'', [%g %g], ...\n', ...
            opts.membraneTauBoundsNs(1), opts.membraneTauBoundsNs(2));
        fprintf('      ''fractionStep'', 0.2, ...\n');
        fprintf('      ''minimumMembraneFraction'', 0.1, ...\n');
        fprintf('      ''slbCountRelTol'', 0.0025);\n');
    end
    fprintf('\n  wrote %s\n', csvFile);
end
