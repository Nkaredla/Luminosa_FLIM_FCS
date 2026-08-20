function calibration = study_flim_component_thresholds(opts)
%STUDY_FLIM_COMPONENT_THRESHOLDS Calibrate the third-component decision rule.
%
% calibration = study_flim_component_thresholds()
% calibration = study_flim_component_thresholds(opts)
%
% Produces a photon-count-dependent threshold for declaring a third lifetime
% component, with a stated false-positive rate.
%
% WHY THIS IS NEEDED
%
% The pipeline gates the third component on a fixed posterior threshold,
% cfg.componentMaps.posteriorThreshold(2), applied identically to native
% pixels and to 4x4 windows carrying ~16x more photons. That is only
% defensible if the posterior probability is a consistent estimator of
% component count - which holds in the M-closed case, where the truth is one
% of the candidate models.
%
% This analysis is M-open: the truth is outside the candidate set, because the
% lifetime grid is discrete (an off-grid lifetime is not representable), the
% SLB amplitude is a hard constraint, MIET lifetimes are continuously
% distributed, and a sliding window mixes lifetimes spatially. Under M-open
% the standard result is that the posterior concentrates on whichever
% candidate is closest in Kullback-Leibler divergence, not on the true
% component count - so the evidence for the richer model grows with photon
% count whether or not a third component exists. Measured on this pipeline,
% the false three-component rate rose from 9% to 54% between 2000 and 50000
% photons at a fixed threshold.
%
% When a Bayesian statistic is used outside the regime where its
% probabilistic reading holds, the way to keep it usable is to characterise
% its operating characteristics and report a decision rule with a known error
% rate. That is what this does: it fixes the false-positive rate and solves
% for the threshold, per photon count.
%
% METHOD
%
%   null        truth = fixed SLB + ONE membrane lifetime (no third component)
%   alternative truth = fixed SLB + TWO membrane lifetimes
%
% The decision statistic is the same one the pipeline gates on:
%   q3 = P(M3) / (P(M2) + P(M3)),  the probability that a third component is
%   present given that any membrane component is.
%
% For each photon count the threshold is the (1 - alpha) quantile of q3 under
% the pooled null, so by construction the false-positive rate is alpha.
% Detection power is then measured by applying that threshold to the
% alternative.
%
% IMPORTANT: the fit settings default to the PRODUCTION values used by
% run_batch_immune_cell_MIET, because a threshold calibrated with a finer grid
% does not transfer to a coarser one. If the pipeline settings change, this
% must be re-run.
%
% Two limits to state when citing the output:
%   1. It bounds the consequences of mis-specification, it does not remove
%      them. The output is an error rate, not a correctness guarantee.
%   2. Every number inherits the simulator's IRF, periodicity and grid. Run
%      with placeTruthOnGrid = true to see how much of the effect is the
%      pipeline's own lifetime discretisation.
%
% opts fields (all optional)
%   photonTotals          default [500 1000 2000 4000 8000 16000 32000 64000]
%   targetFalsePositive   alpha, default 0.05
%   pixelsPerCondition    null/alt pixels per condition, default 2000
%   nullAmplitudes        rows of [slb comp2], default three splits
%   altAmplitudes         rows of [slb comp2 comp3], default three splits
%   nullTau2Range         default [0.5 2.5]
%   altTau2Range          default [0.5 2.0]
%   altTau3Range          default [2.5 5.0]
%   tauSlbNs              default 0.3
%   pulsePeriodNs, dtNs   default 12.5 and 0.16
%   irfFwhmNs             default 0.35
%   membraneTauCount      default 10   (production)
%   fractionStep          default 0.2  (production)
%   minimumMembraneFraction default 0.1 (production)
%   membraneTauBoundsNs   default [] = replicate the pipeline default
%   placeTruthOnGrid      snap true lifetimes onto the fit grid, default false
%   outputDir             default pwd
%   seed                  default 7

    if nargin < 1 || isempty(opts)
        opts = struct();
    end
    opts = fillDefaults(opts);
    rng(opts.seed);

    if exist('flim_bayes_fixed_slb', 'file') ~= 2
        error('study_flim_component_thresholds:MissingDependency', ...
            'flim_bayes_fixed_slb.m must be on the MATLAB path.');
    end
    if ~isfolder(opts.outputDir)
        mkdir(opts.outputDir);
    end

    nBins = round(opts.pulsePeriodNs / opts.dtNs);
    timeNs = (0:nBins-1) * opts.dtNs;
    irf = gaussianIrf(timeNs, opts.irfFwhmNs);
    fitGrid = logspace(log10(opts.membraneTauBoundsNs(1)), ...
        log10(opts.membraneTauBoundsNs(2)), opts.membraneTauCount);

    fprintf('\nstudy_flim_component_thresholds\n');
    fprintf(['  instrument   : period %.3g ns, dt %.3g ns, %d bins, ' ...
        'IRF FWHM %.3g ns\n'], opts.pulsePeriodNs, opts.dtNs, nBins, ...
        opts.irfFwhmNs);
    fprintf(['  fit grid     : %d points over %.3g-%.3g ns (production ' ...
        'settings)\n'], opts.membraneTauCount, ...
        opts.membraneTauBoundsNs(1), opts.membraneTauBoundsNs(2));
    fprintf('  truth on grid: %d\n', opts.placeTruthOnGrid);
    fprintf('  target false-positive rate: %.3g\n', opts.targetFalsePositive);
    fprintf('  pixels per condition: %d\n\n', opts.pixelsPerCondition);

    totals = opts.photonTotals(:).';
    rows = struct([]);
    nullStats = cell(1, numel(totals));
    altStats = cell(1, numel(totals));

    for t = 1:numel(totals)
        photonTotal = totals(t);

        % ---- null: SLB + one membrane component ----
        q3Null = [];
        for r = 1:size(opts.nullAmplitudes, 1)
            fractions = [opts.nullAmplitudes(r, :), 0];
            taus = drawTaus(opts, fitGrid, opts.nullTau2Range, [], ...
                opts.pixelsPerCondition);
            q3 = runCondition(opts, irf, timeNs, nBins, photonTotal, ...
                fractions, taus);
            q3Null = [q3Null; q3(:)]; %#ok<AGROW>
        end

        % ---- alternative: SLB + two membrane components ----
        q3Alt = [];
        for r = 1:size(opts.altAmplitudes, 1)
            fractions = opts.altAmplitudes(r, :);
            taus = drawTaus(opts, fitGrid, opts.altTau2Range, ...
                opts.altTau3Range, opts.pixelsPerCondition);
            q3 = runCondition(opts, irf, timeNs, nBins, photonTotal, ...
                fractions, taus);
            q3Alt = [q3Alt; q3(:)]; %#ok<AGROW>
        end

        nullStats{t} = q3Null;
        altStats{t} = q3Alt;

        % Threshold is the (1-alpha) quantile of the null statistic, so the
        % false-positive rate equals alpha by construction.
        validNull = q3Null(isfinite(q3Null));
        validAlt = q3Alt(isfinite(q3Alt));
        if isempty(validNull)
            threshold = NaN;
        else
            threshold = quantile(validNull, 1 - opts.targetFalsePositive);
        end
        % A threshold at or above 1 means no attainable operating point:
        % the null statistic saturates and specificity cannot be bought.
        saturated = ~isfinite(threshold) || threshold >= 1 - 1e-9;

        if isempty(validAlt) || ~isfinite(threshold)
            power = NaN;
        else
            power = mean(validAlt >= threshold);
        end
        fixedRate = NaN;
        if ~isempty(validNull)
            fixedRate = mean(validNull >= opts.referenceThreshold);
        end
        fixedPower = NaN;
        if ~isempty(validAlt)
            fixedPower = mean(validAlt >= opts.referenceThreshold);
        end

        entry = struct( ...
            'photonTotal', photonTotal, ...
            'targetFalsePositive', opts.targetFalsePositive, ...
            'calibratedThreshold', threshold, ...
            'thresholdSaturated', saturated, ...
            'powerAtCalibrated', power, ...
            'referenceThreshold', opts.referenceThreshold, ...
            'falsePositiveAtReference', fixedRate, ...
            'powerAtReference', fixedPower, ...
            'nullMedianQ3', median(validNull), ...
            'altMedianQ3', median(validAlt), ...
            'nullCount', numel(validNull), ...
            'altCount', numel(validAlt));
        if isempty(rows)
            rows = entry;
        else
            rows(end+1) = entry; %#ok<AGROW>
        end

        fprintf(['  N=%6d  threshold=%.4f%s  power=%.3f  |  at fixed ' ...
            '%.2f: FP=%.3f power=%.3f\n'], photonTotal, threshold, ...
            ternaryText(saturated, ' (SATURATED)', ''), power, ...
            opts.referenceThreshold, fixedRate, fixedPower);
    end

    summary = struct2table(rows);
    calibration = struct('summary', summary, 'opts', opts, ...
        'nullStatistic', {nullStats}, 'altStatistic', {altStats}, ...
        'fitGridNs', fitGrid, 'photonTotals', totals);

    csvFile = fullfile(opts.outputDir, 'flim_component_thresholds.csv');
    matFile = fullfile(opts.outputDir, 'flim_component_thresholds.mat');
    writetable(summary, csvFile);
    save(matFile, 'calibration', '-v7.3');
    figureFile = plotCalibration(calibration);

    fprintf('\n  wrote %s\n  wrote %s\n  wrote %s\n', csvFile, matFile, ...
        figureFile);
    printUsage(calibration);
end

% ------------------------------------------------------------------ options

function opts = fillDefaults(opts)
    defaults = struct( ...
        'photonTotals', [500 1000 2000 4000 8000 16000 32000 64000], ...
        'targetFalsePositive', 0.05, ...
        'referenceThreshold', 0.95, ...
        'pixelsPerCondition', 2000, ...
        'nullAmplitudes', [0.60 0.40; 0.45 0.55; 0.30 0.70], ...
        'altAmplitudes', [0.50 0.30 0.20; 0.34 0.33 0.33; 0.60 0.30 0.10], ...
        'nullTau2Range', [0.5 2.5], ...
        'altTau2Range', [0.5 2.0], ...
        'altTau3Range', [2.5 5.0], ...
        'tauSlbNs', 0.3, ...
        'pulsePeriodNs', 12.5, 'dtNs', 0.16, 'irfFwhmNs', 0.35, ...
        'membraneTauCount', 10, ...
        'fractionStep', 0.2, ...
        'minimumMembraneFraction', 0.1, ...
        'signalGrid', [0.25 0.5 0.75 1], ...
        'membraneTauBoundsNs', [], ...
        'placeTruthOnGrid', false, ...
        'outputDir', pwd, 'seed', 7);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if isempty(opts.membraneTauBoundsNs)
        % Replicate the pipeline default exactly, so a threshold calibrated
        % here is valid for the production configuration.
        tauMinimum = max([1.15 * opts.tauSlbNs, ...
            opts.tauSlbNs + 2 * opts.dtNs, 0.05]);
        tauMaximum = min(0.8 * opts.pulsePeriodNs, ...
            max(5, 4 * opts.tauSlbNs));
        tauMaximum = max(tauMaximum, 1.25 * tauMinimum);
        opts.membraneTauBoundsNs = [tauMinimum tauMaximum];
    end
    if any(abs(sum(opts.nullAmplitudes, 2) - 1) > 1e-9)
        error('study_flim_component_thresholds:NullAmplitudes', ...
            'Each nullAmplitudes row must sum to 1.');
    end
    if any(abs(sum(opts.altAmplitudes, 2) - 1) > 1e-9)
        error('study_flim_component_thresholds:AltAmplitudes', ...
            'Each altAmplitudes row must sum to 1.');
    end
end

% ----------------------------------------------------------------- sampling

function taus = drawTaus(opts, fitGrid, range2, range3, count)
    tau2 = range2(1) + (range2(2) - range2(1)) * rand(count, 1);
    if isempty(range3)
        tau3 = nan(count, 1);
    else
        tau3 = range3(1) + (range3(2) - range3(1)) * rand(count, 1);
    end
    if opts.placeTruthOnGrid
        tau2 = snapToGrid(tau2, fitGrid);
        if ~isempty(range3)
            tau3 = snapToGrid(tau3, fitGrid);
        end
    end
    taus = [tau2, tau3];
end

function value = snapToGrid(value, grid)
    for k = 1:numel(value)
        [~, nearest] = min(abs(grid - value(k)));
        value(k) = grid(nearest);
    end
end

function q3 = runCondition(opts, irf, timeNs, nBins, photonTotal, ...
        fractions, taus)
    nPix = size(taus, 1);
    Y = zeros(nPix, 1, nBins);
    for k = 1:nPix
        componentTaus = [opts.tauSlbNs, taus(k, 1), taus(k, 2)];
        shape = mixtureShape(irf, timeNs, opts.pulsePeriodNs, ...
            componentTaus, fractions);
        Y(k, 1, :) = reshape(multinomialDraw(shape, photonTotal), 1, 1, nBins);
    end

    fitOpts = struct( ...
        'analysisMask', true(nPix, 1), ...
        'minPhotons', 1, ...
        'useGPU', false, ...
        'batchSize', 2048, ...
        'includeBackground', true, ...
        'membraneTauBoundsNs', opts.membraneTauBoundsNs, ...
        'membraneTauCount', opts.membraneTauCount, ...
        'signalGrid', opts.signalGrid, ...
        'fractionStep', opts.fractionStep, ...
        'minimumMembraneFraction', opts.minimumMembraneFraction, ...
        'slbCountRelTol', 0, ...
        'fixedSlbPhotonCount', fractions(1) * photonTotal, ...
        'fixedSlbPhotonCountStd', 0.1 * fractions(1) * photonTotal);

    out = flim_bayes_fixed_slb(Y, irf, opts.pulsePeriodNs, opts.dtNs, ...
        opts.tauSlbNs, fitOpts);

    % The pipeline's own statistic: P(M3) conditional on any membrane
    % component being present. See immune_cell_MIET_sorted_components.
    p2 = double(out.modelProbability(:, :, 2));
    p3 = double(out.modelProbability(:, :, 3));
    q2 = p2 + p3;
    q3 = nan(size(q2));
    usable = isfinite(q2) & q2 > 0;
    q3(usable) = p3(usable) ./ q2(usable);
end

% --------------------------------------------------------------- simulation

function irf = gaussianIrf(timeNs, fwhmNs)
    sigma = fwhmNs / (2 * sqrt(2 * log(2)));
    centre = max(4 * sigma, 3 * (timeNs(2) - timeNs(1)));
    irf = exp(-0.5 * ((timeNs - centre) / sigma).^2);
    irf = irf / sum(irf);
end

function shape = mixtureShape(irf, timeNs, periodNs, taus, fractions)
    shape = zeros(1, numel(timeNs));
    for k = 1:numel(taus)
        if fractions(k) <= 0 || ~isfinite(taus(k)) || taus(k) <= 0
            continue;
        end
        shape = shape + fractions(k) * ...
            periodicDecay(irf, timeNs, periodNs, taus(k));
    end
    shape = shape + 1e-12;
    shape = shape / sum(shape);
end

function pattern = periodicDecay(irf, timeNs, periodNs, tauNs)
    nBins = numel(timeNs);
    decay = zeros(1, nBins);
    for repeat = 0:3
        decay = decay + exp(-(timeNs + repeat * periodNs) / tauNs);
    end
    full = conv(irf, decay);
    pattern = zeros(1, nBins);
    for start = 1:nBins:numel(full)
        stop = min(start + nBins - 1, numel(full));
        span = stop - start + 1;
        pattern(1:span) = pattern(1:span) + full(start:stop);
    end
    pattern = pattern / sum(pattern);
end

function counts = multinomialDraw(shape, photonTotal)
    edges = [0, cumsum(shape)];
    edges(end) = 1;
    counts = histcounts(rand(1, photonTotal), edges);
end

% ------------------------------------------------------------------ reports

function text = ternaryText(condition, a, b)
    if condition
        text = a;
    else
        text = b;
    end
end

function printUsage(calibration)
    summary = calibration.summary;
    opts = calibration.opts;
    fprintf('\n  --- how to use this ---\n');
    fprintf(['  The pipeline applies one fixed threshold at every binning ' ...
        'stage. Replace it\n  with the photon-count-dependent value ' ...
        'below, indexed by the window photon\n  total (expectedPhotonCount ' ...
        'summed over components):\n\n']);
    fprintf('    photons   threshold\n');
    for k = 1:height(summary)
        fprintf('    %7d   %.4f%s\n', summary.photonTotal(k), ...
            summary.calibratedThreshold(k), ...
            ternaryText(summary.thresholdSaturated(k), ...
                '   <- saturated, no usable operating point', ''));
    end
    fprintf(['\n  At the current fixed threshold of %.2f the false-positive ' ...
        'rate ranges from\n  %.3f to %.3f across this photon range, which ' ...
        'is the drift the calibration\n  removes.\n'], ...
        opts.referenceThreshold, ...
        min(summary.falsePositiveAtReference), ...
        max(summary.falsePositiveAtReference));
    if any(summary.thresholdSaturated)
        fprintf(['\n  WARNING: at least one photon count has no attainable ' ...
            'operating point at\n  alpha = %.3g. There, a third component ' ...
            'cannot be declared at that error\n  rate by any threshold, ' ...
            'and the component map should be reported as\n  ' ...
            'indeterminate rather than gated.\n'], opts.targetFalsePositive);
    end
end

% ------------------------------------------------------------------ figures

function name = plotCalibration(calibration)
    summary = calibration.summary;
    opts = calibration.opts;
    totals = calibration.photonTotals;

    h = figure('Color', 'w', 'Visible', 'off', 'Position', [60 60 1320 420]);
    layout = tiledlayout(h, 1, 3, 'Padding', 'compact', ...
        'TileSpacing', 'compact');

    % (a) threshold vs photon count
    ax = nexttile(layout);
    semilogx(ax, summary.photonTotal, summary.calibratedThreshold, ...
        'o-', 'LineWidth', 1.5, 'Color', [0.00 0.45 0.74], ...
        'DisplayName', sprintf('calibrated (alpha=%.2g)', ...
            opts.targetFalsePositive));
    hold(ax, 'on');
    yline(ax, opts.referenceThreshold, 'r--', ...
        'DisplayName', sprintf('current fixed %.2f', opts.referenceThreshold));
    hold(ax, 'off');
    grid(ax, 'on');
    ylim(ax, [0 1.02]);
    xlabel(ax, 'photons per pixel or window');
    ylabel(ax, 'threshold on q3');
    title(ax, 'Calibrated decision threshold', 'FontSize', 9);
    legend(ax, 'Location', 'southeast', 'FontSize', 7, 'Box', 'off');

    % (b) false-positive rate at the fixed threshold
    ax = nexttile(layout);
    semilogx(ax, summary.photonTotal, summary.falsePositiveAtReference, ...
        's-', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10], ...
        'DisplayName', sprintf('at fixed %.2f', opts.referenceThreshold));
    hold(ax, 'on');
    yline(ax, opts.targetFalsePositive, 'k--', ...
        'DisplayName', sprintf('target %.2g', opts.targetFalsePositive));
    hold(ax, 'off');
    grid(ax, 'on');
    xlabel(ax, 'photons per pixel or window');
    ylabel(ax, 'false three-component rate');
    title(ax, 'Why a fixed threshold drifts', 'FontSize', 9);
    legend(ax, 'Location', 'northwest', 'FontSize', 7, 'Box', 'off');

    % (c) power at the calibrated threshold
    ax = nexttile(layout);
    semilogx(ax, summary.photonTotal, summary.powerAtCalibrated, ...
        'd-', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19], ...
        'DisplayName', 'at calibrated threshold');
    hold(ax, 'on');
    semilogx(ax, summary.photonTotal, summary.powerAtReference, ...
        'd:', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5], ...
        'DisplayName', 'at fixed threshold');
    hold(ax, 'off');
    grid(ax, 'on');
    ylim(ax, [0 1.02]);
    xlabel(ax, 'photons per pixel or window');
    ylabel(ax, 'detection power');
    title(ax, 'Cost in sensitivity', 'FontSize', 9);
    legend(ax, 'Location', 'southeast', 'FontSize', 7, 'Box', 'off');

    title(layout, sprintf(['Third-component decision rule | null = SLB + 1 ' ...
        'membrane, alternative = SLB + 2 | truth on grid = %d'], ...
        opts.placeTruthOnGrid), 'FontWeight', 'bold');
    subtitle(layout, sprintf(['production fit settings: %d grid points over ' ...
        '%.2g-%.2g ns, fractionStep %.3g'], opts.membraneTauCount, ...
        opts.membraneTauBoundsNs(1), opts.membraneTauBoundsNs(2), ...
        opts.fractionStep), 'FontSize', 8);

    name = fullfile(opts.outputDir, 'flim_component_thresholds.png');
    exportgraphics(h, name, 'Resolution', 200);
    pdfName = fullfile(opts.outputDir, 'flim_component_thresholds.pdf');
    exportgraphics(h, pdfName, 'ContentType', 'vector');
    close(h);
    if ~isempty(totals)
        % keeps totals referenced for clarity of intent
    end
end
