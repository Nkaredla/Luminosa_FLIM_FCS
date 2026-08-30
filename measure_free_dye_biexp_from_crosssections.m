function out = measure_free_dye_biexp_from_crosssections(dataRoot, opts)
%MEASURE_FREE_DYE_BIEXP_FROM_CROSSSECTIONS Resolve tau_0 into two components.
%
% out = measure_free_dye_biexp_from_crosssections()
% out = measure_free_dye_biexp_from_crosssections(dataRoot, opts)
%
% Re-fits the pooled above-surface decay of every XZ/YZ cross-section in a
% session with ONE and TWO exponentials, and reports the longer component.
%
% WHY THE MONO FIT IS NOT tau_0
%
% estimate_free_dye_lifetime_from_crosssection reports a single-exponential
% tail fit of the region far above the metal, and on this session that gives
% 1.92-2.01 ns. The saved residuals show that fit is wrong in a specific way:
% systematically below the data from about 1.5 to 5 ns and systematically above
% it from 6 to 15 ns. That S-shaped residual is the signature of a mono fit to
% a two-component decay - it is too long early and too short late - so the
% reported number is a PHOTON-WEIGHTED BLEND of a short and a long component,
% not the lifetime of either.
%
% This matters for MIET specifically. A MIET curve approaches tau_0 from below
% and can never invert a lifetime at or above it, so tau_0 sets the ceiling of
% the whole calibration. Using a blended value that is pulled DOWN by a short
% contaminating pool makes the ceiling too low and throws away exactly the
% pixels that sit near it. The longer component is the candidate for the
% unquenched, bilayer-inserted dye that the cell membrane consists of.
%
% WHAT IS AND IS NOT SETTLED BY THIS
%
% Two components are reported with their photon shares. This routine does NOT
% assert which one is the membrane pool - it measures both, sweeps the tail
% start so the reader can see whether each is stable, and reports the deviance
% improvement over the mono fit so the second component has to earn its place.
% A component that moves with the tail start is a fitting artefact regardless
% of how physical it looks.
%
% An IRF-free TAIL fit is used, as in the original routine, so the answer does
% not depend on which IRF is adopted. The cost is sensitivity to where the tail
% starts, which is why the start is swept rather than chosen.
%
% INPUTS
%   dataRoot  session folder, default the RT MemGlow session of 2026-08-13
%   opts .tailStartNs   after the peak, swept (default [0.6 1.0 1.5 2.0 3.0])
%        .outputDir     '' puts results in <dataRoot>\free_dye_above_surface
%        .makeFigure    true
%        .showFigure    false
%
% See also ESTIMATE_FREE_DYE_LIFETIME_FROM_CROSSSECTION, TAILFIT,
% MAKE_RT_MIET_CALIBRATION.

    if nargin < 1 || isempty(dataRoot)
        dataRoot = 'D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_20260813_1';
    end
    if nargin < 2 || isempty(opts); opts = struct(); end
    opts = withDefaults(opts, struct( ...
        'tailStartNs', [0.6 1.0 1.5 2.0 3.0], ...
        'outputDir', '', ...
        'makeFigure', true, ...
        'showFigure', false));
    if isempty(opts.outputDir)
        opts.outputDir = fullfile(dataRoot, 'free_dye_above_surface');
    end
    if ~isfolder(opts.outputDir); mkdir(opts.outputDir); end

    found = dir(fullfile(dataRoot, '**', 'free_dye_above_surface.mat'));
    found = found(~[found.isdir]);
    if isempty(found)
        error('measure_free_dye_biexp_from_crosssections:NoInputs', ...
            ['No free_dye_above_surface.mat under %s. Run ' ...
             'run_free_dye_lifetime_above_surface first.'], dataRoot);
    end

    fprintf('measure_free_dye_biexp_from_crosssections\n');
    fprintf('  session %s\n', dataRoot);
    fprintf('  %d cross-sections\n\n', numel(found));

    records = struct('acquisition', {}, 'plane', {}, 'dtNs', {}, ...
        'decay', {}, 'peakBin', {}, 'monoTauNs', {}, 'sweep', {}, ...
        'tauShortNs', {}, 'tauLongNs', {}, 'longPhotonShare', {}, ...
        'monoDeviance', {}, 'biexpDeviance', {}, 'tauLongSpreadNs', {});

    for index = 1:numel(found)
        matFile = fullfile(found(index).folder, found(index).name);
        loaded = load(matFile, 'out');
        if ~isfield(loaded, 'out'); continue; end
        source = loaded.out;
        if ~isfield(source, 'pooledDecay') || isempty(source.pooledDecay)
            fprintf('  %s: no pooled decay stored, skipped\n', source.acquisition);
            continue;
        end
        record = fitOne(source, opts);
        records(end + 1) = record; %#ok<AGROW>
        fprintf(['  %-18s %s  mono %.3f ns | biexp %.3f + %.3f ns ' ...
                 '(long share %.2f) | deviance %.0f -> %.0f | long spread %.3f ns\n'], ...
            record.acquisition, record.plane, record.monoTauNs, ...
            record.tauShortNs, record.tauLongNs, record.longPhotonShare, ...
            record.monoDeviance, record.biexpDeviance, record.tauLongSpreadNs);
    end

    if isempty(records)
        error('measure_free_dye_biexp_from_crosssections:NoFits', ...
            'No cross-section held a usable pooled decay.');
    end

    out = struct();
    out.dataRoot = dataRoot;
    out.opts = opts;
    out.records = records;
    out.summary = buildTable(records);
    out.monoMedianNs = median([records.monoTauNs]);
    out.tauShortMedianNs = median([records.tauShortNs]);
    out.tauLongMedianNs = median([records.tauLongNs]);
    out.tauLongRangeNs = [min([records.tauLongNs]) max([records.tauLongNs])];
    out.longPhotonShareMedian = median([records.longPhotonShare]);
    out.medianTauLongSpreadNs = median([records.tauLongSpreadNs]);
    out.recommendation = recommend(out);

    csvFile = fullfile(opts.outputDir, 'free_dye_biexp_summary.csv');
    matOut = fullfile(opts.outputDir, 'free_dye_biexp_summary.mat');
    writetable(out.summary, csvFile);
    save(matOut, 'out', '-v7.3');
    out.outputFiles = struct('csv', csvFile, 'mat', matOut, 'figure', '');

    if opts.makeFigure
        figureFile = fullfile(opts.outputDir, 'free_dye_biexp_summary.png');
        try
            writeFigure(records, out, figureFile, opts);
            out.outputFiles.figure = figureFile;
            save(matOut, 'out', '-v7.3');
        catch figureError
            warning('measure_free_dye_biexp_from_crosssections:Figure', ...
                'Numbers were saved, figure failed: %s', figureError.message);
        end
    end

    fprintf('\n--- session ---\n');
    fprintf('mono (blended)     median %.3f ns\n', out.monoMedianNs);
    fprintf('short component    median %.3f ns\n', out.tauShortMedianNs);
    fprintf('LONG component     median %.3f ns, range %.3f-%.3f ns, ', ...
        out.tauLongMedianNs, out.tauLongRangeNs(1), out.tauLongRangeNs(2));
    fprintf('photon share %.2f\n', out.longPhotonShareMedian);
    fprintf('long-component stability across tail starts: %.3f ns\n', ...
        out.medianTauLongSpreadNs);
    fprintf('%s\n', out.recommendation);
    fprintf('saved %s\n', csvFile);
end

% ----------------------------------------------------------------- fitting

function record = fitOne(source, opts)
    decay = double(source.pooledDecay(:));
    dtNs = double(source.dtNs);
    if isfield(source, 'fit') && isfield(source.fit, 'peakBin') && ...
            ~isempty(source.fit.peakBin)
        peakBin = double(source.fit.peakBin);
    else
        [~, peakBin] = max(decay);
    end

    sweep = struct('tailStartNs', {}, 'monoTauNs', {}, 'tauNs', {}, ...
        'photonShare', {}, 'monoDeviance', {}, 'biexpDeviance', {});
    for start = opts.tailStartNs(:).'
        firstBin = peakBin + max(1, round(start / dtNs));
        if firstBin > numel(decay) - 20; continue; end
        tail = decay(firstBin:end);

        [tauMono, ampMono, offMono] = Tailfit(tail, dtNs, 2.0, [], 'mle');
        [tauBi, ampBi, offBi] = Tailfit(tail, dtNs, [0.6 2.5], [], 'mle');

        entry = struct();
        entry.tailStartNs = start;
        entry.monoTauNs = tauMono(1);
        [entry.tauNs, order] = sort(tauBi(:).', 'ascend');
        share = photonShare(tauBi(:).', ampBi(:).');
        entry.photonShare = share(order);
        entry.monoDeviance = poissonDeviance(tail, modelCurve(tail, dtNs, ...
            tauMono, ampMono, offMono));
        entry.biexpDeviance = poissonDeviance(tail, modelCurve(tail, dtNs, ...
            tauBi, ampBi, offBi));
        sweep(end + 1) = entry; %#ok<AGROW>
    end
    if isempty(sweep)
        error('measure_free_dye_biexp_from_crosssections:TailTooShort', ...
            'No tail start left enough bins for %s.', source.acquisition);
    end

    longTaus = arrayfun(@(e) e.tauNs(end), sweep);
    shortTaus = arrayfun(@(e) e.tauNs(1), sweep);
    longShares = arrayfun(@(e) e.photonShare(end), sweep);

    record = struct();
    record.acquisition = source.acquisition;
    record.plane = source.plane;
    record.dtNs = dtNs;
    record.decay = decay;
    record.peakBin = peakBin;
    record.monoTauNs = median(arrayfun(@(e) e.monoTauNs, sweep));
    record.sweep = sweep;
    record.tauShortNs = median(shortTaus);
    record.tauLongNs = median(longTaus);
    record.longPhotonShare = median(longShares);
    record.monoDeviance = median(arrayfun(@(e) e.monoDeviance, sweep));
    record.biexpDeviance = median(arrayfun(@(e) e.biexpDeviance, sweep));
    % Spread across tail starts is the honesty check: a component that is real
    % does not move when the fit window does.
    record.tauLongSpreadNs = max(longTaus) - min(longTaus);
end

function share = photonShare(tau, amp)
% Photons, not amplitudes: a component's photon count is amp*tau, and quoting
% amplitude shares makes a short component look far more abundant than it is.
    weight = amp(:).' .* tau(:).';
    weight = max(weight, 0);
    total = sum(weight);
    if total <= 0
        share = zeros(size(weight));
    else
        share = weight / total;
    end
end

function model = modelCurve(tail, dtNs, tau, amp, offset)
    t = (0:numel(tail) - 1).' * dtNs;
    model = offset + exp(-t ./ tau(:).') * amp(:);
end

function deviance = poissonDeviance(observed, model)
    observed = double(observed(:));
    model = max(double(model(:)), 1e-9);
    positive = observed > 0;
    terms = 2 * (model - observed);
    terms(positive) = terms(positive) + ...
        2 * observed(positive) .* log(observed(positive) ./ model(positive));
    deviance = sum(terms) / max(numel(observed) - 4, 1);
end

function text = recommend(out)
    stable = out.medianTauLongSpreadNs < 0.15;
    improves = out.monoMedianNs > 0 && ...
        median([out.records.biexpDeviance]) < 0.8 * median([out.records.monoDeviance]);
    if stable && improves
        text = sprintf(['=> Two components are supported and the long one is ' ...
            'stable. Use tau_free = %.2f ns for the membrane calibration, not ' ...
            'the blended %.2f ns.'], out.tauLongMedianNs, out.monoMedianNs);
    elseif ~improves
        text = ['=> The second component does not improve the fit enough to ' ...
            'be trusted; keep the mono value.'];
    else
        text = sprintf(['=> The long component moves %.2f ns across tail ' ...
            'starts, so it is not well determined. Treat %.2f ns as an upper ' ...
            'bound rather than a measurement.'], out.medianTauLongSpreadNs, ...
            out.tauLongMedianNs);
    end
end

% ----------------------------------------------------------------- outputs

function summary = buildTable(records)
    summary = table();
    summary.acquisition = string({records.acquisition}.');
    summary.plane = string({records.plane}.');
    summary.monoTauNs = [records.monoTauNs].';
    summary.tauShortNs = [records.tauShortNs].';
    summary.tauLongNs = [records.tauLongNs].';
    summary.longPhotonShare = [records.longPhotonShare].';
    summary.tauLongSpreadNs = [records.tauLongSpreadNs].';
    summary.monoReducedDeviance = [records.monoDeviance].';
    summary.biexpReducedDeviance = [records.biexpDeviance].';
end

function writeFigure(records, out, figureFile, opts)
    visibility = 'off';
    if opts.showFigure; visibility = 'on'; end
    h = figure('Color', 'w', 'Visible', visibility, 'Position', [40 40 1500 820]);
    layout = tiledlayout(h, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % One representative decay with both fits and their residuals.
    record = records(1);
    best = record.sweep(ceil(numel(record.sweep) / 2));
    firstBin = record.peakBin + max(1, round(best.tailStartNs / record.dtNs));
    tail = record.decay(firstBin:end);
    t = (0:numel(tail) - 1).' * record.dtNs;
    [tauMono, ampMono, offMono] = Tailfit(tail, record.dtNs, 2.0, [], 'mle');
    [tauBi, ampBi, offBi] = Tailfit(tail, record.dtNs, [0.6 2.5], [], 'mle');
    monoModel = modelCurve(tail, record.dtNs, tauMono, ampMono, offMono);
    biModel = modelCurve(tail, record.dtNs, tauBi, ampBi, offBi);

    ax = nexttile(layout);
    semilogy(ax, t, max(tail, 0.5), '.', 'Color', [0.35 0.35 0.35]); hold(ax, 'on');
    semilogy(ax, t, monoModel, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6);
    semilogy(ax, t, biModel, '-', 'Color', [0.10 0.35 0.70], 'LineWidth', 1.6);
    xlabel(ax, 'time after tail start [ns]'); ylabel(ax, 'photons per bin');
    legend(ax, {'pooled tail', sprintf('mono %.3f ns', tauMono(1)), ...
        sprintf('biexp %.2f + %.2f ns', min(tauBi), max(tauBi))}, ...
        'Location', 'northeast');
    title(ax, sprintf('%s (%s), tail start %.1f ns', record.acquisition, ...
        record.plane, best.tailStartNs), 'Interpreter', 'none');
    grid(ax, 'on'); box(ax, 'off');

    ax = nexttile(layout);
    plot(ax, t, poissonResidual(tail, monoModel), '-', ...
        'Color', [0.85 0.33 0.10], 'LineWidth', 1); hold(ax, 'on');
    plot(ax, t, poissonResidual(tail, biModel), '-', ...
        'Color', [0.10 0.35 0.70], 'LineWidth', 1);
    yline(ax, 0, 'k-');
    xlabel(ax, 'time after tail start [ns]'); ylabel(ax, 'Poisson residual');
    title(ax, 'the S-shaped mono residual is what the second component removes');
    legend(ax, {'mono', 'biexp'}, 'Location', 'northeast');
    grid(ax, 'on'); box(ax, 'off');

    % Stability of the long component across tail starts, all acquisitions.
    ax = nexttile(layout);
    hold(ax, 'on');
    colours = lines(numel(records));
    for index = 1:numel(records)
        starts = arrayfun(@(e) e.tailStartNs, records(index).sweep);
        longs = arrayfun(@(e) e.tauNs(end), records(index).sweep);
        monos = arrayfun(@(e) e.monoTauNs, records(index).sweep);
        plot(ax, starts, longs, '-o', 'Color', colours(index,:), 'LineWidth', 1.4, ...
            'MarkerFaceColor', colours(index,:), 'MarkerSize', 4);
        plot(ax, starts, monos, '--', 'Color', colours(index,:), 'LineWidth', 1);
    end
    xlabel(ax, 'tail start after the peak [ns]'); ylabel(ax, 'lifetime [ns]');
    title(ax, 'long component (solid) versus mono blend (dashed)');
    grid(ax, 'on'); box(ax, 'off');

    ax = nexttile(layout);
    values = [[records.monoTauNs].'; [records.tauLongNs].'];
    groups = [repmat("mono blend", numel(records), 1); ...
              repmat("long component", numel(records), 1)];
    boxchart(ax, categorical(groups), values);
    ylabel(ax, 'lifetime [ns]');
    title(ax, sprintf('per-acquisition values (n = %d cross-sections)', ...
        numel(records)));
    grid(ax, 'on'); box(ax, 'off');

    title(layout, 'Above-surface pooled decay is biexponential: the mono tau_0 is a blend', ...
        'FontWeight', 'bold');
    subtitle(layout, out.recommendation, 'Interpreter', 'none', 'FontSize', 9);
    exportgraphics(h, figureFile, 'Resolution', 220);
    if ~opts.showFigure; close(h); end
end

function residual = poissonResidual(observed, model)
    observed = double(observed(:));
    model = max(double(model(:)), 1e-9);
    residual = sign(observed - model) .* sqrt(max(2 * ( ...
        model - observed + observed .* log(max(observed, 1e-9) ./ model)), 0));
end

function opts = withDefaults(opts, defaults)
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(opts, names{index}) || isempty(opts.(names{index}))
            opts.(names{index}) = defaults.(names{index});
        end
    end
end
