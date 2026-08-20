function results = study_flim_triexp_recovery(opts)
%STUDY_FLIM_TRIEXP_RECOVERY Can the fixed-SLB method recover a triexponential?
%
% results = study_flim_triexp_recovery()
% results = study_flim_triexp_recovery(opts)
%
% Simulates 2D pixel TCSPC images with a KNOWN triexponential decay and fits
% them with flim_bayes_fixed_slb - the same function the pipeline uses - then
% reports how well the two free lifetimes and the model order are recovered.
%
% Image layout: the second-component lifetime varies along x and the third
% along y, so every pixel is a different (tau2, tau3) pair. The SLB component
% is held at a fixed lifetime, as in the real analysis. The sweep repeats over
% several amplitude splits and several photon budgets.
%
% What is measured is recoverability, not goodness of fit. The truth is known,
% so the output is bias and scatter in tau2 and tau3 plus how often the
% three-component model is actually selected.
%
% Two things to keep in mind when reading the result.
%
%   1. The fixed SLB photon count is supplied EXACTLY here. In the real
%      pipeline it comes from an outside-cell reference and carries its own
%      error, so these numbers are a best case. Set opts.slbCountBiasFraction
%      nonzero to see how fast that assumption matters.
%
%   2. Recovery is bounded by the lifetime grid, not only by photon noise. The
%      grid step is printed with the results: if the median error sits near
%      half the local step, the grid is the limit rather than the method. Note
%      the pipeline default puts tauMinimum at
%      max(1.15*tauSlb, tauSlb + 2*dt, 0.05), which for tauSlb = 0.3 ns and
%      dt = 0.16 ns is 0.62 ns - above the 0.5 ns end of this sweep. The study
%      therefore sets the fit grid explicitly.
%
% opts fields (all optional)
%   tau2Range, tau2Count      second component sweep, default [0.5 2.0], 12
%   tau3Range, tau3Count      third component sweep, default [2.5 5.0], 12
%   tauSlbNs                  fixed short component, default 0.3
%   photonTotals              default [500 2000 8000 20000 50000]
%   amplitudeSets             rows of [slb comp2 comp3] photon fractions
%   pulsePeriodNs, dtNs       instrument, default 12.5 and 0.16
%   irfFwhmNs                 Gaussian IRF width, default 0.35
%   membraneTauCount          fit grid size, default 24
%   membraneTauBoundsNs       fit grid span, default [0.4 5.5]
%   slbCountBiasFraction      mis-specify the fixed SLB count, default 0
%   exampleTaus               [tau2 tau3] rows to draw as example decays,
%                             default [0.8 3.0; 1.5 4.5]
%   examplePhotonTotals       which budgets to draw, default [500 8000 50000]
%   exampleAmplitudeSet       which amplitude row to draw, default 2
%   outputDir                 where PNG/CSV/MAT go, default pwd
%   seed                      default 42

    if nargin < 1 || isempty(opts)
        opts = struct();
    end
    opts = fillStudyDefaults(opts);
    rng(opts.seed);

    if exist('flim_bayes_fixed_slb', 'file') ~= 2
        error('study_flim_triexp_recovery:MissingDependency', ...
            'flim_bayes_fixed_slb.m must be on the MATLAB path.');
    end
    if ~isfolder(opts.outputDir)
        mkdir(opts.outputDir);
    end

    nBins = round(opts.pulsePeriodNs / opts.dtNs);
    timeNs = (0:nBins-1) * opts.dtNs;
    irf = gaussianIrf(timeNs, opts.irfFwhmNs);

    tau2Grid = linspace(opts.tau2Range(1), opts.tau2Range(2), opts.tau2Count);
    tau3Grid = linspace(opts.tau3Range(1), opts.tau3Range(2), opts.tau3Count);
    [tau2True, tau3True] = ndgrid(tau2Grid, tau3Grid);
    nx = opts.tau2Count;
    ny = opts.tau3Count;

    fitGrid = logspace(log10(opts.membraneTauBoundsNs(1)), ...
        log10(opts.membraneTauBoundsNs(2)), opts.membraneTauCount);

    % Snapping the truth onto the fit grid removes discretisation as a source
    % of model mis-specification. Comparing a snapped run against an unsnapped
    % one separates "the grid cannot represent this lifetime" from "the method
    % cannot resolve this lifetime" - identical in the error statistics
    % otherwise.
    if opts.placeTruthOnGrid
        tau2True = snapToGrid(tau2True, fitGrid);
        tau3True = snapToGrid(tau3True, fitGrid);
        tau2Grid = snapToGrid(tau2Grid, fitGrid);
        tau3Grid = snapToGrid(tau3Grid, fitGrid);
    end

    fprintf('\nstudy_flim_triexp_recovery\n');
    fprintf(['  instrument : period %.3g ns, dt %.3g ns, %d bins, ' ...
        'IRF FWHM %.3g ns\n'], opts.pulsePeriodNs, opts.dtNs, nBins, ...
        opts.irfFwhmNs);
    fprintf('  fixed SLB  : %.3g ns', opts.tauSlbNs);
    if opts.slbCountBiasFraction ~= 0
        fprintf(' (count mis-specified by %+.1f%%)', ...
            100 * opts.slbCountBiasFraction);
    end
    fprintf('\n');
    fprintf(['  truth      : tau2 %.2g-%.2g ns (%d steps), ' ...
        'tau3 %.2g-%.2g ns (%d steps)\n'], opts.tau2Range(1), ...
        opts.tau2Range(2), opts.tau2Count, opts.tau3Range(1), ...
        opts.tau3Range(2), opts.tau3Count);
    fprintf(['  fit grid   : %d points over %.2g-%.2g ns; local step ' ...
        '%.3f ns at 1 ns, %.3f ns at 4 ns\n'], opts.membraneTauCount, ...
        opts.membraneTauBoundsNs(1), opts.membraneTauBoundsNs(2), ...
        localGridStep(fitGrid, 1), localGridStep(fitGrid, 4));
    fprintf('\n');

    fprintf('  truth snapped onto fit grid: %d\n', opts.placeTruthOnGrid);
    fprintf('  SLB count prior nodes: %d (0 = hard constraint)\n\n', ...
        opts.slbCountPriorNodes);

    setCount = size(opts.amplitudeSets, 1);
    totalCount = numel(opts.photonTotals);
    rows = struct([]);
    resultGrids = cell(setCount, totalCount);
    exampleStore = cell(setCount, totalCount);

    for setIndex = 1:setCount
        fractions = opts.amplitudeSets(setIndex, :);
        for totalIndex = 1:totalCount
            photonTotal = opts.photonTotals(totalIndex);

            Y = zeros(nx, ny, nBins);
            for ix = 1:nx
                for iy = 1:ny
                    shape = mixtureShape(irf, timeNs, opts.pulsePeriodNs, ...
                        [opts.tauSlbNs, tau2True(ix, iy), ...
                         tau3True(ix, iy)], fractions);
                    Y(ix, iy, :) = reshape( ...
                        multinomialDraw(shape, photonTotal), 1, 1, nBins);
                end
            end

            fitOpts = struct( ...
                'analysisMask', true(nx, ny), ...
                'minPhotons', 1, ...
                'useGPU', false, ...
                'batchSize', 2048, ...
                'includeBackground', true, ...
                'membraneTauBoundsNs', opts.membraneTauBoundsNs, ...
                'membraneTauCount', opts.membraneTauCount, ...
                'signalGrid', [0.25 0.5 0.75 1], ...
                'fractionStep', 0.125, ...
                'minimumMembraneFraction', 0.05, ...
                'slbCountRelTol', 0, ...
                'slbCountPriorNodes', opts.slbCountPriorNodes, ...
                'fixedSlbPhotonCount', fractions(1) * photonTotal * ...
                    (1 + opts.slbCountBiasFraction), ...
                'fixedSlbPhotonCountStd', 0.1 * fractions(1) * photonTotal);

            started = tic;
            out = flim_bayes_fixed_slb(Y, irf, opts.pulsePeriodNs, ...
                opts.dtNs, opts.tauSlbNs, fitOpts);
            elapsed = toc(started);

            % Lifetimes conditional on the three-component model, which is
            % what the pipeline reports for components 2 and 3.
            tau2Fit = double(out.twoMembrane.membraneLifetime1Ns);
            tau3Fit = double(out.twoMembrane.membraneLifetime2Ns);
            modelMap = double(out.completeExponentialCountMAP);

            error2 = tau2Fit - tau2True;
            error3 = tau3Fit - tau3True;

            resultGrids{setIndex, totalIndex} = struct( ...
                'tau2True', tau2True, 'tau3True', tau3True, ...
                'tau2Fit', tau2Fit, 'tau3Fit', tau3Fit, ...
                'error2', error2, 'error3', error3, ...
                'modelMap', modelMap, 'fractions', fractions, ...
                'photonTotal', photonTotal, 'seconds', elapsed);

            % Keep the raw decay and the fitted parameters for a few
            % representative pixels, so the example figure shows real
            % simulated data against the model actually fitted to it.
            exampleStore{setIndex, totalIndex} = captureExamples(opts, ...
                tau2Grid, tau3Grid, Y, out, tau2Fit, tau3Fit, modelMap);

            entry = struct( ...
                'amplitudeSet', setIndex, ...
                'slbFraction', fractions(1), ...
                'comp2Fraction', fractions(2), ...
                'comp3Fraction', fractions(3), ...
                'photonTotal', photonTotal, ...
                'comp2Photons', fractions(2) * photonTotal, ...
                'comp3Photons', fractions(3) * photonTotal, ...
                'medianAbsError2Ns', median(abs(error2(:)), 'omitnan'), ...
                'medianAbsError3Ns', median(abs(error3(:)), 'omitnan'), ...
                'medianBias2Ns', median(error2(:), 'omitnan'), ...
                'medianBias3Ns', median(error3(:), 'omitnan'), ...
                'iqrError2Ns', iqrOmitNan(error2(:)), ...
                'iqrError3Ns', iqrOmitNan(error3(:)), ...
                'correctModelFraction', mean(modelMap(:) == 3), ...
                'gridStepAtTau2Ns', localGridStep(fitGrid, ...
                    mean(opts.tau2Range)), ...
                'gridStepAtTau3Ns', localGridStep(fitGrid, ...
                    mean(opts.tau3Range)), ...
                'seconds', elapsed);
            if isempty(rows)
                rows = entry;
            else
                rows(end+1) = entry; %#ok<AGROW>
            end

            fprintf(['  set %d [%.2f %.2f %.2f]  N=%6d  |dTau2|=%.3f  ' ...
                '|dTau3|=%.3f  M3 chosen %5.1f%%  (%.1f s)\n'], ...
                setIndex, fractions(1), fractions(2), fractions(3), ...
                photonTotal, entry.medianAbsError2Ns, ...
                entry.medianAbsError3Ns, ...
                100 * entry.correctModelFraction, elapsed);
        end
    end

    summary = struct2table(rows);
    results = struct('summary', summary, 'grids', {resultGrids}, ...
        'examples', {exampleStore}, ...
        'opts', opts, 'irf', irf, 'timeNs', timeNs, 'fitGridNs', fitGrid);

    csvFile = fullfile(opts.outputDir, 'flim_triexp_recovery_summary.csv');
    matFile = fullfile(opts.outputDir, 'flim_triexp_recovery.mat');
    writetable(summary, csvFile);
    save(matFile, 'results', '-v7.3');

    figureFiles = plotRecovery(results);
    figureFiles = [figureFiles, plotExampleDecays(results)];
    fprintf('\n  wrote %s\n', csvFile);
    fprintf('  wrote %s\n', matFile);
    for k = 1:numel(figureFiles)
        fprintf('  wrote %s\n', figureFiles{k});
    end
    reportVerdict(results, fitGrid);
end

% ------------------------------------------------------------------ options

function opts = fillStudyDefaults(opts)
    defaults = struct( ...
        'tau2Range', [0.5 2.0], 'tau2Count', 20, ...
        'tau3Range', [2.5 5.0], 'tau3Count', 20, ...
        'tauSlbNs', 0.3, ...
        'photonTotals', [500 2000 8000 20000 50000], ...
        'amplitudeSets', [0.50 0.30 0.20; ...
                          0.34 0.33 0.33; ...
                          0.20 0.30 0.50; ...
                          0.60 0.30 0.10], ...
        'pulsePeriodNs', 12.5, 'dtNs', 0.16, 'irfFwhmNs', 0.35, ...
        'membraneTauCount', 24, 'membraneTauBoundsNs', [0.4 5.5], ...
        'slbCountBiasFraction', 0, ...
        'placeTruthOnGrid', false, ...
        'slbCountPriorNodes', 0, ...
        'exampleTaus', [0.8 3.0; 1.5 4.5], ...
        'examplePhotonTotals', [500 8000 50000], ...
        'exampleAmplitudeSet', 2, ...
        'outputDir', pwd, 'seed', 42);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if any(abs(sum(opts.amplitudeSets, 2) - 1) > 1e-9)
        error('study_flim_triexp_recovery:Amplitudes', ...
            'Each amplitudeSets row must sum to 1 (photon fractions).');
    end
    if opts.membraneTauBoundsNs(1) > opts.tau2Range(1) || ...
            opts.membraneTauBoundsNs(2) < opts.tau3Range(2)
        warning('study_flim_triexp_recovery:GridSpan', ...
            ['The fit grid %.3g-%.3g ns does not cover the truth sweep ' ...
             '%.3g-%.3g ns. Lifetimes outside the grid cannot be recovered ' ...
             'at any photon count.'], opts.membraneTauBoundsNs(1), ...
            opts.membraneTauBoundsNs(2), opts.tau2Range(1), ...
            opts.tau3Range(2));
    end
end

% --------------------------------------------------------------- simulation

function irf = gaussianIrf(timeNs, fwhmNs)
    sigma = fwhmNs / (2 * sqrt(2 * log(2)));
    centre = max(4 * sigma, 3 * (timeNs(2) - timeNs(1)));
    irf = exp(-0.5 * ((timeNs - centre) / sigma).^2);
    irf = irf / sum(irf);
end

function shape = mixtureShape(irf, timeNs, periodNs, taus, fractions)
    % Photon-fraction mixture of periodic single-exponential decays, each
    % convolved with the IRF. Periodicity matters here: with a 12.5 ns period
    % a 5 ns component still carries signal into the next pulse.
    shape = zeros(1, numel(timeNs));
    for k = 1:numel(taus)
        if fractions(k) <= 0
            continue;
        end
        shape = shape + fractions(k) * ...
            periodicDecay(irf, timeNs, periodNs, taus(k));
    end
    shape = shape + 1e-12;      % strictly positive, so sampling edges rise
    shape = shape / sum(shape);
end

function pattern = periodicDecay(irf, timeNs, periodNs, tauNs)
    nBins = numel(timeNs);
    decay = zeros(1, nBins);
    for repeat = 0:3            % include the wrap-around tail
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

% ----------------------------------------------------------------- analysis

function step = localGridStep(fitGrid, tau)
    if numel(fitGrid) < 2
        step = NaN;
        return;
    end
    [~, nearest] = min(abs(fitGrid - tau));
    if nearest == 1
        step = fitGrid(2) - fitGrid(1);
    elseif nearest == numel(fitGrid)
        step = fitGrid(end) - fitGrid(end-1);
    else
        step = (fitGrid(nearest+1) - fitGrid(nearest-1)) / 2;
    end
end

function value = snapToGrid(value, grid)
    shape = size(value);
    flat = value(:);
    for k = 1:numel(flat)
        [~, nearest] = min(abs(grid - flat(k)));
        flat(k) = grid(nearest);
    end
    value = reshape(flat, shape);
end

function value = iqrOmitNan(x)
    x = x(isfinite(x));
    if isempty(x)
        value = NaN;
    else
        value = diff(prctile(x, [25 75]));
    end
end

function reportVerdict(results, fitGrid)
    summary = results.summary;
    opts = results.opts;
    limit2 = 0.5 * localGridStep(fitGrid, mean(opts.tau2Range));
    limit3 = 0.5 * localGridStep(fitGrid, mean(opts.tau3Range));
    fprintf('\n  --- how to read this ---\n');
    for k = 1:height(summary)
        fprintf('  set %d N=%6d : tau2 %s, tau3 %s\n', ...
            summary.amplitudeSet(k), summary.photonTotal(k), ...
            limitedByText(summary.medianAbsError2Ns(k), limit2), ...
            limitedByText(summary.medianAbsError3Ns(k), limit3));
    end
    fprintf(['  grid-limited means the error is already at half the local ' ...
             'grid step\n  (%.3f ns for tau2, %.3f ns for tau3), so a finer ' ...
             'membraneTauCount\n  would be needed to measure the method ' ...
             'any further.\n'], limit2, limit3);
end

function text = limitedByText(errorValue, gridLimit)
    if ~isfinite(errorValue)
        text = 'not recovered';
    elseif errorValue <= gridLimit
        text = 'grid-limited';
    else
        text = 'noise-limited';
    end
end

% ----------------------------------------------------------- example decays

function examples = captureExamples(opts, tau2Grid, tau3Grid, Y, out, ...
        tau2Fit, tau3Fit, modelMap)
    examples = struct([]);
    for e = 1:size(opts.exampleTaus, 1)
        [~, ix] = min(abs(tau2Grid - opts.exampleTaus(e, 1)));
        [~, iy] = min(abs(tau3Grid - opts.exampleTaus(e, 2)));
        entry = struct( ...
            'tau2True', tau2Grid(ix), 'tau3True', tau3Grid(iy), ...
            'counts', reshape(double(Y(ix, iy, :)), 1, []), ...
            'tau2Fit', tau2Fit(ix, iy), 'tau3Fit', tau3Fit(ix, iy), ...
            'slbFraction', ...
                double(out.twoMembrane.fixedSlbPhotonFraction(ix, iy)), ...
            'comp2Fraction', ...
                double(out.twoMembrane.membrane1PhotonFraction(ix, iy)), ...
            'comp3Fraction', ...
                double(out.twoMembrane.membrane2PhotonFraction(ix, iy)), ...
            'backgroundFraction', ...
                double(out.twoMembrane.backgroundFraction(ix, iy)), ...
            'selectedModel', modelMap(ix, iy));
        if isempty(examples)
            examples = entry;
        else
            examples(end+1) = entry; %#ok<AGROW>
        end
    end
end

function files = plotExampleDecays(results)
%PLOTEXAMPLEDECAYS Presentation and paper figure: data, fit and components.
%
% One column per example (tau2, tau3) pair and one pair of rows per photon
% budget: the decay on a log scale with the fitted total and the three fitted
% components, and beneath it the Pearson residual. Saved as PNG for slides
% and PDF for a paper, because the PDF stays vector.
    opts = results.opts;
    setIndex = min(max(1, round(opts.exampleAmplitudeSet)), ...
        size(opts.amplitudeSets, 1));

    [~, levelIndices] = ismember(opts.examplePhotonTotals, opts.photonTotals);
    levelIndices = levelIndices(levelIndices > 0);
    if isempty(levelIndices)
        levelIndices = 1:numel(opts.photonTotals);
    end
    nLevels = numel(levelIndices);
    nExamples = size(opts.exampleTaus, 1);

    timeNs = results.timeNs;
    irf = results.irf;
    period = opts.pulsePeriodNs;

    h = figure('Color', 'w', 'Visible', 'off', ...
        'Position', [60 60 420 * nExamples, 360 * nLevels]);
    % Three rows per photon level: the decay spans the first two so it gets
    % twice the height of its residual strip, which a single-row decay panel
    % is too short for - MATLAB drops all but one decade label.
    layout = tiledlayout(h, 3 * nLevels, nExamples, ...
        'Padding', 'compact', 'TileSpacing', 'tight');

    for levelSlot = 1:nLevels
        totalIndex = levelIndices(levelSlot);
        photonTotal = opts.photonTotals(totalIndex);
        examples = results.examples{setIndex, totalIndex};
        for e = 1:nExamples
            ex = examples(e);
            counts = ex.counts;

            % Rebuild the fitted decay from the recovered parameters, so the
            % curve drawn is the model actually fitted rather than the truth.
            fractions = [ex.slbFraction, ex.comp2Fraction, ex.comp3Fraction];
            taus = [opts.tauSlbNs, ex.tau2Fit, ex.tau3Fit];
            components = zeros(3, numel(timeNs));
            for k = 1:3
                if isfinite(taus(k)) && taus(k) > 0 && ...
                        isfinite(fractions(k)) && fractions(k) > 0
                    components(k, :) = fractions(k) * ...
                        periodicDecay(irf, timeNs, period, taus(k));
                end
            end
            background = zeros(1, numel(timeNs));
            if isfinite(ex.backgroundFraction) && ex.backgroundFraction > 0
                background(:) = ex.backgroundFraction / numel(timeNs);
            end
            model = sum(components, 1) + background;
            scale = sum(model);
            if scale > 0
                components = components / scale;
                model = model / scale;
            end
            expected = photonTotal * model;

            % --- decay panel ---
            ax = nexttile(layout, (3 * levelSlot - 3) * nExamples + e, ...
                [2 1]);
            hold(ax, 'on');
            stairs(ax, timeNs, max(counts, 0.1), 'Color', [0.45 0.45 0.45], ...
                'LineWidth', 0.6, 'DisplayName', 'simulated data');
            plot(ax, timeNs, max(expected, 1e-3), 'k-', 'LineWidth', 1.6, ...
                'DisplayName', 'fitted total');
            colours = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
            names = {sprintf('SLB %.2g ns', opts.tauSlbNs), ...
                sprintf('c2 %.2f ns', ex.tau2Fit), ...
                sprintf('c3 %.2f ns', ex.tau3Fit)};
            for k = 1:3
                plot(ax, timeNs, max(photonTotal * components(k, :), 1e-3), ...
                    '-', 'Color', colours(k, :), 'LineWidth', 1.1, ...
                    'DisplayName', names{k});
            end
            set(ax, 'YScale', 'log');
            upperLimit = max(4, 2 * max(counts));
            ylim(ax, [0.5 upperLimit]);
            % Force one tick per decade; otherwise a short panel silently
            % renders with a single label and the scale is unreadable.
            set(ax, 'YTick', 10 .^ (0:ceil(log10(upperLimit))));
            xlim(ax, [0 timeNs(end)]);
            grid(ax, 'on');
            ax.GridAlpha = 0.12;
            ylabel(ax, 'photons / bin', 'FontSize', 8);
            title(ax, sprintf(['N=%d | true %.2f / %.2f ns | fitted ' ...
                '%.2f / %.2f ns | M%d'], photonTotal, ex.tau2True, ...
                ex.tau3True, ex.tau2Fit, ex.tau3Fit, ex.selectedModel), ...
                'FontSize', 8);
            if levelSlot == 1 && e == 1
                legend(ax, 'Location', 'northeast', 'FontSize', 6, ...
                    'Box', 'off');
            end
            set(ax, 'XTickLabel', []);
            hold(ax, 'off');

            % --- residual panel ---
            ax = nexttile(layout, (3 * levelSlot - 1) * nExamples + e);
            residual = (counts - expected) ./ sqrt(max(expected, 1e-6));
            stem(ax, timeNs, residual, 'Marker', 'none', ...
                'Color', [0.2 0.2 0.2], 'LineWidth', 0.5);
            yline(ax, 0, 'k-');
            yline(ax, 2, 'r:');
            yline(ax, -2, 'r:');
            xlim(ax, [0 timeNs(end)]);
            ylim(ax, [-4 4]);
            set(ax, 'YTick', [-2 0 2]);
            grid(ax, 'on');
            ax.GridAlpha = 0.12;
            ylabel(ax, 'Pearson res.', 'FontSize', 8);
            if levelSlot == nLevels
                xlabel(ax, 'time [ns]', 'FontSize', 8);
            end
        end
    end

    fractionsRow = opts.amplitudeSets(setIndex, :);
    title(layout, sprintf(['Simulated triexponential decays and fixed-SLB ' ...
        'fits | photon fractions %.2f / %.2f / %.2f'], fractionsRow(1), ...
        fractionsRow(2), fractionsRow(3)), 'FontWeight', 'bold');
    subtitle(layout, ['grey = simulated counts, black = fitted total, ' ...
        'coloured = fitted components; dotted lines at +/- 2 sigma'], ...
        'FontSize', 8);

    pngFile = fullfile(opts.outputDir, 'flim_triexp_example_decays.png');
    pdfFile = fullfile(opts.outputDir, 'flim_triexp_example_decays.pdf');
    exportgraphics(h, pngFile, 'Resolution', 300);
    exportgraphics(h, pdfFile, 'ContentType', 'vector');
    close(h);
    files = {pngFile, pdfFile};
end

% ------------------------------------------------------------------ figures

function files = plotRecovery(results)
    opts = results.opts;
    setCount = size(opts.amplitudeSets, 1);
    totalCount = numel(opts.photonTotals);
    files = {};

    for setIndex = 1:setCount
        h = figure('Color', 'w', 'Visible', 'off', ...
            'Position', [60 60 1160 260 * totalCount]);
        layout = tiledlayout(h, totalCount, 3, 'Padding', 'compact', ...
            'TileSpacing', 'compact');
        for totalIndex = 1:totalCount
            g = results.grids{setIndex, totalIndex};

            ax = nexttile(layout);
            drawMap(ax, g.error2, [-0.5 0.5], 'tau2 error [ns]', opts);
            title(ax, sprintf('N=%d | tau2 error, median |.| %.3f ns', ...
                g.photonTotal, median(abs(g.error2(:)), 'omitnan')), ...
                'FontSize', 8);

            ax = nexttile(layout);
            drawMap(ax, g.error3, [-1.5 1.5], 'tau3 error [ns]', opts);
            title(ax, sprintf('tau3 error, median |.| %.3f ns', ...
                median(abs(g.error3(:)), 'omitnan')), 'FontSize', 8);

            ax = nexttile(layout);
            drawMap(ax, g.modelMap, [1 3], 'selected model', opts);
            title(ax, sprintf('model MAP, 3 chosen in %.1f%%', ...
                100 * mean(g.modelMap(:) == 3)), 'FontSize', 8);
        end
        fractions = opts.amplitudeSets(setIndex, :);
        title(layout, sprintf(['Triexponential recovery | photon fractions ' ...
            'SLB %.2f / c2 %.2f / c3 %.2f'], fractions(1), fractions(2), ...
            fractions(3)), 'FontWeight', 'bold');
        subtitle(layout, sprintf(['x = true tau2 %.2g-%.2g ns, ' ...
            'y = true tau3 %.2g-%.2g ns; SLB fixed at %.2g ns'], ...
            opts.tau2Range(1), opts.tau2Range(2), opts.tau3Range(1), ...
            opts.tau3Range(2), opts.tauSlbNs), 'FontSize', 8);

        name = fullfile(opts.outputDir, ...
            sprintf('flim_triexp_recovery_set%d.png', setIndex));
        exportgraphics(h, name, 'Resolution', 200);
        close(h);
        files{end+1} = name; %#ok<AGROW>
    end
end

function drawMap(ax, data, limits, barLabel, opts)
    imagesc(ax, opts.tau2Range, opts.tau3Range, double(data).');
    axis(ax, 'tight');
    set(ax, 'YDir', 'normal');
    colormap(ax, turbo(256));
    clim(ax, limits);
    bar = colorbar(ax);
    bar.Label.String = barLabel;
    xlabel(ax, 'true tau2 [ns]', 'FontSize', 7);
    ylabel(ax, 'true tau3 [ns]', 'FontSize', 7);
end
