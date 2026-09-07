function result = PTU_MitoLysoCytosol_TailFLIM(input, cfg)
%PTU_MITOLYSOCYTOSOL_TAILFLIM Segment B12 compartments and fit only decay tails.
% input is a PTU path or an existing result from the main B12 analysis.
% Each curve is cut independently at its TCSPC peak + cfg.tailOffsetNs (1 ns).
% Only retained photons count: <10 skips fitting; 10-499 allows one exponential;
% >=500 allows one/two model selection. No IRF estimation or Bayesian maps run.
% Fractions/amplitudes describe the retained tail, not the original full decay.
%
% result = PTU_MitoLysoCytosol_TailFLIM('D:\data\RawImage.ptu');
% tailResult = PTU_MitoLysoCytosol_TailFLIM(existingResult, ...
%     struct('saveOutputs', false, 'showFigures', false));

    if nargin < 2 || isempty(cfg)
        cfg = struct();
    end
    helperFolder = fileparts(mfilename('fullpath'));
    parentFolder = fileparts(helperFolder);
    originalPath = path;
    pathCleanup = onCleanup(@() path(originalPath));
    addpath(parentFolder, helperFolder, '-begin');
    if isstruct(input) && isfield(input, 'config')
        inherited = input.config;
        fields = fieldnames(cfg);
        for k = 1:numel(fields)
            inherited.(fields{k}) = cfg.(fields{k});
        end
        cfg = inherited;
    end
    cfg = defaultConfig(cfg);
    if ~isfield(cfg, 'tailOffsetNs') || isempty(cfg.tailOffsetNs)
        cfg.tailOffsetNs = 1;
    end
    validateattributes(cfg.tailOffsetNs, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    if strcmp(cfg.outputStem, 'B12_MitoLysoCytosol_BayesFLIM')
        cfg.outputStem = 'B12_MitoLysoCytosol_TailFLIM';
    end
    cfg.segmentOnly = true;
    if isstruct(input)
        source = input;
    else
        source = PTU_MitoLysoCytosol_BayesFLIM(input, cfg);
    end
    % Carry segmentation and raw histograms, without stale fits or posteriors.
    result = struct();
    carry = {'file', 'head', 'gates', 'intensity', 'masks', 'labels', ...
        'tcspc', 'segmentation', 'regions', 'channelMapping', 'detectorIdsRaw'};
    for k = 1:numel(carry)
        if isfield(source, carry{k})
            result.(carry{k}) = source.(carry{k});
        end
    end
    result.config = cfg;
    result.analysis = 'Independent segmented TCSPC tail fits; no reconvolution';
    dtNs = result.gates.dtNs;
    names = {'globalB12', 'cytosol', 'pooledMitochondria', ...
        'pooledLysosome', 'mitochondria', 'lysosome'};
    curves = {result.tcspc.globalB12, result.tcspc.cytosol, ...
        sum(result.tcspc.mitochondria, 2), sum(result.tcspc.lysosome, 2), ...
        result.tcspc.mitochondria, result.tcspc.lysosome};
    for group = 1:numel(names)
        counts = curves{group};
        [ok, order, tau, tauMean, coeff, fraction, fitted, metrics, messages] = ...
            fitRegions(counts, names{group}, [], dtNs, cfg.maxExp, ...
            cfg.tauSeeds1, cfg.tauSeeds2, cfg.tauBoundsNs, ...
            cfg.includeBackground, cfg.minPhotonsPerRegionFit, cfg.fitMaxIter, ...
            cfg.bicImprovementForTwoExp, cfg.minStateSeparationFraction, ...
            cfg.minStateAmplitudeFraction, cfg.minPhotonsForTwoExp, cfg.tailOffsetNs);
        [~, peakBins] = max(counts, [], 1);
        startBins = min(size(counts, 1) + 1, peakBins + ceil(cfg.tailOffsetNs / dtNs));
        startTimes = (startBins(:) - 0.5) * dtNs;
        startTimes(startBins(:) > size(counts, 1)) = NaN;
        result.fits.(names{group}) = struct('ok', ok, 'modelOrder', order, ...
            'tauNs', tau, 'tauMeanNs', tauMean, ...
            'coefficientsBgComp1Comp2', coeff, 'speciesFraction', fraction, ...
            'fitCounts', fitted, 'metrics', metrics, 'messages', {messages}, ...
            'peakBin', peakBins(:), 'firstFitBin', startBins(:), ...
            'firstFitTimeNs', startTimes, 'totalPhotonCount', sum(counts, 1).', ...
            'retainedPhotonCount', metrics(:, 4));
    end
    result.fits.metricColumns = {'deviance', 'BIC', 'reducedChiSquare', ...
        'photonCount', 'oneExpBIC', 'twoExpBIC', 'oneExpChi2', 'twoExpChi2'};

    % Map connected-region estimates only. Ineligible regions remain NaN/0.
    tauMap = nan(size(result.masks.cell), 'single');
    orderMap = zeros(size(tauMap), 'uint8');
    for group = 2:numel(names)
        if group == 3 || group == 4
            continue; % Pooled fits are diagnostics only.
        end
        f = result.fits.(names{group});
        for region = 1:numel(f.ok)
            if ~f.ok(region)
                continue;
            end
            if group == 2
                mask = result.masks.cytosol;
            else
                mask = result.labels.(names{group}) == region;
            end
            tauMap(mask) = f.tauMeanNs(region);
            orderMap(mask) = uint8(f.modelOrder(region));
        end
    end
    result.flim = struct('tauMeanNs', tauMap, 'modelOrder', orderMap, ...
        'compartment', result.masks.compartmentMap);

    figures = gobjects(0);
    if cfg.showFigures || cfg.saveOutputs
        visibility = 'off';
        if cfg.showFigures
            visibility = 'on';
        end
        figures(1) = showSegmentation(result.intensity.mitochondria640Det2, ...
            result.intensity.b12_560Det1, result.masks.mitochondria, ...
            result.masks.lysosome, result.masks.cellBoundary, ...
            result.masks.compartmentMap, size(curves{5}, 2), size(curves{6}, 2), visibility);
        plotNames = {};
        plotCurves = {};
        plotFits = {};
        plotOrders = [];
        plotTau = zeros(0, 2);
        plotMessages = {};
        for group = [2 5 6]
            f = result.fits.(names{group});
            [~, ranking] = sort(f.totalPhotonCount, 'descend');
            for k = 1:min(4, numel(ranking))
                region = ranking(k);
                plotNames{end+1} = sprintf('%s %d tail', names{group}, region); %#ok<AGROW>
                plotCurves{end+1} = curves{group}(:, region); %#ok<AGROW>
                plotFits{end+1} = f.fitCounts(:, region); %#ok<AGROW>
                plotOrders(end+1, 1) = f.modelOrder(region); %#ok<AGROW>
                plotTau(end+1, :) = f.tauNs(region, :); %#ok<AGROW>
                plotMessages{end+1, 1} = f.messages{region}; %#ok<AGROW>
            end
        end
        figures(2) = showFits(result.tcspc.timeNs, plotNames, plotCurves, ...
            plotFits, plotOrders, plotTau, plotMessages, cfg.maxFitPlots, visibility);
        figures(3) = figure('Name', 'B12 regional tail lifetime', ...
            'Color', 'w', 'Visible', visibility);
        tiledlayout(1, 2);
        nexttile; showMap(tauMap, 'Regional tail lifetime', 'ns', turbo(256));
        nexttile; imagesc(orderMap'); axis image off; caxis([0 2]); colorbar;
        title('Tail model order (0 = not fitted)');
        figures(4) = showLifetimeHistograms(result.fits.mitochondria.tauNs, ...
            result.fits.mitochondria.modelOrder, result.fits.lysosome.tauNs, ...
            result.fits.lysosome.modelOrder, cfg.tauBoundsNs, ...
            cfg.lifetimeHistogramBins, visibility);
    end
    result.savedOutputs = struct();
    if cfg.saveOutputs
        outputDirectory = cfg.outputDir;
        if isempty(outputDirectory)
            outputDirectory = fullfile(fileparts(result.file), 'B12_TailFLIM');
        end
        if ~isfolder(outputDirectory)
            mkdir(outputDirectory);
        end
        matPath = fullfile(outputDirectory, [cfg.outputStem '.mat']);
        suffixes = {'segmentation', 'tail_fits', 'tail_map', 'lifetime_histograms'};
        figurePaths = cell(numel(figures), 1);
        for k = 1:numel(figures)
            figurePaths{k} = fullfile(outputDirectory, ...
                sprintf('%s_%s.png', cfg.outputStem, suffixes{k}));
            exportgraphics(figures(k), figurePaths{k}, 'Resolution', 200);
        end
        result.savedOutputs = struct('directory', outputDirectory, ...
            'mat', matPath, 'figures', {figurePaths});
        save(matPath, 'result', '-v7.3');
    end
    if ~cfg.showFigures
        close(figures);
    end
    fprintf('B12 compartment tail-only fitting complete.\n');
end
