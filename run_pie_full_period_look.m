function out = run_pie_full_period_look(source, opts)
%RUN_PIE_FULL_PERIOD_LOOK What is actually in the decay, both PIE colours at once.
%
% out = run_pie_full_period_look(acquisitionFolderOrPtu)
% out = run_pie_full_period_look(source, opts)
%
% Re-reads the raw PTU to get the UNGATED full-period TCSPC, then fits a
% lifetime distribution for each excitation pulse simultaneously and reports
% every lifetime it finds.
%
% WHY THE RAW FILE HAS TO BE RE-READ
%
% The saved analyses kept only the 640 gate. result.gateDetection.profile spans
% 314 bins and looks like a full period, but it is zero from 0 to 23.7 ns - the
% 485 photons were discarded at read time. Nothing downstream can recover them.
%
% WHAT THIS IS FOR
%
% Deciding, from the data rather than by assumption, which lifetime is the SLB.
% Two things are visible only over the full period:
%
%   - the 485-excited emission, which has not finished decaying when the 640
%     gate opens. Inside that gate it looks like a constant offset, so a gated
%     model absorbs it into the background term and every amplitude is
%     corrupted. Over the full period it has its own rise 25 ns earlier and is
%     separable.
%   - whether the shortest component is a lifetime at all. A component at or
%     below one bin width, arriving with the pulse, is what scattered excitation
%     light and short-lived substrate photoluminescence look like - not the
%     bilayer. The distribution shows it next to the IRF instead of quietly
%     reporting it as the SLB lifetime.
%
% NOTHING IS DECIDED HERE. The output lists every mode with its photon count so
% the choice of which lifetime to call the SLB can be made on the evidence.
%
% opts fields
%   binNs        output bin width (default 0.16, matching the analyses)
%   nNodes       lifetime nodes per pulse (default 60)
%   discreteK    also fit this many discrete components per pulse (default 3)
%   maxPhotons   cap the read for a quick look (default Inf)
%   outputRoot   collect results under one folder at the data root, named by
%                acquisition (default: the parent of the acquisition folder,
%                <dataRoot>/biexp_slb_anchored). Everything for a session then
%                lives in one place instead of being scattered through the data
%                tree.
%   outputDir    override the destination completely

    if nargin < 2 || isempty(opts); opts = struct(); end
    defaults = struct('binNs', 0.16, 'nNodes', 60, 'discreteK', 3, ...
        'maxPhotons', Inf, 'outputDir', '', 'outputRoot', '');
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end

    % ---- locate the PTU and the matching analysis -------------------------
    source = char(source);
    if isfolder(source)
        p = dir(fullfile(source, '**', '*.ptu'));
        if isempty(p)
            error('run_pie_full_period_look:NoPTU', ...
                'No .ptu below %s', source);
        end
        ptuFile = fullfile(p(1).folder, p(1).name);
    else
        ptuFile = source;
    end
    acqFolder = fileparts(ptuFile);
    analysis = dir(fullfile(acqFolder, '**', ...
        'immune_cell_MIET_640nm_red_analysis.mat'));
    if isempty(analysis)
        error('run_pie_full_period_look:NoAnalysis', ...
            'No analysis MAT below %s (needed for the IRF and gate)', acqFolder);
    end
    analysisMat = fullfile(analysis(1).folder, analysis(1).name);

    L = load(analysisMat, 'result');
    result = L.result;
    gateIrf = double(immune_cell_MIET_explorer_field(result, 'irf.curve'));
    gateIrf = max(gateIrf(:), 0);
    gateStart = double(immune_cell_MIET_explorer_field(result, ...
        'channel.gateStartBin'));

    % ---- the full-period decay -------------------------------------------
    ex = extract_full_period_tcspc(ptuFile, struct('binNs', opts.binNs, ...
        'maxPhotons', opts.maxPhotons));
    decay = ex.decay;
    dtNs = ex.dtNs;
    periodNs = ex.periodNs;
    nBin = ex.nBin;

    % ---- place the gate IRF into full-period coordinates ------------------
    % The stored IRF covers the 640 gate only, starting at gateStartBin.
    irfFull = zeros(nBin, 1);
    lo = max(1, round(gateStart));
    hi = min(nBin, lo + numel(gateIrf) - 1);
    irfFull(lo:hi) = gateIrf(1:(hi - lo + 1));
    if sum(irfFull) <= 0
        error('run_pie_full_period_look:BadIrf', ...
            'The gate IRF did not land inside the period.');
    end

    % ---- find both excitation pulses --------------------------------------
    % The 640 pulse is the biggest rise; the 485 pulse is half a period away.
    d = diff([decay(end); decay]);
    [~, redBin] = max(d);
    sepBins = round((periodNs / 2) / dtNs);
    blueBin = mod(redBin - sepBins - 1, nBin) + 1;
    fprintf(['  pulses: 640 at bin %d (%.2f ns), 485 at bin %d (%.2f ns), ' ...
        'separation %.2f ns\n'], redBin, (redBin - 1) * dtNs, blueBin, ...
        (blueBin - 1) * dtNs, sepBins * dtNs);
    % Align the IRF's own peak to the detected 640 pulse.
    [~, irfPeak] = max(irfFull);
    irfFull = circshift(irfFull, redBin - irfPeak);

    % ---- fit -------------------------------------------------------------
    fit = fit_pie_lifetime_distribution(decay, irfFull, dtNs, periodNs, ...
        [redBin, blueBin], struct('nNodes', opts.nNodes, ...
        'pulseNames', {{'640 (red)', '485 (blue)'}}, ...
        'discreteK', opts.discreteK));

    if isempty(opts.outputDir)
        if isempty(opts.outputRoot)
            opts.outputRoot = fullfile(fileparts(acqFolder), ...
                'biexp_slb_anchored');
        end
        [~, acqName] = fileparts(acqFolder);
        opts.outputDir = fullfile(opts.outputRoot, ...
            sprintf('%s_pie_full_period', acqName));
    end
    if ~isfolder(opts.outputDir); mkdir(opts.outputDir); end

    out = fit;
    out.extraction = ex;
    out.ptuFile = ptuFile;
    out.analysisMat = analysisMat;
    out.pulseBins = [redBin, blueBin];
    out.irfFull = irfFull;
    out.outputDir = opts.outputDir;
    save(fullfile(opts.outputDir, 'pie_full_period_fit.mat'), 'out', '-v7.3');
    fprintf('  wrote %s\n', fullfile(opts.outputDir, ...
        'pie_full_period_fit.mat'));
    out.figure = pie_full_period_figure(out, opts.outputDir);
    fprintf('  wrote %s\n', out.figure);
end
