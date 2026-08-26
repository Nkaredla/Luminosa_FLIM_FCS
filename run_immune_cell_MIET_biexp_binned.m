function out = run_immune_cell_MIET_biexp_binned(dataRoot, binSizes, opts)
%RUN_IMMUNE_CELL_MIET_BIEXP_BINNED Soft-SLB biexponential at several bin sizes.
%
% out = run_immune_cell_MIET_biexp_binned()
% out = run_immune_cell_MIET_biexp_binned(dataRoot, binSizes, opts)
%
% Fits every unique acquisition below DATAROOT at each requested spatial binning
% (default [1 2 4], sliding with step 1) using variable projection - BFGS over
% [log tau1, log tau2] with the amplitudes profiled out by whitened Poisson IRLS.
% Results go to <acquisition>/biexp_slb, _bin2, _bin4, so nothing collides.
%
% RESUMABLE. Each completed (acquisition, binSize) pair is appended to
% biexp_binned_done.txt at the data root and skipped on a later call. That
% matters on this machine because MATLAB intermittently aborts AT STARTUP with
%
%     Assertion failed: Graphics system timed out searching for display devices
%
% which is a display-enumeration timeout in MATLAB's graphics init, nothing to do
% with the fitting. Just run this again after such a failure and it continues
% where it stopped. Pass opts.restart = true to ignore the log and refit
% everything.
%
% ONLY THE UNIQUE ACQUISITIONS ARE FITTED. Each acquisition here has two saved
% analyses (immune_cell_MIET and immune_cell_MIET_grid48/cfg_...), and their
% tcspc_pix cubes were verified SHA-256 identical for all six acquisitions - the
% configs differ only in stored metadata. Fitting both would double the runtime
% for identical numbers, so the grid48 copies are skipped. Pass
% opts.includeDuplicates = true to fit them anyway.
%
% Cost, measured: about 2.0 ms/pixel, so roughly 3-4 minutes per acquisition at
% 1x1. Binned stages are slower only because binning lifts more pixels above the
% photon floor.

    if nargin < 1 || isempty(dataRoot)
        dataRoot = ['D:\Luminosa\Data\260813\' ...
            'RT_Jurkat_CD58_memglow_20260813_1'];
    end
    if nargin < 2 || isempty(binSizes); binSizes = [1 2 4]; end
    if nargin < 3 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'restart'); opts.restart = false; end
    if ~isfield(opts, 'includeDuplicates'); opts.includeDuplicates = false; end

    doneFile = fullfile(dataRoot, 'biexp_binned_done.txt');
    done = {};
    if ~opts.restart && isfile(doneFile)
        text = fileread(doneFile);
        done = strtrim(strsplit(text, newline));
        done = done(~cellfun(@isempty, done));
    end

    found = dir(fullfile(dataRoot, '**', ...
        'immune_cell_MIET_640nm_red_analysis.mat'));
    if isempty(found)
        error('run_immune_cell_MIET_biexp_binned:NoAnalyses', ...
            'No analysis MAT below %s', dataRoot);
    end
    paths = {};
    for k = 1:numel(found)
        isDuplicate = contains(found(k).folder, 'grid48');
        if isDuplicate && ~opts.includeDuplicates; continue; end
        paths{end + 1} = fullfile(found(k).folder, found(k).name); %#ok<AGROW>
    end
    paths = sort(paths);

    fprintf(['\nrun_immune_cell_MIET_biexp_binned\n  %d acquisition(s) x ' ...
        '%d binning(s) below\n  %s\n'], numel(paths), numel(binSizes), dataRoot);
    if ~isempty(done)
        fprintf('  resuming: %d pair(s) already done\n', numel(done));
    end

    rows = struct([]);
    for b = binSizes(:)'
        for k = 1:numel(paths)
            tag = sprintf('%s|bin%d', paths{k}, b);
            if any(strcmp(done, tag))
                fprintf('\n### SKIP bin%d (done): %s\n', b, paths{k});
                continue;
            end
            fprintf('\n### bin%d : %s\n', b, paths{k});
            entry = struct('analysisMat', paths{k}, 'binSize', b, ...
                'pixelsFitted', 0, 'tau1MedianNs', NaN, ...
                'tau2MedianNs', NaN, 'tau2Q1Ns', NaN, 'tau2Q3Ns', NaN, ...
                'photonShare2Median', NaN, 'tauMeanMedianNs', NaN, ...
                'reducedDevianceMedian', NaN, 'convergedFraction', NaN, ...
                'status', 'ok');
            try
                o = opts;
                o.binSize = b;
                r = immune_cell_MIET_biexp_slb(paths{k}, o);
                v = isfinite(r.maps.tau1Ns);
                entry.pixelsFitted = nnz(v);
                entry.tau1MedianNs = median(r.maps.tau1Ns(v));
                entry.tau2MedianNs = median(r.maps.tau2Ns(v));
                entry.tau2Q1Ns = quantileLocalBiexp(r.maps.tau2Ns(v), 0.25);
                entry.tau2Q3Ns = quantileLocalBiexp(r.maps.tau2Ns(v), 0.75);
                entry.photonShare2Median = median(r.maps.photonFraction2(v));
                entry.tauMeanMedianNs = median(r.maps.tauMeanNs(v));
                entry.reducedDevianceMedian = median(r.maps.reducedDeviance(v));
                if isfield(r, 'convergedFraction')
                    entry.convergedFraction = r.convergedFraction;
                end
                fid = fopen(doneFile, 'a');
                fprintf(fid, '%s\n', tag);
                fclose(fid);
            catch runError
                fprintf('  FAILED: %s\n', runError.message);
                entry.status = runError.identifier;
            end
            if isempty(rows); rows = entry; else; rows(end + 1) = entry; end %#ok<AGROW>
        end
    end

    if isempty(rows)
        fprintf('\n  nothing to do - every pair was already complete.\n');
        out = struct('summary', table(), 'dataRoot', dataRoot);
        return;
    end

    summary = struct2table(rows);
    csvFile = fullfile(dataRoot, 'biexp_binned_summary.csv');
    writetable(summary, csvFile);

    fprintf('\n=============================================================\n');
    fprintf('BINNED SUMMARY: soft-SLB biexponential, variable projection\n\n');
    fprintf(['  bin   pixels   tau1(ns)  tau2(ns)  tau2 IQR         share2  ' ...
        'meanTau  redDev  acquisition\n']);
    fprintf(['  --------------------------------------------------------' ...
        '--------------------------------\n']);
    for k = 1:height(summary)
        [~, acq] = fileparts(fileparts(fileparts(summary.analysisMat{k})));
        fprintf(['  %3d %8d  %8.4f  %8.4f  [%.3f %.3f]  %6.3f  %7.4f  ' ...
            '%6.3f  %s\n'], summary.binSize(k), summary.pixelsFitted(k), ...
            summary.tau1MedianNs(k), summary.tau2MedianNs(k), ...
            summary.tau2Q1Ns(k), summary.tau2Q3Ns(k), ...
            summary.photonShare2Median(k), summary.tauMeanMedianNs(k), ...
            summary.reducedDevianceMedian(k), acq);
    end
    fprintf('\n  wrote %s\n', csvFile);
    out = struct('summary', summary, 'dataRoot', dataRoot, ...
        'csvFile', csvFile, 'doneFile', doneFile);
end
