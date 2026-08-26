function out = run_immune_cell_MIET_slb_anchored(dataRoot, binSizes, opts)
%RUN_IMMUNE_CELL_MIET_SLB_ANCHORED Measure the bare SLB, then fit the cells to it.
%
% out = run_immune_cell_MIET_slb_anchored()
% out = run_immune_cell_MIET_slb_anchored(dataRoot, binSizes, opts)
%
% Two stages per acquisition:
%
%   1. measure_slb_reference_lifetime on the bare SLB OUTSIDE the cell, which
%      returns a prior centre and width measured from the data. The centre is the
%      SHORT component of a two-component fit of the pooled bare-SLB decay, so
%      the bilayer is separated from the long-lived contamination in that region
%      rather than averaged with it.
%   2. immune_cell_MIET_biexp_slb at each requested binning (default [1 2 4]),
%      with tau1 pinned to that measured centre and a small measured width,
%      over EVERY pixel in the file (opts.pixelMask, default 'valid') rather
%      than only the cell footprint.
%
% WHY FIT THE WHOLE FILE AND NOT JUST THE CELL
%
% Because it makes the anchor self-checking. The bare-SLB region is included in
% the fit, and there the answer is known in advance: tau1 must come back at the
% anchor value and the long component must carry almost no photons. If it does
% not, the anchor is wrong or the two-component model does not hold, and that is
% visible in the same maps rather than needing a separate experiment. Restricting
% the fit to the cell footprint throws that control away.
%
% It also puts cell and bilayer on one continuous set of maps, so the step in
% lifetime at the cell edge can be read directly instead of being compared across
% two separately masked analyses.
%
% WHY ANCHOR RATHER THAN USE THE STORED REFERENCE
%
% The bilayer sits on a spacer above the metal and its height does not change
% because a cell landed on it, so the SLB lifetime under the cell must equal the
% SLB lifetime beside it. That makes the bare-SLB region a direct measurement of
% the prior, and it is a far better anchor than result.slbReference
% .fixedLifetimeNs, which an audit traced to an unseeded multistart in the IRF
% fitter - the same photons returned 0.134897 and 0.344474 ns, and across this
% session the stored value ranged 0.108-0.399 ns. Measured directly on one
% acquisition the bare SLB gives 0.3645 ns, agreeing with that acquisition's
% stored value to 1.2% while exposing the 0.108 and 0.135 ns entries as wrong.
%
% WHY THE WIDTH IS SMALL, AND MEASURED
%
% Because the SLB should be the SAME under the cell as beside it, tau1 needs only
% enough freedom to absorb real bilayer roughness. That roughness is measured
% rather than guessed: at 1x1 the bare-SLB per-pixel IQR is 0.0302 ns, and at 4x4
% (16x the photons) it falls only to 0.0142 ns rather than the 0.0076 that pure
% shot noise would give - so about 0.012 ns of the spread is genuine. The
% suggested sigma comes from that binned spread, around 0.010 ns, five times
% tighter than the 0.05 ns used before.
%
% Pinning tau1 tightly also improves the cell fit: with tau1 effectively fixed
% there is one free nonlinear parameter instead of two, so tau2 - the quantity of
% interest - is better determined.
%
% Results land in <acquisition>/biexp_slb_tau0pXXX (and _binK_tau0pXXX), so the
% anchored fits sit beside the earlier unanchored ones rather than replacing
% them, and the prior used is visible in the path.
%
% RESUMABLE: completed steps are logged to slb_anchored_done.txt at the data
% root, so a rerun after MATLAB's intermittent graphics-init abort continues
% where it stopped.

    if nargin < 1 || isempty(dataRoot)
        dataRoot = ['D:\Luminosa\Data\260813\' ...
            'RT_Jurkat_CD58_memglow_20260813_1'];
    end
    if nargin < 2 || isempty(binSizes); binSizes = [1 2 4]; end
    if nargin < 3 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'restart'); opts.restart = false; end
    if ~isfield(opts, 'includeDuplicates'); opts.includeDuplicates = false; end
    if ~isfield(opts, 'slbOpts'); opts.slbOpts = struct(); end
    % 'valid' is every pixel the detector actually covers. The minPhotons floor
    % inside the fitter then drops the empty corners, so nothing is wasted on
    % pixels with no signal.
    if ~isfield(opts, 'pixelMask'); opts.pixelMask = 'valid'; end

    doneFile = fullfile(dataRoot, 'slb_anchored_done.txt');
    done = {};
    if ~opts.restart && isfile(doneFile)
        done = strtrim(strsplit(fileread(doneFile), newline));
        done = done(~cellfun(@isempty, done));
    end

    found = dir(fullfile(dataRoot, '**', ...
        'immune_cell_MIET_640nm_red_analysis.mat'));
    paths = {};
    for k = 1:numel(found)
        if contains(found(k).folder, 'grid48') && ~opts.includeDuplicates
            continue;
        end
        paths{end + 1} = fullfile(found(k).folder, found(k).name); %#ok<AGROW>
    end
    paths = sort(paths);
    fprintf(['\nrun_immune_cell_MIET_slb_anchored\n  %d acquisition(s), ' ...
        'binnings %s\n  %s\n'], numel(paths), mat2str(binSizes), dataRoot);

    rows = struct([]);
    for k = 1:numel(paths)
        folder = fileparts(paths{k});

        % ---- stage 1: the SLB anchor ------------------------------------
        refFile = fullfile(folder, 'slb_reference', 'slb_reference_fit.mat');
        if isfile(refFile) && ~opts.restart
            L = load(refFile, 'out');
            ref = L.out;
            fprintf('\n### SLB anchor loaded: %.4f +/- %.4f ns  (%s)\n', ...
                ref.suggestedPrior.slbTauNs, ...
                ref.suggestedPrior.slbSigmaNs, folder);
        else
            fprintf('\n### measuring the SLB anchor: %s\n', folder);
            try
                ref = measure_slb_reference_lifetime(paths{k}, opts.slbOpts);
            catch refError
                fprintf('  SLB MEASUREMENT FAILED: %s\n', refError.message);
                continue;
            end
        end
        centre = ref.suggestedPrior.slbTauNs;
        sigma = ref.suggestedPrior.slbSigmaNs;

        % ---- stage 2: the cell fits, anchored --------------------------
        for b = binSizes(:)'
            tag = sprintf('%s|bin%d', paths{k}, b);
            if any(strcmp(done, tag))
                fprintf('### SKIP bin%d (done)\n', b);
                continue;
            end
            fprintf('### bin%d anchored at %.4f +/- %.4f ns\n', b, ...
                centre, sigma);
            entry = struct('analysisMat', paths{k}, 'binSize', b, ...
                'pixelMask', string(opts.pixelMask), ...
                'slbPooledNs', anchorOf(ref), ...
                'slbPriorCentreNs', centre, 'slbPriorSigmaNs', sigma, ...
                'storedReferenceNs', NaN, 'pixelsFitted', 0, ...
                'tau1MedianNs', NaN, 'tau1IqrNs', NaN, ...
                'tau2MedianNs', NaN, 'tau2Q1Ns', NaN, 'tau2Q3Ns', NaN, ...
                'photonShare2Median', NaN, 'tauMeanMedianNs', NaN, ...
                'reducedDevianceMedian', NaN, 'status', 'ok');
            if ~isempty(ref.storedReferenceNs)
                entry.storedReferenceNs = ref.storedReferenceNs(1);
            end
            try
                o = opts;
                o.binSize = b;
                o.slbTauNs = centre;
                o.slbSigmaNs = sigma;
                o.pixelMask = opts.pixelMask;
                r = immune_cell_MIET_biexp_slb(paths{k}, o);
                v = isfinite(r.maps.tau1Ns);
                entry.pixelsFitted = nnz(v);
                entry.tau1MedianNs = median(r.maps.tau1Ns(v));
                entry.tau1IqrNs = ...
                    quantileLocalBiexp(r.maps.tau1Ns(v), 0.75) - ...
                    quantileLocalBiexp(r.maps.tau1Ns(v), 0.25);
                entry.tau2MedianNs = median(r.maps.tau2Ns(v));
                entry.tau2Q1Ns = quantileLocalBiexp(r.maps.tau2Ns(v), 0.25);
                entry.tau2Q3Ns = quantileLocalBiexp(r.maps.tau2Ns(v), 0.75);
                entry.photonShare2Median = median(r.maps.photonFraction2(v));
                entry.tauMeanMedianNs = median(r.maps.tauMeanNs(v));
                entry.reducedDevianceMedian = ...
                    median(r.maps.reducedDeviance(v));
                fid = fopen(doneFile, 'a');
                fprintf(fid, '%s\n', tag);
                fclose(fid);
            catch fitError
                fprintf('  FIT FAILED: %s\n', fitError.message);
                entry.status = fitError.identifier;
            end
            if isempty(rows); rows = entry; else; rows(end + 1) = entry; end %#ok<AGROW>
        end
    end

    if isempty(rows)
        fprintf('\n  nothing to do - every step was already complete.\n');
        out = struct('summary', table(), 'dataRoot', dataRoot);
        return;
    end

    summary = struct2table(rows);
    csvFile = fullfile(dataRoot, 'slb_anchored_summary.csv');
    writetable(summary, csvFile);

    fprintf('\n=============================================================\n');
    fprintf('SLB-ANCHORED SUMMARY\n\n');
    fprintf(['  bin   pixels   SLBpool  prior+/-sig      tau1(IQR)        ' ...
        'tau2(ns)  tau2 IQR        share2  meanTau  acquisition\n']);
    fprintf('  %s\n', repmat('-', 1, 120));
    for k = 1:height(summary)
        [~, acq] = fileparts(fileparts(fileparts(summary.analysisMat{k})));
        fprintf(['  %3d %8d  %7.4f  %.4f+-%.4f  %.4f(%.4f)  %8.4f  ' ...
            '[%.3f %.3f]  %6.3f  %7.4f  %s\n'], ...
            summary.binSize(k), summary.pixelsFitted(k), ...
            summary.slbPooledNs(k), summary.slbPriorCentreNs(k), ...
            summary.slbPriorSigmaNs(k), summary.tau1MedianNs(k), ...
            summary.tau1IqrNs(k), summary.tau2MedianNs(k), ...
            summary.tau2Q1Ns(k), summary.tau2Q3Ns(k), ...
            summary.photonShare2Median(k), summary.tauMeanMedianNs(k), acq);
    end
    fprintf('\n  wrote %s\n', csvFile);
    out = struct('summary', summary, 'dataRoot', dataRoot, ...
        'csvFile', csvFile, 'doneFile', doneFile);
end
