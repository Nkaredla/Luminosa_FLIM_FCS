function out = run_free_dye_lifetime_crosssections(dataRoot, opts)
%RUN_FREE_DYE_LIFETIME_CROSSSECTIONS tau0 from every XZ/YZ scan in a session.
%
% out = run_free_dye_lifetime_crosssections()
% out = run_free_dye_lifetime_crosssections(dataRoot, opts)
%
% Finds every non-XY acquisition below dataRoot, runs
% estimate_free_dye_lifetime_from_crosssection on each, and combines the
% results.
%
% Five cross-sections exist in the 260813 Jurkat/CD58/memglow session and the
% XY analysis never used them - the batch driver catalogues them and skips
% them, because the MIET pipeline only fits XY planes. They are the only
% acquisitions that span height, so they are the only ones that can measure the
% unquenched lifetime.
%
% AGREEMENT ACROSS SCANS IS THE POINT
%
% Each scan gives its own tau0. They were taken at different times (15:48 to
% 16:55) on different cells, so if they agree the number is a property of the
% dye and not of one field of view. Spread across scans is reported next to the
% within-scan spread over tail starts, because those two numbers fail for
% different reasons: within-scan spread means the tail is contaminated, while
% between-scan spread means the sample is heterogeneous or some scans never
% reach unquenched dye.
%
% Scans that do NOT reach a plateau contribute a lower bound only and are
% excluded from the combined value, with that stated rather than averaged in.

    if nargin < 1 || isempty(dataRoot)
        dataRoot = ['D:\Luminosa\Data\260813\' ...
            'RT_Jurkat_CD58_memglow_20260813_1'];
    end
    if nargin < 2 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'outputRoot') || isempty(opts.outputRoot)
        opts.outputRoot = fullfile(dataRoot, 'free_dye_lifetime');
    end
    if ~isfolder(opts.outputRoot); mkdir(opts.outputRoot); end

    files = dir(fullfile(dataRoot, '**', 'RawImage.ptu'));
    if isempty(files)
        error('run_free_dye_lifetime_crosssections:NoFiles', ...
            'No RawImage.ptu below %s', dataRoot);
    end
    fullNames = arrayfun(@(f) fullfile(f.folder, f.name), files, ...
        'UniformOutput', false);
    fullNames = sort(fullNames);

    fprintf(['run_free_dye_lifetime_crosssections: screening %d PTU files ' ...
        'for cross-sections\n'], numel(fullNames));
    crossSections = {};
    planes = {};
    for k = 1:numel(fullNames)
        try
            head = PTU_Read_Head(fullNames{k});
        catch headError
            fprintf('  skip (header failed): %s [%s]\n', fullNames{k}, ...
                headError.identifier);
            continue;
        end
        if ~isfield(head, 'ImgHdr_ScanDirection'); continue; end
        direction = double(head.ImgHdr_ScanDirection);
        switch direction
            case 1; plane = 'YZ';
            case 2; plane = 'XZ';
            otherwise; plane = '';
        end
        if isempty(plane); continue; end
        crossSections{end + 1} = fullNames{k}; %#ok<AGROW>
        planes{end + 1} = plane; %#ok<AGROW>
    end
    fprintf('  found %d cross-section(s)\n', numel(crossSections));
    if isempty(crossSections)
        error('run_free_dye_lifetime_crosssections:NoCrossSections', ...
            'No XZ or YZ acquisitions were found below %s', dataRoot);
    end

    rows = struct([]);
    results = cell(numel(crossSections), 1);
    for k = 1:numel(crossSections)
        [folder, ~] = fileparts(crossSections{k});
        [~, acquisition] = fileparts(folder);
        localOpts = opts;
        localOpts.outputDir = fullfile(opts.outputRoot, acquisition);
        try
            result = estimate_free_dye_lifetime_from_crosssection( ...
                crossSections{k}, localOpts);
            results{k} = result;
            entry = struct('acquisition', acquisition, ...
                'plane', planes{k}, ...
                'axialAxis', result.axial.name, ...
                'axialConcentration', max(result.axial.concentrationLine, ...
                    result.axial.concentrationCol), ...
                'metalEnd', result.metal.name, ...
                'arrivalContrastNs', result.metal.contrastNs, ...
                'bandFwhmPixels', result.bandFwhmPixels, ...
                'internalUmPerPixel', result.internalUmPerPixel, ...
                'plateauReached', result.plateau.ok, ...
                'tau0Ns', result.plateau.tau0Ns, ...
                'tau0SdNs', result.plateau.tau0SdNs, ...
                'plateauOnsetPixels', result.plateau.startPixels, ...
                'plateauPhotons', result.plateau.photons, ...
                'minTauNs', result.plateau.minTauNs, ...
                'maxTauNs', result.plateau.maxTauNs, ...
                'irfCheckTauNs', result.crossCheck.tauNs);
        catch runError
            fprintf('\n  FAILED on %s: %s\n', acquisition, runError.message);
            entry = struct('acquisition', acquisition, 'plane', planes{k}, ...
                'axialAxis', 'failed', 'axialConcentration', NaN, ...
                'metalEnd', '', 'arrivalContrastNs', NaN, ...
                'bandFwhmPixels', NaN, 'internalUmPerPixel', NaN, ...
                'plateauReached', false, 'tau0Ns', NaN, 'tau0SdNs', NaN, ...
                'plateauOnsetPixels', NaN, 'plateauPhotons', 0, ...
                'minTauNs', NaN, 'maxTauNs', NaN, 'irfCheckTauNs', NaN);
        end
        if isempty(rows); rows = entry; else; rows(end + 1) = entry; end %#ok<AGROW>
    end

    summary = struct2table(rows);
    csvFile = fullfile(opts.outputRoot, 'free_dye_lifetime_summary.csv');
    writetable(summary, csvFile);

    fprintf('\n=============================================================\n');
    fprintf('SUMMARY: unquenched lifetime from %d cross-section(s)\n', ...
        numel(crossSections));
    fprintf(['  acquisition          plane  plateau  tau0(ns)   sd    ' ...
        'range seen (ns)   IRF check\n']);
    for k = 1:height(summary)
        if summary.plateauReached(k)
            plateauText = 'yes';
            tauText = sprintf('%7.3f %6.3f', summary.tau0Ns(k), ...
                summary.tau0SdNs(k));
        else
            plateauText = 'NO ';
            tauText = sprintf('%7s %6s', '-', '-');
        end
        if isfinite(summary.irfCheckTauNs(k))
            irfText = sprintf('%.3f', summary.irfCheckTauNs(k));
        else
            irfText = '-';
        end
        fprintf('  %-20s %-5s  %-7s %s   %5.3f - %5.3f     %s\n', ...
            summary.acquisition{k}, summary.plane{k}, plateauText, tauText, ...
            summary.minTauNs(k), summary.maxTauNs(k), irfText);
    end

    reached = summary.plateauReached;
    out = struct('summary', summary, 'results', {results}, ...
        'dataRoot', dataRoot, 'opts', opts, ...
        'combinedTau0Ns', NaN, 'combinedSdNs', NaN, 'lowerBoundNs', NaN);
    if any(reached)
        values = summary.tau0Ns(reached);
        weights = summary.plateauPhotons(reached);
        weights = weights / sum(weights);
        out.combinedTau0Ns = sum(weights .* values);
        out.combinedSdNs = std(values);
        fprintf(['\n  combined tau0 = %.3f ns (photon-weighted over %d ' ...
            'scan(s)); between-scan sd %.3f ns\n'], out.combinedTau0Ns, ...
            nnz(reached), out.combinedSdNs);
        fprintf(['  Between-scan spread is the number to judge this by: ' ...
            'the scans are different\n  cells at different times, so ' ...
            'agreement means tau0 belongs to the dye rather\n  than to one ' ...
            'field of view.\n']);
    end
    if any(~reached)
        out.lowerBoundNs = max(summary.maxTauNs(~reached));
        fprintf(['\n  %d scan(s) never reached a plateau and are EXCLUDED ' ...
            'from the combined value.\n  Their longest observed lifetime, ' ...
            '%.3f ns, is a lower bound only.\n'], nnz(~reached), ...
            out.lowerBoundNs);
    end

    if isfinite(out.combinedTau0Ns)
        fprintf(['\n  What this settles: calibrationCurve.mat tops out at ' ...
            '3.362 ns at its last\n  tabulated height (30 nm) and its most ' ...
            'quenched value is 0.865 ns. A measured\n  tau0 of %.3f ns is ' ...
            'the asymptote that curve must approach, so it decides\n  ' ...
            'whether the curve is merely truncated or built on the wrong ' ...
            'stack.\n'], out.combinedTau0Ns);
    end
    fprintf('\n  wrote %s\n', csvFile);
    matFile = fullfile(opts.outputRoot, 'free_dye_lifetime_summary.mat');
    save(matFile, 'out', '-v7.3');
    fprintf('  wrote %s\n', matFile);
end
