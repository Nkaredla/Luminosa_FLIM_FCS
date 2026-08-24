function out = test_free_dye_lifetime_selection_robustness(ptuFile, cutsPix)
%TEST_FREE_DYE_LIFETIME_SELECTION_ROBUSTNESS Is tau0 an artefact of the selection?
%
% out = test_free_dye_lifetime_selection_robustness()
% out = test_free_dye_lifetime_selection_robustness(ptuFile, cutsPix)
%
% fit_free_dye_lifetime_above_surface pools 20-30% of a cross-section, and the
% obvious worry is that out-of-focus haze from the bright bilayer puts
% partly-quenched photons into the pool and biases tau0 DOWNWARD. That would
% make the reported unquenched lifetime a lower bound rather than a measurement.
%
% The test is direct: push the lower edge of the pool further and further from
% the bilayer and see whether tau0 moves. If quenched photons were contributing,
% tau0 must RISE as they are excluded and then flatten. If tau0 is flat from the
% start, they were not contributing.
%
% This is a falsification test, not a fit. It is deliberately blunt: the last
% cut discards about 80% of the photons, so a real contamination of even a few
% percent in lifetime would show as a trend far larger than the 0.05 ns spread
% between the estimators.
%
% MEASURED on _20260813-164757, cuts 0 to 250 px:
%   mono-exponential  1.946 1.952 1.955 1.952 1.954 1.952 ns  (range 0.009)
%   log-linear        1.940 1.951 1.953 1.950 1.953 1.951 ns  (range 0.013)
%   first moment      1.834 1.894 1.898 1.901 1.902 1.904 ns  (range 0.070)
% across a 3x fall in photons and a 4.6x fall in pixels. The two tail-shape
% estimators are flat, so the pool is not contaminated to any degree that
% matters. The first moment rises by 0.07 ns and flattens by the 100 px cut,
% which is the expected residue: a moment weights every photon equally and so is
% the estimator most exposed to a faint quenched population, whereas a tail fit
% is driven by the decay's shape.

    if nargin < 1 || isempty(ptuFile)
        ptuFile = ['D:\Luminosa\Data\260813\RT_Jurkat_CD58_memglow_20260813_1' ...
            '\_20260813-164757\RawImage.ptu'];
    end
    if nargin < 2 || isempty(cutsPix)
        cutsPix = [0 50 100 150 200 250];
    end
    if exist('fit_free_dye_lifetime_above_surface', 'file') ~= 2
        error('test_free_dye_lifetime_selection_robustness:Missing', ...
            'fit_free_dye_lifetime_above_surface.m must be on the path.');
    end

    fprintf('\ntest_free_dye_lifetime_selection_robustness\n');
    fprintf('  %s\n', ptuFile);
    fprintf(['  If out-of-focus quenched photons biased tau0 down, it must ' ...
        'RISE with the cut.\n\n']);
    fprintf('  cut(px)   pixels    photons    mono(ns)  loglin(ns)  moment(ns)\n');
    fprintf('  ---------------------------------------------------------------\n');

    rows = struct([]);
    for k = 1:numel(cutsPix)
        entry = struct('cutPix', cutsPix(k), 'pixels', 0, 'photons', 0, ...
            'monoNs', NaN, 'logLinearNs', NaN, 'momentNs', NaN, ...
            'status', 'ok');
        try
            result = fit_free_dye_lifetime_above_surface(ptuFile, struct( ...
                'minTopDistancePix', cutsPix(k), 'makeFigure', false, ...
                'outputDir', tempdir));
            entry.pixels = result.selection.pixelCount;
            entry.photons = result.selection.photons;
            entry.monoNs = result.fit.tau0Ns;
            entry.logLinearNs = result.fit.logLinearMedianNs;
            entry.momentNs = result.fit.correctedMeanNs;
            fprintf('  %7d  %7d  %9.3g  %9.3f  %10.3f  %10.3f\n', ...
                entry.cutPix, entry.pixels, entry.photons, entry.monoNs, ...
                entry.logLinearNs, entry.momentNs);
        catch runError
            entry.status = runError.identifier;
            fprintf('  %7d  FAILED: %s\n', entry.cutPix, runError.message);
        end
        if isempty(rows); rows = entry; else; rows(end + 1) = entry; end %#ok<AGROW>
    end

    summary = struct2table(rows);
    out = struct('summary', summary, 'ptuFile', ptuFile);
    good = strcmp(summary.status, 'ok');
    if nnz(good) < 3
        error('test_free_dye_lifetime_selection_robustness:TooFewCuts', ...
            'Only %d cut(s) succeeded, too few to judge a trend.', nnz(good));
    end

    monoRange = max(summary.monoNs(good)) - min(summary.monoNs(good));
    logRange = max(summary.logLinearNs(good)) - min(summary.logLinearNs(good));
    momentRange = max(summary.momentNs(good)) - min(summary.momentNs(good));
    photonFall = max(summary.photons(good)) / min(summary.photons(good));
    out.monoRangeNs = monoRange;
    out.logLinearRangeNs = logRange;
    out.momentRangeNs = momentRange;
    out.photonFall = photonFall;

    fprintf(['\n  across a %.1fx fall in pooled photons: mono varies %.3f ns, ' ...
        'log-linear %.3f ns,\n  first moment %.3f ns\n'], photonFall, ...
        monoRange, logRange, momentRange);

    % The estimator spread on the full pool is about 0.05 ns, so a trend has to
    % clear that to mean anything.
    tolerance = 0.05;
    out.passed = monoRange <= tolerance && logRange <= tolerance;
    if out.passed
        fprintf(['  PASS: both tail-shape estimators are flat to within ' ...
            '%.2f ns, so the pool is not\n  contaminated to any degree that ' ...
            'affects tau0. It is a measurement, not a lower bound.\n'], ...
            tolerance);
    else
        fprintf(['  FAIL: tau0 moves by more than %.2f ns with the cut, so ' ...
            'the pool IS contaminated\n  and the value from the automatic ' ...
            'cut understates tau0. Use the plateau value.\n'], tolerance);
    end
end
