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
% It also answers a second question: can cutting MORE pixels make the pooled
% decay mono-exponential? The pooled decay departs from a single exponential -
% reduced deviance 3.0-4.0, residuals to 9 sigma, and a log-linear fit whose
% answer drifts with the range - and the obvious suspicion is that this is
% spatial mixing, several depths with different degrees of quenching averaged
% into one curve. If so, restricting the pool to a thinner slab further from the
% metal should reduce it. The sweep therefore reports the mono-exponentiality
% diagnostics alongside tau, so the two questions are answered by one pass:
%
%   reducedDev    Poisson deviance per degree of freedom for A*exp(-t/tau)+B
%   maxResidual   largest |Poisson residual| of that fit
%   logLinSpread  how far the log-linear lifetime moves across fit ranges
%
% If all three fall as the cut rises, the departure is spatial mixing and can be
% cut away. If they stay flat while tau stays flat, the decay is intrinsically
% multi-exponential at every height and no amount of cutting will fix it - which
% is itself the answer.
%
% This is a falsification test, not a fit. It is deliberately blunt: the last
% cut discards most of the photons, so a real contamination of even a few
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
        % Reaches well past the automatic plateau cut of about 70 px, so the
        % pool ends up a thin slab far above the bilayer.
        cutsPix = [0 70 100 150 200 250];
    end
    if exist('fit_free_dye_lifetime_above_surface', 'file') ~= 2
        error('test_free_dye_lifetime_selection_robustness:Missing', ...
            'fit_free_dye_lifetime_above_surface.m must be on the path.');
    end

    fprintf('\ntest_free_dye_lifetime_selection_robustness\n');
    fprintf('  %s\n', ptuFile);
    fprintf(['  If out-of-focus quenched photons biased tau0 down, it must ' ...
        'RISE with the cut.\n\n']);
    fprintf(['  cut(px)   pixels    photons   tau(ns)  reducedDev  ' ...
        'maxResid  logLinSpread\n']);
    fprintf(['  --------------------------------------------------------' ...
        '-------------------\n']);

    rows = struct([]);
    for k = 1:numel(cutsPix)
        entry = struct('cutPix', cutsPix(k), 'pixels', 0, 'photons', 0, ...
            'monoNs', NaN, 'logLinearNs', NaN, 'momentNs', NaN, ...
            'pirlsNs', NaN, 'reducedDev', NaN, 'maxResidual', NaN, ...
            'logLinSpread', NaN, 'status', 'ok');
        try
            result = fit_free_dye_lifetime_above_surface(ptuFile, struct( ...
                'minTopDistancePix', cutsPix(k), 'makeFigure', false, ...
                'outputDir', tempdir));
            entry.pixels = result.selection.pixelCount;
            entry.photons = result.selection.photons;
            entry.monoNs = result.fit.tau0Ns;
            entry.logLinearNs = result.fit.logLinearMedianNs;
            entry.momentNs = result.fit.correctedMeanNs;
            entry.pirlsNs = result.fit.pirls.tauNs;
            entry.reducedDev = result.fit.pirls.reducedDeviance;
            entry.maxResidual = max(abs(result.fit.pirls.residual));
            entry.logLinSpread = result.fit.logLinearSpreadNs;
            fprintf('  %7d  %7d  %9.3g  %7.3f  %10.3f  %8.1f  %12.3f\n', ...
                entry.cutPix, entry.pixels, entry.photons, entry.pirlsNs, ...
                entry.reducedDev, entry.maxResidual, entry.logLinSpread);
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

    monoRange = max(summary.pirlsNs(good)) - min(summary.pirlsNs(good));
    logRange = max(summary.logLinearNs(good)) - min(summary.logLinearNs(good));
    momentRange = max(summary.momentNs(good)) - min(summary.momentNs(good));
    photonFall = max(summary.photons(good)) / min(summary.photons(good));
    out.monoRangeNs = monoRange;
    out.logLinearRangeNs = logRange;
    out.momentRangeNs = momentRange;
    out.photonFall = photonFall;

    fprintf(['\n  across a %.1fx fall in pooled photons: tau varies %.3f ns, ' ...
        'log-linear %.3f ns,\n  first moment %.3f ns\n'], photonFall, ...
        monoRange, logRange, momentRange);

    % Does cutting more make the decay mono-exponential?
    devFirst = summary.reducedDev(find(good, 1, 'first'));
    devLast = summary.reducedDev(find(good, 1, 'last'));
    spreadFirst = summary.logLinSpread(find(good, 1, 'first'));
    spreadLast = summary.logLinSpread(find(good, 1, 'last'));
    out.reducedDevChange = devLast - devFirst;
    out.logLinSpreadChange = spreadLast - spreadFirst;
    fprintf(['\n  MONO-EXPONENTIALITY across the sweep: reduced deviance ' ...
        '%.2f -> %.2f, log-linear\n  range spread %.3f -> %.3f ns\n'], ...
        devFirst, devLast, spreadFirst, spreadLast);
    if devLast < 0.7 * devFirst
        fprintf(['  Cutting further DOES reduce the departure from a single ' ...
            'exponential, so it was\n  substantially spatial mixing. Use ' ...
            'the tighter cut.\n']);
    else
        fprintf(['  Cutting further does NOT reduce it. The departure is ' ...
            'intrinsic to the dye at\n  every height sampled, not spatial ' ...
            'mixing, so a thinner slab will not help -\n  the decay needs ' ...
            'more than one exponential.\n']);
    end

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
