function test_B12_tail_fitting
%TEST_B12_TAIL_FITTING Synthetic regression checks; no PTU data required.
    folder = fileparts(mfilename('fullpath'));
    oldPath = path;
    cleanup = onCleanup(@() path(oldPath));
    addpath(folder, fileparts(folder), '-begin');
    cfg = defaultConfig(struct());
    assert(cfg.minPhotonsPerRegionFit == 10 && cfg.minPhotonsForTwoExp == 500);
    dt = 0.125;
    irf = [1; zeros(79, 1)];
    args = {irf, dt, 2, 1.5, [0.4 2.5], [0.05 10], false, 10, 150, 2, 0.1, 0.02, 500};
    % Integer photon totals test eligibility, independently of selected BIC.
    for photons = [0 9 10 499 500]
        curve = photonHistogram(photons, 80);
        [ok, order, ~, ~, ~, ~, fitted, metrics] = selectModel(curve, args{:});
        assert(metrics(4) == photons);
        if photons < 10
            assert(~ok && order == 0 && all(isnan(fitted)));
        elseif photons < 500
            assert(ok && order == 1 && isinf(metrics(6)));
        else
            assert(ok && isfinite(metrics(6)));
        end
    end
    % A bright peak must not admit biexponential fitting of a photon-poor tail.
    for photons = [9 10 499 500]
        curve = zeros(80, 1);
        curve(3) = 10000;
        curve(11:end) = photonHistogram(photons, 70);
        [ok, order, ~, ~, ~, ~, fitted, metrics] = selectModel(curve, args{:}, 1);
        assert(metrics(4) == photons && all(isnan(fitted(1:10))));
        assert(ok == (photons >= 10));
        if photons >= 10 && photons < 500
            assert(order == 1 && isinf(metrics(6)));
        elseif photons == 500
            assert(isfinite(metrics(6)));
        end
    end
    % Known tail lifetime and independently shifted peaks, with noninteger
    % offset/bin ratio (ceil(1/0.128) = 8 bins = 1.024 ns).
    dt = 0.128;
    curve = [zeros(4, 1); 1000 * exp(-(0:75)' * dt / 1.8)];
    args{2} = dt;
    args{3} = 1;
    [ok, ~, tau, ~, ~, ~, fitted] = selectModel(curve, args{:}, 1);
    assert(ok && abs(tau(1) - 1.8) < 0.01);
    assert(all(isnan(fitted(1:12))) && all(isfinite(fitted(13:end))));
    late = zeros(80, 1); late(end) = 1000;
    [ok, order, ~, ~, ~, ~, fitted, metrics] = selectModel(late, args{:}, 1);
    assert(~ok && order == 0 && metrics(4) == 0 && all(isnan(fitted)));
    [ok, order, ~, ~, ~, ~, fitted] = fitRegions(zeros(80, 0), 'empty', args{:}, 1);
    assert(isempty(ok) && isempty(order) && isequal(size(fitted), [80 0]));

    % Public entry point on presegmented data: no stale Bayesian output,
    % skipped maps remain NaN, and each region retains its own cutoff.
    source = struct('file', fullfile(tempdir, 'synthetic.ptu'));
    source.gates = struct('dtNs', dt);
    source.masks = struct('cell', true(2), 'cytosol', logical([0 0; 0 1]), ...
        'compartmentMap', uint8([1 2; 1 3]));
    source.labels = struct('mitochondria', [1 0; 2 0], 'lysosome', [0 1; 0 0]);
    poor = zeros(80, 1); poor(2) = 1000; poor(10) = 9;
    shifted = [zeros(3, 1); curve(1:end-3)];
    source.tcspc = struct('globalB12', curve + shifted + poor, ...
        'cytosol', curve, 'mitochondria', [curve poor], 'lysosome', shifted, ...
        'timeNs', ((0:79)' + 0.5) * dt);
    source.bayes = struct('stale', true);
    result = PTU_MitoLysoCytosol_TailFLIM(source, ...
        struct('showFigures', false, 'saveOutputs', false, 'maxExp', 1));
    assert(~isfield(result, 'bayes'));
    assert(isequal(result.fits.mitochondria.ok, [true; false]));
    assert(result.fits.mitochondria.firstFitBin(1) == 13);
    assert(result.fits.lysosome.firstFitBin == 16);
    assert(isnan(result.flim.tauMeanNs(2, 1)) && result.flim.modelOrder(2, 1) == 0);
    assert(abs(result.flim.tauMeanNs(1, 1) - 1.8) < 0.01);
    fprintf('PASS: B12 photon thresholds, independent tails, lifetime recovery, and maps.\n');
end

function counts = photonHistogram(photons, bins)
    weights = exp(-(0:bins-1)' / 8) + 0.4 * exp(-(0:bins-1)' / 24);
    counts = floor(photons * weights / sum(weights));
    remaining = photons - sum(counts);
    counts(1:remaining) = counts(1:remaining) + 1;
end
