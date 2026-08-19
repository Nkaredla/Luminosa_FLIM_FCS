function test_flim_bayes_fixed_slb_gpu_grouping()
%TEST_FLIM_BAYES_FIXED_SLB_GPU_GROUPING Pick a safe slbCountRelTol.
%
% Covers what test_immune_cell_MIET_sliding_synthetic does not: the fixed-SLB
% evaluator itself.
%
% Two findings this script exists to record and re-check:
%
%  1. The GPU path is SLOWER here (measured 0.64x on a Quadro T1000 Max-Q).
%     The per-batch likelihood matmul is too small to amortise the transfers,
%     and its single precision perturbs the model argmax. useGPU stays false.
%
%  2. Relative-tolerance grouping of photon totals is the real speedup, but
%     it is an approximation. At slbCountRelTol = 0.02 modelProbability moved
%     by 0.07, which is far too much. This script sweeps the tolerance and
%     reports accuracy against the exact reference so the largest safe value
%     can be chosen on evidence.
%
% Acceptance used below: max |dP(model)| <= 0.01 AND model-selection
% agreement >= 99.5%. The script asserts that the value configured in
% run_batch_immune_cell_MIET would pass, and names the best safe tolerance.

    rng(7);

    TOLERANCES = [0.0005 0.001 0.0025 0.005 0.01 0.02];
    MAX_DP = 0.01;         % model-probability budget
    MIN_AGREEMENT = 0.995; % model-selection budget

    % ---- synthetic instrument ----
    dtNs = 0.16;
    pulsePeriodNs = 50;
    nBins = round(pulsePeriodNs / dtNs);
    timeNs = (0:nBins-1) * dtNs;

    irf = exp(-0.5 * ((timeNs - 2) / 0.25).^2);
    irf = irf / sum(irf);

    tauSlbNs = 0.5;
    tauMembraneNs = 1.4;
    decay = @(tau) conv(irf, exp(-timeNs / tau), 'full');
    slbShape = decay(tauSlbNs); slbShape = slbShape(1:nBins);
    memShape = decay(tauMembraneNs); memShape = memShape(1:nBins);
    slbShape = slbShape / sum(slbShape);
    memShape = memShape / sum(memShape);

    % ---- synthetic pixels ----
    % Totals span a 4x4-like range and are nearly all distinct, which mirrors
    % real binned data: the group count is set by the number of distinct
    % totals, so grouping has a lot to collapse.
    nPix = 600;
    totals = round(linspace(120, 4000, nPix));
    slbFractionTrue = 0.35;

    shape = slbFractionTrue * slbShape + (1 - slbFractionTrue) * memShape;
    shape = shape + 1e-12;   % keeps histcounts edges strictly increasing
    shape = shape / sum(shape);
    edges = [0, cumsum(shape(:).')];
    edges(end) = 1;
    assert(all(diff(edges) > 0), 'Sampling edges are not strictly increasing.');

    Y = zeros(nPix, 1, nBins);
    for k = 1:nPix
        Y(k, 1, :) = reshape(histcounts(rand(1, totals(k)), edges), 1, 1, nBins);
    end

    baseOpts = struct( ...
        'analysisMask', true(nPix, 1), ...
        'minPhotons', 10, ...
        'batchSize', 2048, ...
        'includeBackground', true, ...
        'signalGrid', [0.25 0.5 0.75 1], ...
        'membraneTauCount', 10, ...
        'fractionStep', 0.2, ...
        'minimumMembraneFraction', 0.1, ...
        'fixedSlbPhotonCount', slbFractionTrue * mean(totals), ...
        'fixedSlbPhotonCountStd', 0.1 * slbFractionTrue * mean(totals));
    run = @(useGPU, relTol) timeRun(baseOpts, useGPU, relTol, Y, irf, ...
        pulsePeriodNs, dtNs, tauSlbNs);

    % ---- exact reference ----
    [reference, tRef] = run(false, 0);
    exactGroups = numel(unique(totals));
    fprintf('\nreference: CPU exact grouping, %.2f s, %d groups (%d pixels)\n', ...
        tRef, exactGroups, nPix);

    % ---- sweep ----
    fprintf('\n relTol    groups  speedup   max|dP|   max|dTau|   agreement  verdict\n');
    fprintf(  ' ------    ------  -------   -------   ---------   ---------  -------\n');
    safe = NaN(1, numel(TOLERANCES));
    for k = 1:numel(TOLERANCES)
        relTol = TOLERANCES(k);
        [out, elapsed] = run(false, relTol);
        groups = countGroups(totals, relTol);
        dP = maxAbsDiff(reference, out, 'modelProbability');
        dTau = maxAbsDiff(reference, out, 'tauMeanArithmetic');
        agree = modelAgreement(reference, out);
        ok = dP <= MAX_DP && agree >= MIN_AGREEMENT;
        clipOk = abs( ...
            reference.fixedSlbPhotonConstraint.clippedPixelFraction - ...
            out.fixedSlbPhotonConstraint.clippedPixelFraction) < 1e-12;
        assert(clipOk, ...
            'slbCountRelTol=%g changed clippedPixelFraction; must stay exact.', ...
            relTol);
        if ok; safe(k) = relTol; end
        fprintf(' %6.4f    %6d  %6.2fx   %7.4f   %9.3g   %8.4f%%  %s\n', ...
            relTol, groups, tRef/elapsed, dP, dTau, 100*agree, ...
            ternary(ok, 'ok', 'TOO COARSE'));
    end

    best = max(safe(~isnan(safe)));
    fprintf('\nclipped-pixel fraction stayed exact at every tolerance: OK\n');
    if isempty(best)
        fprintf(['NO tested tolerance met the budget (|dP|<=%.3g, ' ...
                 'agreement>=%.2f%%). Use slbCountRelTol = 0.\n'], ...
                 MAX_DP, 100*MIN_AGREEMENT);
    else
        fprintf('largest tolerance meeting the budget: slbCountRelTol = %g\n', best);
    end

    % ---- confirm what run_batch_immune_cell_MIET is actually configured with ----
    configured = configuredRelTol();
    if isnan(configured)
        fprintf('could not read slbCountRelTol from run_batch_immune_cell_MIET\n');
    else
        fprintf('run_batch_immune_cell_MIET currently uses slbCountRelTol = %g\n', ...
            configured);
        assert(configured == 0 || (~isempty(best) && configured <= best), ...
            ['run_batch_immune_cell_MIET uses slbCountRelTol=%g, which did ' ...
             'not meet the accuracy budget. Lower it to %g or 0.'], ...
            configured, ternary(isempty(best), 0, best));
    end

    % ---- GPU, informational only ----
    if gpuDeviceCount > 0
        [gpuOut, tGpu] = run(true, 0);
        fprintf(['\nGPU (exact grouping): %.2f s = %.2fx vs CPU, model ' ...
                 'agreement %.4f%%, max|dP| %.4f\n'], ...
            tGpu, tRef/tGpu, 100*modelAgreement(reference, gpuOut), ...
            maxAbsDiff(reference, gpuOut, 'modelProbability'));
        fprintf('  -> confirms useGPU should stay false for this hardware.\n');
    end

    fprintf('\ntest_flim_bayes_fixed_slb_gpu_grouping: PASS\n');
end

function n = countGroups(totals, relTol)
    % Mirrors the binning in evaluate_fixed_count_grid.
    if relTol <= 0
        n = numel(unique(totals));
        return;
    end
    positive = totals > 0;
    key = -ones(1, numel(totals));
    key(positive) = floor(log(totals(positive)) / log1p(relTol));
    n = numel(unique(key));
end

function value = configuredRelTol()
    value = NaN;
    here = fileparts(mfilename('fullpath'));
    scriptFile = fullfile(here, 'run_batch_immune_cell_MIET.m');
    if ~isfile(scriptFile); return; end
    text = fileread(scriptFile);
    token = regexp(text, '''slbCountRelTol''\s*,\s*([0-9.eE+-]+)', 'tokens', 'once');
    if ~isempty(token); value = str2double(token{1}); end
end

function out = ternary(condition, a, b)
    if condition; out = a; else; out = b; end
end

function [out, elapsed] = timeRun(baseOpts, useGPU, relTol, Y, irf, ...
        pulsePeriodNs, dtNs, tauSlbNs)
    opts = baseOpts;
    opts.useGPU = useGPU;
    opts.slbCountRelTol = relTol;
    tic;
    out = flim_bayes_fixed_slb(Y, irf, pulsePeriodNs, dtNs, tauSlbNs, opts);
    elapsed = toc;
end

function v = pick(out, name)
    assert(isfield(out, name), ...
        'Expected output field "%s" is absent; update the test.', name);
    v = double(out.(name)(:));
    assert(~isempty(v), 'Output field "%s" is empty.', name);
end

function d = maxAbsDiff(x, y, name)
    a = pick(x, name); b = pick(y, name);
    assert(numel(a) == numel(b), 'Field "%s" size mismatch.', name);
    good = isfinite(a) & isfinite(b);
    assert(any(good), 'Field "%s" has no finite values.', name);
    d = max(abs(a(good) - b(good)));
end

function frac = modelAgreement(x, y)
    a = pick(x, 'completeExponentialCountMAP');
    b = pick(y, 'completeExponentialCountMAP');
    assert(numel(a) == numel(b), 'Model map sizes differ.');
    good = a > 0 & b > 0;
    assert(any(good), 'No pixels were assigned a model.');
    frac = mean(a(good) == b(good));
end
