function test_flim_bayes_fixed_slb_two_state
%TEST_FLIM_BAYES_FIXED_SLB_TWO_STATE Verify the maxMembraneStates ceiling.
%
% Checks that opts.maxMembraneStates = 1 truncates the model comparison to
% the fixed SLB and fixed SLB + one free membrane lifetime, without
% perturbing anything the biexponential model already decided:
%
%   1. the triexponential grid is never built (its state count is NaN)
%   2. P(M3) is exactly zero and no pixel selects M3
%   3. the model-averaged maps stay finite and the probabilities sum to one
%   4. the M2-conditional lifetime map is bit-identical to the three-model run
%   5. the sorted component-3 display layer is empty while component 2 is not
%
% Run it after any change to the model set, the grids, or the model average.

    rng(7);
    dtNs = 0.16;
    nt = 128;
    pulsePeriodNs = nt * dtNs;
    t = (0:nt-1)' * dtNs;

    irf = exp(-((t - 0.5) .^ 2) / (2 * 0.12 ^ 2));
    irf = irf / sum(irf);

    tauSlb = 0.55;
    pSlb = decayPattern(irf, t, tauSlb, nt);
    p1 = decayPattern(irf, t, 1.40, nt);
    p2 = decayPattern(irf, t, 3.10, nt);

    % Three pixel families: pure SLB, SLB + one membrane state, SLB + two.
    nx = 6; ny = 3;
    photons = 4000;
    weights = [0.95 0 0; 0.55 0.40 0; 0.40 0.30 0.25];
    Ypix = zeros(nx, ny, nt);
    for ix = 1:nx
        for iy = 1:ny
            shape = weights(iy, 1) * pSlb + weights(iy, 2) * p1 + ...
                weights(iy, 3) * p2 + 0.05 / nt;
            shape = shape / sum(shape);
            Ypix(ix, iy, :) = reshape(poissrnd(photons * shape), 1, 1, nt);
        end
    end

    opts = struct('minPhotons', 50, 'useGPU', false, 'batchSize', 64, ...
        'membraneTauCount', 10, 'fractionStep', 0.2, ...
        'minimumMembraneFraction', 0.1, 'signalGrid', [0.5 0.75 1]);

    full3 = flim_bayes_fixed_slb(Ypix, irf, pulsePeriodNs, dtNs, tauSlb, opts);
    opts2 = opts;
    opts2.maxMembraneStates = 1;
    two = flim_bayes_fixed_slb(Ypix, irf, pulsePeriodNs, dtNs, tauSlb, opts2);

    assert(full3.maxMembraneStates == 2, 'the default must stay at two states');
    assert(two.maxMembraneStates == 1, 'the ceiling was not applied');
    assert(isequal(two.evaluatedModels, [true true false]), ...
        'M3 must be marked as not evaluated');
    assert(isnan(two.gridInfo.modelStateCount(3)), ...
        'the M3 grid must never be built');
    assert(all(isfinite(full3.gridInfo.modelStateCount)), ...
        'the three-model run must still build every grid');
    assert(two.modelPrior(3) == 0, 'M3 must carry zero prior mass');

    valid = two.validMask;
    assert(all(valid(:)), 'every synthetic pixel should be fitted');
    assert(all(two.probabilityTriexponential(valid) == 0), ...
        'M3 must carry exactly zero posterior probability');
    assert(all(two.completeExponentialCountMAP(valid) <= 2), ...
        'no pixel may select M3');
    assert(all(isfinite(two.tauMeanArithmetic(valid))), ...
        'the model-averaged lifetime must stay finite');
    assert(all(isfinite(two.fixedSlbFraction(valid))), ...
        'the model-averaged SLB fraction must stay finite');
    assert(all(two.membrane2PhotonFraction(valid) == 0), ...
        'the component-3 photon fraction must be identically zero');
    assert(max(abs(sum(two.modelProbability, 3) - 1), [], 'all') < 1e-4, ...
        'model probabilities must still sum to one');

    % The M2 grid and likelihood are untouched, so its conditional maps must
    % not shift when the model set around them shrinks.
    delta = abs(two.oneMembrane.membraneLifetime1Ns(valid) - ...
        full3.oneMembrane.membraneLifetime1Ns(valid));
    assert(max(delta) == 0, ...
        'the M2-conditional fit must not depend on whether M3 was evaluated');

    % Where the three-model run already rejects M3, the model-averaged answer
    % must be unchanged too.
    weakM3 = full3.probabilityTriexponential < 1e-3 & valid;
    if any(weakM3(:))
        d = abs(two.tauMeanArithmetic(weakM3) - full3.tauMeanArithmetic(weakM3));
        assert(max(d) < 5e-3, ...
            'restricting M3 changed a pixel where M3 had already lost');
    end

    layers = immune_cell_MIET_sorted_components(two, (1:(nx * ny))', ...
        [nx ny], true(nx, ny), sum(Ypix, 3), ...
        struct('posteriorThreshold', [0.8 0.95], 'minExpectedPhotons', [10 10]));
    assert(nnz(isfinite(layers.display.thirdLifetimeNs)) == 0, ...
        'no component-3 pixel may be displayed');
    assert(nnz(isfinite(layers.display.secondLifetimeNs)) > 0, ...
        'component 2 should still be recovered');

    fprintf('test_flim_bayes_fixed_slb_two_state: PASS\n');
    fprintf('  M2 grid states           : %d\n', two.gridInfo.modelStateCount(2));
    fprintf('  M3 grid states skipped   : %d\n', ...
        full3.gridInfo.modelStateCount(3));
    fprintf('  component-2 display pix  : %d of %d\n', ...
        nnz(isfinite(layers.display.secondLifetimeNs)), nx * ny);
end

function p = decayPattern(irf, t, tauNs, nt)
    d = exp(-t / tauNs);
    p = conv(irf, d);
    p = p(1:nt);
    p = p / sum(p);
end
