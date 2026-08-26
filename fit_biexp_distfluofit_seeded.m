function out = fit_biexp_distfluofit_seeded(decay, irf, dtNs, periodNs, opts)
%FIT_BIEXP_DISTFLUOFIT_SEEDED Biexponential seeded from a lifetime distribution.
%
% out = fit_biexp_distfluofit_seeded(decay, irf, dtNs, periodNs)
% out = fit_biexp_distfluofit_seeded(decay, irf, dtNs, periodNs, opts)
%
% Fits  y(t) = B + a1*P(t;tau1) + a2*P(t;tau2)  with NO prior on either
% lifetime, and returns them sorted so tau1 <= tau2. Used to measure the bare-SLB
% reference, where the SLB lifetime is taken as the SHORTER component.
%
% WHY SEED FROM A LIFETIME DISTRIBUTION
%
% A two-component fit of a decay with a small long-lived tail is easy to start
% badly: seed both lifetimes near the dominant short component and the search can
% sit there, reporting a near-degenerate pair and missing the tail entirely. A
% distribution fit finds where the components are without being told how many
% there are, so its modes are good starting values.
%
% WHAT DistFluofit ACTUALLY SAYS HERE, ONCE ITS UNITS ARE RIGHT
%
% Its tau output is a RATE unless it plots - the inversion sits inside the
% "if bild" block - so distfluofit_peaks inverts it. Done correctly, the grid
% spans dt to p (0.16 to 50 ns here) and the whole distribution piles onto the
% SHORTEST node, 0.16 ns.
%
% That is not a failure. It is the grid saying the lifetime is at or below one
% bin width, which is right: the two-component fit puts the SLB at 0.077 ns, and
% DistFluofit cannot represent anything below dt. Its own source carries the
% remedy commented out one line above the live grid,
% %tau = (1/dt/10)./exp((0:N)/N*log(p/dt/10)), which shifts the grid down a
% decade for exactly this case.
%
% So a pile-up on the LOWEST node is treated as information - seed below it -
% rather than as a degenerate result to discard. A pile-up on the highest node,
% or modes outside a plausible range, still falls back to a coarse
% two-dimensional scan of the pooled decay, which is data-driven and cheap for
% one decay.
%
% This matters for the bare SLB specifically. On acquisition 155036 the pooled
% bare-SLB decay is dominated by a ~0.08 ns component but carries a few percent
% of photons at ~1.3 ns; a mono fit of the same data returns 0.0915 ns with a
% reduced deviance of 3924 and a background three times the observed pre-pulse
% level, because it has nowhere to put that tail.
%
% WHY THE FIT IS UNPENALISED
%
% This IS the measurement of the SLB lifetime, so constraining it to a prior
% would be circular. The prior belongs downstream, in the cell fit, centred on
% what comes out of here.
%
% opts fields
%   nDistNodes   lifetime nodes for DistFluofit (default 100)
%   shiftBins    colour-shift search bounds in channels (default [-3 3])
%   tauSeedsNs   override the seeds entirely, as [short long]
%   gtol         outer gradient tolerance (default 1e-4, tighter than the
%                per-pixel fits because this is one pooled decay)
%
% out fields
%   tau1Ns, tau2Ns        sorted lifetimes
%   slbTauNs              = tau1Ns, the value to anchor the cell fit on
%   amplitude             [a1 a2] as PHOTON counts (the patterns are unit-sum)
%   photonFraction        [f1 f2]
%   background            counts per bin
%   deviance, reducedDeviance
%   model, residual
%   distTauNs, distWeights, distSeedsNs   what DistFluofit reported
%   monoTauNs, monoReducedDeviance        a mono fit, for comparison

    if nargin < 5 || isempty(opts); opts = struct(); end
    defaults = struct('nDistNodes', 100, 'shiftBins', [-3 3], ...
        'tauSeedsNs', [], 'gtol', 1e-4);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end

    decay = double(decay(:));
    nBin = numel(decay);
    irf = double(irf(:));
    irf = max(irf, 0);
    if sum(irf) > 0; irf = irf / sum(irf); end
    basis = biexp_slb_basis(irf, dtNs, periodNs, nBin);

    % ---- stage 1: where are the components? ------------------------------
    distTau = [];
    distWeights = [];
    distRejected = false;
    seeds = opts.tauSeedsNs;
    if isempty(seeds)
        try
            % distfluofit_peaks takes the RAW distribution and reduces it here,
            % because DistFluofit's own peak finder cannot handle a peak at the
            % first node - which is exactly where the SLB component sits.
            % Its tau output is already in ns; do NOT multiply by dt.
            [peakTau, peakW] = distfluofit_peaks(irf, decay, periodNs, ...
                dtNs, opts.shiftBins, opts.nDistNodes);
            distTau = peakTau;
            distWeights = peakW;
            good = isfinite(distTau) & distTau > 0;
            distTau = distTau(good);
            distWeights = distWeights(good);
            atShortEdge = ~isempty(distTau) && ...
                abs(min(distTau) - dtNs) < 0.05 * dtNs;
            if numel(distTau) >= 2
                seeds = sort(distTau(1:2));
                if atShortEdge
                    % The short mode is at the grid floor, so the true value is
                    % below it. Start half a bin down.
                    seeds(1) = dtNs / 2;
                end
            elseif numel(distTau) == 1 && atShortEdge
                fprintf(['      DistFluofit puts all weight on its LOWEST ' ...
                    'node (%.4f ns = dt), i.e.\n      the lifetime is below ' ...
                    'its grid; seeding below it\n'], distTau(1));
                seeds = [dtNs / 2, 1.5];
            elseif numel(distTau) == 1
                seeds = sort([distTau(1), distTau(1) * 6]);
            else
                distRejected = true;
            end
        catch distError
            fprintf(['      DistFluofit seeding failed (%s); using a coarse ' ...
                'scan instead\n'], distError.message);
        end
    end

    % Coarse two-dimensional scan: still data-driven, and cheap for one decay.
    if isempty(seeds) || numel(seeds) < 2 || any(~isfinite(seeds))
        g1 = logspace(log10(0.03), log10(1.0), 14);
        g2 = logspace(log10(0.3), log10(8.0), 16);
        bestDev = inf;
        for a = g1
            for b = g2
                if b <= a * 1.3; continue; end
                Pa = biexp_slb_pattern_batch(basis, a);
                Pb = biexp_slb_pattern_batch(basis, b);
                bta = poisson_nnls_irls_batch(decay, Pa, Pb, ...
                    struct('maxIter', 30, 'tol', 1e-10));
                mdl = max(bta(1) + bta(2) * Pa + bta(3) * Pb, 1e-12);
                dv = poisson_nnls_em_deviance(decay, mdl);
                if dv < bestDev
                    bestDev = dv;
                    seeds = [a b];
                end
            end
        end
        fprintf('      coarse scan seeds [%.4f %.4f]\n', seeds(1), seeds(2));
    end
    seeds = reshape(sort(seeds(:)), 1, []);
    if seeds(2) <= seeds(1) * 1.2; seeds(2) = seeds(1) * 4; end

    % ---- stage 2: the biexponential, unpenalised -------------------------
    % A huge prior width makes the penalty vanish, so biexp_slb_bfgs_batch
    % performs a free two-lifetime fit. Reusing it keeps one implementation of
    % the variable-projection machinery rather than a second, divergent copy.
    fitOpts = struct('slbTauNs', seeds(1), 'slbSigmaNs', 1e6, ...
        'tau2SeedNs', seeds(2), 'gtol', opts.gtol, 'maxIter', 200, ...
        'innerSolver', 'irls');
    vp = biexp_slb_bfgs_batch(decay, basis, fitOpts);

    taus = [vp.tau1Ns, vp.tau2Ns];
    amps = [vp.beta(2), vp.beta(3)];
    [taus, order] = sort(taus);
    taus = reshape(taus, 1, []);
    amps = amps(order);

    model = max(vp.beta(1) + amps(1) * ...
        biexp_slb_pattern_batch(basis, taus(1)) + amps(2) * ...
        biexp_slb_pattern_batch(basis, taus(2)), 1e-12);
    dev = poisson_nnls_em_deviance(decay, model);

    % ---- a mono fit for comparison ---------------------------------------
    monoGrid = logspace(log10(0.03), log10(8), 120);
    mono = fit_mono_poisson_batch(decay, basis, monoGrid, struct());

    out = struct();
    out.tau1Ns = taus(1);
    out.tau2Ns = taus(2);
    out.slbTauNs = taus(1);
    out.amplitude = amps;
    out.photonFraction = amps / max(sum(amps), eps);
    out.background = vp.beta(1);
    out.deviance = dev;
    out.reducedDeviance = dev / max(nBin - 5, 1);
    out.model = model;
    out.residual = (decay - model) ./ sqrt(max(model, 1));
    out.photons = sum(decay);
    out.distTauNs = distTau;
    out.distWeights = distWeights;
    out.distSeedsNs = seeds;
    out.distRejected = distRejected;
    out.monoTauNs = mono.tauNs;
    out.monoReducedDeviance = mono.reducedDeviance;
    out.converged = vp.converged;

    fprintf('      DistFluofit peaks: %s ns -> seeds [%.4f %.4f]\n', ...
        mat2str(round(distTau, 4)), seeds(1), seeds(2));
    fprintf(['      biexp: tau = %.4f / %.4f ns, photon share %.3f / %.3f, ' ...
        'B %.1f/bin, redDev %.1f\n'], taus(1), taus(2), ...
        out.photonFraction(1), out.photonFraction(2), out.background, ...
        out.reducedDeviance);
    fprintf('      mono for comparison: %.4f ns, redDev %.1f\n', ...
        mono.tauNs, mono.reducedDeviance);
end
