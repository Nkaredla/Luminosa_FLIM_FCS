function fit = fit_biexp_slb_soft(decay, irf, dtNs, periodNs, opts)
%FIT_BIEXP_SLB_SOFT Biexponential with a SOFT-FIXED SLB component, plus background.
%
% fit = fit_biexp_slb_soft(decay, irf, dtNs, periodNs)
% fit = fit_biexp_slb_soft(decay, irf, dtNs, periodNs, opts)
%
% Fits
%
%     y(t) = B + a1 * P(t; tau1) + a2 * P(t; tau2)
%
% where P(t; tau) is a periodic exponential convolved with the IRF and
% normalised to unit sum, tau1 is the SLB lifetime held SOFTLY near a reference
% value, tau2 is free, and B is a constant background. All three of B, a1 and a2
% are non-negative and are solved together by PIRLSnonneg.
%
% WHY TWO COMPONENTS AND NOT THREE
%
% Because three cannot be distinguished from two in this data. Measured, not
% assumed: the per-pixel explorer shows the third component contributing 0-10%
% amplitude with an IDENTICAL maximum residual to the biexponential, at 1x1, at
% 5x5 with 45k photons, and at 15x15 with 214k photons. Earlier simulation work
% in this project found the same thing from the other direction - the false
% three-component rate RISES with photon count, from 9% to 54%, because the
% likelihood cannot separate a real third exponential from two on-grid ones
% approximating an off-grid one.
%
% So the honest model here has two components, and this function is the
% replacement for the three-model comparison rather than an addition to it.
%
% WHY THE SLB IS SOFT-FIXED AND NOT HARD-FIXED
%
% Hard-fixing tau1 makes the SLB reference an assumption that no pixel can
% contradict, and that reference is not solid: in the saved analyses the SLB
% reference fit itself reports a reduced chi-square of about 8400, so its
% 0.355 ns is a value from a poor fit rather than a measured constant. A soft
% constraint keeps the benefit - the short component cannot wander off to
% absorb the IRF or the background - while letting a pixel that genuinely
% disagrees say so, and the disagreement is reported in sigma so it is visible.
%
% The penalty is Gaussian in tau1:
%
%     objective = PoissonDeviance + ((tau1 - slbTauNs) / slbSigmaNs)^2
%
% which is the log-posterior for a Gaussian prior on tau1, so slbSigmaNs is a
% genuine prior width rather than a tuning knob. At the default 0.05 ns a pixel
% must gain 4 units of deviance to move tau1 by one sigma.
%
% opts fields
%   slbTauNs      centre of the prior, ns (default 0.3549, the saved reference)
%   slbSigmaNs    prior width, ns (default 0.05)
%   tau2BoundsNs  search bounds for the free component (default [0.5 8])
%   firstBin      first bin to fit (default 1: the IRF is reconvolved, so the
%                 rise carries information and is not discarded)
%   restarts      seeds for tau2 (default 4, log-spaced across the bounds)
%
% OUTPUT fit
%   tau1Ns, tau2Ns      fitted lifetimes
%   slbPullSigma        (tau1 - slbTauNs) / slbSigmaNs, signed
%   amplitude           [a1 a2], species fractions summing to 1
%   photonFraction      [f1 f2], the shares of the fitted photons
%   background          B, counts per bin
%   backgroundFraction  B * nBins / total photons
%   model, residual     the fitted curve and its Poisson residual
%   deviance, reducedDeviance, chiSquare
%   photons, ok, message

    if nargin < 5 || isempty(opts); opts = struct(); end
    defaults = struct('slbTauNs', 0.3549, 'slbSigmaNs', 0.05, ...
        'tau2BoundsNs', [0.5 8], 'firstBin', 1, 'restarts', 4, ...
        'minPhotons', 40);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if exist('PIRLSnonneg', 'file') ~= 2
        error('fit_biexp_slb_soft:NoPIRLS', ...
            'PIRLSnonneg.m must be on the MATLAB path.');
    end

    decay = double(decay(:));
    nBin = numel(decay);
    irf = double(irf(:));
    if numel(irf) ~= nBin
        irf = interp1(linspace(0, 1, numel(irf))', irf, ...
            linspace(0, 1, nBin)', 'linear', 0);
    end
    irf = max(irf, 0);
    if sum(irf) > 0; irf = irf / sum(irf); end

    first = max(1, min(round(opts.firstBin), nBin - 8));
    y = decay(first:end);
    photons = sum(y);
    fit = struct('tau1Ns', NaN, 'tau2Ns', NaN, 'slbPullSigma', NaN, ...
        'amplitude', [NaN NaN], 'photonFraction', [NaN NaN], ...
        'background', NaN, 'backgroundFraction', NaN, ...
        'model', [], 'residual', [], 'deviance', NaN, ...
        'reducedDeviance', NaN, 'chiSquare', NaN, 'photons', photons, ...
        'firstBin', first, 'ok', false, 'message', '');
    if photons < opts.minPhotons
        fit.message = sprintf('only %.0f photons', photons);
        return;
    end

    % Two nonlinear parameters only; the three amplitudes are linear and come
    % from PIRLSnonneg at every trial pair.
    objective = @(p) biexp_slb_objective(p, y, irf, dtNs, periodNs, ...
        first, nBin, opts);
    seeds = logspace(log10(opts.tau2BoundsNs(1) * 1.2), ...
        log10(opts.tau2BoundsNs(2) * 0.6), max(1, opts.restarts));
    best = []; bestValue = inf;
    options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-4, ...
        'MaxFunEvals', 600, 'MaxIter', 400);
    for s = 1:numel(seeds)
        start = [log(opts.slbTauNs), log(seeds(s))];
        [candidate, value] = fminsearch(objective, start, options);
        if value < bestValue; bestValue = value; best = candidate; end
    end
    if isempty(best)
        fit.message = 'the search returned nothing';
        return;
    end

    tau1 = exp(best(1));
    tau2 = exp(best(2));
    [~, beta, model] = biexp_slb_objective(best, y, irf, dtNs, periodNs, ...
        first, nBin, opts);

    amplitude = beta(2:3)';
    share = amplitude / max(sum(amplitude), eps);
    % Photon share weights the amplitude by the lifetime, which is what the
    % brightness of each population actually contributes to the decay.
    weighted = amplitude .* [tau1 tau2];
    positive = model > 0;
    dof = max(1, nnz(positive) - 5);   % 2 lifetimes + 2 amplitudes + background

    fit.tau1Ns = tau1;
    fit.tau2Ns = tau2;
    fit.slbPullSigma = (tau1 - opts.slbTauNs) / opts.slbSigmaNs;
    fit.amplitude = share;
    fit.photonFraction = weighted / max(sum(weighted), eps);
    fit.background = beta(1);
    fit.backgroundFraction = beta(1) * numel(y) / max(photons, 1);
    fit.model = model;
    fit.residual = (y - model) ./ sqrt(max(model, 1));
    fit.deviance = biexp_slb_deviance(y, model);
    fit.reducedDeviance = fit.deviance / dof;
    fit.chiSquare = sum((y(positive) - model(positive)) .^ 2 ./ ...
        model(positive)) / dof;
    fit.ok = isfinite(tau1) && isfinite(tau2) && tau1 > 0 && tau2 > 0;
    if ~fit.ok; fit.message = 'non-physical lifetime'; end
end
