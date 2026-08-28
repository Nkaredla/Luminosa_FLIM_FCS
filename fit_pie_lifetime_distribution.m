function out = fit_pie_lifetime_distribution(decay, irfFull, dtNs, periodNs, ...
        pulseBins, opts)
%FIT_PIE_LIFETIME_DISTRIBUTION All lifetimes in a PIE period, split by pulse.
%
% out = fit_pie_lifetime_distribution(decay, irfFull, dtNs, periodNs, pulseBins)
% out = fit_pie_lifetime_distribution(..., opts)
%
% Fits the FULL-PERIOD decay with a lifetime DISTRIBUTION for each excitation
% pulse simultaneously:
%
%     y(t) = D + sum_i a_i P(t; tau_i | pulse 1) + sum_j b_j P(t; tau_j | pulse 2)
%
% All amplitudes non-negative, D a genuine constant dark term. Nothing is
% assumed about how many components there are or which pulse they belong to -
% the fit decides, and every lifetime is reported.
%
% THIS EXISTS TO ANSWER ONE QUESTION: WHAT IS ACTUALLY IN THE DECAY?
%
% Fitting the 640 gate alone forces a choice that the data does not support. The
% 485-excited emission has not finished decaying when that gate opens - there is
% a flat pedestal of ~70-100 counts/bin before the 640 pulse, and in SLB-only
% pixels the tail falls from ~300 to ~100 counts/bin over the last 9 ns, which
% is a component near 8 ns, not an offset. A gated model with a constant
% background has nowhere to put that, so the background term absorbs it and
% every amplitude is corrupted - and with unit-sum patterns those amplitudes ARE
% the per-pixel photon counts.
%
% Over the full period the two colours separate because each has its own rise:
% a component excited at 485 rises at the 485 pulse, 25 ns before the 640 one.
% That is the information a single gate throws away.
%
% A DISTRIBUTION RATHER THAN A FIXED NUMBER OF EXPONENTIALS
%
% Because the question is what is there, not how well two components fit. Fixing
% the count first is what produced the earlier trouble: three components could
% not be distinguished from two, yet the third was still being reported. A
% non-negative fit over a dense lifetime grid puts weight only where the data
% needs it, so the number of components is an output.
%
% Note this is a regularisation-free NNLS over a fine grid, so it is not a
% smooth distribution - NNLS concentrates weight onto a few nodes. Read the
% MODES, not the width; the widths are an artefact of the basis.
%
% opts fields
%   tauRangeNs   grid span (default [dtNs/4, periodNs/2])
%   nNodes       lifetime nodes per pulse (default 60)
%   pulseNames   labels for reporting (default {'pulse1','pulse2'})
%   discreteK    also fit this many DISCRETE components per pulse, seeded from
%                the distribution modes (default 0, meaning skip)
%
% out fields
%   tauNs                 the lifetime grid
%   weightsByPulse        nPulse-by-nNodes photon counts
%   photonFractionByPulse same, normalised over everything
%   background            counts per bin
%   model, residual, deviance, reducedDeviance
%   modes                 struct array of the peaks found, per pulse
%   discrete              the discrete fit, when asked for

    if nargin < 6 || isempty(opts); opts = struct(); end
    defaults = struct('tauRangeNs', [], 'nNodes', 60, 'pulseNames', [], ...
        'discreteK', 0);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if isempty(opts.tauRangeNs)
        opts.tauRangeNs = [dtNs / 4, periodNs / 2];
    end
    if isempty(opts.pulseNames)
        opts.pulseNames = arrayfun(@(k) sprintf('pulse%d', k), ...
            1:numel(pulseBins), 'UniformOutput', false);
    end

    decay = double(decay(:));
    basis = pie_two_pulse_basis(irfFull, dtNs, periodNs, pulseBins);
    nPulse = numel(basis.shifts);
    tauNs = logspace(log10(opts.tauRangeNs(1)), log10(opts.tauRangeNs(2)), ...
        opts.nNodes);

    % One column per (lifetime, pulse), plus the constant dark term.
    cols = cell(1, nPulse);
    totals = cell(1, nPulse);
    for p = 1:nPulse
        [cols{p}, totals{p}] = pie_pattern_columns(basis, tauNs, p);
    end
    design = [ones(basis.nBin, 1), cols{:}];

    fprintf(['fit_pie_lifetime_distribution\n  %d bins, dt %.4f ns, period ' ...
        '%.2f ns, pulses at bins %s\n  %d lifetime node(s) per pulse over ' ...
        '[%.4f %.4f] ns, %d column(s)\n'], basis.nBin, dtNs, periodNs, ...
        mat2str(basis.pulseBins), opts.nNodes, tauNs(1), tauNs(end), ...
        size(design, 2));

    [beta, dev] = poisson_nnls_whitened(design, decay);
    background = beta(1);
    weights = reshape(beta(2:end), opts.nNodes, nPulse)';

    model = max(design * beta, 1e-12);
    dof = max(1, basis.nBin - nnz(beta) - 1);

    out = struct();
    out.tauNs = tauNs;
    out.weightsByPulse = weights;
    out.photonFractionByPulse = weights / max(sum(weights(:)), eps);
    out.background = background;
    out.backgroundPhotons = background * basis.nBin;
    out.model = model;
    out.residual = (decay - model) ./ sqrt(max(model, 1));
    out.deviance = dev;
    out.reducedDeviance = dev / dof;
    out.decay = decay;
    out.basis = basis;
    out.patternTotals = totals;
    out.pulseNames = opts.pulseNames;
    out.photons = sum(decay);

    % ---- report every mode, per pulse ------------------------------------
    fprintf(['  background %.2f counts/bin (%.3g photons, %.2f%% of the ' ...
        'total)\n  reduced deviance %.1f, max |residual| %.1f sigma\n'], ...
        background, out.backgroundPhotons, ...
        100 * out.backgroundPhotons / max(out.photons, 1), ...
        out.reducedDeviance, max(abs(out.residual)));
    modes = struct('pulse', {}, 'tauNs', {}, 'photons', {}, 'fraction', {});
    for p = 1:nPulse
        w = weights(p, :);
        share = sum(w) / max(sum(weights(:)), eps);
        fprintf('  %s: %.3g photons (%.1f%% of all fitted)\n', ...
            opts.pulseNames{p}, sum(w), 100 * share);
        active = find(w > 0.005 * max(max(weights), [], 'all'));
        for a = active
            isPeak = (a == 1 || w(a) >= w(a - 1)) && ...
                (a == numel(w) || w(a) >= w(a + 1));
            if ~isPeak; continue; end
            fprintf('      tau %8.4f ns   %10.4g photons   %6.2f%%\n', ...
                tauNs(a), w(a), 100 * w(a) / max(sum(weights(:)), eps));
            modes(end + 1) = struct('pulse', opts.pulseNames{p}, ...
                'tauNs', tauNs(a), 'photons', w(a), ...
                'fraction', w(a) / max(sum(weights(:)), eps)); %#ok<AGROW>
        end
    end
    out.modes = modes;

    % ---- optional discrete fit, seeded from the modes --------------------
    if opts.discreteK > 0
        out.discrete = fit_pie_discrete(decay, basis, modes, tauNs, ...
            opts.discreteK, opts.pulseNames);
    else
        out.discrete = [];
    end
end
