function out = fit_pie_discrete(decay, basis, modes, tauGrid, K, pulseNames)
%FIT_PIE_DISCRETE Discrete multi-exponential fit of a PIE period, both pulses.
%
% out = fit_pie_discrete(decay, basis, modes, tauGrid, K, pulseNames)
%
% Fits K components per excitation pulse plus a constant dark term, with the
% lifetimes free and seeded from the distribution modes. Amplitudes come from a
% whitened non-negative Poisson solve at every trial point, so only the
% lifetimes are searched (variable projection).
%
% Seeded from the distribution rather than from constants because that is the
% whole reason the distribution is fitted first: a two-colour decay has
% components 25 ns apart in excitation time and several decades apart in
% amplitude, and a search started at arbitrary values can settle on a
% near-degenerate pair and miss a real component entirely.

    nPulse = numel(basis.shifts);
    if nargin < 6 || isempty(pulseNames)
        pulseNames = arrayfun(@(k) sprintf('pulse%d', k), 1:nPulse, ...
            'UniformOutput', false);
    end

    % Seeds: the K strongest modes of each pulse, spread if there are too few.
    seeds = zeros(nPulse, K);
    for p = 1:nPulse
        mine = modes(strcmp({modes.pulse}, pulseNames{p}));
        if isempty(mine)
            seeds(p, :) = logspace(log10(tauGrid(1) * 2), ...
                log10(tauGrid(end) / 4), K);
        else
            [~, order] = sort([mine.photons], 'descend');
            take = [mine(order).tauNs];
            if numel(take) >= K
                seeds(p, :) = sort(take(1:K));
            else
                pad = logspace(log10(max(take) * 2), ...
                    log10(tauGrid(end) / 4), K - numel(take));
                seeds(p, :) = sort([take, pad]);
            end
        end
    end

    x0 = log(reshape(seeds', 1, []));
    objective = @(x) pieDiscreteObjective(x, decay, basis, nPulse, K);
    options = optimset('Display', 'off', 'TolX', 1e-5, 'TolFun', 1e-5, ...
        'MaxFunEvals', 4000, 'MaxIter', 3000);
    [xBest, devBest] = fminsearch(objective, x0, options);

    [~, beta, model, design] = pieDiscreteObjective(xBest, decay, basis, ...
        nPulse, K);
    taus = reshape(exp(xBest), K, nPulse)';
    amps = reshape(beta(2:end), K, nPulse)';
    dof = max(1, basis.nBin - numel(beta) - numel(xBest));

    out = struct();
    out.tauNs = taus;
    out.photons = amps;
    out.photonFraction = amps / max(sum(amps(:)), eps);
    out.background = beta(1);
    out.backgroundPhotons = beta(1) * basis.nBin;
    out.model = model;
    out.residual = (decay - model) ./ sqrt(max(model, 1));
    out.deviance = devBest;
    out.reducedDeviance = devBest / dof;
    out.seedsNs = seeds;
    out.design = design;
    out.pulseNames = pulseNames;

    fprintf('  DISCRETE FIT, %d component(s) per pulse\n', K);
    for p = 1:nPulse
        for k = 1:K
            fprintf('      %-8s tau %8.4f ns  %10.4g photons  %6.2f%%\n', ...
                pulseNames{p}, taus(p, k), amps(p, k), ...
                100 * out.photonFraction(p, k));
        end
    end
    fprintf(['      background %.2f/bin | reduced deviance %.1f | ' ...
        'max |residual| %.1f sigma\n'], out.background, ...
        out.reducedDeviance, max(abs(out.residual)));
end
