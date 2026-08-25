function [value, beta, model] = biexp_slb_objective(logTau, y, irf, dtNs, ...
        periodNs, firstBin, nBin, opts)
%BIEXP_SLB_OBJECTIVE Penalised Poisson deviance for the soft-SLB biexponential.
%
% [value, beta, model] = biexp_slb_objective(logTau, y, irf, dtNs, periodNs, ...
%                                            firstBin, nBin, opts)
%
% Variable projection: the background and the two amplitudes enter linearly, so
% only the two log-lifetimes are searched and PIRLSnonneg supplies
% [B; a1; a2] >= 0 at each trial pair.
%
% PIRLSnonneg rather than plain non-negative least squares because the counts
% span several decades. Unweighted least squares would let the bright bins near
% the peak dominate and leave the background essentially unconstrained, which is
% exactly the failure that produced a fitted background twice the measured
% pedestal earlier in this project.
%
% The value returned is the Poisson deviance PLUS the Gaussian prior penalty on
% tau1, so minimising it is maximising the posterior under a Gaussian prior of
% width opts.slbSigmaNs centred on opts.slbTauNs.

    tau1 = exp(logTau(1));
    tau2 = exp(logTau(2));
    bounds = opts.tau2BoundsNs;
    % tau2 must stay above tau1, otherwise the two components can swap roles and
    % the prior silently ends up constraining the LONG component instead.
    if ~isfinite(tau1) || ~isfinite(tau2) || tau1 <= 0.02 || tau1 > 3 || ...
            tau2 <= max(tau1, bounds(1)) || tau2 > bounds(2)
        value = 1e12;
        beta = [0; 0; 0];
        model = ones(numel(y), 1);
        return;
    end

    p1 = biexp_slb_pattern(irf, dtNs, periodNs, tau1, nBin);
    p2 = biexp_slb_pattern(irf, dtNs, periodNs, tau2, nBin);
    design = [ones(numel(y), 1), p1(firstBin:end), p2(firstBin:end)];
    beta = PIRLSnonneg(design, y);
    model = max(design * beta, 1e-12);

    penalty = ((tau1 - opts.slbTauNs) / max(opts.slbSigmaNs, eps)) ^ 2;
    value = biexp_slb_deviance(y, model) + penalty;
    if ~isfinite(value); value = 1e12; end
end
