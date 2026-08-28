function [P, total, dP] = biexp_slb_pattern_batch(basis, tau)
%BIEXP_SLB_PATTERN_BATCH Unit-sum patterns, and their tau-derivatives, for many taus.
%
% [P, total]     = biexp_slb_pattern_batch(basis, tau)
% [P, total, dP] = biexp_slb_pattern_batch(basis, tau)
%
% BASIS comes from biexp_slb_basis. TAU is 1-by-P. P is nBin-by-P and matches
% biexp_slb_pattern(irf, dt, period, tau(k), nBin) column by column. TOTAL is
% 1-by-P, the PRE-normalisation sum, which is what turns a fitted amplitude into
% a pre-exponential species weight. DP is nBin-by-P, d(pattern)/d(tau).
%
% THE DERIVATIVE
%
% With u(tau) = C * exp(-t/tau) the un-normalised pattern and S = sum(u),
%
%     P  = u / S
%     dP = u' / S - u * S' / S^2,       S' = sum(u')
%
% and u' = C * ((t / tau^2) .* exp(-t / tau)), since d/dtau exp(-t/tau) is
% (t/tau^2) exp(-t/tau). Both u and u' come from the SAME shared circulant, so
% the derivative costs one extra matrix product rather than a new build - which
% is why supplying an analytic gradient to the outer optimizer is cheap.
%
% The four-period factor cancels in P (it is a scalar multiple removed by the
% unit-sum normalisation) but not in TOTAL, so it is applied only there.

    % Match the basis's class rather than forcing double, so a single-precision
    % or gpuArray basis stays that way instead of being silently promoted.
    tau = cast(reshape(tau, 1, []), 'like', basis.C);
    if any(tau <= 0)
        error('biexp_slb_pattern_batch:NonPositiveTau', ...
            'All lifetimes must be positive; min is %.4g.', min(tau));
    end
    t = basis.timeNs;

    E = exp(-t ./ tau);              % nBin-by-P
    U = basis.C * E;
    S = sum(U, 1);                   % 1-by-P
    S(S <= 0) = eps;
    P = U ./ S;

    if nargout > 1
        total = S .* basis.periodFactor(tau);
    end

    if nargout > 2
        dE = (t ./ (tau .^ 2)) .* E;
        dU = basis.C * dE;
        dS = sum(dU, 1);
        dP = dU ./ S - U .* (dS ./ (S .^ 2));
    end
end
