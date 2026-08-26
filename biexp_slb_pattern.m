function [pattern, total] = biexp_slb_pattern(irf, dtNs, periodNs, tauNs, nBin)
%BIEXP_SLB_PATTERN Unit-sum IRF-convolved periodic exponential.
%
% Four laser periods are summed before convolution so incomplete decay from
% previous pulses is represented, and the convolution tail is folded back into
% one period, which is what a TCSPC histogram physically does.
%
% Unit-sum normalisation matters for interpretation: with it, the fitted
% amplitude of a component EQUALS its photon count (verified to 2e-16 by
% recovering 200 and 400 from a decay synthesised with those counts), so
% amplitudes are directly comparable between components of different lifetime.
%
% The second output TOTAL is the sum BEFORE that normalisation. It is what
% converts a fitted amplitude back into a pre-exponential (species) weight:
% species_i = a_i / total_i. Using tau in place of total is only asymptotically
% correct - at tau = 0.058 ns with dt = 0.16 ns, tau/dt is 0.34 times the true
% total, so that shortcut errs by a factor of ~2.4 at short lifetimes.

    timeNs = (0:nBin - 1)' * dtNs;
    decay = zeros(nBin, 1);
    for repeat = 0:3
        decay = decay + exp(-(timeNs + repeat * periodNs) / tauNs);
    end
    full = conv(irf(:), decay);
    pattern = zeros(nBin, 1);
    for start = 1:nBin:numel(full)
        stop = min(start + nBin - 1, numel(full));
        span = stop - start + 1;
        pattern(1:span) = pattern(1:span) + full(start:stop);
    end
    total = sum(pattern);
    if total > 0; pattern = pattern / total; end
    % TOTAL is returned un-normalised on purpose - see the header.
end
