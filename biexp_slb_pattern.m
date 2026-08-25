function pattern = biexp_slb_pattern(irf, dtNs, periodNs, tauNs, nBin)
%BIEXP_SLB_PATTERN Unit-sum IRF-convolved periodic exponential.
%
% Four laser periods are summed before convolution so incomplete decay from
% previous pulses is represented, and the convolution tail is folded back into
% one period, which is what a TCSPC histogram physically does.
%
% Unit-sum normalisation matters for interpretation: with it, the fitted
% amplitude of a component is proportional to its PHOTON count, so amplitude
% ratios are directly comparable between components of different lifetime.

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
end
