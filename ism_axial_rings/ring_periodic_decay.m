function pattern = ring_periodic_decay(irf, timeNs, periodNs, tauNs)
%RING_PERIODIC_DECAY IRF-convolved exponential folded into one laser period.
%
% pattern = ring_periodic_decay(irf, timeNs, periodNs, tauNs)
%
% Returns a unit-sum decay pattern on the same bins as timeNs. Four laser
% periods are summed before convolution so that incomplete decay from previous
% pulses is represented, and the convolution tail is folded back into the
% period, which is what a TCSPC histogram physically does.

    nBins = numel(timeNs);
    decay = zeros(nBins, 1);
    for repeat = 0:3
        decay = decay + exp(-(timeNs(:) + repeat * periodNs) / tauNs);
    end
    full = conv(irf(:), decay);
    pattern = zeros(nBins, 1);
    for start = 1:nBins:numel(full)
        stop = min(start + nBins - 1, numel(full));
        span = stop - start + 1;
        pattern(1:span) = pattern(1:span) + full(start:stop);
    end
    total = sum(pattern);
    if total > 0; pattern = pattern / total; end
end
