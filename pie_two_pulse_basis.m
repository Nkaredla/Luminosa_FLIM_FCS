function basis = pie_two_pulse_basis(irfFull, dtNs, periodNs, pulseBins)
%PIE_TWO_PULSE_BASIS Pattern builder for a PIE period with two excitation pulses.
%
% basis = pie_two_pulse_basis(irfFull, dtNs, periodNs, pulseBins)
%
% IRFFULL is the instrument response over ONE FULL PERIOD, positioned at the
% reference pulse. PULSEBINS lists the bin index of each excitation pulse; the
% FIRST entry is the reference the IRF is aligned to.
%
% WHY A SINGLE SHIFTED PATTERN IS ENOUGH
%
% Over one full laser period the histogram is circular: a photon arriving late
% wraps to the start of the next period, which is what a TCSPC histogram
% physically records. A component excited at pulse B is therefore exactly the
% same shape as one excited at pulse A, circularly shifted by the pulse
% separation. So one pattern per lifetime is built and then shifted, rather than
% convolving separately per pulse.
%
% This also means the pattern is a fixed circulant acting on a geometric vector -
% the same structure the per-pixel fits use - so a whole lifetime grid costs one
% matrix product.
%
% WHY THE FULL PERIOD, NOT A GATE
%
% Because the two colours cannot be separated inside one gate. The 485-excited
% emission is still decaying when the 640 gate opens, and inside that gate it
% looks like an offset, so a gated model with a constant background silently
% absorbs it. Over the full period the 485 component has its own rise and its
% own decay and is identifiable.
%
% FIELDS
%   C           nBin-by-nBin circulant of irfFull
%   shifts      bin offset of each pulse from the reference
%   nBin, dtNs, periodNs, timeNs, pulseBins

    irfFull = double(irfFull(:));
    irfFull = max(irfFull, 0);
    if sum(irfFull) > 0; irfFull = irfFull / sum(irfFull); end
    nBin = numel(irfFull);

    C = zeros(nBin, nBin);
    for j = 1:nBin
        C(:, j) = circshift(irfFull, j - 1);
    end

    pulseBins = round(pulseBins(:)');
    basis = struct();
    basis.C = C;
    basis.nBin = nBin;
    basis.dtNs = dtNs;
    basis.periodNs = periodNs;
    basis.timeNs = (0:nBin - 1)' * dtNs;
    basis.pulseBins = pulseBins;
    basis.shifts = pulseBins - pulseBins(1);
    basis.irfFull = irfFull;
end
