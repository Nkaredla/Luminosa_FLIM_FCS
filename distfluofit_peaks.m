function [peakTauNs, peakWeight, tauNs, weight] = distfluofit_peaks(irf, ...
        decay, periodNs, dtNs, shiftBins, nNodes)
%DISTFLUOFIT_PEAKS Lifetime-distribution peaks, safe at the ends of the grid.
%
% [peakTauNs, peakWeight, tauNs, weight] = distfluofit_peaks(irf, decay, ...
%       periodNs, dtNs, shiftBins, nNodes)
%
% Runs DistFluofit to get a lifetime DISTRIBUTION and returns its modes, sorted
% by weight. Used to seed a two-component fit so the search starts where the
% components actually are.
%
% WHY NOT JUST CALL DistFluofit WITH flag > 0
%
% Its own peak extraction cannot handle a peak at the FIRST node. It locates
% peaks by rising edges of cx > 0.2*max(cx):
%
%     t1 = t(tmp(2:end) > tmp(1:end-1)) + 1;
%     ...
%     if t1(1) > t2(1)
%
% If the distribution is already above threshold at index 1 there is no rising
% edge, t1 comes back empty, and t1(1) errors with "Index exceeds the number of
% array elements". That is not a corner case here: the bare SLB's dominant
% component sits at the short end of the lifetime grid, so it fires every time.
% DistFluofit is left unmodified - other code depends on it - and the
% distribution is taken raw (flag = 0) and reduced here instead.
%
% THE UNITS TRAP. DistFluofit's tau output is a RATE, not a lifetime, unless you
% asked it to plot. Line 51 builds
%
%     tau = (1/dt)./exp((0:N)/N*log(p/dt))       % spans 1/dt down to 1/p
%
% and the inversion  tau = 1./tau  is at line 103, INSIDE the "if bild" plotting
% block. Called with bild = 0 - which is what any programmatic use does - the
% array comes back as rates in 1/ns. Inverting here gives lifetimes spanning
% dt to p, which is the sensible range: one bin width up to one full period.
%
% Reading those rates as lifetimes inverts the whole interpretation: a
% distribution piled up at the SHORTEST lifetime looks like one piled up at the
% LONGEST, and a perfectly sensible answer looks degenerate.
%
% THE GRID CANNOT GO BELOW ONE BIN. Its short end is dt by construction, so a
% lifetime shorter than a bin can only be reported as "the lowest node". That is
% not a failure, it is the grid's resolution limit, and it happens here: the SLB
% is near 0.08 ns against dt = 0.16 ns. DistFluofit's own source carries a
% commented-out alternative one line above the live one,
%
%     %tau = (1/dt/10)./exp((0:N)/N*log(p/dt/10));
%
% which shifts the grid down a decade for exactly this situation.
%
% A NOTE ON PERIOD. p sets only the grid's upper end. The incomplete-decay
% correction in the basis uses dt*n - the stored WINDOW - not p:
%
%     Convol(irf, exp(-t*dt*tau)) ./ (1 - exp(-dt*n*tau))
%
% so with a PIE half-window (25 ns stored, 50 ns period) the two disagree. For
% tau below a couple of ns both factors are 1 to within 1e-5, so it does not
% matter for the SLB or the membrane; it matters only at the long end.

    if nargin < 5 || isempty(shiftBins); shiftBins = [-3 3]; end
    if nargin < 6 || isempty(nNodes); nNodes = 100; end

    irf = double(irf(:));
    decay = double(decay(:));

    [cx, rate] = DistFluofit(irf, decay, periodNs, dtNs, shiftBins, 0, 0, ...
        nNodes);

    cx = cx(:)';
    % Rates -> lifetimes. See the header: the inversion inside DistFluofit only
    % runs when it plots.
    tau = 1 ./ rate(:)';
    m = min(numel(cx), numel(tau));
    weight = max(cx(1:m), 0);
    tauNs = tau(1:m);
    [tauNs, order] = sort(tauNs);
    weight = weight(order);
    if sum(weight) > 0
        weight = weight / sum(weight);
    end

    % Local maxima INCLUDING both ends, which is the whole point.
    isPeak = false(1, m);
    for k = 1:m
        left = weight(max(k - 1, 1));
        right = weight(min(k + 1, m));
        isPeak(k) = weight(k) >= left && weight(k) >= right && weight(k) > 0;
    end
    % Drop peaks that are negligible next to the largest.
    isPeak = isPeak & weight > 0.05 * max(weight);

    peakTauNs = tauNs(isPeak);
    peakWeight = weight(isPeak);
    % A peak on the lowest node means the true lifetime is at or below the grid's
    % short end (one bin). Report it as such rather than as a measurement AT that
    % node, so a caller can seed below it instead of being anchored to dt.
    if ~isempty(peakTauNs) && abs(peakTauNs(1) - tauNs(1)) < 1e-12
        peakTauNs(1) = tauNs(1);
    end
    [peakWeight, order] = sort(peakWeight, 'descend');
    peakTauNs = peakTauNs(order);

    % If the distribution is unimodal there is nothing to seed a two-component
    % fit with, so split it at its weighted median and take a centroid either
    % side. That gives two distinct, data-driven starting values rather than a
    % pair of arbitrary constants.
    if numel(peakTauNs) < 2 && sum(weight) > 0
        cumw = cumsum(weight);
        split = find(cumw >= 0.5, 1, 'first');
        split = max(2, min(m - 1, split));
        lowW = weight(1:split);
        highW = weight(split + 1:end);
        low = sum(tauNs(1:split) .* lowW) / max(sum(lowW), eps);
        high = sum(tauNs(split + 1:end) .* highW) / max(sum(highW), eps);
        if isfinite(low) && isfinite(high) && high > low
            peakTauNs = [low high];
            peakWeight = [sum(lowW) sum(highW)];
        end
    end
end
