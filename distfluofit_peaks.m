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
% A UNITS TRAP WORTH NAMING. DistFluofit's tau output is already in NANOSECONDS,
% not channels: line 51 builds a rate as (1/dt)./exp(...) and line 103 inverts
% it. Multiplying by dt "to convert" makes the seeds about 6x too small at this
% bin width.

    if nargin < 5 || isempty(shiftBins); shiftBins = [-3 3]; end
    if nargin < 6 || isempty(nNodes); nNodes = 100; end

    irf = double(irf(:));
    decay = double(decay(:));

    [cx, tau] = DistFluofit(irf, decay, periodNs, dtNs, shiftBins, 0, 0, nNodes);

    cx = cx(:)';
    tau = tau(:)';
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
