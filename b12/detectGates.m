function [gateStarts, gateStops, gateLength, profile, smoothProfile, threshold, ...
        baseline, riseBins, peakBins, peakHeights, rawStarts, rawStops, rawScores] = ...
        detectGates(histogram, laserCount, thresholdFraction, preBins, ...
        minimumSeparation, weakPeakFraction)
%DETECTGATES Detect PIE windows from the 15%% TCSPC rising edge.
% Gates begin preBins before each independently refined rising edge. This is
% the same convention used by AutoDetectTimeGates in Luminosa_GUI.

    profile = max(mean(double(histogram), 2), 0);
    if numel(profile) < 2 || all(profile <= 0)
        error('No TCSPC signal is available for automatic gate detection.');
    end

    smoothBins = max(3, min(21, 2 * floor(numel(profile) / 200) + 1));
    smoothProfile = conv(profile, ones(smoothBins, 1) / smoothBins, 'same');
    sortedProfile = sort(smoothProfile);
    baseline = mean(sortedProfile(1:max(1, round(0.2 * numel(sortedProfile)))));
    threshold = baseline + thresholdFraction * (max(smoothProfile) - baseline);

    % First try connected above-threshold pulse windows.
    edges = diff([false; smoothProfile > threshold; false]);
    rawStarts = find(edges == 1);
    rawStops = find(edges == -1) - 1;
    rawScores = zeros(numel(rawStarts), 1);
    for k = 1:numel(rawStarts)
        rawScores(k) = sum(smoothProfile(rawStarts(k):rawStops(k)));
    end

    % If thresholding merged or missed pulses, select separated local maxima.
    if numel(rawStarts) < laserCount
        isPeak = false(size(smoothProfile));
        isPeak(2:end-1) = smoothProfile(2:end-1) >= smoothProfile(1:end-2) & ...
            smoothProfile(2:end-1) >= smoothProfile(3:end);
        candidates = find(isPeak & smoothProfile > 0);
        [~, order] = sort(smoothProfile(candidates), 'descend');
        candidates = candidates(order);
        selected = zeros(0, 1);
        for k = 1:numel(candidates)
            if isempty(selected) || all(abs(candidates(k) - selected) >= minimumSeparation)
                selected(end + 1, 1) = candidates(k); %#ok<AGROW>
            end
            if numel(selected) == laserCount
                break;
            end
        end
        selected = sort(selected);
        halfWidth = max(5, floor(minimumSeparation / 2));
        rawStarts = max(1, selected - halfWidth);
        rawStops = min(numel(smoothProfile), selected + halfWidth);
        rawScores = smoothProfile(selected);
    end
    if numel(rawStarts) < laserCount
        error('Auto gate detection found only %d of %d laser pulses.', numel(rawStarts), laserCount);
    end

    % Keep the strongest requested pulses, then restore chronological order.
    [~, scoreOrder] = sort(rawScores, 'descend');
    keep = scoreOrder(1:laserCount);
    [rawStarts, chronologicalOrder] = sort(rawStarts(keep));
    rawStops = rawStops(keep(chronologicalOrder));
    rawScores = rawScores(keep(chronologicalOrder));

    % Find each peak and its own threshold crossing on the rising edge.
    peakBins = zeros(laserCount, 1);
    peakHeights = zeros(laserCount, 1);
    riseBins = zeros(laserCount, 1);
    for k = 1:laserCount
        [peakHeights(k), relativePeak] = max(smoothProfile(rawStarts(k):rawStops(k)));
        peakBins(k) = rawStarts(k) + relativePeak - 1;
        if k == 1
            leftBound = 1;
        else
            leftBound = max(1, floor((peakBins(k - 1) + peakBins(k)) / 2));
        end
        localThreshold = baseline + thresholdFraction * (peakHeights(k) - baseline);
        lastBelow = find(smoothProfile(leftBound:peakBins(k)) <= localThreshold, 1, 'last');
        if isempty(lastBelow)
            riseBins(k) = rawStarts(k);
        else
            riseBins(k) = min(peakBins(k), leftBound + lastBelow);
        end
    end

    gateStarts = max(1, riseBins - preBins);
    nextStart = [gateStarts(2:end) - 1; numel(smoothProfile)];
    gateLength = max(1, min(nextStart - gateStarts + 1));
    gateStops = gateStarts + gateLength - 1;

    minimumPeak = baseline + weakPeakFraction * (max(peakHeights) - baseline);
    if min(peakHeights) < minimumPeak
        error('At least one detected laser peak is too weak for reliable gating.');
    end
end
