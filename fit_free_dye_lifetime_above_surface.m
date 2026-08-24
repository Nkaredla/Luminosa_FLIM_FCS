function out = fit_free_dye_lifetime_above_surface(source, opts)
%FIT_FREE_DYE_LIFETIME_ABOVE_SURFACE Unquenched lifetime from pixels off the metal.
%
% out = fit_free_dye_lifetime_above_surface(ptuFileOrFolder)
% out = fit_free_dye_lifetime_above_surface(source, opts)
%
% Measures the dye's unquenched lifetime tau0 from an XZ or YZ cross-section by
% pooling the TCSPC of every pixel that sits away from the metal surface and
% carries real signal, then fitting one decay.
%
% No axial calibration is used or needed. The surface is found from the data,
% pixels beyond it are selected by count rate, and the fit is of the pooled
% decay - so nothing in the answer depends on the micron size of an axial pixel.
%
% HOW THE SURFACE IS FOUND
%
% Quenching is strongest at the metal, so the surface is the narrow band of
% LOWEST mean arrival time. It is located per lateral column rather than once
% for the whole image, because the coverslip is never perfectly level across a
% 600-pixel scan and a single global row would cut into the bilayer on one side
% and leave a quenched margin on the other.
%
% Intensity is deliberately not used to find it. The bilayer is bright, but so
% is a cell membrane, and dye density varies; the arrival-time minimum marks the
% metal by physics rather than by brightness.
%
% WHICH SIDE IS "ABOVE"
%
% Also taken from the data: the side of the band where the mean arrival time
% RISES is the side away from the metal, because that is what unquenching does.
% Photon counts on each side are printed as a cross-check, since below the metal
% there should be little but reflection.
%
% THE RATE FILTER
%
% Beyond the surface band most pixels are empty medium, and including them would
% pull the pooled decay toward background. A pixel is kept only if its count
% exceeds a threshold derived from the image itself: the background level is
% taken as a high quantile of the emptiest part of the field, and the cut is set
% a stated number of Poisson sigmas above it. The count of pixels and photons
% passing is reported, and so is the fraction of the pooled photons that the
% background could account for - because if that fraction is large, tau0 is
% biased no matter how good the fit looks.
%
% THE DETECTOR IS HARDWARE-GATED, AND THAT SETS THE FIT WINDOW
%
% Inspecting the pooled decay directly settles what the window has to be. In
% scan 154813 it runs 40128 counts at 2.44 ns, 17372 at 4.88 ns, then 11.6 at
% 7.32 ns - a factor of 1500 in 2.4 ns. No fluorescence decay falls that fast;
% that is a gate closing. A second, much weaker PIE pulse then appears near
% 9.9 ns at about 2% of the first peak.
%
% So each PIE window is only about 9 ns long and the usable decay is roughly
% the first 6-7 ns after the rise. The window is therefore ended at the FIRST
% point after the peak where the decay falls to the pedestal, which is the gate
% edge or true extinction, whichever comes first.
%
% Two earlier attempts got this wrong in opposite directions and both are worth
% recording. Ending the window "where signal stops being significant" ran it to
% ~30 ns, swallowing the gate edge, the second pulse and its decay. Before that,
% a second-pulse test requiring a rise above 10% of the first peak missed the
% real second pulse entirely, because it is only 2% of it.
%
% Photons beyond 15 ns are 0.39% of the total and beyond 25 ns 0.17%, so there
% is NO long-lived component here - an earlier reading of a 19% "excess" as a
% long component came from equating a mono-exponential's photon total with
% peak*tau/dt, which ignores both IRF convolution and the gate.
%
% THREE ESTIMATES ARE REPORTED, BECAUSE ONE NUMBER WOULD MISLEAD
%
% The pooled decay is NOT mono-exponential, and a single fitted lifetime hides
% that. Comparing the observed photon total against what a mono-exponential of
% the fitted lifetime would contain leaves about a fifth unaccounted for, and
% signal stays significant out to 27-31 ns where a 2.1 ns exponential would be
% at 0.3 counts against a 21-count threshold. So:
%
%   1. FIRST MOMENT, raw and TRUNCATION-CORRECTED. The raw first moment of the
%      background-subtracted decay from the peak is only the lifetime when the
%      whole decay is observed. Here the gate cuts at about 2.3 lifetimes, and
%      truncation biases it badly low: for a single exponential
%          <t> = tau - T*exp(-T/tau) / (1 - exp(-T/tau))
%      which for T = 6.15 ns turns a true 2.74 ns into an observed 2.01 ns. The
%      corrected value inverts that expression. An earlier version reported the
%      RAW moment as the headline "model-free" number, which understated tau0
%      by about 0.8 ns - a first moment is not model-free in a truncated window.
%   2. MONO-EXPONENTIAL TAIL FIT, swept over tail starts. With the window
%      correctly ended at the gate the sweep is stable to about 0.14 ns, so
%      this is now a real estimate rather than only a diagnostic.
%   3. LIFETIME DISTRIBUTION (DistTailfit), which infers the components without
%      being told how many there are, and accepts the pulse period so
%      incomplete decay from the previous pulse is modelled rather than
%      absorbed.
%   4. LOG-LINEAR REGRESSION over several ranges. Weighted least squares on
%      log(counts - pedestal), which is transparent enough to check by hand and
%      shares no code with the other two fitters. Its drift ACROSS ranges is
%      the most honest statement of how single-exponential the decay really is.
%
% A two-exponential Tailfit was tried and collapsed to 0%/100% on all five
% scans, so it is not reported; the distribution fit replaces it.
%
% THE FIT
%
% Primary is an IRF-free Poisson-MLE tail fit, with the tail start SWEPT. A
% lifetime that is stable as the start moves past the IRF is real; one that
% drifts means the tail is still contaminated. That stability check is the whole
% reason to prefer a tail fit here, and it is reported rather than hidden.
%
% An IRF-deconvolved Fluofit runs as a cross-check, seeded at a FIXED lifetime
% rather than at the tail-fit answer. An earlier version seeded it with the
% tail-fit value and got back exactly that value, to four significant figures,
% on all five scans - the optimiser had not moved, so the "cross-check" was
% only echoing its own input. It now reports whether the fit actually moved,
% and a check that did not move is reported as inconclusive rather than as
% agreement.
%
% OUTPUTS PER FILE
%   <acq>/free_dye_above_surface/mean_flim.png    intensity, mean arrival time,
%                                                 and the selected mask
%   <acq>/free_dye_above_surface/fitted_tcspc.png pooled decay, fit, residuals
%   <acq>/free_dye_above_surface/*.mat, *.csv
%
% opts
%   piePulseIndex      PIE window (default 2, the 640 nm one)
%   maxNgate           TCSPC bins after rebinning (default 512)
%   surfaceMarginPix   pixels to skip beyond the surface band (default 6)
%   bandWidthPix       assumed half-width of the low-lifetime band, used only
%                      as a floor for the margin (default 2)
%   minPixelCounts     absolute floor on counts per pixel (default 5)
%   backgroundSigmas   Poisson sigmas above background for the rate filter
%                      (default 5)
%   tailStartNs        swept tail starts (default [0.6 1.0 1.5 2.0 3.0])
%   irfSource          '' to search the session, a MAT path, or 'none'
%   outputDir, makeFigure

    if nargin < 2 || isempty(opts); opts = struct(); end
    defaults = struct( ...
        'piePulseIndex', 2, 'maxNgate', 512, 'tcspcBinNs', 0.05, ...
        'surfaceMarginPix', 6, 'bandWidthPix', 2, ...
        'minPixelCounts', 5, 'backgroundSigmas', 5, ...
        'minColumnCounts', 20, ...
        'tailStartNs', [0.4 0.8 1.2 1.8 2.5], ...
        'irfSource', '', 'outputDir', '', 'makeFigure', true, ...
        'photonsPerChunk', 5e6, 'photonBlockSize', 2e6);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    for required = {'PTU_Read_Head', 'PTU_FLIM_GPU', 'Tailfit'}
        if exist(required{1}, 'file') ~= 2
            error('fit_free_dye_lifetime_above_surface:Missing', ...
                '%s.m must be on the MATLAB path.', required{1});
        end
    end

    ptuFile = resolvePtu(source);
    [acqFolder, ~] = fileparts(ptuFile);
    [~, acqName] = fileparts(acqFolder);
    if isempty(opts.outputDir)
        opts.outputDir = fullfile(acqFolder, 'free_dye_above_surface');
    end
    if ~isfolder(opts.outputDir); mkdir(opts.outputDir); end

    head = PTU_Read_Head(ptuFile);
    scanDirection = fieldOr(head, 'ImgHdr_ScanDirection', NaN);
    switch scanDirection
        case 1; plane = 'YZ';
        case 2; plane = 'XZ';
        case 0; plane = 'XY';
        otherwise; plane = 'unknown';
    end
    if scanDirection == 0
        error('fit_free_dye_lifetime_above_surface:NotACrossSection', ...
            ['This is an XY scan (ImgHdr_ScanDirection = 0), which has no ' ...
             'axial axis and therefore no pixels away from the surface.']);
    end
    pixX = round(fieldOr(head, 'ImgHdr_PixX', NaN));
    pixY = round(fieldOr(head, 'ImgHdr_PixY', NaN));
    syncRate = fieldOr(head, 'TTResult_SyncRate', NaN);
    periodNs = 1e9 / syncRate;
    dwellMs = fieldOr(head, 'ImgHdr_TimePerPixel', NaN);

    fprintf('\n=============================================================\n');
    fprintf('fit_free_dye_lifetime_above_surface  |  %s  (%s)\n', acqName, plane);
    fprintf('  %s\n', ptuFile);
    fprintf('  %d x %d pixels, period %.2f ns, dwell %.4f ms/pixel\n', ...
        pixX, pixY, periodNs, dwellMs);

    % ---- photons -----------------------------------------------------------
    ptu = PTU_FLIM_GPU(ptuFile, struct( ...
        'photonsPerChunk', opts.photonsPerChunk, 'computePerFrame', false, ...
        'storeTcspcPix', false, 'storePhotonLists', true, ...
        'storeTimeCell', false, 'showWaitbar', false, ...
        'minLifetimeBin_ns', opts.tcspcBinNs, 'maxNgate', opts.maxNgate, ...
        'computeGlobalMaps', false, 'storePhotonSync', false, ...
        'storePhotonFrame', false, 'lowMemoryPhotonLists', true, ...
        'photonBlockSize', opts.photonBlockSize));
    nGate = double(ptu.Ngate);
    dtNs = periodNs / nGate;
    fprintf('  %d photons, %d TCSPC bins of %.4f ns\n', ...
        numel(ptu.im_tcspc), nGate, dtNs);

    % ---- PIE gate ----------------------------------------------------------
    % im_tcspc is ALREADY a 1-based bin index in this reader - see
    % immune_cell_MIET_reassigned_sliding_tcspc, which tests
    % "timeBin >= 1 & timeBin <= gateLength" and offsets with (timeBin - 1).
    % Adding one here pushed the maximum to Ngate+1 and accumarray refused it.
    rawBin = double(ptu.im_tcspc(:));
    binLow = min(rawBin); binHigh = max(rawBin);
    fprintf('  im_tcspc spans bins %g to %g against Ngate = %d\n', ...
        binLow, binHigh, nGate);
    if binHigh > nGate || binLow < 1
        error('fit_free_dye_lifetime_above_surface:BinRange', ...
            ['im_tcspc spans %g to %g, outside the 1..%d range implied by ' ...
             'Ngate. The bin convention is not what this code assumes.'], ...
            binLow, binHigh, nGate);
    end
    globalCurve = accumarray(rawBin, 1, [nGate, 1]);
    gate = pickGate(globalCurve, dtNs);
    gateLength = gate.length;
    fprintf(['  gate: peak at bin %d (%.2f ns); rise rolled to bin 1; ' ...
        'period window %d bins = %.2f ns\n'], gate.peakBin, ...
        (gate.peakBin - 1) * dtNs, gateLength, gateLength * dtNs);
    fprintf('      %s\n', gate.note);

    % The excitation pulse sits near the end of the period, so the decay wraps
    % round to low bin numbers. Rolling the histogram so the rise is at bin 1
    % makes the decay contiguous; without it the tail is cut at the period
    % boundary and the fit loses the very bins that constrain a long lifetime.
    bin = mod(rawBin - gate.riseBin, nGate) + 1;
    row = double(ptu.im_line(:));
    col = double(ptu.im_col(:));
    keep = bin >= 1 & bin <= gateLength & ...
        row >= 1 & row <= pixY & col >= 1 & col <= pixX;
    bin = bin(keep);
    row = row(keep);
    col = col(keep);
    fprintf('  %d photons inside the gate and image bounds\n', numel(bin));

    counts = accumarray([row, col], 1, [pixY, pixX]);
    binSum = accumarray([row, col], bin, [pixY, pixX]);
    meanArrivalNs = nan(pixY, pixX);
    hasSignal = counts >= max(1, opts.minPixelCounts);
    meanArrivalNs(hasSignal) = (binSum(hasSignal) ./ counts(hasSignal) - 1) * dtNs;

    % ---- which axis is axial, and which side is away from the metal --------
    axial = pickAxialAxis(counts, binSum, dtNs);
    fprintf(['\n  axial axis: %s  (mean-arrival range along rows %.3f ns, ' ...
        'along cols %.3f ns)\n'], axial.name, axial.rowRangeNs, ...
        axial.colRangeNs);
    if axial.margin < 1.5
        warning('fit_free_dye_lifetime_above_surface:AmbiguousAxis', ...
            ['The two axes show similar arrival-time structure (ratio ' ...
             '%.2f); check mean_flim.png before trusting the surface ' ...
             'line.'], axial.margin);
    end
    % Work in [axial, lateral] orientation from here on.
    if axial.isRow
        countsAx = counts; arrivalAx = meanArrivalNs;
    else
        countsAx = counts'; arrivalAx = meanArrivalNs';
    end

    [surfaceIndex, direction] = locateSurface(countsAx, arrivalAx, opts);
    fprintf(['  surface: median index %.1f of %d; "away from metal" is the ' ...
        '%s direction\n'], median(surfaceIndex, 'omitnan'), ...
        size(countsAx, 1), direction.name);
    fprintf(['      mean arrival %.3f ns in the band, %.3f ns just beyond it ' ...
        'on the chosen side,\n      %.3f ns on the other side; photons ' ...
        '%.3g beyond vs %.3g on the other side\n'], ...
        direction.bandArrivalNs, direction.chosenArrivalNs, ...
        direction.otherArrivalNs, direction.chosenPhotons, ...
        direction.otherPhotons);
    if direction.arrivalGainNs < 0.05
        warning('fit_free_dye_lifetime_above_surface:WeakGradient', ...
            ['The arrival time rises by only %.3f ns away from the band, so ' ...
             'the metal side is not clearly resolved.'], ...
            direction.arrivalGainNs);
    end

    % ---- rate filter -------------------------------------------------------
    selection = selectAbovePixels(countsAx, surfaceIndex, direction, opts);
    fprintf(['\n  rate filter: background %.2f counts/pixel, cut at %.2f ' ...
        '(%g sigma), floor %d\n'], selection.backgroundCounts, ...
        selection.threshold, opts.backgroundSigmas, opts.minPixelCounts);
    fprintf(['      %d of %d pixels beyond the surface pass (%.2f%% of the ' ...
        'frame), holding %.4g photons\n'], selection.pixelCount, ...
        selection.candidateCount, ...
        100 * selection.pixelCount / numel(countsAx), selection.photons);
    fprintf('      background could account for %.2f%% of those photons\n', ...
        100 * selection.backgroundFraction);
    if selection.pixelCount < 20 || selection.photons < 2000
        error('fit_free_dye_lifetime_above_surface:TooFewPixels', ...
            ['Only %d pixels and %.4g photons survived the rate filter. ' ...
             'Lower backgroundSigmas or minPixelCounts.'], ...
            selection.pixelCount, selection.photons);
    end
    if selection.backgroundFraction > 0.2
        warning('fit_free_dye_lifetime_above_surface:BackgroundHeavy', ...
            ['Background may be %.0f%% of the pooled photons; tau0 will be ' ...
             'biased. Raise backgroundSigmas.'], ...
            100 * selection.backgroundFraction);
    end

    % ---- pool the TCSPC and fit -------------------------------------------
    if axial.isRow
        selectedMaskImage = selection.mask;
        pixelIsSelected = selection.mask(sub2ind(size(selection.mask), row, col));
    else
        selectedMaskImage = selection.mask';
        pixelIsSelected = selection.mask(sub2ind(size(selection.mask), col, row));
    end
    pooled = accumarray(bin(pixelIsSelected), 1, [gateLength, 1]);
    pedestal = measurePedestal(pooled);
    fprintf('\n  pooled decay: %.4g photons in %d bins, peak %.4g\n', ...
        sum(pooled), gateLength, max(pooled));
    fprintf(['      pedestal %.0f counts/bin; decay returns to it at bin ' ...
        '%d (%.2f ns after the rise),\n      so the PIE window is %.2f ns ' ...
        'long and the fit uses only that\n'], pedestal.level, ...
        pedestal.signalEnd, pedestal.signalEnd * dtNs, ...
        pedestal.signalEnd * dtNs);
    if isfinite(pedestal.secondPulseBin)
        fprintf(['      second PIE pulse at bin %d (%.2f ns), peak %.0f ' ...
            'counts = %.1f%% of the first;\n      it is EXCLUDED from the ' ...
            'fit\n'], pedestal.secondPulseBin, ...
            pedestal.secondPulseBin * dtNs, pedestal.secondPulseCounts, ...
            100 * pedestal.secondPulseCounts / max(max(pooled), 1));
    else
        fprintf('      no second PIE pulse detected in this period\n');
    end
    usableNs = (pedestal.signalEnd - pedestal.peakBin) * dtNs;
    fprintf('      usable decay after the peak: %.2f ns\n', usableNs);
    if usableNs < 3
        warning('fit_free_dye_lifetime_above_surface:ShortSignal', ...
            ['Only %.2f ns of decay is available after the peak, so a ' ...
             'multi-ns lifetime is poorly constrained.'], usableNs);
    end

    fit = fitPooledDecay(pooled, dtNs, periodNs, opts, ptuFile, pedestal);
    fit.pedestal = pedestal;

    fprintf('\n  RESULT for %s\n', acqName);
    fprintf(['      [1] first moment %.3f ns raw -> %.3f ns corrected for ' ...
        'the %.2f ns window\n'], fit.meanLifetimeNs, fit.correctedMeanNs, ...
        fit.momentWindowNs);
    fprintf('      [2] mono-exponential tail fit      = %.3f ns   (spread %.3f ns)\n', ...
        fit.tau0Ns, fit.tailSpreadNs);
    fprintf(['      [4] log-linear regression          = %.3f ns   ' ...
        '(range spread %.3f ns: %s)\n'], fit.logLinearMedianNs, ...
        fit.logLinearSpreadNs, mat2str(round(fit.logLinearTauNs, 3)));
    for k = 1:numel(opts.tailStartNs)
        fprintf('          tail from +%.1f ns : %-7s (%d bins, %.4g photons)\n', ...
            opts.tailStartNs(k), fmt(fit.tauPerTailStart(k)), ...
            fit.tailBins(k), fit.tailPhotons(k));
    end
    if fit.distOk
        fprintf('      [3] lifetime distribution (DistTailfit), period-corrected:\n');
        for k = 1:numel(fit.distTauNs)
            fprintf(['              tau = %6.3f ns   amplitude %5.1f%%   ' ...
                'photons %5.1f%%\n'], fit.distTauNs(k), ...
                100 * fit.distAmpShare(k), 100 * fit.distPhotonShare(k));
        end
        fprintf('              intensity-weighted mean of those: %.3f ns\n', ...
            fit.distMeanNs);
    else
        fprintf('      [3] lifetime distribution unavailable: %s\n', ...
            fit.distReason);
    end
    if fit.irfOk
        fprintf(['      IRF-deconvolved cross-check: %s ns from a %.1f ns ' ...
            'seed  (IRF from %s)\n'], fmt(fit.irfTauNs(1)), ...
            fit.irfSeedNs, fit.irfSourceShort);
        fprintf(['          TREAT WITH CAUTION: the IRF is borrowed from ' ...
            'another acquisition and\n          resampled onto this ' ...
            'binning, so any timing shift between them goes\n          ' ...
            'straight into this number. Across the five scans it scatters ' ...
            'threefold.\n']);
    else
        fprintf('      IRF cross-check unavailable: %s\n', fit.irfReason);
    end

    out = struct('ptuFile', ptuFile, 'acquisition', acqName, 'plane', plane, ...
        'pixX', pixX, 'pixY', pixY, 'periodNs', periodNs, 'dtNs', dtNs, ...
        'dwellMs', dwellMs, 'gate', gate, 'axial', axial, ...
        'counts', counts, 'meanArrivalNs', meanArrivalNs, ...
        'surfaceIndex', surfaceIndex, 'direction', direction, ...
        'selection', selection, 'selectedMask', selectedMaskImage, ...
        'pooledDecay', pooled, 'fit', fit, 'opts', opts);

    matFile = fullfile(opts.outputDir, 'free_dye_above_surface.mat');
    save(matFile, 'out', '-v7.3');
    fprintf('\n  wrote %s\n', matFile);
    if opts.makeFigure
        out.flimFigure = plotMeanFlim(out);
        out.decayFigure = plotFittedDecay(out);
        fprintf('  wrote %s\n', out.flimFigure);
        fprintf('  wrote %s\n', out.decayFigure);
    end
end

% ================================================================== plumbing

function ptuFile = resolvePtu(source)
    source = char(source);
    if isfolder(source)
        candidate = fullfile(source, 'RawImage.ptu');
        if isfile(candidate); ptuFile = candidate; return; end
        found = dir(fullfile(source, '**', 'RawImage.ptu'));
        if isempty(found)
            error('fit_free_dye_lifetime_above_surface:NoPtu', ...
                'No RawImage.ptu below %s', source);
        end
        ptuFile = fullfile(found(1).folder, found(1).name);
    elseif isfile(source)
        ptuFile = source;
    else
        error('fit_free_dye_lifetime_above_surface:NoSource', ...
            'Not a file or folder: %s', source);
    end
end

function value = fieldOr(s, name, fallback)
    if isfield(s, name) && isnumeric(s.(name)) && isscalar(s.(name)) && ...
            isfinite(s.(name))
        value = double(s.(name));
    else
        value = fallback;
    end
end

function text = fmt(value)
    if isnan(value); text = '-'; else; text = sprintf('%.3f', value); end
end

function gate = pickGate(globalCurve, dtNs)
%PICKGATE Roll the period so the excitation rise sits at bin 1.
%
% The pulse in these acquisitions lands near the END of the 50 ns period, so a
% naive contiguous-region gate stops at the period boundary and discards the
% part of the decay that has wrapped round to low bin numbers. For a 2.7 ns
% lifetime with a 6.5 ns window that throws away roughly a tenth of the
% photons, all of it from the late tail - exactly the bins that constrain a
% long lifetime.
%
% Rolling the histogram instead keeps the whole decay contiguous. The window
% then runs the full period unless a SECOND pulse is present, in which case it
% stops just before it so the two PIE channels are not mixed.
    curve = double(globalCurve(:));
    nGate = numel(curve);
    width = max(3, round(nGate / 100));
    smooth = conv(curve, ones(width, 1) / width, 'same');
    [peakHeight, peakBin] = max(smooth);
    baseline = median(smooth);

    % Walk back from the peak to where the signal leaves the baseline; that is
    % the rise. Circular, because the rise can precede bin 1.
    threshold = baseline + 0.05 * (peakHeight - baseline);
    riseBin = peakBin;
    for step = 1:nGate
        candidate = mod(peakBin - step - 1, nGate) + 1;
        if smooth(candidate) <= threshold
            riseBin = candidate;
            break;
        end
    end
    % Keep a few bins of pre-rise baseline so the fit can see where zero is.
    preBins = max(2, round(0.3 / dtNs));
    riseBin = mod(riseBin - preBins - 1, nGate) + 1;

    rolled = circshift(smooth, -(riseBin - 1));

    % A second PIE pulse must be a genuine LOCAL MAXIMUM - the decay has to
    % fall to a minimum and then rise again. An earlier version simply took
    % the largest value beyond a fixed offset, which is always the left edge
    % of a monotonically decaying segment: it "found" a second pulse a few ns
    % after the first on every scan and clipped the window to 4 ns, discarding
    % most of the tail. Requiring a preceding trough is what distinguishes a
    % real pulse from the decay itself.
    firstPeak = mod(peakBin - riseBin, nGate) + 1;
    searchFrom = min(nGate - 2, firstPeak + max(5, round(1 / dtNs)));
    gateLength = nGate;
    note = 'single pulse: the window spans the whole period';
    if searchFrom < nGate - 4
        segment = rolled(searchFrom:end);
        [troughHeight, troughOffset] = min(segment);
        troughBin = searchFrom + troughOffset - 1;
        if troughBin < nGate - 2
            afterTrough = rolled(troughBin:end);
            [riseHeight, riseOffset] = max(afterTrough);
            isRealPulse = riseHeight > max(3 * troughHeight, ...
                baseline + 0.10 * (peakHeight - baseline)) && ...
                riseOffset > 2;
            if isRealPulse
                gateLength = troughBin - 1;
                note = sprintf(['second pulse at bin %d (trough at %d); the ' ...
                    'window stops at the trough'], ...
                    troughBin + riseOffset - 1, troughBin);
            end
        end
    end
    if gateLength < round(8 / dtNs)
        error('fit_free_dye_lifetime_above_surface:GateTooShort', ...
            ['The gate came out %d bins (%.2f ns), too short to fit a ' ...
             'multi-ns lifetime. Check the decay figure and the pulse ' ...
             'detection.'], gateLength, gateLength * dtNs);
    end
    gate = struct('method', 'peak-aligned circular roll', ...
        'peakBin', peakBin, 'riseBin', riseBin, 'length', gateLength, ...
        'preBins', preBins, 'note', note);
    gate.photons = sum(curve);
end


function axial = pickAxialAxis(counts, binSum, dtNs)
%PICKAXIALAXIS The axial axis is the one carrying the quenching gradient.
%
% An earlier version compared how CONCENTRATED the intensity was along each
% axis, on the argument that a flat bilayer is thin. On this data that
% discriminator failed: with a cell in the field and out-of-focus haze the two
% scores came out at 0.254 versus 0.254 on one scan, and the choice then
% flipped between acquisitions. Intensity concentration is simply not specific
% to the axial direction.
%
% The arrival time is. Quenching makes the lifetime rise with distance from the
% metal, so the axial axis carries a large monotonic spread in mean arrival
% time while the lateral axis carries almost none. Comparing the RANGE of the
% marginal mean-arrival profile is therefore a test of the actual physics
% rather than of image contrast.
    rowRange = marginalArrivalRange(counts, binSum, dtNs, 2);
    colRange = marginalArrivalRange(counts, binSum, dtNs, 1);
    axial = struct('rowRangeNs', rowRange, 'colRangeNs', colRange);
    axial.isRow = rowRange > colRange;
    if axial.isRow
        axial.name = 'im_line (image rows)';
        axial.margin = rowRange / max(colRange, eps);
    else
        axial.name = 'im_col (image columns)';
        axial.margin = colRange / max(rowRange, eps);
    end
end

function range = marginalArrivalRange(counts, binSum, dtNs, sumDim)
%MARGINALARRIVALRANGE Spread of photon-weighted mean arrival along one axis.
% Slices holding less than a thousandth of the signal are ignored, so a few
% nearly empty edge rows cannot manufacture a large range out of noise.
    sliceCounts = sum(counts, sumDim);
    sliceBinSum = sum(binSum, sumDim);
    sliceCounts = sliceCounts(:); sliceBinSum = sliceBinSum(:);
    enough = sliceCounts >= max(50, 0.001 * sum(sliceCounts));
    if nnz(enough) < 5; range = 0; return; end
    arrival = (sliceBinSum(enough) ./ sliceCounts(enough) - 1) * dtNs;
    smooth = conv(arrival, ones(5, 1) / 5, 'same');
    inner = smooth(3:end - 2);
    if isempty(inner); range = 0; return; end
    range = max(inner) - min(inner);
end

% =========================================================== surface finding

function [surfaceIndex, direction] = locateSurface(countsAx, arrivalAx, opts)
%LOCATESURFACE Per-column axial index of the low-arrival-time band.
% Per column rather than one global row: over a 600-pixel scan the coverslip is
% never perfectly level, and a single row would cut into the bilayer at one end
% while leaving a quenched margin at the other.
    [nAxial, nLateral] = size(countsAx);
    surfaceIndex = nan(1, nLateral);
    for c = 1:nLateral
        column = countsAx(:, c);
        if sum(column) < opts.minColumnCounts; continue; end
        arrival = arrivalAx(:, c);
        valid = isfinite(arrival) & column >= max(1, opts.minPixelCounts);
        if nnz(valid) < 3; continue; end
        % Weight by counts so a single dim pixel with a freak short arrival
        % time cannot define the surface.
        weighted = arrival;
        weighted(~valid) = inf;
        smooth = smoothIgnoringInf(weighted, 3);
        [~, index] = min(smooth);
        surfaceIndex(c) = index;
    end
    if all(isnan(surfaceIndex))
        error('fit_free_dye_lifetime_above_surface:NoSurface', ...
            'No lateral column had enough signal to locate the surface.');
    end
    % Fill gaps and de-noise: the surface is a physical plane, so it must be
    % smooth across the field even where individual columns are dim.
    lateral = 1:nLateral;
    good = isfinite(surfaceIndex);
    surfaceIndex = interp1(lateral(good), surfaceIndex(good), lateral, ...
        'linear', 'extrap');
    surfaceIndex = movmedian(surfaceIndex, max(3, round(0.05 * nLateral)));
    % Extrapolation at the field edges can put the line outside the image, and
    % every later step indexes with it. Clamp so a dim edge column cannot throw
    % an out-of-bounds error deep inside the selection.
    margin = max(opts.surfaceMarginPix, 3 * opts.bandWidthPix);
    guard = margin + opts.bandWidthPix + 1;
    surfaceIndex(~isfinite(surfaceIndex)) = round(nAxial / 2);
    surfaceIndex = min(max(surfaceIndex, guard), nAxial - guard);
    if nAxial < 4 * guard
        error('fit_free_dye_lifetime_above_surface:AxisTooShort', ...
            ['The axial axis is only %d pixels, too short for a %d-pixel ' ...
             'margin on both sides.'], nAxial, guard);
    end
    lowSide = 0; highSide = 0; lowPhotons = 0; highPhotons = 0;
    bandArrival = [];
    for c = 1:nLateral
        s = round(surfaceIndex(c));
        band = max(1, s - opts.bandWidthPix):min(nAxial, s + opts.bandWidthPix);
        bandArrival = [bandArrival; arrivalAx(band, c)]; %#ok<AGROW>
        low = 1:max(1, s - margin);
        high = min(nAxial, s + margin):nAxial;
        lowSide = lowSide + sum(arrivalAx(low, c) .* countsAx(low, c), 'omitnan');
        highSide = highSide + sum(arrivalAx(high, c) .* countsAx(high, c), 'omitnan');
        lowPhotons = lowPhotons + sum(countsAx(low, c));
        highPhotons = highPhotons + sum(countsAx(high, c));
    end
    lowMean = lowSide / max(lowPhotons, 1);
    highMean = highSide / max(highPhotons, 1);
    direction = struct('bandArrivalNs', mean(bandArrival, 'omitnan'), ...
        'marginPix', margin);
    % Unquenching raises the arrival time, so the side where it rises is away
    % from the metal.
    if highMean >= lowMean
        direction.sign = +1;
        direction.name = 'increasing index';
        direction.chosenArrivalNs = highMean;
        direction.otherArrivalNs = lowMean;
        direction.chosenPhotons = highPhotons;
        direction.otherPhotons = lowPhotons;
    else
        direction.sign = -1;
        direction.name = 'decreasing index';
        direction.chosenArrivalNs = lowMean;
        direction.otherArrivalNs = highMean;
        direction.chosenPhotons = lowPhotons;
        direction.otherPhotons = highPhotons;
    end
    direction.arrivalGainNs = direction.chosenArrivalNs - ...
        direction.bandArrivalNs;
end

function smooth = smoothIgnoringInf(values, width)
    finiteMask = isfinite(values);
    filled = values;
    filled(~finiteMask) = 0;
    kernel = ones(width, 1);
    numerator = conv(filled, kernel, 'same');
    denominator = conv(double(finiteMask), kernel, 'same');
    smooth = inf(size(values));
    ok = denominator > 0;
    smooth(ok) = numerator(ok) ./ denominator(ok);
end

% ============================================================= rate filtering

function selection = selectAbovePixels(countsAx, surfaceIndex, direction, opts)
%SELECTABOVEPIXELS Pixels beyond the surface band whose count clears background.
    [nAxial, nLateral] = size(countsAx);
    candidate = false(nAxial, nLateral);
    axialIndex = (1:nAxial)';
    for c = 1:nLateral
        offset = direction.sign * (axialIndex - surfaceIndex(c));
        candidate(:, c) = offset >= direction.marginPix;
    end

    % Background: the emptiest tenth of the candidate region, so it is measured
    % where there is no sample rather than assumed.
    candidateCounts = countsAx(candidate);
    if isempty(candidateCounts)
        error('fit_free_dye_lifetime_above_surface:NoCandidates', ...
            'No pixels lie beyond the surface margin on the chosen side.');
    end
    sorted = sort(candidateCounts);
    background = mean(sorted(1:max(1, round(0.10 * numel(sorted)))));
    threshold = max(opts.minPixelCounts, ...
        background + opts.backgroundSigmas * sqrt(max(background, 1)));

    mask = candidate & countsAx >= threshold;
    photons = sum(countsAx(mask));
    selection = struct('mask', mask, 'candidate', candidate, ...
        'candidateCount', nnz(candidate), 'pixelCount', nnz(mask), ...
        'photons', photons, 'backgroundCounts', background, ...
        'threshold', threshold, ...
        'backgroundFraction', background * nnz(mask) / max(photons, 1));
end

% ==================================================================== fitting

function tau = correctTruncatedMoment(observedMeanNs, windowNs)
%CORRECTTRUNCATEDMOMENT Invert <t> = tau - T*e^(-T/tau)/(1-e^(-T/tau)).
% Observing a single exponential for only T ns after the peak shortens its
% first moment; this recovers the tau that would give the observed moment.
% Monotonic in tau, so bisection is safe and needs no toolbox.
    tau = NaN;
    if ~isfinite(observedMeanNs) || observedMeanNs <= 0 || ...
            ~isfinite(windowNs) || windowNs <= 0
        return;
    end
    model = @(x) x - windowNs * exp(-windowNs / x) / ...
        max(1 - exp(-windowNs / x), eps);
    low = 0.05; high = 50;
    if model(low) > observedMeanNs || model(high) < observedMeanNs
        return;   % outside the invertible range; leave NaN rather than guess
    end
    for iteration = 1:200
        mid = 0.5 * (low + high);
        if model(mid) < observedMeanNs; low = mid; else; high = mid; end
    end
    tau = 0.5 * (low + high);
end

function pedestal = measurePedestal(pooled)
%MEASUREPEDESTAL Background level, and where the FIRST PIE window ends.
%
% The level comes from the most-decayed bins at the end of the rolled period,
% as a median so a stray bright bin cannot inflate it.
%
% The window end is the FIRST bin after the peak where the smoothed decay falls
% to the pedestal, not the last bin where it is still detectable. That
% distinction is the whole point: the detector is hardware-gated, so after the
% gate closes the histogram returns to background and a second, weak PIE pulse
% follows. Taking the LAST detectable bin ran the window past both of them and
% put the gate edge and a second decay inside the fit.
    pooled = double(pooled(:));
    nBin = numel(pooled);
    tailStart = max(1, round(0.90 * nBin));
    level = median(pooled(tailStart:end));
    smooth = conv(pooled, ones(5, 1) / 5, 'same');
    [~, peakBin] = max(smooth);

    % Fall back to the pedestal: within a couple of Poisson sigmas of it the
    % decay is over, whether because the dye has finished emitting or because
    % the gate shut.
    closed = level + 2 * sqrt(max(level, 1));
    signalEnd = nBin;
    for b = peakBin + 2:nBin
        if smooth(b) <= closed
            signalEnd = b;
            break;
        end
    end

    % Report the next pulse if there is one, so the PIE structure is visible
    % rather than inferred. A PIE partner can be a small fraction of the first
    % peak - here about 2% - so the test is against the PEDESTAL, not against
    % the first peak.
    secondPulseBin = NaN;
    secondPulseCounts = NaN;
    if signalEnd < nBin - 4
        rest = smooth(signalEnd:end);
        [height, offset] = max(rest);
        if height > level + 10 * sqrt(max(level, 1)) && offset > 2
            secondPulseBin = signalEnd + offset - 1;
            secondPulseCounts = height;
        end
    end

    pedestal = struct('level', level, 'signalEnd', signalEnd, ...
        'peakBin', peakBin, 'closedLevel', closed, ...
        'secondPulseBin', secondPulseBin, ...
        'secondPulseCounts', secondPulseCounts, ...
        'fraction', level * nBin / max(sum(pooled), 1));
end

function fit = fitPooledDecay(pooled, dtNs, periodNs, opts, ptuFile, pedestal)
%FITPOOLEDDECAY Swept-tail-start MLE tail fit, plus an IRF-deconvolved check.
    pooled = double(pooled(:));
    % Everything beyond signalEnd is pedestal; including it drags a
    % mono-exponential fit downwards and eventually fits pure noise.
    pooled = pooled(1:pedestal.signalEnd);
    gateLength = numel(pooled);
    [~, peakBin] = max(pooled);
    starts = opts.tailStartNs(:)';
    tauPerStart = nan(1, numel(starts));
    tailBins = zeros(1, numel(starts));
    tailPhotons = zeros(1, numel(starts));
    for k = 1:numel(starts)
        first = peakBin + max(1, round(starts(k) / dtNs));
        if first >= numel(pooled) - 8; continue; end
        tail = pooled(first:end);
        tailBins(k) = numel(tail);
        tailPhotons(k) = sum(tail);
        if sum(tail) < 200; continue; end
        try
            tau = Tailfit(tail, dtNs, 3.0, [0.1 15], 'mle', false);
            tauPerStart(k) = tau(1);
        catch
            tauPerStart(k) = NaN;
        end
    end
    % [1] Model-free mean lifetime: the first moment of the pedestal-subtracted
    % decay from the peak onward. Equals the lifetime exactly for a single
    % exponential and the intensity-weighted mean for a mixture, with no fit
    % involved and so nothing that can converge badly.
    fromPeak = pooled(peakBin:end) - pedestal.level;
    fromPeak = max(fromPeak, 0);
    lagNs = (0:numel(fromPeak) - 1)' * dtNs;
    signalPhotons = sum(fromPeak);
    if signalPhotons > 0
        meanLifetimeNs = sum(lagNs .* fromPeak) / signalPhotons;
    else
        meanLifetimeNs = NaN;
    end
    % Undo the truncation bias. The window covers only T ns after the peak, so
    % the observed moment is systematically short; inverting the single
    % exponential relation recovers the lifetime that would have produced it.
    windowNs = numel(fromPeak) * dtNs;
    correctedMeanNs = correctTruncatedMoment(meanLifetimeNs, windowNs);

    fit = struct('peakBin', peakBin, 'tauPerTailStart', tauPerStart, ...
        'tailBins', tailBins, 'tailPhotons', tailPhotons, ...
        'meanLifetimeNs', meanLifetimeNs, ...
        'correctedMeanNs', correctedMeanNs, ...
        'momentWindowNs', windowNs, ...
        'meanLifetimePhotons', signalPhotons, ...
        'logLinearTauNs', [], 'logLinearRanges', [], ...
        'logLinearMedianNs', NaN, 'logLinearSpreadNs', NaN, ...
        'distOk', false, 'distTauNs', [], 'distAmpShare', [], ...
        'distPhotonShare', [], 'distMeanNs', NaN, 'distReason', '', ...
        'tau0Ns', median(tauPerStart, 'omitnan'), ...
        'tailSpreadNs', max(tauPerStart) - min(tauPerStart), ...
        'irfOk', false, 'irfTauNs', NaN, 'irfAmp', NaN, ...
        'irfMoved', false, 'irfSeedNs', NaN, 'gateLength', gateLength, ...
        'irfReason', '', 'irfSource', '', 'irfSourceShort', '', ...
        'irf', [], 'model', [], 'modelFirstBin', NaN);

    % [4] Log-linear regression over several ranges, weighted by counts so
    % bright early bins are not swamped by noisy late ones. Sharing no code
    % with Tailfit or DistTailfit makes this a genuine cross-check, and the
    % spread across ranges measures departure from a single exponential.
    rangeStarts = round([0.4 0.6 1.0 1.5] / dtNs);
    logTau = nan(1, numel(rangeStarts));
    for k = 1:numel(rangeStarts)
        i0 = peakBin + max(1, rangeStarts(k));
        if i0 >= numel(pooled) - 8; continue; end
        y = pooled(i0:end) - pedestal.level;
        t = (i0:numel(pooled))' * dtNs;
        usable = y > max(50, 5 * sqrt(max(pedestal.level, 1)));
        if nnz(usable) < 8; continue; end
        yy = y(usable); tt = t(usable);
        weights = yy;                       % Poisson: var(log y) ~ 1/y
        sw = sum(weights);
        meanT = sum(weights .* tt) / sw;
        meanL = sum(weights .* log(yy)) / sw;
        slope = sum(weights .* (tt - meanT) .* (log(yy) - meanL)) / ...
            max(sum(weights .* (tt - meanT) .^ 2), eps);
        if slope < 0; logTau(k) = -1 / slope; end
    end
    fit.logLinearTauNs = logTau;
    fit.logLinearRanges = rangeStarts * dtNs;
    fit.logLinearMedianNs = median(logTau, 'omitnan');
    if any(isfinite(logTau))
        fit.logLinearSpreadNs = max(logTau) - min(logTau);
    end

    % [3] Lifetime distribution. The number of components is inferred rather
    % than assumed, and the pulse period is passed so incomplete decay from the
    % previous pulse is modelled instead of being absorbed into an offset or
    % into a spuriously long component.
    first = peakBin + max(1, round(min(starts) / dtNs));
    if exist('DistTailfit', 'file') ~= 2
        fit.distReason = 'DistTailfit.m is not on the path';
    elseif first >= numel(pooled) - 12 || sum(pooled(first:end)) < 1000
        fit.distReason = 'too few tail photons';
    else
        try
            % N rates must not exceed what a short tail can support. With the
            % window now ending at the gate there are only ~60 tail bins, and
            % asking for 200 rates made the solve fail outright on one scan.
            tail = pooled(first:end);
            nRate = max(10, min(120, floor(numel(tail) / 2)));
            [amplitude, rate] = DistTailfit(tail, dtNs, 1, false, nRate, ...
                min(20, 4 * numel(tail) * dtNs), periodNs);
            amplitude = double(amplitude(:)');
            rate = double(rate(:)');
            keep = isfinite(rate) & rate > 0 & isfinite(amplitude) & ...
                amplitude > 0;
            amplitude = amplitude(keep); rate = rate(keep);
            if isempty(amplitude)
                fit.distReason = 'the distribution fit returned no components';
            else
                tau = 1 ./ rate;
                [tau, order] = sort(tau);
                amplitude = amplitude(order);
                % Amplitude share is what fraction of molecules; photon share
                % weights by lifetime, which is what a brightness-weighted
                % reading of the decay actually reflects.
                photonWeight = amplitude .* tau;
                fit.distTauNs = tau;
                fit.distAmpShare = amplitude / sum(amplitude);
                fit.distPhotonShare = photonWeight / sum(photonWeight);
                fit.distMeanNs = sum(photonWeight .* tau) / sum(photonWeight);
                fit.distOk = true;
            end
        catch distError
            fit.distReason = distError.message;
        end
    end

    if strcmpi(char(opts.irfSource), 'none')
        fit.irfReason = 'disabled'; return;
    end
    if exist('Fluofit', 'file') ~= 2
        fit.irfReason = 'Fluofit.m not on the path'; return;
    end
    matFile = char(opts.irfSource);
    if isempty(matFile); matFile = findSessionIrf(ptuFile); end
    if isempty(matFile) || ~isfile(matFile)
        fit.irfReason = 'no XY analysis MAT with a stored IRF in this session';
        return;
    end
    loaded = load(matFile, 'result');
    if ~isfield(loaded, 'result') || ~isfield(loaded.result, 'irf') || ...
            ~isfield(loaded.result.irf, 'curve')
        fit.irfReason = 'the analysis MAT has no stored IRF'; return;
    end
    irf = double(loaded.result.irf.curve(:));
    irf = resampleTo(irf, gateLength);
    irf = max(irf, 0);
    if sum(irf) <= 0; fit.irfReason = 'the stored IRF is empty'; return; end
    irf = irf / sum(irf);
    fit.irf = irf;
    fit.irfSource = matFile;
    [~, shortName] = fileparts(fileparts(fileparts(matFile)));
    fit.irfSourceShort = shortName;
    % Seed at a FIXED value, not at the tail-fit answer. Seeding with the
    % tail-fit result returned exactly that result on all five scans, which
    % means the optimiser had not moved and the check was echoing its input.
    seedNs = 2.5;
    fit.irfSeedNs = seedNs;
    try
        [~, ~, amp, tau] = Fluofit(irf, pooled, periodNs, dtNs, ...
            seedNs, [], 0, 'mle', false);
        fit.irfTauNs = double(tau(:)');
        fit.irfAmp = double(amp(:)');
        fit.irfMoved = abs(fit.irfTauNs(1) - seedNs) > 1e-4;
        fit.irfOk = fit.irfMoved;
        if ~fit.irfMoved
            fit.irfReason = sprintf(['Fluofit returned its %.3f ns seed ' ...
                'unchanged, so the check is inconclusive'], seedNs);
        end
    catch fitError
        fit.irfReason = fitError.message;
    end
end

function file = findSessionIrf(ptuFile)
    sessionRoot = fileparts(fileparts(ptuFile));
    found = dir(fullfile(sessionRoot, '**', ...
        'immune_cell_MIET_640nm_red_analysis.mat'));
    file = '';
    if isempty(found); return; end
    [~, order] = sort([found.datenum], 'descend');
    file = fullfile(found(order(1)).folder, found(order(1)).name);
end

function resampled = resampleTo(curve, targetLength)
    curve = double(curve(:));
    if numel(curve) == targetLength; resampled = curve; return; end
    resampled = interp1(linspace(0, 1, numel(curve))', curve, ...
        linspace(0, 1, targetLength)', 'linear', 0);
end

% ==================================================================== figures

function name = plotMeanFlim(out)
    counts = out.counts;
    arrival = out.meanArrivalNs;
    mask = out.selectedMask;
    if ~out.axial.isRow
        counts = counts'; arrival = arrival'; mask = mask';
    end

    h = figure('Color', 'w', 'Visible', 'off', 'Position', [40 40 1500 460]);
    layout = tiledlayout(h, 1, 3, 'Padding', 'compact', ...
        'TileSpacing', 'compact');

    ax = nexttile(layout);
    imagesc(ax, log10(double(counts) + 1));
    axis(ax, 'image'); colormap(ax, gray); colorbar(ax);
    hold(ax, 'on');
    plot(ax, 1:numel(out.surfaceIndex), out.surfaceIndex, 'r-', 'LineWidth', 1);
    xlabel(ax, 'lateral (pixels)'); ylabel(ax, 'axial (pixels)');
    title(ax, sprintf('%s intensity, log_{10} photons, surface in red', ...
        out.plane));

    ax = nexttile(layout);
    handle = imagesc(ax, arrival);
    set(handle, 'AlphaData', ~isnan(arrival));
    axis(ax, 'image'); colormap(ax, parula); colorbar(ax);
    % Setting c.Label.String triggers a colorbar MarkedDirty listener that
    % errors on an invisible figure under software OpenGL, so the quantity is
    % named in the axes title instead.
    finiteValues = arrival(isfinite(arrival));
    if ~isempty(finiteValues)
        % Local quantile rather than prctile: the Statistics Toolbox is not
        % assumed anywhere else in this folder and should not be assumed here.
        caxis(ax, [quantileLocal(finiteValues, 0.02) ...
            quantileLocal(finiteValues, 0.98)]);
    end
    hold(ax, 'on');
    plot(ax, 1:numel(out.surfaceIndex), out.surfaceIndex, 'r-', 'LineWidth', 1);
    xlabel(ax, 'lateral (pixels)');
    title(ax, 'mean FLIM: arrival time in ns, not IRF-corrected');

    ax = nexttile(layout);
    overlay = zeros([size(counts) 3]);
    grey = double(counts);
    if max(grey(:)) > 0; grey = grey / max(grey(:)); end
    grey = grey .^ 0.4;
    overlay(:, :, 1) = grey; overlay(:, :, 2) = grey; overlay(:, :, 3) = grey;
    green = overlay(:, :, 2); green(mask) = 1; overlay(:, :, 2) = green;
    red = overlay(:, :, 1); red(mask) = 0.15; overlay(:, :, 1) = red;
    blue = overlay(:, :, 3); blue(mask) = 0.15; overlay(:, :, 3) = blue;
    image(ax, overlay); axis(ax, 'image');
    xlabel(ax, 'lateral (pixels)');
    title(ax, sprintf('%d selected pixels, %.3g photons', ...
        out.selection.pixelCount, out.selection.photons));

    name = fullfile(out.opts.outputDir, 'mean_flim.png');
    exportgraphics(h, name, 'Resolution', 200);
    close(h);
end

function value = quantileLocal(values, level)
    values = sort(double(values(:)));
    if isempty(values); value = NaN; return; end
    if numel(values) == 1; value = values; return; end
    position = min(max(level * numel(values) + 0.5, 1), numel(values));
    low = floor(position); high = ceil(position);
    if low == high
        value = values(low);
    else
        value = values(low) + (position - low) * (values(high) - values(low));
    end
end

function name = plotFittedDecay(out)
    decay = double(out.pooledDecay);
    timeNs = (0:numel(decay) - 1)' * out.dtNs;
    fit = out.fit;

    h = figure('Color', 'w', 'Visible', 'off', 'Position', [40 40 900 620]);
    layout = tiledlayout(h, 4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax = nexttile(layout, 1, [3 1]);
    semilogy(ax, timeNs, max(decay, 0.5), 'k.', 'MarkerSize', 5, ...
        'DisplayName', sprintf('pooled data (%.3g photons)', sum(decay)));
    hold(ax, 'on');
    % The pedestal and the fit boundary are drawn, because the fit deliberately
    % stops where the decay stops beating the background and a figure spanning
    % the whole period would imply otherwise.
    modelCurve = [];
    fitLast = numel(decay);
    if isfield(fit, 'pedestal') && ~isempty(fit.pedestal)
        fitLast = min(numel(decay), fit.pedestal.signalEnd);
        yline(ax, max(fit.pedestal.level, 0.5), '--', 'Color', [0.4 0.4 0.4], ...
            'DisplayName', sprintf('pedestal %.0f counts/bin', ...
            fit.pedestal.level));
        if fitLast < numel(decay)
            xline(ax, timeNs(fitLast), '-.', 'Color', [0.2 0.5 0.2], ...
                'DisplayName', 'fit window ends');
        end
    end
    if isfinite(fit.tau0Ns)
        first = fit.peakBin + max(1, round(min(out.opts.tailStartNs) / out.dtNs));
        first = min(first, fitLast - 2);
        span = first:fitLast;
        tailTime = timeNs(span) - timeNs(first);
        headBins = first:min(first + 3, fitLast);
        scale = sum(decay(headBins)) / max(numel(headBins), 1);
        modelCurve = scale * exp(-tailTime / fit.tau0Ns) + ...
            max(fit.pedestal.level, 0);
        semilogy(ax, timeNs(span), max(modelCurve, 0.5), 'r-', ...
            'LineWidth', 1.6, 'DisplayName', ...
            sprintf('tail fit \\tau = %.3f ns + pedestal', fit.tau0Ns));
        xline(ax, timeNs(first), ':k', 'tail start');
    end
    grid(ax, 'on');
    ylabel(ax, 'photons per bin');
    legend(ax, 'Location', 'northeast');
    titleText = sprintf('%s (%s): pooled decay above the surface', ...
        out.acquisition, out.plane);
    if fit.irfOk
        titleText = sprintf('%s   |   IRF check %.3f ns', titleText, ...
            fit.irfTauNs(1));
    end
    title(ax, titleText);

    ax = nexttile(layout, 4, [1 1]);
    if ~isempty(modelCurve)
        first = fit.peakBin + max(1, round(min(out.opts.tailStartNs) / out.dtNs));
        first = min(first, fitLast - 2);
        span = first:fitLast;
        observed = decay(span);
        expected = max(modelCurve, 1e-9);
        residual = (observed - expected) ./ sqrt(max(expected, 1));
        plot(ax, timeNs(span), residual, 'b-', 'LineWidth', 1);
        yline(ax, 0, 'k-');
        ylabel(ax, 'Poisson residual');
    end
    xlabel(ax, 'time (ns)');
    grid(ax, 'on');

    name = fullfile(out.opts.outputDir, 'fitted_tcspc.png');
    exportgraphics(h, name, 'Resolution', 200);
    close(h);
end
