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
% THE BIN WIDTH COMES FROM THE READER, NEVER FROM period/nGate
%
% PTU_FLIM_GPU rebins the native 20 ps channels by an integer divisor and then
% CAPS the bin count: chDiv = round(minLifetimeBin_ns / nativeRes) = 3, so the
% stored bin width is 0.060 ns, and Ngate = min(maxNgate, ceil(period/0.060)+1)
% = min(512, 835) = 512. The 512 stored bins therefore tile only 30.72 ns of the
% 50 ns period, and PTU_FLIM_GPU discards every photon after 30.72 ns
% (photonValid = ... & tmptcspc < Ngate*chDiv).
%
% Deriving dtNs as periodNs/nGate assumed the bins tile the period. They do not,
% and that made every nanosecond in this program 1.6276 times too large. The
% reader exports the truth as ptu.Resolution_ns and also overwrites
% head.MeasDesc_Resolution with it; this code now reads that and refuses to run
% if the stored array does not cover the period.
%
% The error was self-concealing in two ways worth recording. A range check on
% im_tcspc can never fire, because the reader already clamps indices into
% [1, Ngate] - so it gave false reassurance. And the wrong width made the two
% PIE peaks appear 41 ns apart instead of 25 ns, which is what led to
% "a hardware gate at 7.2 ns" and "a second 640 pulse at 2% of the peak". Both
% were artefacts: the first is the end of the stored array, the second is the
% 485 nm peak reached by wrapping past it.
%
% The check that settles the bin width is the PIE geometry itself. The two peaks
% sit 421 stored bins apart; at 0.060 ns that is 25.3 ns, exactly half the
% period and exactly the header's PIETimeGate of 1250 native bins. No other bin
% width makes the interleaving physical.
%
% WHICH PIE WINDOW - TAKEN FROM THE HEADER, NOT GUESSED
%
% The file records LaserCount = 2, LaserWL = [485 640], PIENumPIEWindows = 2 and
% PIETimeGate = [0 1250] native bins. At 20 ps native resolution that is windows
% of 0-25 ns for the 485 nm laser and 25-50 ns for the 640 nm laser. The window
% is therefore selected by matching the requested excitation wavelength against
% LaserWL, and the choice is checked against where the photons actually are.
%
% An earlier version detected peaks in the decay instead. It happened to land on
% the right window, but it also mistook the 485 window's pulse for a "second
% 640 pulse" and mistook the 25 ns window boundary for a hardware gate.
%
% One consequence is unavoidable and must be stated: the 640 pulse fires at
% 44.0 ns, only 6.0 ns before its window closes at 50 ns. At a ~2.7 ns lifetime
% that is 2.2 lifetimes of usable decay, so the tail is short whatever the fit
% does. Photons arriving after 50 ns fall in the next period's 485 window and
% are not recoverable for the 640 decay.
%
% HOW THE SURFACE IS FOUND - ONE FLAT ROW
%
% The supported bilayer is FLAT and it is the brightest, most quenched feature
% in the frame, so the surface is a single axial row: the minimum of the
% marginal photon-weighted mean-arrival profile.
%
% A per-column minimum was tried first and it failed visibly. Where a cell sits
% on the bilayer the per-column minimum climbs the cell, so the "surface" line
% jumped ~130 pixels upward across the cell footprint and dragged the selection
% with it. A flat plane must be fitted as a flat plane.
%
% Intensity does not CHOOSE the row - the arrival-time minimum does, because a
% cell membrane is bright too. But intensity does restrict WHERE to look, and
% that turned out to be essential. A cross-section is half sample and half
% opaque substrate, and in the dark half the few hundred photons per row give
% arbitrary arrival times; the unrestricted minimum landed there on four of five
% scans, putting the "surface" 100-150 rows into empty space. The search is
% therefore confined to rows carrying a stated fraction of the brightest row's
% signal, which excludes the empty half without deciding anything about which
% bright row is the bilayer.
%
% WHICH SIDE IS "ABOVE", AND HOW FAR UP TO GO
%
% Primarily the side holding MORE PHOTONS. Above the metal is sample; below it
% is substrate, so the asymmetry is large and unambiguous. The arrival-time
% gradient is kept as a cross-check and a disagreement between the two is
% warned about rather than silently resolved - relying on the gradient alone put
% the sample below the surface on two scans once the dark half was in play. But taking
% everything on that side is wrong: just above the bilayer the dye is still
% partly quenched, and pooling that with free dye higher up gives a decay that
% is a mixture of lifetimes, which showed as strongly curved residuals reaching
% 20 sigma. The selection therefore starts only where the arrival time has
% PLATEAUED - the top part of the section - and that distance is measured from
% the arrival profile rather than assumed.
%
% HOW FAR ABOVE THE PLATEAU, AND WHY 110 PIXELS
%
% The arrival-time plateau alone puts the pool about 70 px above the bilayer, and
% there the pooled decay still departs from a single exponential. Sweeping the
% cut with test_free_dye_lifetime_selection_robustness shows the departure is
% largely SPATIAL MIXING - partly-quenched dye nearer the bilayer - and that it
% can be cut away:
%
%   cut(px)   photons   tau(ns)  reducedDev  maxResid  logLinSpread(ns)
%       0     3.24e6     1.908      8.68      22.1        0.052
%      70     2.71e6     1.934      3.40       8.5        0.031
%     180     1.86e6     1.939      2.59       5.8        0.028
%     300     2.57e5     1.935      1.24       4.2        0.031
%
% tau is flat at 1.934-1.939 for every cut past 70 px, so the lifetime does not
% depend on this choice; only the goodness of fit does.
%
% The gain is real and not merely lost statistical power. A fixed RELATIVE model
% error makes the reduced-deviance excess scale with photon count, and the
% observed excess is only 35-40% of that prediction at every cut, so the relative
% misfit itself falls by a factor of about 2.5 to 3. The log-linear range spread
% - a lifetime drift, hence nearly photon-independent - halves, 0.052 to
% 0.028 ns.
%
% But it SATURATES: the genuine fraction is 37% at 70 px, 36% at 180 and 40% at
% 300, i.e. flat, and the log-linear spread reaches its floor of 0.028 ns by
% 180 px. Cut 300's reduced deviance of 1.24 looks best only because it holds
% 12.6 times fewer photons.
%
% extraTopMarginPix nevertheless DEFAULTS TO 0, i.e. the plateau alone. Running
% with 110 px was tried on all five scans and does what the sweep promises -
% reduced deviance 3.0-4.0 down to 1.9-3.1, max residual 8.5-9.2 down to
% 6.0-7.0 sigma, and a fitted background within about 2 counts/bin of the
% measured pedestal instead of roughly twice it - but it shrinks the pool from
% 20-36% of the frame to 11-18%, discarding 30-45% of the photons. That was
% judged too expensive for the gain.
%
% The decisive point is that it does not matter for the answer: tau0 is
% 1.89-1.99 ns at the plateau cut and 1.90-2.01 ns with the extra margin, and
% the per-scan lifetime moves by less than 0.03 ns. Only the goodness of fit
% changes. Set extraTopMarginPix = 110 when the fit quality matters more than
% the photon budget.
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
%   2. MONO-EXPONENTIAL PLUS FITTED BACKGROUND, by variable projection. The
%      model is
%          y(t) = A * exp(-t/tau) + B
%      with A and B BOTH non-negative and solved together by PIRLSnonneg -
%      Poisson iteratively reweighted least squares - inside a one-dimensional
%      search over tau that minimises the Poisson deviance. This is the primary
%      estimator.
%
%      An earlier version instead fitted tau with Tailfit and then, for display,
%      scaled an exponential to the first few bins and added the MEDIAN pedestal
%      on top. That is not a fit: the amplitude and the background were never
%      solved against the data, and the residuals showed it, running +12 sigma
%      at the start of the tail, down to -4 in the middle and back to +3 at the
%      end. Fitting the background as a free non-negative parameter removes
%      that structure, and because the offset now absorbs the flat part the fit
%      can use every bin out to the end of the PIE window instead of stopping
%      where signal meets background.
%
%      The swept Tailfit is retained as a diagnostic only, since its drift
%      across tail starts still reports departures from a single exponential.
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
%   excitationNm       excitation wavelength whose PIE window to use
%                      (default 640, matched against the header's LaserWL)
%   minTopDistancePix  how far above the surface to start; [] (default) finds
%                      the arrival-time plateau from the data
%   extraTopMarginPix  pixels added BEYOND the arrival-time plateau before the
%                      pool starts. DEFAULT 0: the plateau alone. Setting it to
%                      110 measurably improves the fit but costs 30-45% of the
%                      photons, which was judged too expensive - see the note
%                      on spatial mixing below, which records both sides so the
%                      trade is one number away.
%   surfaceBrightFraction  only rows holding at least this fraction of the
%                      brightest row's photons are considered for the bilayer
%                      (default 0.25). This excludes the dark half of a
%                      cross-section, where arrival times are pure noise.
%   maxNgate           cap on stored TCSPC bins (default 1024). It MUST be
%                      large enough for the stored bins to span the whole laser
%                      period, or the reader silently discards the late part of
%                      the period; at 0.060 ns bins and a 50 ns period that
%                      needs 835, and the previous default of 512 did not.
%   minPixelCounts     absolute floor on counts per pixel (default 5)
%   backgroundSigmas   Poisson sigmas above background for the rate filter
%                      (default 5)
%   tailStartNs        swept tail starts (default [0.6 1.0 1.5 2.0 3.0])
%   irfSource          '' to search the session, a MAT path, or 'none'
%   outputDir, makeFigure

    if nargin < 2 || isempty(opts); opts = struct(); end
    defaults = struct( ...
        'excitationNm', 640, 'maxNgate', 1024, 'tcspcBinNs', 0.05, ...
        'minTopDistancePix', [], 'plateauFraction', 0.90, ...
        'extraTopMarginPix', 0, ...
        'surfaceBrightFraction', 0.25, 'minSurfaceRowCounts', 500, ...
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
    dtNs = readerBinWidthNs(ptu);
    storedSpanNs = nGate * dtNs;
    fprintf('  %d photons, %d TCSPC bins of %.4f ns = %.2f ns stored\n', ...
        numel(ptu.im_tcspc), nGate, dtNs, storedSpanNs);
    if storedSpanNs < 0.99 * periodNs
        error('fit_free_dye_lifetime_above_surface:TruncatedPeriod', ...
            ['The reader stored only %.2f ns of the %.2f ns period ' ...
             '(Ngate = %d, capped by maxNgate = %d; dt = %.4f ns). Photons ' ...
             'from %.2f to %.2f ns were DISCARDED, the histogram cannot be ' ...
             'wrapped, and the two PIE windows cannot be separated. Raise ' ...
             'maxNgate to at least %d.'], storedSpanNs, periodNs, nGate, ...
            opts.maxNgate, dtNs, storedSpanNs, periodNs, ...
            ceil(periodNs / dtNs) + 1);
    end

    % ---- PIE gate ----------------------------------------------------------
    % im_tcspc is ALREADY a 1-based bin index in this reader - see
    % immune_cell_MIET_reassigned_sliding_tcspc, which tests
    % "timeBin >= 1 & timeBin <= gateLength" and offsets with (timeBin - 1).
    % Adding one here pushed the maximum to Ngate+1 and accumarray refused it.
    rawBin = double(ptu.im_tcspc(:));
    if isfield(ptu, 'im_chan')
        detectorIds = unique(double(ptu.im_chan(:)));
        fprintf(['  detectors: %d channel(s) present (ids %d-%d), ALL summed ' ...
            '- no channel filter is applied\n'], numel(detectorIds), ...
            min(detectorIds), max(detectorIds));
    end
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
    gate = pieWindowFromHeader(head, globalCurve, nGate, periodNs, ...
        opts.excitationNm);
    gateLength = gate.length;
    fprintf(['  PIE: %d window(s) in the header at %s ns for lasers %s nm; ' ...
        'using the %g nm window\n'], gate.windowCount, ...
        mat2str(round(gate.windowStartNs, 2)), mat2str(gate.laserNm), ...
        gate.selectedNm);
    fprintf(['      window spans %.2f-%.2f ns; pulse peak at %.2f ns, so ' ...
        '%.2f ns of decay remains\n      before the window closes\n'], ...
        gate.windowStartNs(gate.index), gate.windowStopNs, ...
        (gate.peakBin - 1) * dtNs, gate.usableNs);
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

    [surfaceIndex, direction] = locateFlatSurface(countsAx, arrivalAx, opts);
    fprintf(['  surface: FLAT at axial row %d of %d (the bilayer), chosen ' ...
        'among %d rows holding\n      at least %.0f%% of the brightest ' ...
        'row; "away from metal" is the %s direction\n'], surfaceIndex, ...
        size(countsAx, 1), direction.candidateRows, ...
        100 * opts.surfaceBrightFraction, direction.name);
    fprintf(['      photons %.3g on the chosen side vs %.3g on the other; ' ...
        'arrival %.3f -> %.3f ns\n'], direction.chosenPhotons, ...
        direction.otherPhotons, direction.bandArrivalNs, ...
        direction.chosenArrivalNs);
    if direction.disagree
        warning('fit_free_dye_lifetime_above_surface:DirectionDisagreement', ...
            ['Photon count and arrival-time gradient disagree about which ' ...
             'side is away from the metal. Photon count was used. Inspect ' ...
             'mean_flim.png.']);
    end
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
    profile = arrivalVersusDistance(countsAx, arrivalAx, surfaceIndex, ...
        direction);
    top = choosePlateauDistance(profile, opts);
    top.plateauOnlyPix = top.distancePix;
    if isempty(opts.minTopDistancePix)
        top.distancePix = top.distancePix + opts.extraTopMarginPix;
    end
    fprintf(['\n  arrival time versus distance above the bilayer (photon-' ...
        'weighted):\n']);
    for d = [0 20 40 60 80 100 150 200 250 300]
        if d <= max(profile.distance)
            fprintf('        +%3d px : %.3f ns  (%.4g photons)\n', d, ...
                interp1(profile.distance, profile.arrivalNs, d), ...
                interp1(profile.distance, profile.photons, d));
        end
    end
    fprintf(['      rises from %.3f ns at the bilayer to %.3f ns; %.0f%% of ' ...
        'that rise is reached\n      by +%d px\n'], profile.arrivalNs(1), ...
        top.plateauArrivalNs, 100 * opts.plateauFraction, ...
        top.plateauOnlyPix);
    fprintf(['      pool starts at +%d px = plateau + %d px extra margin, ' ...
        'which is where the\n      departure from a single exponential ' ...
        'stops improving (see the help)\n'], top.distancePix, ...
        top.distancePix - top.plateauOnlyPix);
    if ~isempty(opts.minTopDistancePix)
        fprintf('      (overridden by opts.minTopDistancePix = %d)\n', ...
            top.distancePix);
    end

    selection = selectTopPixels(countsAx, surfaceIndex, direction, ...
        top.distancePix, opts);
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
    fprintf(['      [2] A*exp(-t/tau)+B, PIRLSnonneg   = %.3f ns   ' ...
        '<- PRIMARY\n'], fit.pirls.tauNs);
    fprintf(['              amplitude %.4g, fitted background %.2f ' ...
        'counts/bin (pedestal said %.0f)\n'], fit.pirls.amplitude, ...
        fit.pirls.offset, fit.pedestal.level);
    backgroundBins = fit.pirls.lastBin - fit.pedestal.signalEnd;
    fprintf(['              %d bins fitted to %.2f ns, of which %d are past ' ...
        'the decay and pin B\n'], fit.pirls.nBin, ...
        fit.pirls.lastBin * dtNs, max(backgroundBins, 0));
    fprintf(['              reduced deviance %.3f, max |residual| %.1f ' ...
        'sigma\n'], fit.pirls.reducedDeviance, ...
        max(abs(fit.pirls.residual)));
    fprintf('      [2b] Tailfit sweep (diagnostic)    = %.3f ns   (spread %.3f ns)\n', ...
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
        fprintf(['              photon-weighted mean %.3f ns is THE ' ...
            'estimate; the components\n              above are rate-grid ' ...
            'points, not separate populations\n'], fit.distMeanNs);
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
        'arrivalProfile', profile, 'topChoice', top, ...
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

function dtNs = readerBinWidthNs(ptu)
%READERBINWIDTHNS The width of one STORED TCSPC bin, from the reader itself.
%
% Never period/nGate: the reader caps the bin count, so the stored bins need not
% tile the period, and assuming they do scaled every lifetime in this program by
% 1.6276. PTU_FLIM_GPU exports the true value as Resolution_ns and also
% overwrites head.MeasDesc_Resolution with the coarse width, so either is
% authoritative and both are preferred to any arithmetic on the period.
    dtNs = [];
    if isfield(ptu, 'Resolution_ns') && isscalar(ptu.Resolution_ns) && ...
            isfinite(ptu.Resolution_ns) && ptu.Resolution_ns > 0
        dtNs = double(ptu.Resolution_ns);
        return;
    end
    if isfield(ptu, 'head') && isfield(ptu.head, 'MeasDesc_Resolution')
        candidate = 1e9 * double(ptu.head.MeasDesc_Resolution);
        if isfinite(candidate) && candidate > 0
            dtNs = candidate;
            return;
        end
    end
    if isempty(dtNs)
        error('fit_free_dye_lifetime_above_surface:NoBinWidth', ...
            ['The reader returned neither Resolution_ns nor a usable ' ...
             'head.MeasDesc_Resolution, so the stored bin width is unknown. ' ...
             'It must not be inferred from the period.']);
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

function gate = pieWindowFromHeader(head, globalCurve, nGate, periodNs, wantedNm)
%PIEWINDOWFROMHEADER Select the PIE window of the requested laser.
%
% PIETimeGate gives each window's start in native TCSPC bins and LaserWL gives
% the wavelength per laser, index for index. The window is chosen by matching
% the requested wavelength, then the choice is CHECKED against where the
% photons actually are - if the global peak is not inside the chosen window,
% something about the index correspondence is wrong and that is worth a warning
% rather than a silent wrong answer.
    curve = double(globalCurve(:));
    dtNs = periodNs / nGate;
    if ~isfield(head, 'PIETimeGate') || isempty(head.PIETimeGate)
        error('fit_free_dye_lifetime_above_surface:NoPieHeader', ...
            ['The PTU header has no PIETimeGate, so the excitation window ' ...
             'cannot be identified. Detecting peaks instead was tried and ' ...
             'misread the window boundary as a hardware gate.']);
    end
    gatesNative = double(head.PIETimeGate(:));
    windowCount = numel(gatesNative);
    if isfield(head, 'PIENumPIEWindows') && ~isempty(head.PIENumPIEWindows)
        windowCount = min(windowCount, round(double(head.PIENumPIEWindows(1))));
    end
    gatesNative = gatesNative(1:windowCount);
    nativeDtNs = 1e9 * double(head.MeasDesc_Resolution);

    windowStartNs = gatesNative(:)' * nativeDtNs;
    startBin = floor(windowStartNs / dtNs) + 1;
    stopBin = [startBin(2:end) - 1, nGate];

    laserNm = [];
    if isfield(head, 'LaserWL') && ~isempty(head.LaserWL)
        laserNm = double(head.LaserWL(:))';
    end
    if numel(laserNm) >= windowCount
        [~, index] = min(abs(laserNm(1:windowCount) - wantedNm));
    else
        index = windowCount;      % last window by convention
        laserNm = nan(1, windowCount);
    end

    [~, peakBin] = max(curve);
    if peakBin < startBin(index) || peakBin > stopBin(index)
        warning('fit_free_dye_lifetime_above_surface:PieWindowMismatch', ...
            ['The brightest bin (%d) lies outside the %g nm window ' ...
             '(bins %d-%d). The laser-to-window index correspondence may ' ...
             'not hold for this file.'], peakBin, laserNm(index), ...
            startBin(index), stopBin(index));
    end

    % Roll so the rise of THIS window's pulse sits at bin 1, and stop at the
    % window's own end - photons past it belong to the other laser.
    inWindow = false(nGate, 1);
    inWindow(startBin(index):stopBin(index)) = true;
    masked = curve; masked(~inWindow) = 0;
    smooth = conv(masked, ones(5, 1) / 5, 'same');
    [peakHeight, localPeak] = max(smooth);
    baseline = median(curve(inWindow));
    threshold = baseline + 0.05 * (peakHeight - baseline);
    riseBin = localPeak;
    while riseBin > startBin(index) && smooth(riseBin) > threshold
        riseBin = riseBin - 1;
    end
    preBins = max(2, round(0.3 / dtNs));
    riseBin = max(startBin(index), riseBin - preBins);

    gate = struct('method', 'PIE window from the PTU header', ...
        'windowCount', windowCount, 'index', index, ...
        'laserNm', laserNm(1:windowCount), 'selectedNm', laserNm(index), ...
        'windowStartNs', windowStartNs, ...
        'windowStopNs', stopBin(index) * dtNs, ...
        'startBin', startBin(index), 'stopBin', stopBin(index), ...
        'riseBin', riseBin, 'peakBin', localPeak - riseBin + 1, ...
        'length', stopBin(index) - riseBin + 1, ...
        'usableNs', (stopBin(index) - localPeak) * dtNs, ...
        'photons', sum(curve(startBin(index):stopBin(index))));
    gate.note = sprintf(['%.2f ns of decay after the peak = %.1f lifetimes ' ...
        'at 2.7 ns'], gate.usableNs, gate.usableNs / 2.7);
end

function gate = pickGate(globalCurve, dtNs) %#ok<DEFNU>
% Retained for reference only; pieWindowFromHeader supersedes it.
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

function [surfaceIndex, direction] = locateFlatSurface(countsAx, arrivalAx, opts)
%LOCATEFLATSURFACE One axial row for the supported bilayer, plus the "up" side.
%
% The bilayer is a flat plane, so it is ONE row: the minimum of the marginal
% photon-weighted mean-arrival profile, searched only among rows that carry real
% signal.
%
% Two failures shaped this. A per-column minimum climbed the cell where a cell
% sat on the bilayer, moving the "surface" ~130 rows and dragging the selection
% with it - a flat plane must be fitted as flat. Then an unrestricted flat
% search put the row in the DARK half of the section, where a few hundred
% photons per row give meaningless arrival times. Intensity does not pick the
% row, but it must bound the search.
    [nAxial, ~] = size(countsAx);
    sliceCounts = sum(countsAx, 2);
    sliceWeighted = sum(arrivalAx .* countsAx, 2, 'omitnan');

    bright = sliceCounts >= opts.surfaceBrightFraction * max(sliceCounts) & ...
        sliceCounts >= opts.minSurfaceRowCounts;
    if nnz(bright) < 3
        error('fit_free_dye_lifetime_above_surface:NoBrightRows', ...
            ['Only %d axial row(s) hold %.0f%% of the brightest row and at ' ...
             'least %d photons, which is too few to locate the bilayer.'], ...
            nnz(bright), 100 * opts.surfaceBrightFraction, ...
            opts.minSurfaceRowCounts);
    end
    profile = inf(nAxial, 1);
    profile(bright) = sliceWeighted(bright) ./ sliceCounts(bright);
    smooth = smoothIgnoringInf(profile, 5);
    smooth(~bright) = inf;
    [~, surfaceIndex] = min(smooth);

    % Which side is away from the metal. Photon count decides: above the metal
    % is sample, below is substrate. The arrival gradient only cross-checks.
    lowRows = 1:max(1, surfaceIndex - 1);
    highRows = min(nAxial, surfaceIndex + 1):nAxial;
    lowPhotons = sum(sliceCounts(lowRows));
    highPhotons = sum(sliceCounts(highRows));
    span = max(20, round(0.05 * nAxial));
    lowNear = max(1, surfaceIndex - span):max(1, surfaceIndex - 1);
    highNear = min(nAxial, surfaceIndex + 1):min(nAxial, surfaceIndex + span);
    lowArrival = weightedArrival(arrivalAx, countsAx, lowNear);
    highArrival = weightedArrival(arrivalAx, countsAx, highNear);

    byPhotons = sign(highPhotons - lowPhotons);
    if byPhotons == 0; byPhotons = 1; end
    byArrival = sign(highArrival - lowArrival);

    direction = struct('sign', byPhotons, ...
        'candidateRows', nnz(bright), ...
        'bandArrivalNs', smooth(surfaceIndex), ...
        'disagree', byArrival ~= 0 && byArrival ~= byPhotons);
    if byPhotons > 0
        direction.name = 'increasing index';
        direction.chosenArrivalNs = highArrival;
        direction.otherArrivalNs = lowArrival;
        direction.chosenPhotons = highPhotons;
        direction.otherPhotons = lowPhotons;
    else
        direction.name = 'decreasing index';
        direction.chosenArrivalNs = lowArrival;
        direction.otherArrivalNs = highArrival;
        direction.chosenPhotons = lowPhotons;
        direction.otherPhotons = highPhotons;
    end
    direction.arrivalGainNs = direction.chosenArrivalNs - ...
        direction.bandArrivalNs;
end

function value = weightedArrival(arrivalAx, countsAx, rows)
    if isempty(rows); value = NaN; return; end
    w = countsAx(rows, :);
    a = arrivalAx(rows, :);
    good = isfinite(a) & w > 0;
    if ~any(good(:)); value = NaN; return; end
    value = sum(a(good) .* w(good)) / sum(w(good));
end

function profile = arrivalVersusDistance(countsAx, arrivalAx, surfaceIndex, ...
        direction)
%ARRIVALVERSUSDISTANCE Photon-weighted mean arrival as a function of height.
% This is the curve that decides where the quenched zone ends, so it is
% returned and printed rather than being consumed silently inside a threshold.
    nAxial = size(countsAx, 1);
    axialIndex = (1:nAxial)';
    offset = direction.sign * (axialIndex - surfaceIndex);
    maxOffset = max(offset);
    distance = (0:maxOffset)';
    arrivalNs = nan(size(distance));
    photons = zeros(size(distance));
    for k = 1:numel(distance)
        rows = axialIndex(offset == distance(k));
        if isempty(rows); continue; end
        photons(k) = sum(sum(countsAx(rows, :)));
        arrivalNs(k) = weightedArrival(arrivalAx, countsAx, rows);
    end
    keep = isfinite(arrivalNs) & photons > 0;
    profile = struct('distance', distance(keep), ...
        'arrivalNs', arrivalNs(keep), 'photons', photons(keep));
end

function top = choosePlateauDistance(profile, opts)
%CHOOSEPLATEAUDISTANCE Where the arrival time has stopped rising.
% Pooling everything above the bilayer mixed partly-quenched dye with free dye
% and produced residuals reaching 20 sigma. The pool therefore starts only
% where the arrival time has reached the requested fraction of its total rise.
    arrival = profile.arrivalNs;
    smooth = conv(arrival, ones(9, 1) / 9, 'same');
    base = smooth(1);
    plateau = quantileLocal(smooth, 0.95);
    target = base + opts.plateauFraction * (plateau - base);
    index = find(smooth >= target, 1, 'first');
    if isempty(index); index = numel(smooth); end
    top = struct('distancePix', profile.distance(index), ...
        'plateauArrivalNs', plateau, 'targetArrivalNs', target);
    if ~isempty(opts.minTopDistancePix)
        top.distancePix = opts.minTopDistancePix;
    end
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

function selection = selectTopPixels(countsAx, surfaceIndex, direction, ...
        topDistancePix, opts)
%SELECTTOPPIXELS Pixels at least topDistancePix above the bilayer, with signal.
    [nAxial, nLateral] = size(countsAx);
    axialIndex = (1:nAxial)';
    offset = direction.sign * (axialIndex - surfaceIndex);
    candidate = repmat(offset >= topDistancePix, 1, nLateral);

    candidateCounts = countsAx(candidate);
    if isempty(candidateCounts)
        error('fit_free_dye_lifetime_above_surface:NoCandidates', ...
            ['No pixels lie %d or more pixels above the bilayer on the ' ...
             'chosen side.'], topDistancePix);
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
        'threshold', threshold, 'topDistancePix', topDistancePix, ...
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

function fit = fitMonoPlusOffset(decay, dtNs, firstBin, lastBin, periodNs)
%FITMONOPLUSOFFSET  y = A*exp(-t/tau) + B, with A,B >= 0 and tau searched.
%
% Variable projection: A and B enter linearly, so only tau is searched and the
% two amplitudes come from PIRLSnonneg at each trial tau. PIRLSnonneg is the
% right inner solver here rather than plain least squares because the data are
% photon counts spanning four decades - unweighted least squares would let the
% bright early bins dominate completely and leave the background essentially
% unconstrained, which is how the previous version ended up with a background
% taken from a median instead of from the fit.
%
% The background is a FITTED parameter, not a pre-subtracted constant. That lets
% the fit run to the end of the PIE window, using the flat region to pin B
% rather than discarding it.
    decay = double(decay(:));
    lastBin = min(lastBin, numel(decay));
    y = decay(firstBin:lastBin);
    t = (0:numel(y) - 1)' * dtNs;

    objective = @(logTau) monoDeviance(exp(logTau), t, y, periodNs);
    logTauHat = fminbnd(objective, log(0.05), log(20), ...
        optimset('TolX', 1e-5, 'Display', 'off'));
    tau = exp(logTauHat);
    [deviance, beta, model] = monoDeviance(tau, t, y, periodNs);

    % Two free amplitudes plus one lifetime.
    dof = max(1, numel(y) - 3);
    fit = struct('tauNs', tau, 'amplitude', beta(2), 'offset', beta(1), ...
        'deviance', deviance, 'reducedDeviance', deviance / dof, ...
        'firstBin', firstBin, 'lastBin', lastBin, 'nBin', numel(y), ...
        'timeNs', t, 'data', y, 'model', model, ...
        'residual', (y - model) ./ sqrt(max(model, 1)));
end

function [deviance, beta, model] = monoDeviance(tau, t, y, periodNs)
    if ~isfinite(tau) || tau <= 0
        deviance = 1e12; beta = [0; 0]; model = ones(size(y)); return;
    end
    column = exp(-t / tau);
    if isfinite(periodNs) && periodNs > 0
        % Incomplete decay from previous pulses. Negligible at tau << period,
        % but it costs nothing and keeps the model honest for longer lifetimes.
        column = column / max(1 - exp(-periodNs / tau), eps);
    end
    design = [ones(numel(t), 1), column];
    beta = PIRLSnonneg(design, y);
    model = max(design * beta, 1e-12);
    good = y > 0;
    deviance = 2 * (sum(model - y) + ...
        sum(y(good) .* log(y(good) ./ model(good))));
    if ~isfinite(deviance); deviance = 1e12; end
end

function fit = fitPooledDecay(pooled, dtNs, periodNs, opts, ptuFile, pedestal)
%FITPOOLEDDECAY Swept-tail-start MLE tail fit, plus an IRF-deconvolved check.
    fullDecay = double(pooled(:));
    % Two windows, deliberately different, because the two fits need different
    % things from the late bins.
    %
    % The PIRLS fit gets the WHOLE window to the end of the PIE gate. Its
    % background B is a fitted non-negative parameter, and the flat bins beyond
    % the decay are precisely what pins it - discarding them leaves B
    % constrained only by the bins where signal and background overlap, which is
    % the weakest possible place to determine it.
    %
    % The swept Tailfit gets the truncated window. It has no proper background
    % model, so trailing pedestal bins drag a mono-exponential downwards and
    % eventually it fits pure noise. That was measured: the sweep fell
    % 2.23 -> 1.99 -> 1.55 ns as the window filled with pedestal.
    pooled = fullDecay(1:pedestal.signalEnd);
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

    % [2] The primary estimator: A*exp(-t/tau) + B with both amplitudes from
    % PIRLSnonneg. Runs to the end of the window because the fitted offset
    % absorbs the background instead of the window having to avoid it.
    pirlsFirst = peakBin + max(1, round(0.4 / dtNs));
    if pirlsFirst < numel(fullDecay) - 8
        pirls = fitMonoPlusOffset(fullDecay, dtNs, pirlsFirst, ...
            numel(fullDecay), periodNs);
    else
        pirls = struct('tauNs', NaN, 'amplitude', NaN, 'offset', NaN, ...
            'reducedDeviance', NaN, 'firstBin', NaN, 'lastBin', NaN, ...
            'nBin', 0, 'timeNs', [], 'data', [], 'model', [], 'residual', []);
    end

    fit = struct('pirls', pirls, ...
        'peakBin', peakBin, 'tauPerTailStart', tauPerStart, ...
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
            taumax = min(20, 4 * numel(tail) * dtNs);
            try
                [amplitude, rate] = DistTailfit(tail, dtNs, 1, false, ...
                    nRate, taumax, periodNs);
            catch
                % The clustered branch is fragile on some window lengths;
                % the unclustered distribution is still usable, so fall back
                % to it rather than losing the estimate entirely.
                [amplitude, rate] = DistTailfit(tail, dtNs, 0, false, ...
                    nRate, taumax, periodNs);
            end
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
    yline(ax, out.surfaceIndex, 'r-', 'LineWidth', 1.2);
    poolStartRow = out.surfaceIndex + out.direction.sign * ...
        out.selection.topDistancePix;
    yline(ax, poolStartRow, 'c--', 'LineWidth', 1);
    xlabel(ax, 'lateral (pixels)'); ylabel(ax, 'axial (pixels)');
    title(ax, sprintf(['%s intensity, log_{10}; bilayer red, pool starts ' ...
        'cyan'], out.plane));

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
    yline(ax, out.surfaceIndex, 'r-', 'LineWidth', 1.2);
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
    % The curve drawn IS the fitted model - amplitude, background and lifetime
    % all from PIRLSnonneg. The previous version drew an exponential scaled to
    % the first few bins with a median pedestal added, which is why the plotted
    % curve and the residuals disagreed with each other.
    pirls = fit.pirls;
    haveFit = isstruct(pirls) && isfinite(pirls.tauNs) && ~isempty(pirls.model);
    if haveFit
        yline(ax, max(pirls.offset, 0.5), '--', 'Color', [0.4 0.4 0.4], ...
            'DisplayName', sprintf('fitted background %.1f counts/bin', ...
            pirls.offset));
        span = pirls.firstBin:pirls.lastBin;
        semilogy(ax, timeNs(span), max(pirls.model, 0.5), 'r-', ...
            'LineWidth', 1.6, 'DisplayName', sprintf( ...
            'A e^{-t/\\tau}+B, \\tau = %.3f ns', pirls.tauNs));
        xline(ax, timeNs(pirls.firstBin), ':k', 'fit start');
    end
    grid(ax, 'on');
    ylabel(ax, 'photons per bin');
    legend(ax, 'Location', 'northeast');
    % The borrowed-IRF cross-check is deliberately NOT in the title: it uses an
    % IRF from another acquisition resampled onto this binning and scatters
    % threefold between scans, so putting it beside the fitted lifetime invites
    % it to be read as a second measurement.
    title(ax, sprintf(['%s (%s): pooled decay above the surface, ' ...
        '\\tau = %.3f ns'], out.acquisition, out.plane, fit.pirls.tauNs));

    ax = nexttile(layout, 4, [1 1]);
    if haveFit
        span = pirls.firstBin:pirls.lastBin;
        plot(ax, timeNs(span), pirls.residual, 'b-', 'LineWidth', 1);
        yline(ax, 0, 'k-');
        yline(ax, 3, ':', 'Color', [0.6 0.6 0.6]);
        yline(ax, -3, ':', 'Color', [0.6 0.6 0.6]);
        ylabel(ax, 'Poisson residual');
        ylim(ax, [-6 6]);
    end
    xlabel(ax, 'time (ns)');
    grid(ax, 'on');

    name = fullfile(out.opts.outputDir, 'fitted_tcspc.png');
    exportgraphics(h, name, 'Resolution', 200);
    close(h);
end
