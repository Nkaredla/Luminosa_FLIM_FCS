function out = estimate_free_dye_lifetime_from_crosssection(source, opts)
%ESTIMATE_FREE_DYE_LIFETIME_FROM_CROSSSECTION Unquenched lifetime from XZ/YZ scans.
%
% out = estimate_free_dye_lifetime_from_crosssection(ptuFileOrFolder)
% out = estimate_free_dye_lifetime_from_crosssection(source, opts)
%
% WHY THIS EXISTS
%
% The MIET analysis needs the dye's UNQUENCHED lifetime tau0. Everything the
% XY acquisitions measure is quenched by the metal to some degree, so tau0
% cannot be read off them. The XZ and YZ cross-sections span height, so far
% from the metal the quenching vanishes and the lifetime there IS tau0.
%
% This also matters because the calibration currently on disk
% (calibrationCurve.mat) is inconsistent with the measured pools at both ends:
% it spans only 0.1-30 nm, its most-quenched value is 0.865 ns while the SLB
% component is measured at ~0.34 ns, and it maps the 3.30 ns component to
% 28.5 nm rather than the hundreds of nanometres the internalised-dye reading
% implies. A measured tau0 is the single number most likely to settle which of
% those is wrong, because it pins the asymptote the whole curve must approach.
%
% As a by-product the routine returns lifetime versus axial position, which is
% an EXPERIMENTAL quenching curve for this exact sample and stack - no model
% assumptions about metal thickness, spacer or quantum yield.
%
% THREE THINGS ARE DETERMINED FROM THE DATA, NOT ASSUMED
%
% 1. WHICH IMAGE AXIS IS AXIAL. The PTU header cannot say: PixX and PixY are
%    both 602, so the reader's row/column disambiguation is degenerate for
%    these scans. Instead, the supported bilayer is a flat plane, so it is
%    sharply confined along the axial axis and extended along the lateral one.
%    The axis whose marginal intensity profile is CONCENTRATED is the axial
%    axis. Both are printed with the deciding statistic.
%
% 2. WHICH END IS THE METAL. Quenching shortens the lifetime, so the metal is
%    at whichever end the mean photon arrival time is shortest. Using the
%    lifetime rather than the intensity matters: intensity is confounded by dye
%    density, whereas the arrival-time gradient points at the metal by
%    construction. The axis is then flipped so the index increases AWAY from
%    the metal.
%
% 3. WHERE THE LIFETIME PLATEAUS. tau0 is read off the plateau, not at an
%    arbitrary "far" slice, so the answer does not depend on the axial
%    calibration being known - see the warning below.
%
% THE AXIAL CALIBRATION IS AMBIGUOUS AND IS NOT GUESSED
%
% ImgHdr_PixResol (0.0241 um here) is the lateral pixel size.
% ImgHdr_ZCalibration is 200000 against ImgHdr_XCalibration 62075, a ratio of
% 3.22, but the units of those constants are not documented in the file and
% the two candidate conversions differ by a factor of 10 (0.0075 um/pixel
% versus 0.078 um/pixel, i.e. a 4.5 um or a 47 um sweep over 602 pixels). This
% routine therefore reports distances in PIXELS as the primary axis, gives
% both candidate micron scales side by side, and adds a third, internal
% estimate from the apparent thickness of the bilayer plane - which is a real
% 2D object, so its axial extent is the axial PSF width and that is known
% independently (~0.5-0.6 um FWHM at NA 1.45, 690 nm emission).
%
% tau0 itself does not depend on any of this. Only the quenching curve's
% distance axis does.
%
% METHOD
%
% The primary fit is an IRF-FREE TAIL FIT (Tailfit, Poisson MLE). For a
% lifetime of a few ns with a sub-ns IRF this is the robust choice, and it
% removes any dependence on which IRF is adopted. Because a tail fit's answer
% can drift with where the tail starts, the start is SWEPT and the spread
% reported - a stable tau0 across tail starts is the evidence that the value
% is real. An IRF-deconvolved Fluofit using an IRF from a nearby XY
% acquisition of the same session is run as a cross-check when one is
% available.
%
% INPUTS
%   source   a RawImage.ptu path, or an acquisition folder containing one
%   opts     .piePulseIndex     PIE window to use (default 2, the 640 nm one)
%            .maxNgate          TCSPC bins after rebinning (default 512)
%            .axialBinPixels    slices per fitted bin (default 8)
%            .minPhotonsPerBin  skip bins below this (default 2000)
%            .tailStartNs       tail starts after the peak, swept
%                               (default [0.6 1.0 1.5 2.0])
%            .plateauTolerance  fraction of the total rise within which the
%                               curve counts as flat (default 0.05)
%            .irfSource         path to an XY analysis MAT for the cross-check,
%                               or '' to search the session, or 'none'
%            .outputDir, .makeFigure

    if nargin < 2 || isempty(opts); opts = struct(); end
    defaults = struct( ...
        'piePulseIndex', 2, 'maxNgate', 512, 'tcspcBinNs', 0.05, ...
        'axialBinPixels', 8, 'minPhotonsPerBin', 2000, ...
        'tailStartNs', [0.6 1.0 1.5 2.0], 'plateauTolerance', 0.05, ...
        'irfSource', '', 'outputDir', '', 'makeFigure', true, ...
        'photonsPerChunk', 5e6, 'photonBlockSize', 2e6, ...
        'axialPsfFwhmUm', 0.55);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    for required = {'PTU_Read_Head', 'PTU_FLIM_GPU', 'Tailfit'}
        if exist(required{1}, 'file') ~= 2
            error('estimate_free_dye_lifetime_from_crosssection:Missing', ...
                '%s.m must be on the MATLAB path.', required{1});
        end
    end

    ptuFile = resolvePtuFile(source);
    if isempty(opts.outputDir)
        opts.outputDir = fullfile(fileparts(ptuFile), 'free_dye_lifetime');
    end
    if ~isfolder(opts.outputDir); mkdir(opts.outputDir); end

    head = PTU_Read_Head(ptuFile);
    geometry = crossSectionGeometry(head);
    fprintf('\n=============================================================\n');
    fprintf('estimate_free_dye_lifetime_from_crosssection\n');
    fprintf('  %s\n', ptuFile);
    fprintf('  scan plane %s (ImgHdr_ScanDirection = %g), %d x %d pixels\n', ...
        geometry.plane, geometry.rawScanDirection, geometry.pixX, geometry.pixY);
    fprintf('  lateral pixel %.5f um; period %.3f ns; native resolution %.1f ps\n', ...
        geometry.lateralUm, geometry.periodNs, 1e3 * geometry.nativeDtNs);
    fprintf('  axial scale candidates: %.5f or %.5f um/pixel (see help)\n', ...
        geometry.axialUmCandidateA, geometry.axialUmCandidateB);

    % ---- read photons ------------------------------------------------------
    readOpts = struct('photonsPerChunk', opts.photonsPerChunk, ...
        'computePerFrame', false, 'storeTcspcPix', false, ...
        'storePhotonLists', true, 'storeTimeCell', false, ...
        'showWaitbar', false, ...
        'minLifetimeBin_ns', opts.tcspcBinNs, 'maxNgate', opts.maxNgate, ...
        'computeGlobalMaps', false, 'storePhotonSync', false, ...
        'storePhotonFrame', false, 'lowMemoryPhotonLists', true, ...
        'photonBlockSize', opts.photonBlockSize);
    ptu = PTU_FLIM_GPU(ptuFile, readOpts);
    for field = {'im_tcspc', 'im_line', 'im_col', 'Ngate'}
        if ~isfield(ptu, field{1})
            error('estimate_free_dye_lifetime_from_crosssection:Reader', ...
                'The reader did not return %s.', field{1});
        end
    end
    nGate = double(ptu.Ngate);
    dtNs = geometry.periodNs / nGate;
    fprintf('  read %d photons into %d TCSPC bins of %.4f ns\n', ...
        numel(ptu.im_tcspc), nGate, dtNs);

    % ---- PIE gate ----------------------------------------------------------
    % im_tcspc is already a 1-based bin index in this reader; the previous
    % "+1 then drop the first bin" shifted every arrival time by one bin.
    rawBin = double(ptu.im_tcspc(:));
    if max(rawBin) > nGate || min(rawBin) < 1
        error('estimate_free_dye_lifetime_from_crosssection:BinRange', ...
            'im_tcspc spans %g to %g, outside 1..%d.', ...
            min(rawBin), max(rawBin), nGate);
    end
    globalCurve = accumarray(rawBin, 1, [nGate, 1]);
    gate = selectPieGate(globalCurve, opts.piePulseIndex, dtNs);
    fprintf('  gate detection: %s\n', gate.method);
    fprintf(['  PIE: %d window(s) detected; using window %d, bins %d-%d ' ...
        '(%.2f-%.2f ns), %.3g photons (%.1f%% of all)\n'], ...
        gate.windowCount, gate.index, gate.startBin, gate.stopBin, ...
        (gate.startBin - 1) * dtNs, (gate.stopBin - 1) * dtNs, ...
        gate.photons, 100 * gate.photons / max(sum(globalCurve), 1));

    inGate = rawBin >= gate.startBin & rawBin <= gate.stopBin;
    lineIndex = double(ptu.im_line(:));
    colIndex = double(ptu.im_col(:));
    timeBin = rawBin - gate.startBin + 1;
    gateLength = gate.stopBin - gate.startBin + 1;
    valid = inGate & lineIndex >= 1 & lineIndex <= geometry.pixY & ...
        colIndex >= 1 & colIndex <= geometry.pixX;
    lineIndex = lineIndex(valid);
    colIndex = colIndex(valid);
    timeBin = timeBin(valid);
    fprintf('  %d photons inside the gate and the image bounds\n', numel(timeBin));

    intensity = accumarray([lineIndex, colIndex], 1, ...
        [geometry.pixY, geometry.pixX]);

    % ---- 1. which axis is axial -------------------------------------------
    axial = chooseAxialAxis(intensity);
    fprintf(['\n  [1] axial axis: %s. Concentration statistic (fraction of ' ...
        'signal in the\n      narrowest 10%% of the axis) is %.3f along ' ...
        'im_line and %.3f along im_col;\n      the flat bilayer plane makes ' ...
        'the axial axis the concentrated one.\n'], ...
        axial.name, axial.concentrationLine, axial.concentrationCol);
    if axial.margin < 1.5
        warning(['estimate_free_dye_lifetime_from_crosssection:' ...
            'AmbiguousAxis'], ...
            ['The two axes are similarly concentrated (ratio %.2f). The ' ...
             'axial assignment may be wrong; inspect the intensity figure.'], ...
            axial.margin);
    end
    if axial.isLine
        axialIndex = lineIndex; lateralIndex = colIndex;
        nAxial = geometry.pixY; nLateral = geometry.pixX;
    else
        axialIndex = colIndex; lateralIndex = lineIndex;
        nAxial = geometry.pixX; nLateral = geometry.pixY;
    end

    % ---- 2. which end is the metal ----------------------------------------
    arrival = accumarray(axialIndex, timeBin, [nAxial, 1]);
    perSlice = accumarray(axialIndex, 1, [nAxial, 1]);
    meanArrivalNs = nan(nAxial, 1);
    enough = perSlice >= max(200, 0.02 * max(perSlice));
    meanArrivalNs(enough) = (arrival(enough) ./ perSlice(enough) - 1) * dtNs;
    metal = chooseMetalEnd(meanArrivalNs, perSlice);
    fprintf(['\n  [2] metal end: %s. Mean arrival time %.3f ns over the ' ...
        'first quarter of the\n      populated range and %.3f ns over the ' ...
        'last; the metal is the short end.\n'], ...
        metal.name, metal.leadingMeanNs, metal.trailingMeanNs);
    if metal.contrastNs < 0.05
        warning(['estimate_free_dye_lifetime_from_crosssection:' ...
            'NoArrivalGradient'], ...
            ['The arrival-time contrast across the axis is only %.3f ns, ' ...
             'so the metal end is not clearly identified.'], metal.contrastNs);
    end
    if metal.flip
        axialIndex = nAxial + 1 - axialIndex;
        perSlice = flipud(perSlice);
        meanArrivalNs = flipud(meanArrivalNs);
        intensity = flipAlong(intensity, axial.isLine);
    end

    % Origin: the brightest slice is the bilayer plane, which sits ON the
    % metal, so distances are measured from there rather than from pixel 1.
    [~, originSlice] = max(smoothProfile(perSlice, 5));
    distancePixels = (1:nAxial)' - originSlice;
    fprintf(['      bilayer plane (brightest slice) at index %d; distances ' ...
        'run from there\n'], originSlice);
    bandFwhm = profileFwhm(smoothProfile(perSlice, 5), originSlice);
    internalUmPerPixel = NaN;
    if isfinite(bandFwhm) && bandFwhm > 0
        internalUmPerPixel = opts.axialPsfFwhmUm / bandFwhm;
        fprintf(['      that plane is %.1f pixels FWHM; if its true axial ' ...
            'width is the %.2f um\n      axial PSF, the axial step is ' ...
            '%.5f um/pixel (third, internal estimate)\n'], ...
            bandFwhm, opts.axialPsfFwhmUm, internalUmPerPixel);
    end

    % ---- 3. lifetime versus axial distance --------------------------------
    binEdges = 1:opts.axialBinPixels:(nAxial + 1);
    nBin = numel(binEdges) - 1;
    binAssignment = min(nBin, floor((axialIndex - 1) / opts.axialBinPixels) + 1);
    decays = accumarray([binAssignment, timeBin], 1, [nBin, gateLength]);
    binPhotons = sum(decays, 2);
    binCentrePixels = nan(nBin, 1);
    for b = 1:nBin
        binCentrePixels(b) = mean(distancePixels( ...
            binEdges(b):min(binEdges(b + 1) - 1, nAxial)));
    end

    fprintf(['\n  [3] tail fit per axial bin (%d pixels per bin, Poisson ' ...
        'MLE, no IRF)\n'], opts.axialBinPixels);
    fprintf(['      dist(px)  photons   tau1(ns) over tail starts %s   ' ...
        'spread\n'], mat2str(opts.tailStartNs));
    profile = struct('distancePixels', {}, 'photons', {}, 'tauNs', {}, ...
        'tauPerTailStart', {}, 'spreadNs', {}, 'tauBiNs', {}, ...
        'ampBi', {}, 'decay', {});
    peakBin = argmaxLocal(sum(decays, 1));
    for b = 1:nBin
        if binPhotons(b) < opts.minPhotonsPerBin; continue; end
        decay = decays(b, :)';
        [tauPerStart, tauBi, ampBi] = tailFitSweep(decay, dtNs, ...
            peakBin, opts.tailStartNs);
        entry = struct('distancePixels', binCentrePixels(b), ...
            'photons', binPhotons(b), ...
            'tauNs', median(tauPerStart, 'omitnan'), ...
            'tauPerTailStart', tauPerStart, ...
            'spreadNs', max(tauPerStart) - min(tauPerStart), ...
            'tauBiNs', tauBi, 'ampBi', ampBi, 'decay', decay);
        profile(end + 1) = entry; %#ok<AGROW>
        fprintf('      %8.1f  %7.0f   %s   %.3f\n', ...
            entry.distancePixels, entry.photons, ...
            sprintf('%6.3f ', tauPerStart), entry.spreadNs);
    end
    if isempty(profile)
        error('estimate_free_dye_lifetime_from_crosssection:NoBins', ...
            ['No axial bin reached %d photons. Lower minPhotonsPerBin or ' ...
             'raise axialBinPixels.'], opts.minPhotonsPerBin);
    end

    % ---- 4. the plateau = tau0 --------------------------------------------
    plateau = findPlateau([profile.distancePixels]', [profile.tauNs]', ...
        [profile.photons]', opts.plateauTolerance);
    fprintf('\n  [4] unquenched lifetime from the plateau\n');
    fprintf(['      lifetime rises from %.3f ns at the metal to %.3f ns ' ...
        'far away\n'], plateau.minTauNs, plateau.maxTauNs);
    if plateau.ok
        fprintf(['      plateau starts at %.1f px and covers %d bin(s); ' ...
            'tau0 = %.3f ns (sd %.3f)\n'], ...
            plateau.startPixels, plateau.binCount, plateau.tau0Ns, ...
            plateau.tau0SdNs);
        for label = {'A', 'B', 'internal'}
            switch label{1}
                case 'A'; scale = geometry.axialUmCandidateA;
                case 'B'; scale = geometry.axialUmCandidateB;
                otherwise; scale = internalUmPerPixel;
            end
            if isfinite(scale)
                fprintf('      plateau onset = %.3f um under scale %s\n', ...
                    plateau.startPixels * scale, label{1});
            end
        end
    else
        fprintf(['      NO PLATEAU: the lifetime is still changing at the ' ...
            'far edge, so this scan\n      does not reach unquenched dye. ' ...
            'The largest value seen, %.3f ns, is a LOWER\n      BOUND on ' ...
            'tau0, not a measurement of it.\n'], plateau.maxTauNs);
    end

    % ---- 5. IRF-deconvolved cross-check ------------------------------------
    crossCheck = irfCrossCheck(profile, plateau, dtNs, geometry.periodNs, ...
        opts, ptuFile, gateLength);

    % The per-bin decays were only needed for the cross-check above; keeping
    % them would put a gateLength-by-nBin array in every saved MAT for no
    % further use.
    profile = rmfield(profile, 'decay');

    out = struct('ptuFile', ptuFile, 'geometry', geometry, 'gate', gate, ...
        'axial', axial, 'metal', metal, 'originSlice', originSlice, ...
        'bandFwhmPixels', bandFwhm, 'internalUmPerPixel', internalUmPerPixel, ...
        'distancePixels', distancePixels, 'slicePhotons', perSlice, ...
        'meanArrivalNs', meanArrivalNs, 'intensity', intensity, ...
        'profile', profile, 'plateau', plateau, 'crossCheck', crossCheck, ...
        'dtNs', dtNs, 'opts', opts, 'nLateral', nLateral);

    matFile = fullfile(opts.outputDir, 'free_dye_lifetime.mat');
    save(matFile, 'out', '-v7.3');
    fprintf('\n  wrote %s\n', matFile);
    writeProfileCsv(out, fullfile(opts.outputDir, 'free_dye_lifetime.csv'));
    if opts.makeFigure
        out.figure = plotCrossSection(out);
        fprintf('  wrote %s\n', out.figure);
    end
end

% ===================================================================== setup

function ptuFile = resolvePtuFile(source)
    source = char(source);
    if isfolder(source)
        candidate = fullfile(source, 'RawImage.ptu');
        if ~isfile(candidate)
            found = dir(fullfile(source, '**', 'RawImage.ptu'));
            if isempty(found)
                error(['estimate_free_dye_lifetime_from_crosssection:' ...
                    'NoPtu'], 'No RawImage.ptu below %s', source);
            end
            candidate = fullfile(found(1).folder, found(1).name);
        end
        ptuFile = candidate;
    elseif isfile(source)
        ptuFile = source;
    else
        error('estimate_free_dye_lifetime_from_crosssection:NoSource', ...
            'Not a file or folder: %s', source);
    end
end

function geometry = crossSectionGeometry(head)
    geometry = struct();
    geometry.rawScanDirection = fieldOr(head, 'ImgHdr_ScanDirection', NaN);
    switch geometry.rawScanDirection
        case 0; geometry.plane = 'XY';
        case 1; geometry.plane = 'YZ';
        case 2; geometry.plane = 'XZ';
        otherwise; geometry.plane = 'unknown';
    end
    if geometry.rawScanDirection == 0
        error(['estimate_free_dye_lifetime_from_crosssection:' ...
            'NotACrossSection'], ...
            ['ImgHdr_ScanDirection = 0 means an XY scan, which has no axial ' ...
             'axis. This routine needs an XZ or YZ cross-section.']);
    end
    geometry.pixX = round(fieldOr(head, 'ImgHdr_PixX', NaN));
    geometry.pixY = round(fieldOr(head, 'ImgHdr_PixY', NaN));
    geometry.lateralUm = fieldOr(head, 'ImgHdr_PixResol', NaN);
    syncRate = fieldOr(head, 'TTResult_SyncRate', NaN);
    if isfinite(syncRate) && syncRate > 0
        geometry.periodNs = 1e9 / syncRate;
    else
        geometry.periodNs = 1e9 * fieldOr(head, ...
            'MeasDesc_GlobalResolution', NaN);
    end
    geometry.nativeDtNs = 1e9 * fieldOr(head, 'MeasDesc_Resolution', NaN);
    xCal = fieldOr(head, 'ImgHdr_XCalibration', NaN);
    zCal = fieldOr(head, 'ImgHdr_ZCalibration', NaN);
    % The units of these constants are not documented in the file, and the two
    % readings differ by more than a factor of ten, so both are carried and
    % neither is used silently. See the function help.
    geometry.xCalibration = xCal;
    geometry.zCalibration = zCal;
    geometry.axialUmCandidateA = geometry.lateralUm * xCal / zCal;
    geometry.axialUmCandidateB = geometry.lateralUm * zCal / xCal;
end

function value = fieldOr(s, name, fallback)
    if isfield(s, name) && isnumeric(s.(name)) && isscalar(s.(name)) && ...
            isfinite(s.(name))
        value = double(s.(name));
    else
        value = fallback;
    end
end

% ============================================================ gate selection

function gate = selectPieGate(globalCurve, wantedIndex, dtNs)
%SELECTPIEGATE Split the period at its intensity minima and take one window.
% AutodetectTimeGates is used when available because it is the routine the rest
% of the project relies on; otherwise a local equivalent runs. Either way the
% chosen window is printed so it can be checked against the decay figure.
    globalCurve = double(globalCurve(:));
    nGate = numel(globalCurve);
    starts = []; lengths = [];
    method = 'AutodetectTimeGates';
    if exist('AutodetectTimeGates', 'file') == 2
        try
            [starts, lengths] = AutodetectTimeGates(globalCurve, wantedIndex);
            starts = round(double(starts(:)));
            lengths = round(double(lengths(:)));
        catch gateError
            % AutodetectTimeGates calls mCluster, which is absent from this
            % checkout, so the fallback is the EXPECTED path here rather than
            % an anomaly. Which path ran is printed so the gate stays checkable.
            method = ['local fallback (AutodetectTimeGates failed: ' ...
                gateError.identifier ')'];
            starts = []; lengths = [];
        end
    else
        method = 'local fallback (AutodetectTimeGates not on the path)';
    end
    if isempty(starts) || any(~isfinite(starts)) || any(starts < 1)
        if strcmp(method, 'AutodetectTimeGates')
            method = ['local fallback (AutodetectTimeGates returned ' ...
                'nothing usable)'];
        end
        [starts, lengths] = fallbackGates(globalCurve);
    end
    keep = starts >= 1 & starts <= nGate & lengths > 0;
    starts = starts(keep); lengths = lengths(keep);
    if isempty(starts)
        starts = 1; lengths = nGate;
    end
    index = min(max(1, wantedIndex), numel(starts));
    gate = struct('windowCount', numel(starts), 'index', index, ...
        'startBin', starts(index), ...
        'stopBin', min(nGate, starts(index) + lengths(index) - 1), ...
        'allStarts', starts', 'allLengths', lengths', 'dtNs', dtNs, ...
        'method', method);
    gate.photons = sum(globalCurve(gate.startBin:gate.stopBin));
end

function [starts, lengths] = fallbackGates(curve)
    smooth = smoothProfile(curve, max(3, round(numel(curve) / 100)));
    baseline = median(smooth);
    above = smooth > baseline + 0.15 * (max(smooth) - baseline);
    edges = diff([false; above(:); false]);
    starts = find(edges == 1);
    stops = find(edges == -1) - 1;
    lengths = stops - starts + 1;
end

% ============================================================ axis decisions

function axial = chooseAxialAxis(intensity)
%CHOOSEAXIALAXIS The flat bilayer is confined along the axial axis.
% Statistic: the largest fraction of the total signal that any contiguous 10%
% of the axis holds. A thin plane gives a large value; an extended lateral
% direction gives something near 0.1.
    lineProfile = sum(intensity, 2);
    colProfile = sum(intensity, 1)';
    axial = struct();
    axial.concentrationLine = peakWindowFraction(lineProfile, 0.10);
    axial.concentrationCol = peakWindowFraction(colProfile, 0.10);
    axial.isLine = axial.concentrationLine > axial.concentrationCol;
    if axial.isLine
        axial.name = 'im_line (image rows)';
        axial.margin = axial.concentrationLine / max(axial.concentrationCol, eps);
    else
        axial.name = 'im_col (image columns)';
        axial.margin = axial.concentrationCol / max(axial.concentrationLine, eps);
    end
    axial.lineProfile = lineProfile;
    axial.colProfile = colProfile;
end

function fraction = peakWindowFraction(profile, windowFraction)
    profile = double(profile(:));
    total = sum(profile);
    if total <= 0; fraction = 0; return; end
    width = max(1, round(windowFraction * numel(profile)));
    running = conv(profile, ones(width, 1), 'same');
    fraction = max(running) / total;
end

function metal = chooseMetalEnd(meanArrivalNs, perSlice)
%CHOOSEMETALEND Quenching shortens the lifetime, so the metal is the short end.
% Intensity is deliberately not used: it is confounded by dye density, whereas
% the arrival-time gradient points at the metal by construction.
    populated = find(perSlice >= max(200, 0.02 * max(perSlice)));
    metal = struct('leadingMeanNs', NaN, 'trailingMeanNs', NaN, ...
        'contrastNs', 0, 'flip', false, 'name', 'undetermined');
    if numel(populated) < 8; return; end
    span = numel(populated);
    lead = populated(1:max(1, round(0.25 * span)));
    trail = populated(end - max(1, round(0.25 * span)) + 1:end);
    metal.leadingMeanNs = mean(meanArrivalNs(lead), 'omitnan');
    metal.trailingMeanNs = mean(meanArrivalNs(trail), 'omitnan');
    metal.contrastNs = abs(metal.trailingMeanNs - metal.leadingMeanNs);
    if metal.leadingMeanNs <= metal.trailingMeanNs
        metal.flip = false;
        metal.name = 'low index (no flip needed)';
    else
        metal.flip = true;
        metal.name = 'high index (axis flipped so distance grows upward)';
    end
end

function image = flipAlong(image, alongRows)
    if alongRows; image = flipud(image); else; image = fliplr(image); end
end

% =================================================================== fitting

function [tauPerStart, tauBi, ampBi] = tailFitSweep(decay, dtNs, peakBin, ...
        tailStartsNs)
%TAILFITSWEEP Mono-exponential tail fit from several tail starts, plus a 2-exp.
% Sweeping the tail start is the whole robustness argument for an IRF-free
% fit: if the answer is stable as the start moves past the IRF, the IRF is no
% longer contributing and the lifetime is real. A drifting answer means the
% tail is contaminated and the number should not be trusted.
    decay = double(decay(:));
    tauPerStart = nan(1, numel(tailStartsNs));
    for k = 1:numel(tailStartsNs)
        first = peakBin + max(1, round(tailStartsNs(k) / dtNs));
        if first >= numel(decay) - 8; continue; end
        tail = decay(first:end);
        if sum(tail) < 200; continue; end
        try
            tau = Tailfit(tail, dtNs, 2.0, [0.1 12], 'mle', false);
            tauPerStart(k) = tau(1);
        catch
            tauPerStart(k) = NaN;
        end
    end
    tauBi = [NaN NaN]; ampBi = [NaN NaN];
    first = peakBin + max(1, round(min(tailStartsNs) / dtNs));
    if first < numel(decay) - 12
        tail = decay(first:end);
        if sum(tail) >= 1000
            try
                [tau, amp] = Tailfit(tail, dtNs, [0.6 3.0], [0.1 12], ...
                    'mle', false);
                tau = tau(:)'; amp = amp(:)';
                n = min(2, numel(tau));
                tauBi(1:n) = tau(1:n);
                ampBi(1:n) = amp(1:n) / max(sum(amp(1:n)), eps);
            catch
                % leave as NaN
            end
        end
    end
end

function plateau = findPlateau(distance, tau, photons, tolerance)
%FINDPLATEAU The furthest run of bins whose lifetime is flat within tolerance.
% tau0 is read off a plateau rather than at the furthest bin, so the answer
% does not require knowing the axial calibration, and so a scan that never
% reaches unquenched dye reports a LOWER BOUND instead of a wrong number.
    plateau = struct('ok', false, 'tau0Ns', NaN, 'tau0SdNs', NaN, ...
        'startPixels', NaN, 'binCount', 0, ...
        'minTauNs', min(tau), 'maxTauNs', max(tau), 'photons', 0);
    positive = distance > 0;
    if nnz(positive) < 3; return; end
    d = distance(positive); t = tau(positive); p = photons(positive);
    [d, order] = sort(d); t = t(order); p = p(order);
    rise = max(t) - min(t);
    if rise <= 0; return; end
    band = tolerance * rise;
    reference = t(end);
    % Walk inward from the far edge while the lifetime stays within the band.
    first = numel(t);
    while first > 1 && abs(t(first - 1) - reference) <= band
        first = first - 1;
    end
    if numel(t) - first + 1 < 2; return; end
    weights = p(first:end) / sum(p(first:end));
    plateau.ok = true;
    plateau.tau0Ns = sum(weights .* t(first:end));
    plateau.tau0SdNs = std(t(first:end));
    plateau.startPixels = d(first);
    plateau.binCount = numel(t) - first + 1;
    plateau.photons = sum(p(first:end));
end

function crossCheck = irfCrossCheck(profile, plateau, dtNs, periodNs, opts, ...
        ptuFile, gateLength)
%IRFCROSSCHECK IRF-deconvolved fit of the plateau decay, if an IRF is available.
% The IRF is taken from an XY acquisition of the same session rather than from
% this scan, which avoids estimating the IRF from the same photons being
% fitted. The tail fit remains the primary number precisely because it needs
% none of this.
    crossCheck = struct('available', false, 'reason', '', 'tauNs', NaN, ...
        'tauSdNs', NaN, 'photons', 0, 'irfSource', '');
    if strcmpi(char(opts.irfSource), 'none')
        crossCheck.reason = 'disabled by opts.irfSource';
        fprintf('\n  [5] IRF cross-check skipped (disabled)\n');
        return;
    end
    if ~plateau.ok
        crossCheck.reason = 'no plateau to cross-check';
        fprintf('\n  [5] IRF cross-check skipped (no plateau)\n');
        return;
    end
    matFile = char(opts.irfSource);
    if isempty(matFile)
        matFile = findSessionIrf(ptuFile);
    end
    if isempty(matFile) || ~isfile(matFile)
        crossCheck.reason = 'no XY analysis MAT found in the session';
        fprintf('\n  [5] IRF cross-check skipped (%s)\n', crossCheck.reason);
        return;
    end
    if exist('Fluofit', 'file') ~= 2
        crossCheck.reason = 'Fluofit.m is not on the path';
        fprintf('\n  [5] IRF cross-check skipped (%s)\n', crossCheck.reason);
        return;
    end
    loaded = load(matFile, 'result');
    if ~isfield(loaded, 'result') || ~isfield(loaded.result, 'irf') || ...
            ~isfield(loaded.result.irf, 'curve')
        crossCheck.reason = 'the analysis MAT has no stored IRF';
        fprintf('\n  [5] IRF cross-check skipped (%s)\n', crossCheck.reason);
        return;
    end
    irf = double(loaded.result.irf.curve(:));
    % The stored IRF belongs to the pipeline's own binning; resample it onto
    % this scan's bins rather than assuming they match.
    irf = resampleCurve(irf, gateLength);
    irf = max(irf, 0); irf = irf / max(sum(irf), eps);
    far = [profile.distancePixels] >= plateau.startPixels;
    crossCheck.irfSource = matFile;
    try
        decay = plateauDecay(profile, far, gateLength);
        % Fluofit returns [c, offset, A, tau, ...]; plotFlag off, no initial
        % guess routine, single exponential seeded at the tail-fit answer.
        [~, ~, ~, tau, ~, dtau] = Fluofit(irf, decay, periodNs, dtNs, ...
            plateau.tau0Ns, [], 0, 'mle', false);
        tau = double(tau(:)');
        crossCheck.tauNs = tau(1);
        if ~isempty(dtau); crossCheck.tauSdNs = double(dtau(1)); end
        crossCheck.photons = sum(decay);
        crossCheck.available = true;
    catch fitError
        crossCheck.reason = fitError.message;
    end
    if crossCheck.available
        fprintf(['\n  [5] IRF cross-check: tau = %.3f ns (tail fit gave ' ...
            '%.3f ns)\n      IRF from %s\n'], crossCheck.tauNs, ...
            plateau.tau0Ns, matFile);
    else
        fprintf('\n  [5] IRF cross-check failed: %s\n', crossCheck.reason);
    end
end

function decay = plateauDecay(profile, mask, gateLength)
    decay = zeros(gateLength, 1);
    found = false;
    for k = find(mask)
        if isfield(profile(k), 'decay') && ~isempty(profile(k).decay)
            decay = decay + double(profile(k).decay(:));
            found = true;
        end
    end
    if ~found
        error(['estimate_free_dye_lifetime_from_crosssection:' ...
            'NoStoredDecay'], ...
            'Per-bin decays were not retained, so the IRF fit cannot run.');
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

function resampled = resampleCurve(curve, targetLength)
    curve = double(curve(:));
    if numel(curve) == targetLength; resampled = curve; return; end
    source = linspace(0, 1, numel(curve))';
    target = linspace(0, 1, targetLength)';
    resampled = interp1(source, curve, target, 'linear', 0);
end

% =================================================================== helpers

function smooth = smoothProfile(profile, width)
    profile = double(profile(:));
    width = max(1, round(width));
    if width <= 1; smooth = profile; return; end
    smooth = conv(profile, ones(width, 1) / width, 'same');
end

function index = argmaxLocal(values)
    [~, index] = max(double(values(:)));
end

function fwhm = profileFwhm(profile, peakIndex)
    profile = double(profile(:));
    baseline = median(profile(profile <= median(profile)));
    peak = profile(peakIndex);
    if ~isfinite(peak) || peak <= baseline; fwhm = NaN; return; end
    half = baseline + 0.5 * (peak - baseline);
    left = peakIndex;
    while left > 1 && profile(left) > half; left = left - 1; end
    right = peakIndex;
    while right < numel(profile) && profile(right) > half; right = right + 1; end
    fwhm = right - left;
end

function writeProfileCsv(out, csvFile)
    profile = out.profile;
    if isempty(profile); return; end
    rows = struct('distancePixels', {}, 'distanceUmCandidateA', {}, ...
        'distanceUmInternal', {}, 'photons', {}, 'tauNs', {}, ...
        'tauSpreadNs', {}, 'tauShortNs', {}, 'tauLongNs', {}, ...
        'ampShort', {});
    for k = 1:numel(profile)
        rows(k) = struct( ...
            'distancePixels', profile(k).distancePixels, ...
            'distanceUmCandidateA', profile(k).distancePixels * ...
                out.geometry.axialUmCandidateA, ...
            'distanceUmInternal', profile(k).distancePixels * ...
                out.internalUmPerPixel, ...
            'photons', profile(k).photons, ...
            'tauNs', profile(k).tauNs, ...
            'tauSpreadNs', profile(k).spreadNs, ...
            'tauShortNs', profile(k).tauBiNs(1), ...
            'tauLongNs', profile(k).tauBiNs(2), ...
            'ampShort', profile(k).ampBi(1));
    end
    writetable(struct2table(rows), csvFile);
    fprintf('  wrote %s\n', csvFile);
end

function name = plotCrossSection(out)
    h = figure('Color', 'w', 'Visible', 'off', 'Position', [40 40 1500 420]);
    layout = tiledlayout(h, 1, 3, 'Padding', 'compact', ...
        'TileSpacing', 'compact');

    ax = nexttile(layout);
    image = out.intensity;
    if ~out.axial.isLine; image = image'; end
    imagesc(ax, log10(double(image) + 1));
    axis(ax, 'image'); colormap(ax, gray); colorbar(ax);
    xlabel(ax, 'lateral (pixels)'); ylabel(ax, 'axial (pixels, metal at top)');
    title(ax, sprintf('%s intensity, log10', out.geometry.plane));

    ax = nexttile(layout);
    yyaxis(ax, 'left');
    plot(ax, out.distancePixels, out.slicePhotons, 'LineWidth', 1.4);
    ylabel(ax, 'photons per axial slice');
    yyaxis(ax, 'right');
    plot(ax, out.distancePixels, out.meanArrivalNs, 'LineWidth', 1.4);
    ylabel(ax, 'mean arrival time (ns)');
    xlabel(ax, 'distance from the bilayer plane (pixels)');
    grid(ax, 'on');
    title(ax, 'intensity and arrival time versus height');

    ax = nexttile(layout);
    hold(ax, 'on');
    distance = [out.profile.distancePixels];
    tau = [out.profile.tauNs];
    spread = [out.profile.spreadNs];
    errorbar(ax, distance, tau, spread / 2, '-o', 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'auto', 'DisplayName', 'tail fit (mono)');
    if out.plateau.ok
        yline(ax, out.plateau.tau0Ns, '--k', ...
            sprintf('\\tau_0 = %.3f ns', out.plateau.tau0Ns));
        xline(ax, out.plateau.startPixels, ':k', 'plateau onset');
    end
    xlabel(ax, 'distance from the bilayer plane (pixels)');
    ylabel(ax, 'lifetime (ns)');
    grid(ax, 'on');
    title(ax, 'quenching curve, measured');
    legend(ax, 'Location', 'southeast');

    name = fullfile(out.opts.outputDir, 'free_dye_lifetime.png');
    exportgraphics(h, name, 'Resolution', 200);
    close(h);
end
