function out = extract_full_period_tcspc(ptuFile, opts)
%EXTRACT_FULL_PERIOD_TCSPC Ungated TCSPC over the WHOLE laser period.
%
% out = extract_full_period_tcspc(ptuFile)
% out = extract_full_period_tcspc(ptuFile, opts)
%
% Histograms every photon's microtime over the full period, summed over
% detectors, with NO PIE gating. Returns both pulses in one decay.
%
% WHY THIS IS NEEDED
%
% The saved analyses keep only the 640 gate. result.gateDetection.profile looks
% like a full-period histogram - it is 314 bins spanning 50.24 ns - but it is
% zero from 0 to 23.7 ns: only bins 148-303 carry counts. The 485 photons were
% discarded at read time (piePulseIndex = 2, recordedExcitationNm = 640), so
% nothing downstream can recover them and the raw PTU has to be re-read.
%
% WHY IT MATTERS
%
% The 485-excited emission has not finished decaying when the 640 gate opens.
% Inside the 640 gate there is a flat pedestal of ~70-100 counts/bin BEFORE the
% 640 pulse, and in SLB-only pixels the tail falls from ~300 to ~100 counts/bin
% across the last 9 ns, implying a component near 8 ns. A model with a CONSTANT
% background cannot represent a decaying contaminant, so that term absorbs it -
% which corrupts every amplitude, and the amplitudes are the SLB photon counts.
%
% Fitting both pulses together is the way to separate them, and that needs the
% whole period.
%
% opts fields
%   binNs          output bin width (default 0.16, matching the analyses)
%   photonsPerChunk records per read (default 5e6)
%   maxPhotons     stop after this many (default Inf; set it for a quick look)
%   channels       detector channels to include (default all found)
%
% out fields
%   decay          nBin-by-1 counts over one full period
%   dtNs, periodNs, nBin, tAxisNs
%   photons, channelsUsed, perChannel
%   rawResolutionNs, rawBinsPerPeriod

    if nargin < 2 || isempty(opts); opts = struct(); end
    defaults = struct('binNs', 0.16, 'photonsPerChunk', 5e6, ...
        'maxPhotons', Inf, 'channels', []);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if ~isfile(ptuFile)
        error('extract_full_period_tcspc:MissingFile', ...
            'PTU not found: %s', ptuFile);
    end

    head = PTU_Read_Head(ptuFile);
    rawRes = double(head.MeasDesc_Resolution) * 1e9;     % ns per raw bin
    syncRate = double(head.TTResult_SyncRate);
    periodNs = 1e9 / syncRate;
    % Bin width must be an exact multiple of the hardware resolution, otherwise
    % the rebinning smears photons across bins. This is the same provenance rule
    % the rest of the pipeline uses: never derive dt as period/nBin.
    chDiv = max(1, round(opts.binNs / rawRes));
    dtNs = chDiv * rawRes;
    nBin = ceil(periodNs / dtNs);
    fprintf(['extract_full_period_tcspc\n  %s\n  raw resolution %.4f ns, ' ...
        'sync %.1f Hz -> period %.3f ns\n  binning %d raw bins -> dt %.4f ns, ' ...
        '%d bins per period\n'], ptuFile, rawRes, syncRate, periodNs, ...
        chDiv, dtNs, nBin);

    total = double(head.TTResult_NumberOfRecords);
    decay = zeros(nBin, 1);
    perChannel = containers.Map('KeyType', 'double', 'ValueType', 'any');
    photons = 0;
    cnt = 0;
    started = tic;
    while cnt < total && photons < opts.maxPhotons
        [~, tcspc, chan, markers, num] = PTU_Read(ptuFile, ...
            [cnt + 1, opts.photonsPerChunk], head);
        if isempty(num) || num == 0; break; end
        cnt = cnt + num;
        n = min([numel(tcspc), numel(chan), numel(markers)]);
        if n <= 0; continue; end
        tcspc = double(tcspc(1:n));
        chan = double(chan(1:n));
        markers = double(markers(1:n));
        keep = markers == 0;               % real photons, not line/frame marks
        if ~isempty(opts.channels)
            keep = keep & ismember(chan, opts.channels);
        end
        if ~any(keep); continue; end
        tc = tcspc(keep);
        ch = chan(keep);
        % PTU microtimes are 0-based raw bins.
        bin = floor(tc / chDiv) + 1;
        bin = bin(bin >= 1 & bin <= nBin);
        decay = decay + accumarray(bin, 1, [nBin 1]);
        photons = photons + numel(bin);
        for c = unique(ch)'
            b = floor(tc(ch == c) / chDiv) + 1;
            b = b(b >= 1 & b <= nBin);
            if isKey(perChannel, c)
                perChannel(c) = perChannel(c) + accumarray(b, 1, [nBin 1]);
            else
                perChannel(c) = accumarray(b, 1, [nBin 1]);
            end
        end
        fprintf('    %.3g of %.3g records, %.3g photons (%.0f s)\n', ...
            cnt, total, photons, toc(started));
    end

    out = struct();
    out.decay = decay;
    out.dtNs = dtNs;
    out.periodNs = periodNs;
    out.nBin = nBin;
    out.tAxisNs = (0:nBin - 1)' * dtNs;
    out.photons = photons;
    out.channelsUsed = cell2mat(keys(perChannel));
    out.perChannel = perChannel;
    out.rawResolutionNs = rawRes;
    out.rawBinsPerPeriod = periodNs / rawRes;
    out.ptuFile = ptuFile;
    out.elapsedSeconds = toc(started);
    fprintf('  %.4g photons over %d channel(s) in %.0f s\n', photons, ...
        numel(out.channelsUsed), out.elapsedSeconds);
end
