function [data, tag, life_imm, life_imm_c] = Luminosa_Process_scan(name, num_PIE, timegate, flagint, cutoff, flagfit, pic)
%LUMINOSA_PROCESS_SCAN Local pixel-wise MIET lifetime extraction for Luminosa.
% This replaces the upstream per-pixel PixelTau path with an IRF-aware
% DistFluofit-style fit solved through the local PIRLS GPU/CPU solver.

    estimatorLabel = 'Luminosa_DistFluofit_extension_PIRLSGPU_v1';

    if nargin < 7 || isempty(pic)
        pic = 0;
    end
    if nargin < 6 || isempty(flagfit)
        flagfit = 1;
    end
    if nargin < 5 || isempty(cutoff)
        cutoff = 0.5;
    end
    if nargin < 4 || isempty(flagint) || flagint == 0
        flagint = 0;
    else
        flagint = 1;
    end
    if nargin < 3 || isempty(timegate)
        timegate = [];
    end
    if nargin < 2 || isempty(num_PIE)
        num_PIE = 1;
    end

    [head, im_sync, im_tcspc, im_chan, im_line, im_col] = readScanData(name); %#ok<ASGLU>

    dind = unique(im_chan);
    maxch = numel(dind);
    maxres = max([head.Resolution]);
    Resolution = max([maxres 0.064]);
    chDiv = Resolution / maxres;
    im_tcspc = ceil(double(im_tcspc) ./ chDiv);
    NgateRaw = double(max(im_tcspc));
    nx = double(head.ImgHdr.PixX);
    ny = double(head.ImgHdr.PixY);

    tcspcdata = zeros(maxch, NgateRaw);
    indel = false(maxch, 1);
    for ch = 1:maxch
        tcspcdata(ch, :) = hist(double(im_tcspc(im_chan == dind(ch))), 1:NgateRaw);
        indel(ch) = sum(tcspcdata(ch, :)) < nx * ny;
    end
    tcspcdata(indel, :) = [];
    dind(indel) = [];
    maxch_n = numel(dind);
    tcspcdata = tcspcdata.';

    data = struct();
    data.t = (1:NgateRaw) .* Resolution;
    data.tcspcdata = tcspcdata;

    if isempty(timegate)
        timegate = DetectTimeGates(tcspcdata, num_PIE, Resolution);
        if isempty(timegate)
            Ngate = NgateRaw;
        else
            Ngate = 1 + timegate(1, 2) + timegate(1, 4) - timegate(1, 1);
        end
    else
        Ngate = 1 + timegate(1, 2) + timegate(1, 4) - timegate(1, 1);
    end
    Ngate = max(1, round(double(Ngate)));

    tcspc = zeros(Ngate, maxch_n, num_PIE);
    planePhoton = cell(maxch_n * num_PIE, 1);
    planeLocalPix = cell(maxch_n * num_PIE, 1);
    planeMicro = cell(maxch_n * num_PIE, 1);

    for pulseIdx = 1:num_PIE
        for chIdx = 1:maxch_n
            planeIdx = (pulseIdx - 1) * maxch_n + chIdx;
            gateRow = timegate(planeIdx, :);
            [mask, micro] = buildPlanePhotonMask(im_tcspc, im_chan, dind(chIdx), gateRow);
            planePhoton{planeIdx} = mask;
            planeMicro{planeIdx} = micro(:);
            planeLocalPix{planeIdx} = double(im_line(mask)) + (double(im_col(mask)) - 1) * ny;
            tcspc(:, chIdx, pulseIdx) = buildPlaneTcspcFromGate(tcspcdata(:, chIdx), gateRow, Ngate);
        end
    end

    data.t = (1:Ngate) .* Resolution;
    data.timegate = timegate;
    data.tcspc = tcspc;

    pos = zeros(maxch_n * num_PIE, 1);
    for pulseIdx = 1:num_PIE
        for chIdx = 1:maxch_n
            planeIdx = (pulseIdx - 1) * maxch_n + chIdx;
            [~, pos(planeIdx)] = max(tcspc(:, chIdx, pulseIdx));
        end
    end
    shift = ceil(cutoff / Resolution);

    tag = zeros(ny, nx, maxch_n, num_PIE);
    life_imm = nan(ny, nx, maxch_n, num_PIE);
    life_imm_c = nan(ny, nx, maxch_n, num_PIE);
    fitInfo = cell(maxch_n * num_PIE, 1);
    useGPU = localGpuAvailable();
    pulsePeriodNs = 1e9 / head.SyncRate;

    for pulseIdx = 1:num_PIE
        for chIdx = 1:maxch_n
            planeIdx = (pulseIdx - 1) * maxch_n + chIdx;
            localPix = planeLocalPix{planeIdx};
            micro = planeMicro{planeIdx};
            if isempty(localPix)
                continue;
            end

            cutoffBin = min(Ngate, max(1, pos(planeIdx) + shift));
            if flagint
                planeTagMask = micro > cutoffBin;
            else
                planeTagMask = true(size(micro));
            end
            planeTag = accumarray(localPix(planeTagMask), 1, [nx * ny, 1], @sum, 0);
            tag(:, :, chIdx, pulseIdx) = reshape(planeTag, [ny, nx]);

            if flagfit
                [planeLife, planeInfo] = fitPlaneLifetimeDistFluo( ...
                    planeTag, localPix, micro, squeeze(tcspc(:, chIdx, pulseIdx)), ...
                    head, pulsePeriodNs, Resolution, useGPU);
            else
                [planeLife, planeInfo] = meanArrivalLifetimePlane( ...
                    planeTag, localPix, micro, cutoffBin, Resolution);
            end

            life_imm(:, :, chIdx, pulseIdx) = reshape(planeLife, [ny, nx]);
            life_imm_c(:, :, chIdx, pulseIdx) = reshape(planeLife, [ny, nx]);
            fitInfo{planeIdx} = planeInfo;
        end
    end

    data.tag = tag;
    data.life_imm = life_imm;
    data.life_imm_c = life_imm_c;
    data.pos = pos;
    data.shift = shift;
    data.cutoff = cutoff;
    data.intensityImageFlag = flagint;
    data.flagfit = flagfit;
    data.fitInfo = fitInfo;
    data.lifetimeEstimator = estimatorLabel;
    data.usedGPU = any(cellfun(@(s) isstruct(s) && isfield(s, 'usedGPU') && s.usedGPU, fitInfo));

    if pic
        showLifetimePreview(data, tag, life_imm, head);
    end

    save([name(1:end-4) '_PS'], 'data', 'life_imm', 'life_imm_c', 'tag');
end

function [head, im_sync, im_tcspc, im_chan, im_line, im_col] = readScanData(name)
    if exist([name(1:end-4) '_Core_Scan.mat'], 'file')
        load([name(1:end-4) '_Core_Scan.mat'], 'head', 'im_sync', 'im_tcspc', 'im_chan', 'im_line', 'im_col');
        return;
    end

    if strcmpi(name(end-2:end), 'ht3')
        disp('Read tcspc-data from ht3-file...');
        [head, im_sync, im_tcspc, im_chan, im_line, im_col] = Core_ScanRead(name);
        save([name(1:end-4) '_Core_Scan.mat'], 'head', 'im_sync', 'im_tcspc', 'im_chan', 'im_line', 'im_col');
        return;
    end

    if strcmpi(name(end-2:end), 'ptu')
        disp('Read tcspc-data from ptu-file...');
        [head, im_sync, im_tcspc, im_chan, im_line, im_col] = PTU_ScanRead(name);
        head.ImgHdr.PixX = head.ImgHdr_PixX;
        head.ImgHdr.PixY = head.ImgHdr_PixY;
        head.ImgHdr.X0 = head.ImgHdr_X0;
        head.ImgHdr.Y0 = head.ImgHdr_Y0;
        head.ImgHdr.PixelSize = head.ImgHdr_PixResol;
        head.ImgHdr.PixelTime = head.ImgHdr_DwellTime;
        head.Resolution = 1e9 * head.MeasDesc_Resolution;
        head.SyncRate = 1 ./ head.MeasDesc_GlobalResolution;
        save([name(1:end-4) '_Core_Scan.mat'], 'head', 'im_sync', 'im_tcspc', 'im_chan', 'im_line', 'im_col');
        return;
    end

    error('You have to give an ht3- or ptu-file.');
end

function [mask, micro] = buildPlanePhotonMask(im_tcspc, im_chan, channelId, gateRow)
    chanMask = (im_chan == channelId);
    if gateRow(3) == 0
        mask = chanMask & im_tcspc >= gateRow(1) & im_tcspc <= gateRow(2);
        raw = double(im_tcspc(mask));
        micro = raw - gateRow(1) + 1;
        return;
    end

    firstMask = chanMask & im_tcspc >= gateRow(1) & im_tcspc <= gateRow(2);
    lastMask = chanMask & im_tcspc >= gateRow(3) & im_tcspc <= gateRow(4);
    mask = firstMask | lastMask;
    raw = double(im_tcspc(mask));
    firstRegion = raw >= gateRow(1);
    micro = zeros(size(raw));
    micro(firstRegion) = raw(firstRegion) - gateRow(1) + 1;
    micro(~firstRegion) = raw(~firstRegion) + gateRow(2) - gateRow(1) + 1;
end

function decay = buildPlaneTcspcFromGate(tcspcCol, gateRow, Ngate)
    if gateRow(3) == 0
        decay = tcspcCol(gateRow(1):gateRow(2));
    else
        decay = [tcspcCol(gateRow(1):gateRow(2)); tcspcCol(gateRow(3):gateRow(4))];
    end
    decay = decay(:);
    if numel(decay) < Ngate
        decay(end+1:Ngate, 1) = 0;
    elseif numel(decay) > Ngate
        decay = decay(1:Ngate);
    end
end

function [planeLife, info] = fitPlaneLifetimeDistFluo(planeTag, localPix, micro, globalDecay, head, pulsePeriodNs, dtNs, useGPU)
    minPhotons = 100;
    nLocal = numel(planeTag);
    planeLife = nan(nLocal, 1);
    info = struct('method', 'distfluofit_pirls', 'usedGPU', false, 'globalShiftBins', 0, ...
        'minPhotons', minPhotons, 'nPixelsFit', 0);

    validLocalPix = find(planeTag > minPhotons);
    if isempty(validLocalPix)
        return;
    end

    mapLocalToValid = zeros(nLocal, 1, 'uint32');
    mapLocalToValid(validLocalPix) = uint32(1:numel(validLocalPix));
    validPos = double(mapLocalToValid(localPix));
    keep = validPos > 0 & micro >= 1 & micro <= numel(globalDecay);
    if ~any(keep)
        return;
    end

    Y = accumarray([validPos(keep) micro(keep)], 1, [numel(validLocalPix), numel(globalDecay)], @sum, 0);
    Y = double(Y.');

    globalDecay = max(double(globalDecay(:)), 0);
    if ~any(globalDecay > 0)
        return;
    end

    irf = estimatePlaneIrf(head, globalDecay);
    [tauMean, fitInfo] = fitLifetimeBatchDistFluo(Y, irf, globalDecay, pulsePeriodNs, dtNs, useGPU);
    planeLife(validLocalPix) = tauMean(:);

    info.usedGPU = fitInfo.usedGPU;
    info.globalShiftBins = fitInfo.globalShiftBins;
    info.nPixelsFit = numel(validLocalPix);
    info.tauBasis = fitInfo.tauBasis;
end

function irf = estimatePlaneIrf(head, globalDecay)
    irf = [];
    try
        irfRaw = Calc_mIRF(head, reshape(globalDecay, 1, []));
        if ndims(irfRaw) >= 3
            irf = squeeze(irfRaw(1, :, 1));
        else
            irf = squeeze(irfRaw(1, :));
        end
    catch
    end

    irf = max(double(irf(:)), 0);
    if isempty(irf) || ~any(irf > 0)
        irf = zeros(numel(globalDecay), 1);
        [~, peakIdx] = max(globalDecay);
        irf(max(1, peakIdx)) = 1;
    end
end

function [tauMean, info] = fitLifetimeBatchDistFluo(Y, irf, globalDecay, pulsePeriodNs, dtNs, useGPU)
    nBasis = 100;
    shiftInt = -5:5;
    coarseErr = nan(size(shiftInt));

    for idx = 1:numel(shiftInt)
        [Mfit, Mfull, indFit, ~] = buildDistFluoModel(irf, pulsePeriodNs, dtNs, shiftInt(idx), nBasis, pulsePeriodNs);
        beta = solvePirlsNonnegBatch(Mfit, globalDecay(indFit), useGPU);
        z = Mfull * beta;
        coarseErr(idx) = sum((z(indFit) - globalDecay(indFit)).^2 ./ max(abs(z(indFit)), eps)) ./ max(numel(indFit), 1);
    end

    shiftFine = shiftInt;
    if all(isfinite(coarseErr))
        tmpShift = shiftInt(1):0.1:shiftInt(end);
        tmpErr = interp1(shiftInt, coarseErr, tmpShift, 'linear');
        if any(isfinite(tmpErr))
            shiftFine = tmpShift;
            coarseErr = tmpErr;
        end
    end

    validErr = isfinite(coarseErr);
    if ~any(validErr)
        bestShift = 0;
    else
        shiftCandidates = shiftFine(validErr);
        errCandidates = coarseErr(validErr);
        [~, bestIdx] = min(errCandidates);
        bestShift = shiftCandidates(bestIdx);
    end

    [Mfit, ~, indFit, tauBasis] = buildDistFluoModel(irf, pulsePeriodNs, dtNs, bestShift, nBasis, pulsePeriodNs);
    [beta, solverInfo] = solvePirlsNonnegBatch(Mfit, Y(indFit, :), useGPU);

    betaComp = max(double(beta(2:end, :)), 0);
    tauBasis = double(tauBasis(:));
    weightSum = sum(betaComp, 1);
    tauMean = sum(betaComp .* tauBasis, 1) ./ max(weightSum, eps);
    tauMean(weightSum <= 0) = NaN;
    tauMean = tauMean(:);

    info = struct();
    info.usedGPU = solverInfo.usedGPU;
    info.globalShiftBins = bestShift;
    info.tauBasis = tauBasis;
end

function [Mfit, Mfull, indFit, tauBasis] = buildDistFluoModel(irf, pulsePeriodNs, dtNs, shiftBins, nBasis, tauMaxNs)
    if nargin < 6 || isempty(tauMaxNs)
        tauMaxNs = pulsePeriodNs;
    end

    irf = double(irf(:));
    n = numel(irf);
    nPulse = max(1, floor(pulsePeriodNs / dtNs));
    tp = dtNs * (1:nPulse)';
    tauBasis = (1 / dtNs) ./ exp((0:nBasis) / nBasis * log(tauMaxNs / dtNs));

    M0 = [ones(n, 1) Convol(irf, exp(-tp * tauBasis))];
    M0 = M0 ./ max(sum(M0, 1), eps);

    t = (1:n)';
    shFloor = floor(shiftBins);
    shCeil = ceil(shiftBins);
    Mfull = (1 - shiftBins + shFloor) * M0(rem(rem(t - shFloor - 1, n) + n, n) + 1, :) + ...
        (shiftBins - shFloor) * M0(rem(rem(t - shCeil - 1, n) + n, n) + 1, :);

    c = ceil(abs(shiftBins)) * sign(shiftBins);
    indFit = max(1, 1 + c):min(n, n + c);
    Mfit = Mfull(indFit, :);
    tauBasis = tauBasis(:);
end

function [beta, info] = solvePirlsNonnegBatch(M, Y, useGPU)
    if isvector(Y)
        Y = double(Y(:));
    else
        Y = double(Y);
    end
    M = double(M);

    if size(Y, 1) ~= size(M, 1)
        error('Design matrix and data must have the same number of rows.');
    end

    info = struct('usedGPU', false, 'method', 'cpu_projected_gradient');

    if useGPU
        try
            [betaGpu, infoGpu] = gpuProjectedGradientPirls(M, Y, 10, 25, 4096);
            beta = double(betaGpu);
            info.usedGPU = isfield(infoGpu, 'usedGPU') && infoGpu.usedGPU;
            info.method = infoGpu.method;
            return;
        catch
        end
    end

    beta = cpuProjectedGradientPirls(M, Y, 10, 50);
end

function [beta, info] = gpuProjectedGradientPirls(M, Y, maxPirlsIter, maxPgIter, batchSize)
    [nSamples, nBasis] = size(M);
    [nSamplesY, nPixels] = size(Y);
    if nSamplesY ~= nSamples
        error('M and Y must have the same number of rows.');
    end

    Msingle = single(M);
    Ysingle = single(Y);
    beta = zeros(nBasis, nPixels, 'single');
    info = struct('usedGPU', false, 'method', 'gpu_projected_gradient');

    gpuDevice();
    Mg = gpuArray(Msingle);
    Ikg = eye(nBasis, 'single', 'gpuArray');
    tiny = gpuArray(single(max(1e-6, 0.1 / max(1, nSamples))));
    MG = reshape(Mg, nSamples, nBasis, 1);

    for startIdx = 1:batchSize:nPixels
        pixelIdx = startIdx:min(startIdx + batchSize - 1, nPixels);
        batchCount = numel(pixelIdx);
        Yg = gpuArray(Ysingle(:, pixelIdx));

        MtM = Mg' * Mg;
        Mty = Mg' * Yg;
        betaGpu = max((MtM + 1e-6 * Ikg) \ Mty, 0);

        for iterIdx = 1:maxPirlsIter
            mu = max(Mg * betaGpu, tiny);
            W = 1 ./ mu;
            W3 = reshape(W, nSamples, 1, batchCount);
            MW = MG .* W3;
            Aw = pagemtimes(permute(MG, [2 1 3]), MW);
            YW = reshape(W .* Yg, nSamples, 1, batchCount);
            bw = pagemtimes(permute(MG, [2 1 3]), YW);
            beta3 = reshape(betaGpu, nBasis, 1, batchCount);
            lipschitz = max(sum(abs(Aw), 2), 1e-6);
            lipschitz = max(lipschitz, [], 1);

            for pgIdx = 1:maxPgIter
                grad = pagemtimes(Aw, beta3) - bw;
                beta3 = max(beta3 - grad ./ lipschitz, 0);
            end

            betaGpu = reshape(beta3, nBasis, batchCount);
        end

        beta(:, pixelIdx) = gather(betaGpu);
    end

    info.usedGPU = true;
    info.method = 'gpu_batch_pirls_projected_gradient';
end

function betaAll = cpuProjectedGradientPirls(M, Y, maxPirlsIter, maxPgIter)
    [~, nBasis] = size(M);
    nPixels = size(Y, 2);
    tiny = max(1e-6, 0.1 / max(1, size(M, 1)));
    reg = 1e-6 * eye(nBasis);
    MtM = M' * M + reg;
    Mty = M' * Y;
    betaAll = max(MtM \ Mty, 0);

    for pixIdx = 1:nPixels
        beta = betaAll(:, pixIdx);
        y = Y(:, pixIdx);
        for iterIdx = 1:maxPirlsIter
            mu = max(M * beta, tiny);
            w = 1 ./ mu;
            Aw = M' * (M .* w) + reg;
            bw = M' * (w .* y);
            betaNew = projectedGradientNnls(Aw, bw, maxPgIter, beta);
            if sum((betaNew - beta) .^ 2) < 1e-10
                beta = betaNew;
                break;
            end
            beta = betaNew;
        end
        betaAll(:, pixIdx) = beta;
    end
end

function beta = projectedGradientNnls(A, b, maxIter, beta)
    if nargin < 4 || isempty(beta)
        beta = max(A \ b, 0);
    else
        beta = max(double(beta(:)), 0);
    end

    L = max(sum(abs(A), 2));
    if ~isfinite(L) || L <= 0
        L = 1;
    end

    for iterIdx = 1:maxIter
        grad = A * beta - b;
        betaNew = max(beta - grad ./ L, 0);
        if sum((betaNew - beta) .^ 2) < 1e-10
            beta = betaNew;
            break;
        end
        beta = betaNew;
    end
end

function [planeLife, info] = meanArrivalLifetimePlane(planeTag, localPix, micro, cutoffBin, dtNs)
    planeLife = nan(numel(planeTag), 1);
    info = struct('method', 'cutoff_mean_arrival', 'usedGPU', false, 'globalShiftBins', 0, ...
        'minPhotons', 100, 'nPixelsFit', 0, 'tauBasis', []);

    keep = micro > cutoffBin;
    if ~any(keep)
        return;
    end

    counts = accumarray(localPix(keep), 1, [numel(planeTag), 1], @sum, 0);
    weighted = accumarray(localPix(keep), micro(keep) * dtNs, [numel(planeTag), 1], @sum, 0);
    mask = counts > 0;
    planeLife(mask) = weighted(mask) ./ counts(mask) - cutoffBin * dtNs;
    info.nPixelsFit = nnz(mask);
end

function showLifetimePreview(data, tag, life_imm, head)
    pixel = head.ImgHdr.PixelSize;
    x0 = head.ImgHdr.X0;
    y0 = head.ImgHdr.Y0;

    if size(tag, 3) > 1
        for ch = 1:size(tag, 3)
            figure;
            subplot(1, 2, 1);
            imagesc(x0 + pixel .* (1:size(tag, 2)), y0 + pixel .* (1:size(tag, 1)), squeeze(tag(:, :, ch)));
            axis equal tight;
            title(sprintf('Intensity detector %d', ch));
            colorbar;

            subplot(1, 2, 2);
            imagesc(x0 + pixel .* (1:size(life_imm, 2)), y0 + pixel .* (1:size(life_imm, 1)), squeeze(life_imm(:, :, ch)));
            axis equal tight;
            title(sprintf('Lifetime detector %d', ch));
            colorbar;
        end
    else
        figure;
        subplot(1, 2, 1);
        imagesc(x0 + pixel .* (1:size(tag, 2)), y0 + pixel .* (1:size(tag, 1)), tag);
        axis equal tight;
        title('Intensity image');
        colorbar;

        subplot(1, 2, 2);
        imagesc(x0 + pixel .* (1:size(life_imm, 2)), y0 + pixel .* (1:size(life_imm, 1)), life_imm);
        axis equal tight;
        title(sprintf('Lifetime image (%s)', data.lifetimeEstimator), 'Interpreter', 'none');
        colorbar;
    end
end

function tf = localGpuAvailable()
    tf = false;
    try
        tf = exist('gpuDeviceCount', 'file') == 2 && gpuDeviceCount > 0;
    catch
        tf = false;
    end
end
