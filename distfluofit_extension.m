function out = distfluofit_extension(tcspc_pix, irf, pulsePeriodNs, dtNs, opts)
% DISTFLUOFIT_EXTENSION
% DistFluofit-style distributed lifetime basis followed by per-pixel PIRLS.
%
% INPUTS
%   tcspc_pix      : [nx x ny x t] or [nx x ny x t x nCh] TCSPC cube
%   irf            : [t x 1] IRF on the same grid
%   pulsePeriodNs  : repetition period in ns
%   dtNs           : TCSPC bin width in ns
%   opts           : optional struct with fields
%       .useGPU      : true/false, default auto
%       .batchSize   : pixel chunk size, default 1024
%       .mode        : 'PIRLS' | 'Nonneg' | 'Default', default 'PIRLS'
%       .nTau        : DistFluofit N parameter, default 100
%       .shiftBounds : [min max] shift bounds in bins, default [-10 10]
%
% OUTPUT
%   out struct with global DistFluofit fit, basis, coefficients,
%   per-pixel mean lifetime, and reconstruction.

    if nargin < 5 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'useGPU') || isempty(opts.useGPU)
        opts.useGPU = (gpuDeviceCount > 0);
    end
    if ~isfield(opts, 'batchSize') || isempty(opts.batchSize)
        opts.batchSize = 1024;
    end
    if ~isfield(opts, 'mode') || isempty(opts.mode)
        opts.mode = 'PIRLS';
    end
    if ~isfield(opts, 'nTau') || isempty(opts.nTau)
        opts.nTau = 100;
    end
    if ~isfield(opts, 'shiftBounds') || isempty(opts.shiftBounds)
        opts.shiftBounds = [-10 10];
    end

    Ypix = single(tcspc_pix);
    if ndims(Ypix) == 4
        Ypix = sum(Ypix, 4);
    end
    if ndims(Ypix) ~= 3
        error('tcspc_pix must be a 3D or 4D TCSPC cube.');
    end

    [nx, ny, nt] = size(Ypix);
    irf = max(double(irf(:)), 0);
    if numel(irf) ~= nt
        error('IRF length (%d) must match the TCSPC time axis (%d).', numel(irf), nt);
    end
    if any(irf > 0)
        irf = irf ./ max(sum(irf), eps);
    else
        irf = zeros(nt, 1);
        irf(1) = 1;
    end

    % Use the sampled TCSPC grid as the modeled period to avoid one-bin
    % mismatches between PTU metadata and the stored histogram cube.
    modelPeriodNs = nt * dtNs;

    globalDecay = squeeze(sum(sum(Ypix, 1), 2));
    globalDecay = double(globalDecay(:));

    [cxGlob, tauBasisInv, offset, cshift, globalRecon, tAxis, err, Mfit, Mfull, indFit, globalFitInfo, shiftScan] = ...
        solve_global_distfluofit(irf, globalDecay, pulsePeriodNs, modelPeriodNs, dtNs, opts);

    [Amp, Z, pmInfo] = fit_pixels_distfluofit(Ypix, Mfit, Mfull, indFit, opts);

    speciesIdx = 2:size(Amp, 3);
    speciesSum = sum(Amp(:, :, speciesIdx), 3);
    speciesSumSafe = max(speciesSum, eps('single'));

    AmpFrac = zeros(size(Amp), 'single');
    if ~isempty(speciesIdx)
        AmpFrac(:, :, speciesIdx) = Amp(:, :, speciesIdx) ./ speciesSumSafe;
    end

    tauBasisInv = double(tauBasisInv(:));
    tauBasisNs = 1 ./ max(tauBasisInv, eps);
    tauMean = nan(nx, ny, 'single');
    if ~isempty(speciesIdx)
        tauMean = zeros(nx, ny, 'single');
        for k = 1:numel(speciesIdx)
            tauMean = tauMean + AmpFrac(:, :, speciesIdx(k)) * single(tauBasisNs(k));
        end
        tauMean(speciesSum <= 0) = NaN;
    end

    out = struct();
    out.globalDecay = globalDecay;
    out.irf = irf;
    out.pulsePeriodNs = pulsePeriodNs;
    out.modelPeriodNs = modelPeriodNs;
    out.dtNs = dtNs;
    out.shiftBounds = opts.shiftBounds;
    out.nTau = opts.nTau;
    out.mode = opts.mode;
    out.batchSize = opts.batchSize;

    out.cshift = cshift;
    out.offset = offset;
    out.err = err;
    out.tAxis = double(tAxis(:));
    out.tAxisNs = dtNs * double(tAxis(:));
    out.globalRecon = globalRecon(:);
    out.globalDistribution = cxGlob(:);
    out.tauBasis = tauBasisNs(:);
    out.tauBasisInv = tauBasisInv(:);
    out.patterns = Mfull;
    out.patternsFit = Mfit;
    out.fitIndices = indFit(:);

    out.Amp = Amp;
    out.AmpFrac = AmpFrac;
    out.recon = Z;
    out.tauMeanArithmetic = tauMean;
    out.intensity = sum(Ypix, 3);
    out.patternMatchInfo = pmInfo;
    out.globalFitInfo = globalFitInfo;
    out.shiftScan = shiftScan;
    out.speciesIdx = speciesIdx;
end

function [cxGlob, tauBasisInv, offset, cshift, globalRecon, tAxis, err, Mfit, Mfull, indFit, fitInfo, shiftScan] = ...
        solve_global_distfluofit(irf, globalDecay, pulsePeriodNs, modelPeriodNs, dtNs, opts)
    shiftMin = floor(opts.shiftBounds(1));
    shiftMax = ceil(opts.shiftBounds(end));
    shiftInt = shiftMin:shiftMax;
    if isempty(shiftInt)
        shiftInt = round(mean(opts.shiftBounds));
    end

    coarseErr = nan(size(shiftInt));
    for idx = 1:numel(shiftInt)
        [MfitTmp, MfullTmp, indFitTmp, ~] = build_distfluofit_model( ...
            irf, pulsePeriodNs, modelPeriodNs, dtNs, shiftInt(idx), opts.nTau);
        [betaTmp, ~] = solve_distfluofit_batch(MfitTmp, globalDecay(indFitTmp), opts.mode, opts.useGPU, opts.batchSize);
        zTmp = MfullTmp * betaTmp;
        coarseErr(idx) = sum((zTmp(indFitTmp) - globalDecay(indFitTmp)).^2 ./ max(abs(zTmp(indFitTmp)), eps)) ./ ...
            max(numel(indFitTmp), 1);
    end

    shiftFine = shiftInt;
    errFine = coarseErr;
    if numel(shiftInt) > 1 && all(isfinite(coarseErr))
        tmpShift = shiftInt(1):0.1:shiftInt(end);
        tmpErr = interp1(shiftInt, coarseErr, tmpShift, 'linear');
        if any(isfinite(tmpErr))
            shiftFine = tmpShift;
            errFine = tmpErr;
        end
    end

    validErr = isfinite(errFine);
    if ~any(validErr)
        cshift = 0;
    else
        shiftCandidates = shiftFine(validErr);
        errCandidates = errFine(validErr);
        [~, bestIdx] = min(errCandidates);
        cshift = shiftCandidates(bestIdx);
    end

    [Mfit, Mfull, indFit, tauBasisInv] = build_distfluofit_model( ...
        irf, pulsePeriodNs, modelPeriodNs, dtNs, cshift, opts.nTau);
    [betaGlob, fitInfo] = solve_distfluofit_batch(Mfit, globalDecay(indFit), opts.mode, opts.useGPU, opts.batchSize);

    globalRecon = Mfull * betaGlob;
    err = sum((globalRecon(indFit) - globalDecay(indFit)).^2 ./ max(abs(globalRecon(indFit)), eps)) ./ ...
        max(numel(indFit), 1);
    offset = double(betaGlob(1));
    cxGlob = double(betaGlob(2:end));
    tAxis = (1:numel(irf))';

    shiftScan = struct();
    shiftScan.coarseShiftBins = shiftInt(:);
    shiftScan.coarseErr = coarseErr(:);
    shiftScan.fineShiftBins = shiftFine(:);
    shiftScan.fineErr = errFine(:);
end

function [Amp, Z, info] = fit_pixels_distfluofit(Ypix, Mfit, Mfull, indFit, opts)
    [nx, ny, nt] = size(Ypix);
    Y = reshape(permute(Ypix(:, :, indFit), [3 1 2]), numel(indFit), []);
    [Cflat, info] = solve_distfluofit_batch(Mfit, Y, opts.mode, opts.useGPU, opts.batchSize);
    Zflat = Mfull * Cflat;

    Amp = reshape(Cflat.', nx, ny, size(Mfull, 2));
    Z = ipermute(reshape(Zflat, nt, nx, ny), [3 1 2]);
end

function [Beta, info] = solve_distfluofit_batch(M, Y, mode, useGPU, batchSize)
    if nargin < 3 || isempty(mode)
        mode = 'PIRLS';
    end
    if nargin < 4 || isempty(useGPU)
        useGPU = (gpuDeviceCount > 0);
    end
    if nargin < 5 || isempty(batchSize)
        batchSize = 1024;
    end

    if isvector(Y)
        Y = double(Y(:));
    else
        Y = double(Y);
    end
    M = double(M);

    switch mode
        case 'Default'
            Beta = single(M \ Y);
            info = struct('method', 'least squares', 'usedGPU', false);

        case 'Nonneg'
            nPix = size(Y, 2);
            Beta = zeros(size(M, 2), nPix, 'single');
            for idx = 1:nPix
                Beta(:, idx) = single(lsqnonneg(M, Y(:, idx)));
            end
            info = struct('method', 'pixelwise lsqnonneg', 'usedGPU', false);

        case 'PIRLS'
            [Beta, info] = PIRLSnonneg_batch_gpu_matlab(M, Y, 10, 25, batchSize, useGPU);

        otherwise
            error('Unknown DistFluofit solver mode: %s', mode);
    end
end

function [Mfit, Mfull, indFit, tauBasisInv] = build_distfluofit_model(irf, pulsePeriodNs, modelPeriodNs, dtNs, shiftBins, nBasis)
    n = numel(irf);
    periodBins = max(round(modelPeriodNs / dtNs), 1);
    if periodBins ~= n
        periodBins = n;
    end
    tauMaxBins = max(pulsePeriodNs / dtNs, 1 + eps);
    tp = dtNs * (1:periodBins)';
    t = (1:n)';

    tauBasisInv = (1 / dtNs) ./ exp((0:nBasis) / nBasis * log(tauMaxBins));
    M0 = [ones(size(t)), Convol(irf, exp(-tp * tauBasisInv))];
    colSum = sum(M0, 1);
    colSum(colSum <= 0) = 1;
    M0 = M0 ./ (ones(n, 1) * colSum);

    shFloor = floor(shiftBins);
    shCeil = ceil(shiftBins);
    Mfull = (1 - shiftBins + shFloor) * M0(rem(rem(t - shFloor - 1, n) + n, n) + 1, :) + ...
        (shiftBins - shFloor) * M0(rem(rem(t - shCeil - 1, n) + n, n) + 1, :);

    c = ceil(abs(shiftBins)) * sign(shiftBins);
    indFit = max(1, 1 + c):min(n, n + c);
    Mfit = Mfull(indFit, :);
end
