function [Beta, info] = PIRLSnonneg_batch_gpu_matlab(M, Y, maxPirlsIter, maxNnlsIter, batchSize, useGPU)
% PIRLSNONNEG_BATCH_GPU_MATLAB
% GPU-friendly batch PIRLS for non-negative coefficients.
%
% This is a MATLAB/gpuArray implementation aimed at the same workflow as the
% Python PIRLS_nonneg_kernel, but the inner NNLS step is solved with a
% projected-gradient update that is much easier to vectorize on the GPU.
%
% INPUTS
%   M            : [nSamples x nBasis] design / pattern matrix
%   Y            : [nSamples x nPixels] data matrix, one column per pixel
%   maxPirlsIter : number of PIRLS outer iterations (default 10)
%   maxNnlsIter  : number of projected-gradient inner iterations (default 25)
%   batchSize    : number of pixels processed per GPU batch (default 4000)
%   useGPU       : true/false, auto-detect if omitted
%
% OUTPUTS
%   Beta         : [nBasis x nPixels] non-negative coefficients
%   info         : struct with diagnostics
%
% NOTES
%   1) This is not a bitwise clone of the PyCUDA code. The PyCUDA kernel uses
%      per-pixel coordinate-descent NNLS. Here the inner NNLS is solved with a
%      projected-gradient pagewise update so it can run efficiently with
%      gpuArray + pagemtimes.
%   2) For 3-4 FLIM basis functions this usually behaves well.
%   3) Pixel chunking is essential; otherwise pagewise weighted systems can
%      consume too much GPU memory.
%
% Narain's FLIM workflow: use together with PatternMatchIm_matlab.m

    if nargin < 3 || isempty(maxPirlsIter), maxPirlsIter = 10; end
    if nargin < 4 || isempty(maxNnlsIter),  maxNnlsIter  = 25; end
    if nargin < 5 || isempty(batchSize),    batchSize    = 4000; end
    if nargin < 6 || isempty(useGPU)
        useGPU = (gpuDeviceCount > 0);
    end

    [nSamples, nBasis] = size(M);
    [nSamplesY, nPixels] = size(Y);
    assert(nSamplesY == nSamples, 'M and Y must have the same number of rows.');

    % Keep GPU math in single precision unless you have a strong reason not to.
    M = single(M);
    Y = single(Y);

    Beta = zeros(nBasis, nPixels, 'single');
    info = struct();
    info.usedGPU = false;
    info.method  = 'CPU exact fallback';
    info.maxPirlsIter = maxPirlsIter;
    info.maxNnlsIter  = maxNnlsIter;
    info.batchSize    = batchSize;
    info.requestedBatchSize = batchSize;
    info.gpuBatchSize = [];
    info.gpuFailureMessage = '';
    info.gpuProcessedPixels = 0;

    if useGPU
        try
            g = gpuDevice;
            Mg = gpuArray(M);
            MtM = Mg' * Mg;
            MG = reshape(Mg, nSamples, nBasis, 1);
            Ikg = eye(nBasis, 'single', 'gpuArray');
            tiny = gpuArray(single(max(1e-6, 0.1 / max(1, nSamples))));
            gpuBatchSize = estimateSafeGpuBatchSize(g, nSamples, nBasis, batchSize, nPixels);
            info.gpuBatchSize = gpuBatchSize;

            i0 = 1;
            while i0 <= nPixels
                batchUsed = min(gpuBatchSize, nPixels - i0 + 1);
                batchDone = false;
                while ~batchDone
                    idx = i0:min(i0 + batchUsed - 1, nPixels);
                    try
                        Beta_g = runGpuPirlsBatch(Mg, MG, MtM, Ikg, Y(:, idx), maxPirlsIter, maxNnlsIter, tiny);
                        Beta(:, idx) = gather(Beta_g);
                        info.gpuProcessedPixels = idx(end);
                        i0 = idx(end) + 1;
                        batchDone = true;
                    catch MEbatch
                        if isGpuOutOfMemoryError(MEbatch) && batchUsed > 1
                            batchUsed = max(1, floor(batchUsed / 2));
                            gpuBatchSize = min(gpuBatchSize, batchUsed);
                            info.gpuBatchSize = gpuBatchSize;
                        else
                            rethrow(MEbatch);
                        end
                    end
                end
            end

            info.usedGPU = true;
            info.method  = 'GPU PIRLS + projected-gradient NNLS';
            return;

        catch ME
            info.usedGPU = info.gpuProcessedPixels > 0;
            info.gpuFailureMessage = ME.message;
            if info.usedGPU
                info.method = 'Hybrid GPU PIRLS + CPU exact fallback';
                warning('PIRLSnonneg_batch_gpu_matlab:GPUFallback', ...
                    'GPU PIRLS path failed after %d pixel(s), switching the remainder to CPU exact fallback.\n%s', ...
                    info.gpuProcessedPixels, ME.message);
            else
                info.method = 'CPU exact fallback';
                warning('PIRLSnonneg_batch_gpu_matlab:GPUFallback', ...
                    'GPU PIRLS path failed, switching to CPU exact fallback.\n%s', ME.message);
            end
        end
    end

    % ------------------------------------------------------------------
    % CPU fallback: closer to your original Python CPU routine
    % ------------------------------------------------------------------
    Md = double(M);
    Yd = double(Y);
    tiny = max(1e-6, 0.1 / max(1, nSamples));

    cpuStart = max(1, info.gpuProcessedPixels + 1);
    for i = cpuStart:nPixels
        yi = Yd(:, i);
        beta = lsqnonneg(Md, yi);

        for it = 1:maxPirlsIter
            w = 1 ./ max(Md * beta, tiny);
            Mw = Md .* w;
            Aw = Md' * Mw;
            bw = Md' * (w .* yi);
            betaNew = lsqnonneg(Aw, bw);

            if norm(betaNew - beta) < 1e-10
                beta = betaNew;
                break;
            end
            beta = betaNew;
        end

        Beta(:, i) = single(beta);
    end
end

function Beta_g = runGpuPirlsBatch(Mg, MG, MtM, Ikg, Ybatch, maxPirlsIter, maxNnlsIter, tiny)
    [nSamples, nBasis] = size(Mg);
    Yg = gpuArray(Ybatch);
    B = size(Ybatch, 2);

    % Unweighted initialization.
    Mty = Mg' * Yg;
    Beta_g = max((MtM + 1e-6 * Ikg) \ Mty, 0);

    for it = 1:maxPirlsIter
        Mu = max(Mg * Beta_g, tiny);
        W  = 1 ./ Mu;                               % [nSamples x B]

        % Weighted normal equations for each pixel/page:
        % Aw(:,:,p) = M' * diag(w(:,p)) * M
        % bw(:,:,p) = M' * diag(w(:,p)) * y(:,p)
        W3 = reshape(W, nSamples, 1, B);
        MW = MG .* W3;
        Aw = pagemtimes(permute(MG, [2 1 3]), MW);  % [nBasis x nBasis x B]

        YW = reshape(W .* Yg, nSamples, 1, B);
        bw = pagemtimes(permute(MG, [2 1 3]), YW);  % [nBasis x 1 x B]

        beta3 = reshape(Beta_g, nBasis, 1, B);

        % Pagewise projected-gradient NNLS.
        % L is a safe Lipschitz bound from the max absolute row sum.
        L = max(sum(abs(Aw), 2), 1e-6);             % [nBasis x 1 x B]
        L = max(L, [], 1);                          % [1 x 1 x B]

        for k = 1:maxNnlsIter
            grad  = pagemtimes(Aw, beta3) - bw;
            beta3 = max(beta3 - grad ./ L, 0);
        end

        Beta_g = reshape(beta3, nBasis, B);
    end
end

function batchSizeOut = estimateSafeGpuBatchSize(g, nSamples, nBasis, requestedBatchSize, nPixels)
    batchSizeOut = max(1, min(round(double(requestedBatchSize)), round(double(nPixels))));
    availBytes = NaN;
    try
        availBytes = double(g.AvailableMemory);
    catch
    end
    if ~isfinite(availBytes) || availBytes <= 0
        return;
    end

    bytesPerSingle = 4;
    fixedBytes = bytesPerSingle * max(1, 5 * nSamples * nBasis + nBasis * nBasis);
    perPixelBytes = bytesPerSingle * max(1, ...
        (nSamples * nBasis) + (nBasis * nBasis) + (6 * nSamples) + (6 * nBasis));
    usableBudget = 0.35 * availBytes - fixedBytes;
    if usableBudget <= perPixelBytes
        batchSizeOut = 1;
        return;
    end

    batchSizeOut = min(batchSizeOut, max(1, floor(usableBudget / perPixelBytes)));
end

function tf = isGpuOutOfMemoryError(ME)
    id = '';
    try
        id = lower(char(ME.identifier));
    catch
    end
    msg = '';
    try
        msg = lower(char(ME.message));
    catch
    end
    tf = contains(id, 'outofmemory') || contains(id, 'memory') || ...
         contains(msg, 'out of memory') || contains(msg, 'insufficient memory') || ...
         contains(msg, 'cuda_error_out_of_memory') || contains(msg, 'device memory');
end
