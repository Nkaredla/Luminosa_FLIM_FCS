function out = runBlocks(cubeFlat, idx, basis, tauGrid, blockSize, nBin)
%RUNBLOCKS Mono-exponential Poisson fit of a pixel list, in blocks.
%
% Split out as its own file because the SLB measurement needs the same
% block loop three times (native, binned, and any future stage) and this
% project keeps one function per file rather than using local subfunctions.

    n = numel(idx);
    tauNs = nan(1, n);
    background = nan(1, n);
    amplitude = nan(1, n);
    deviance = nan(1, n);
    photons = nan(1, n);
    nBlock = ceil(n / blockSize);
    for k = 1:nBlock
        lo = (k - 1) * blockSize + 1;
        hi = min(n, k * blockSize);
        Y = double(cubeFlat(idx(lo:hi), :))';
        f = fit_mono_poisson_batch(Y, basis, tauGrid, struct());
        tauNs(lo:hi) = f.tauNs;
        background(lo:hi) = f.background;
        amplitude(lo:hi) = f.amplitude;
        deviance(lo:hi) = f.deviance;
        photons(lo:hi) = f.photons;
        fprintf('      block %d/%d\n', k, nBlock);
    end
    out = struct('tauNs', tauNs, 'background', background, ...
        'amplitude', amplitude, 'deviance', deviance, ...
        'reducedDeviance', deviance / max(nBin - 3, 1), 'photons', photons);
end
