function [C, Z, info] = PatternMatchIm_matlab(y, M, mode, useGPU, batchSize)
% PATTERNMATCHIM_MATLAB
% MATLAB version of PatternMatchIm for FLIM cubes.
%
% INPUTS
%   y         : [nx x ny x t] TCSPC cube
%   M         : [t x nBasis] pattern matrix
%   mode      : 'Default' | 'Nonneg' | 'PIRLS'
%   useGPU    : logical, used for PIRLS mode
%   batchSize : pixel chunk size for PIRLS GPU mode
%
% OUTPUTS
%   C         : [nx x ny x nBasis] coefficient maps
%   Z         : [nx x ny x t] reconstructed cube
%   info      : diagnostics struct

    if nargin < 3 || isempty(mode), mode = 'PIRLS'; end
    if nargin < 4 || isempty(useGPU)
        useGPU = (gpuDeviceCount > 0);
    end
    if nargin < 5 || isempty(batchSize), batchSize = 4000; end

    [nx, ny, t] = size(y);
    [t2, nBasis] = size(M);
    assert(t2 == t, 'time-axis of y and M must match.');

    y = single(y);
    M = single(M);

    % Normalize columns exactly as in the Python workflow.
    M = M ./ max(sum(M, 1), eps('single'));

    % Reshape cube to [t x npix]. MATLAB is column-major, so permute first.
    Y = reshape(permute(y, [3 1 2]), t, []);

    switch mode
        case 'Default'
            Cflat = M \ Y;
            info = struct('method', 'least squares', 'usedGPU', false);

        case 'Nonneg'
            npix = size(Y, 2);
            Cflat = zeros(nBasis, npix, 'single');
            for i = 1:npix
                Cflat(:, i) = single(lsqnonneg(double(M), double(Y(:, i))));
            end
            info = struct('method', 'pixelwise lsqnonneg', 'usedGPU', false);

        case 'PIRLS'
            [Cflat, info] = PIRLSnonneg_batch_gpu_matlab(M, Y, 10, 25, batchSize, useGPU);

        otherwise
            error('Unknown mode: %s', mode);
    end

    Zflat = M * Cflat;

    C = reshape(Cflat.', nx, ny, nBasis);
    Z = ipermute(reshape(Zflat, t, nx, ny), [3 1 2]);
end
