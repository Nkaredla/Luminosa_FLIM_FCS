function out = GlobalMultiExpPatternMatchFromTCSPC(tcspc_pix, irf, pulsePeriodNs, dtNs, tau0, opts)
% GLOBALMULTIEXPPATTERNMATCHFROMTCSPC
% General global multi-exponential FLIM/TCSPC pipeline:
%   1) fit the whole-image decay with N exponentials (N = numel(tau0))
%   2) build a basis / pattern matrix from the fitted decay components
%   3) pattern-match each pixel with nonnegative amplitudes
%
% This version is adapted to J. Enderlein style Fluofit with signature:
%   [c, offset, A, tau, dc, dtau, irs, zz, t, chi] = Fluofit(irf, y, p, dt, tau, lim, init)
%
% REQUIRED EXTERNALS
%   Fluofit.m
%   PatternMatchIm_matlab.m
%
% INPUTS
%   tcspc_pix      : [nx x ny x t] or [nx x ny x t x nPIE]
%   irf            : [t x 1] IRF
%   pulsePeriodNs  : repetition / PIE period in ns
%   dtNs           : TCSPC bin width in ns
%   tau0           : initial lifetime guesses, e.g. [0.4 2.0] or [0.3 1.5 4.2]
%   opts           : optional struct with fields
%       .useGPU            : true/false (default auto)
%       .batchSize         : default 4000
%       .mode              : 'PIRLS' | 'Nonneg' | 'Default' (default 'PIRLS')
%       .includeBackground : true/false (default true)
%       .backgroundPattern : [] or explicit [t x 1] background basis
%       .pieIndex          : for 4D tcspc_pix, default 1
%       .normalizePatterns : true/false (default true)
%       .sortLifetimes     : true/false (default true)
%       .fluoFitLim        : [] or [min_tau ... max_tau] in ns
%       .fluoFitInit       : passed to Fluofit init argument (default 0)
%
% OUTPUT
%   out struct containing global fit, basis, amplitudes, fractions, mean lifetimes,
%   and reconstruction.

    if nargin < 6 || isempty(opts), opts = struct(); end

    if ~isfield(opts, 'useGPU') || isempty(opts.useGPU)
        opts.useGPU = (gpuDeviceCount > 0);
    end
    if ~isfield(opts, 'batchSize') || isempty(opts.batchSize)
        opts.batchSize = 4000;
    end
    if ~isfield(opts, 'mode') || isempty(opts.mode)
        opts.mode = 'PIRLS';
    end
    if ~isfield(opts, 'includeBackground') || isempty(opts.includeBackground)
        opts.includeBackground = true;
    end
    if ~isfield(opts, 'backgroundPattern')
        opts.backgroundPattern = [];
    end
    if ~isfield(opts, 'pieIndex') || isempty(opts.pieIndex)
        opts.pieIndex = 1;
    end
    if ~isfield(opts, 'normalizePatterns') || isempty(opts.normalizePatterns)
        opts.normalizePatterns = true;
    end
    if ~isfield(opts, 'sortLifetimes') || isempty(opts.sortLifetimes)
        opts.sortLifetimes = true;
    end
    if ~isfield(opts, 'fluoFitLim')
        opts.fluoFitLim = [];
    end
    if ~isfield(opts, 'fluoFitInit') || isempty(opts.fluoFitInit)
        opts.fluoFitInit = 0; % use provided tau0 by default
    end

    % Select TCSPC cube.
    if ndims(tcspc_pix) == 4
        Ypix = single(tcspc_pix(:, :, :, opts.pieIndex));
    else
        Ypix = single(tcspc_pix);
    end

    [nx, ny, nt] = size(Ypix);
    globalDecay = squeeze(sum(sum(Ypix, 1), 2));
    globalDecay = double(globalDecay(:));
    irf = double(irf(:));
    tau0 = double(tau0(:)).';
    nExp = numel(tau0);

    if numel(irf) ~= nt
        error('IRF length (%d) must match time bins in tcspc_pix (%d).', numel(irf), nt);
    end
    if nExp < 1
        error('tau0 must contain at least one initial lifetime guess.');
    end

    % Correct call for Enderlein Fluofit.
    [cshift, offset, Aglob, taufit, dc, dtau, irfShifted, patterns0, tAxis, chi] = ...
        Fluofit(irf, globalDecay, pulsePeriodNs, dtNs, tau0, opts.fluoFitLim, opts.fluoFitInit);

    taufit = double(taufit(:));
    Aglob = double(Aglob(:));
    patterns0 = double(patterns0);
    tAxis = double(tAxis(:));

    % Fluofit returns zz with columns: [background, exp1, exp2, ...].
    if size(patterns0, 2) < nExp + 1
        error(['Fluofit returned zz with %d columns, but expected at least %d ' ...
               '(background + %d exponentials).'], size(patterns0, 2), nExp + 1, nExp);
    end

    bgPatternFromFit = patterns0(:, 1);
    expPatterns = patterns0(:, 2:(nExp + 1));

    if opts.sortLifetimes
        [taufit, ord] = sort(taufit, 'ascend');
        expPatterns = expPatterns(:, ord);
        if numel(Aglob) >= nExp
            Aglob = Aglob(ord);
        end
    else
        ord = 1:nExp;
    end

    % Build basis matrix.
    M = expPatterns;
    basisNames = arrayfun(@(k) sprintf('exp_%d', k), 1:nExp, 'UniformOutput', false);

    if opts.includeBackground
        if isempty(opts.backgroundPattern)
            bg = bgPatternFromFit;
            if ~any(bg > 0)
                bg = ones(nt, 1);
            end
        else
            bg = double(opts.backgroundPattern(:));
            if numel(bg) ~= nt
                error('backgroundPattern must have length %d.', nt);
            end
        end
        M = [bg, M];
        basisNames = [{'background'}, basisNames];
    end

    if opts.normalizePatterns
        colsum = sum(M, 1);
        colsum(colsum <= 0) = 1;
        M = M ./ colsum;
    end

    % Pixelwise pattern matching.
    [Amp, Z, pmInfo] = PatternMatchIm_matlab(Ypix, M, opts.mode, opts.useGPU, opts.batchSize);

    % Species fractions excluding the background channel.
    if opts.includeBackground
        speciesIdx = 2:size(Amp, 3);
    else
        speciesIdx = 1:size(Amp, 3);
    end

    speciesSum = sum(Amp(:, :, speciesIdx), 3);
    speciesSumSafe = max(speciesSum, eps('single'));

    AmpFrac = zeros(size(Amp), 'single');
    AmpFrac(:, :, speciesIdx) = Amp(:, :, speciesIdx) ./ speciesSumSafe;

    tauMean = zeros(nx, ny, 'single');
    tauRate = zeros(nx, ny, 'single');
    for k = 1:nExp
        idx = speciesIdx(k);
        tauMean = tauMean + AmpFrac(:, :, idx) * single(taufit(k));
        tauRate = tauRate + AmpFrac(:, :, idx) * single(1 / taufit(k));
    end
    tauHarmonic = 1 ./ max(tauRate, eps('single'));

    intensity = sum(Ypix, 3);

    out = struct();
    out.globalDecay = globalDecay;
    out.irf = irf;
    out.pulsePeriodNs = pulsePeriodNs;
    out.dtNs = dtNs;
    out.tau0 = tau0;
    out.nExp = nExp;

    out.taufit = taufit;
    out.lifetimeOrder = ord;
    out.Aglob = Aglob;
    out.cshift = cshift;
    out.dc = dc;
    out.dtau = dtau;
    out.zfit = sum(patterns0, 2);
    out.patterns = M;
    out.patternsRaw = patterns0;
    out.offset = offset;
    out.irfShifted = irfShifted;
    out.tAxis = tAxis;
    out.chi = chi;
    out.basisNames = basisNames;

    out.Amp = Amp;
    out.AmpFrac = AmpFrac;
    out.recon = Z;
    out.tauMeanArithmetic = tauMean;
    out.tauMeanHarmonic = tauHarmonic;
    out.intensity = intensity;
    out.patternMatchInfo = pmInfo;
    out.speciesIdx = speciesIdx;
end
