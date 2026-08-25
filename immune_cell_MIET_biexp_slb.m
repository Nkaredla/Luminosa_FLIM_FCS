function out = immune_cell_MIET_biexp_slb(source, opts)
%IMMUNE_CELL_MIET_BIEXP_SLB Per-pixel biexponential with a soft-fixed SLB.
%
% out = immune_cell_MIET_biexp_slb(analysisMatOrFolder)
% out = immune_cell_MIET_biexp_slb(source, opts)
%
% Fits every selected pixel with
%
%     y(t) = B + a1 * P(t; tau1) + a2 * P(t; tau2)
%
% tau1 soft-constrained to the SLB reference, tau2 free, and B, a1, a2 all
% non-negative from PIRLSnonneg. Writes maps, a figure and a summary into a NEW
% folder beside the source analysis, leaving the three-model results untouched.
%
% WHY THIS REPLACES THE THREE-MODEL FIT
%
% Two components cannot be told from three in this data. From the per-pixel
% explorer: the third component takes 0-10% amplitude and leaves the maximum
% residual UNCHANGED, at 1x1, at 5x5 with 45k photons, and at 15x15 with 214k
% photons. From simulation earlier in the project: the false three-component
% rate RISES with photon count, 9% to 54%, because the likelihood cannot
% separate a real third exponential from two approximating an off-grid one.
% Fitting three and reporting the third is therefore reporting noise.
%
% HOW IT IS MADE FAST ENOUGH
%
% Around 100k pixels need fitting, and a two-dimensional nonlinear search per
% pixel is far too slow. Two stages instead:
%
%   1. GRID RANKING, vectorised over all pixels at once. For each (tau1, tau2)
%      pair on a grid the design matrix is the SAME for every pixel, so the
%      linear solve is one small matrix product against the whole pixel block.
%      Ranking uses chi-square rather than the Poisson deviance purely for
%      speed - it only has to order the pairs, and it avoids 3e9 logarithms.
%      Negative amplitudes are clamped to zero before the residual so the
%      ranking always reflects a feasible model.
%
%   2. PIRLSnonneg REFINEMENT, once per pixel at its best pair. This is where
%      the reported amplitudes and background come from, so the numbers that
%      matter are Poisson-weighted and properly non-negative. Plain least
%      squares would let the bright bins near the peak dominate and leave the
%      background barely constrained - the failure that earlier produced a
%      fitted background twice the measured pedestal.
%
% opts fields
%   slbTauNs        prior centre; default from the source analysis
%   slbSigmaNs      prior width, ns (default 0.05)
%   tau2BoundsNs    bounds for the free component (default [0.5 8])
%   tau2GridCount   grid points for tau2 (default 40, log-spaced)
%   tau1GridCount   grid points for tau1 across +/-2 sigma (default 5)
%   binSize         k for k-by-k spatial summing (default 1). The output folder
%                   name carries it, so different binnings never overwrite.
%   pixelMask       'cellFootprint' (default), 'valid', 'all', or a logical map
%   minPhotons      skip pixels below this (default 200)
%   blockSize       pixels per vectorised block (default 20000)
%   outputDir       default <acquisition>/biexp_slb[_binK]
%   makeFigure      default true

    if nargin < 2 || isempty(opts); opts = struct(); end
    defaults = struct('slbTauNs', [], 'slbSigmaNs', 0.05, ...
        'tau2BoundsNs', [0.5 8], 'tau2GridCount', 40, 'tau1GridCount', 5, ...
        'binSize', 1, 'pixelMask', 'cellFootprint', 'minPhotons', 200, ...
        'blockSize', 20000, 'outputDir', '', 'makeFigure', true);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    for required = {'PIRLSnonneg', 'biexp_slb_pattern', 'biexp_slb_deviance'}
        if exist(required{1}, 'file') ~= 2
            error('immune_cell_MIET_biexp_slb:Missing', ...
                '%s.m must be on the MATLAB path.', required{1});
        end
    end

    analysisMat = immune_cell_MIET_biexp_resolve(source);
    fprintf('\nimmune_cell_MIET_biexp_slb\n  %s\n', analysisMat);
    stored = whos('-file', analysisMat);
    names = {stored.name};
    if ~ismember('result', names)
        error('immune_cell_MIET_biexp_slb:NoResult', ...
            'No "result" variable in %s', analysisMat);
    end
    cubeVariable = '';
    for candidate = {'tcspc_pix', 'tcspcPix'}
        if ismember(candidate{1}, names); cubeVariable = candidate{1}; break; end
    end
    if isempty(cubeVariable)
        error('immune_cell_MIET_biexp_slb:NoCube', ...
            ['No per-pixel TCSPC cube in this MAT, so nothing can be ' ...
             'refitted. Re-run the pipeline with cfg.saveTcspcPix = true.']);
    end

    loaded = load(analysisMat, 'result');
    result = loaded.result;
    dtNs = immune_cell_MIET_explorer_field(result, 'bayesian.compact.dtNs');
    periodNs = immune_cell_MIET_explorer_field(result, ...
        'bayesian.compact.pulsePeriodNs');
    irf = immune_cell_MIET_explorer_field(result, 'irf.curve');
    if isempty(dtNs) || isempty(periodNs) || isempty(irf)
        error('immune_cell_MIET_biexp_slb:NoTiming', ...
            'The result lacks dtNs, pulsePeriodNs or the IRF.');
    end
    dtNs = double(dtNs); periodNs = double(periodNs);
    irf = double(irf(:)); irf = max(irf, 0);
    if sum(irf) > 0; irf = irf / sum(irf); end

    if isempty(opts.slbTauNs)
        candidate = immune_cell_MIET_explorer_field(result, ...
            'bayesian.compact.fixedSlbLifetimeNs');
        if isempty(candidate)
            candidate = immune_cell_MIET_explorer_field(result, ...
                'slbReference.fixedLifetimeNs');
        end
        if isempty(candidate)
            error('immune_cell_MIET_biexp_slb:NoSlbReference', ...
                ['The result has no SLB reference lifetime; pass ' ...
                 'opts.slbTauNs explicitly.']);
        end
        opts.slbTauNs = double(candidate(1));
    end
    slbChi = immune_cell_MIET_explorer_field(result, ...
        'slbReference.reducedChiSquare');
    fprintf(['  dt %.4f ns, period %.2f ns, %d bins; SLB prior %.4f ' ...
        '+/- %.3f ns\n'], dtNs, periodNs, numel(irf), opts.slbTauNs, ...
        opts.slbSigmaNs);
    if ~isempty(slbChi) && isfinite(slbChi(1))
        fprintf(['      note: the SLB reference fit itself has reduced ' ...
            'chi-square %.4g, which is why\n      the constraint is SOFT ' ...
            'rather than fixed\n'], slbChi(1));
    end

    % ---- pixels ---------------------------------------------------------
    cube = load(analysisMat, cubeVariable);
    cube = cube.(cubeVariable);
    [nRow, nCol, nBin] = size(cube);
    mask = immune_cell_MIET_biexp_mask(result, opts.pixelMask, nRow, nCol);
    if opts.binSize > 1
        kernel = ones(opts.binSize);
        cube = convn(double(cube), kernel, 'same');
        mask = mask & (convn(double(mask), kernel, 'same') >= ...
            opts.binSize ^ 2 - 0.5);
        fprintf('  spatial binning %dx%d applied to the cube\n', ...
            opts.binSize, opts.binSize);
    end
    intensity = sum(double(cube), 3);
    mask = mask & intensity >= opts.minPhotons;
    pixelIndex = find(mask);
    fprintf(['  %d pixel(s) selected (%s, >= %d photons), median %.0f ' ...
        'photons\n'], numel(pixelIndex), maskLabel(opts.pixelMask), ...
        opts.minPhotons, median(intensity(pixelIndex)));
    if isempty(pixelIndex)
        error('immune_cell_MIET_biexp_slb:NoPixels', ...
            'No pixel met the mask and the %d-photon floor.', opts.minPhotons);
    end

    % ---- lifetime grid --------------------------------------------------
    tau1Grid = opts.slbTauNs + opts.slbSigmaNs * ...
        linspace(-2, 2, max(1, opts.tau1GridCount));
    tau1Grid = tau1Grid(tau1Grid > 0.02);
    tau2Grid = logspace(log10(opts.tau2BoundsNs(1)), ...
        log10(opts.tau2BoundsNs(2)), max(2, opts.tau2GridCount));
    fprintf('  grid: %d tau1 x %d tau2 = %d pairs\n', numel(tau1Grid), ...
        numel(tau2Grid), numel(tau1Grid) * numel(tau2Grid));

    out = immune_cell_MIET_biexp_run(cube, mask, pixelIndex, intensity, ...
        irf, dtNs, periodNs, nBin, tau1Grid, tau2Grid, opts);
    out.analysisMat = analysisMat;
    out.imageSize = [nRow nCol];
    out.opts = opts;
    out.slbReferenceReducedChiSquare = slbChi;

    % ---- save -----------------------------------------------------------
    if isempty(opts.outputDir)
        folder = fileparts(analysisMat);
        name = 'biexp_slb';
        if opts.binSize > 1
            name = sprintf('biexp_slb_bin%d', opts.binSize);
        end
        opts.outputDir = fullfile(folder, name);
    end
    if ~isfolder(opts.outputDir); mkdir(opts.outputDir); end
    out.opts.outputDir = opts.outputDir;
    matFile = fullfile(opts.outputDir, 'biexp_slb_maps.mat');
    save(matFile, 'out', '-v7.3');
    fprintf('\n  wrote %s\n', matFile);
    immune_cell_MIET_biexp_report(out, opts.outputDir);
    if opts.makeFigure
        out.figure = immune_cell_MIET_biexp_figure(out, opts.outputDir);
        fprintf('  wrote %s\n', out.figure);
    end
end
