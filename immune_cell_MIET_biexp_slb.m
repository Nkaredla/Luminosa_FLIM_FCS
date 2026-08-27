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
%   shortlist       how many (tau1,tau2) pairs per pixel survive the cheap
%                   pre-rank and get scored by the exact Poisson deviance
%                   (default 8). The cheap score is biased, so it is used only
%                   to shortlist; the reported lifetime always comes from the
%                   deviance. Raise this if the run warns that winners are
%                   landing at the shortlist edge.
%   refineTau2      parabolic sub-node refinement of tau2 (default true).
%                   Without it tau2 can only take grid values, which is why six
%                   acquisitions previously all reported 2.07237064230134.
%   outputDir       default <acquisition>/biexp_slb[_binK]
%   method          'vp' (default) for variable projection - BFGS over
%                   [log tau1, log tau2] with the amplitudes profiled out by
%                   Poisson maximum likelihood - or 'grid' for the older grid
%                   search. 'vp' returns CONTINUOUS lifetimes, reaches a
%                   strictly lower objective on 100% of pixels measured, and
%                   runs at 2.8 ms/pixel against the grid's 6.2.
%   gtol            outer gradient tolerance for 'vp' (default 1e-3)
%   innerSolver     'irls' (default) or 'em'. Whitened IRLS reaches the same
%                   amplitudes as EM (median deviance difference 0 over 3000
%                   real pixels, bound active on half of them) about 4x faster.
%   displayMask     mask applied to the FIGURES only, not to the fit (default
%                   'cellFootprint'). Everything is fitted so the bare SLB can
%                   check the anchor; only the cell is drawn, so the colour
%                   scales are set by the cell rather than flattened by a region
%                   whose answer is already known.
%   fixSlbTau       hold tau1 at slbTauNs and fit only its amplitude
%                   (default false)
%   makeFigure      default true

    if nargin < 2 || isempty(opts); opts = struct(); end
    defaults = struct('slbTauNs', [], 'slbSigmaNs', 0.05, ...
        'tau2BoundsNs', [0.5 8], 'tau2GridCount', 40, 'tau1GridCount', 5, ...
        'binSize', 1, 'pixelMask', 'cellFootprint', 'minPhotons', 200, ...
        'blockSize', 20000, 'shortlist', 8, 'refineTau2', true, ...
        'method', 'vp', 'gtol', 1e-3, 'tau2SeedNs', 2.0, ...
        'displayMask', 'cellFootprint', 'fixSlbTau', false, ...
        'innerSolver', 'irls', ...
        'irls', struct('maxIter', 60, 'tol', 1e-12, 'maxHalvings', 12), ...
        'em', struct('maxIter', 2000, 'tol', 1e-12, 'checkEvery', 10), ...
        'outputDir', '', 'makeFigure', true);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    for required = {'poisson_nnls_em', 'poisson_nnls_em_deviance', ...
            'biexp_slb_pattern', 'biexp_slb_deviance', ...
            'biexp_slb_basis', 'biexp_slb_pattern_batch', ...
            'poisson_nnls_em_batch', 'biexp_slb_profiled_batch', ...
            'biexp_slb_bfgs_batch', 'poisson_nnls_irls_batch', ...
            'poisson_nnls3_exact'}
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

    % A prior supplied by the caller is a deliberate session-wide choice, so it
    % gets its own output folder - otherwise a pinned run would overwrite the
    % per-analysis one and the two could not be compared.
    priorWasPinned = ~isempty(opts.slbTauNs);
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
    % The display region is independent of what gets fitted.
    if isempty(opts.displayMask)
        displayMask = true(nRow, nCol);
    else
        displayMask = immune_cell_MIET_biexp_mask(result, opts.displayMask, ...
            nRow, nCol);
    end
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
    if numel(tau1Grid) < opts.tau1GridCount
        fprintf(['  NOTE: %d of %d tau1 node(s) fell below the positivity ' ...
            'floor and were\n        dropped, so the grid reaches only ' ...
            '%+.2f sigma below the prior centre.\n'], ...
            opts.tau1GridCount - numel(tau1Grid), opts.tau1GridCount, ...
            (min(tau1Grid) - opts.slbTauNs) / max(opts.slbSigmaNs, eps));
    end
    fprintf('  grid: %d tau1 x %d tau2 = %d pairs\n', numel(tau1Grid), ...
        numel(tau2Grid), numel(tau1Grid) * numel(tau2Grid));

    if strcmpi(opts.method, 'vp')
        out = immune_cell_MIET_biexp_vp_run(cube, mask, pixelIndex, ...
            intensity, irf, dtNs, periodNs, nBin, opts);
    else
        out = immune_cell_MIET_biexp_run(cube, mask, pixelIndex, ...
            intensity, irf, dtNs, periodNs, nBin, tau1Grid, tau2Grid, opts);
    end
    if opts.binSize > 1
        % Match the erosion the binned fit applies, so the drawn region is one
        % the fit actually covers.
        kern = ones(opts.binSize);
        displayMask = displayMask & (convn(double(displayMask), kern, ...
            'same') >= opts.binSize ^ 2 - 0.5);
    end
    out.maps.displayMask = displayMask;
    out.analysisMat = analysisMat;
    out.imageSize = [nRow nCol];
    out.opts = opts;
    out.slbReferenceReducedChiSquare = slbChi;

    % ---- save -----------------------------------------------------------
    if isempty(opts.outputDir)
        folder = fileparts(analysisMat);
        name = 'biexp_slb';
        if opts.binSize > 1
            name = sprintf('%s_bin%d', name, opts.binSize);
        end
        if strcmpi(opts.method, 'grid')
            % 'vp' is the default and keeps the plain name; a grid run is the
            % deliberate exception and says so in the path.
            name = sprintf('%s_grid', name);
        end
        % A fit over something other than the cell footprint gets its own
        % folder, so whole-file and footprint-only results never overwrite one
        % another and the path says which is which.
        if ~(ischar(opts.pixelMask) || isstring(opts.pixelMask)) || ...
                ~strcmpi(char(opts.pixelMask), 'cellFootprint')
            if ischar(opts.pixelMask) || isstring(opts.pixelMask)
                name = sprintf('%s_%s', name, char(opts.pixelMask));
            else
                name = sprintf('%s_mask', name);
            end
        end
        if priorWasPinned
            % e.g. biexp_slb_tau0p350 - the choice is visible in the path.
            name = sprintf('%s_tau%s', name, strrep(sprintf('%.3f', ...
                opts.slbTauNs), '.', 'p'));
        end
        opts.outputDir = fullfile(folder, name);
    end
    if ~isfolder(opts.outputDir); mkdir(opts.outputDir); end
    out.opts.outputDir = opts.outputDir;
    matFile = fullfile(opts.outputDir, 'biexp_slb_maps.mat');
    save(matFile, 'out', '-v7.3');
    fprintf('\n  wrote %s\n', matFile);
    immune_cell_MIET_biexp_report(out, opts.outputDir);
    % Companion file in the shape the existing explorer GUI expects, so the
    % biexp maps can be browsed with the same tool. The TCSPC cube is
    % referenced, not copied.
    try
        out.explorerFile = export_biexp_for_explorer(matFile, analysisMat);
    catch exportError
        fprintf('  WARNING: explorer export failed (%s)
', exportError.message);
    end
    % Figures are NOT allowed to fail the run. The fit is already saved by this
    % point, and a plotting error used to propagate out of here and discard
    % minutes of completed fitting - which is how a MATLAB colorbar listener bug
    % came to throw away a whole acquisition. Report and carry on instead; the
    % figures can always be regenerated from the saved MAT.
    if opts.makeFigure
        try
            out.figure = immune_cell_MIET_biexp_figure(out, opts.outputDir);
            fprintf('  wrote %s\n', out.figure);
        catch figError
            fprintf('  WARNING: overview figure failed (%s); the fit is saved.\n', ...
                figError.message);
        end
        try
            out.componentFigure = ...
                immune_cell_MIET_biexp_component_figure(out, opts.outputDir);
            fprintf('  wrote %s\n', out.componentFigure);
        catch figError
            fprintf(['  WARNING: component figure failed (%s); the fit is ' ...
                'saved.\n'], figError.message);
        end
    end
end
