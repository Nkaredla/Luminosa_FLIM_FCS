function out = measure_slb_reference_lifetime(source, opts)
%MEASURE_SLB_REFERENCE_LIFETIME Fit the BARE SLB outside the cell.
%
% out = measure_slb_reference_lifetime(analysisMatOrFolder)
% out = measure_slb_reference_lifetime(source, opts)
%
% Determines the SLB reference lifetime from the segmentation's slbReference
% pixels - the supported lipid bilayer OUTSIDE the cell.
%
% THE ANCHOR IS THE SHORT COMPONENT OF A BIEXPONENTIAL, NOT A MONO FIT
%
% A single exponential does not describe the pooled bare SLB. On 155036 it gives
% reduced deviance 3924, and it does so by inflating the background to three
% times the observed pre-pulse level in order to absorb a long-lived tail it
% cannot represent. Adding a second component drops the reduced deviance to 411
% and shows what the tail is: a few percent of the photons at around 1.3 ns,
% which is residual contamination rather than bilayer.
%
% So the bare-SLB decay is fitted with TWO components, seeded from a DistFluofit
% lifetime distribution so the search starts where the components actually are,
% and the SHORTER lifetime is taken as the SLB. That separates the bilayer from
% whatever else is in the region instead of averaging the two together, which is
% what the mono fit was doing.
%
% The per-pixel passes still use a mono fit: at a few hundred photons per pixel
% a two-component fit is not supportable, and their purpose is only to measure
% how much the bilayer VARIES across the field, not its absolute value.
%
% WHY THIS MEASUREMENT IS THE RIGHT ANCHOR
%
% The bilayer sits on a spacer above the metal, and that geometry does not change
% because a cell landed on top of it. So the SLB lifetime under the cell should
% equal the SLB lifetime beside the cell, up to whatever real bilayer roughness
% exists. That makes the bare-SLB region a direct, assumption-free measurement of
% the quantity the two-component cell fit puts a prior on.
%
% It replaces a bad anchor. The prior centre used until now came from
% result.slbReference.fixedLifetimeNs, which an audit traced to an UNSEEDED
% multistart in the IRF fitter: the same photons returned 0.134897 ns and
% 0.344474 ns in two analyses of one acquisition, and across the session the
% value ranged 0.108-0.399 ns. Anchoring on a number that moves 2.6x on
% identical data cannot be right; anchoring on the bare SLB can.
%
% WHAT IT REPORTS, AND WHY BOTH BINNINGS
%
%   pooled      one lifetime from all bare-SLB photons - the centre to use.
%   per pixel   at 1x1 and at 4x4. The 1x1 spread is dominated by shot noise and
%               says little about real variation; the 4x4 spread, with 16x the
%               photons, is a far better estimate of how much the bilayer
%               ACTUALLY varies across the field, and that is what should set
%               the prior width. Reporting only one of the two would confuse
%               shot noise with structure.
%
% AN SLB-ONLY MASK, BUILT FROM THE MEAN FLIM MAP
%
% The stored slbReference mask is GEOMETRIC - it is the region outside the cell -
% so it still admits anything bright and long-lived that happens to sit there:
% debris, a second cell's skirt, membrane fragments. Since the SLB is the
% shortest-lived thing in the field, those contaminants can be removed by their
% arrival time alone. The candidate region is therefore tightened by keeping only
% pixels whose mean arrival time falls in the lower opts.meanTauQuantile of the
% candidates' own distribution (default 0.6), optionally also below an absolute
% opts.maxMeanTauNs.
%
% The threshold is a QUANTILE OF THE CANDIDATES rather than an absolute number so
% it adapts to each acquisition instead of hard-coding a value that the audit
% already showed varies across this session.
%
% This mask is used ONLY to measure the SLB reference. It is not used for, and
% does not touch, the cell segmentation.
%
% out.suggestedPrior gives (centre, sigma) to pass as opts.slbTauNs and
% opts.slbSigmaNs for the cell fit.

    if nargin < 2 || isempty(opts); opts = struct(); end
    defaults = struct('maskName', 'slbReference', 'minPhotons', 0, ...
        'edgeMarginPix', 6, 'tauRangeNs', [0.05 2.0], 'tauNodes', 60, ...
        'binForVariation', 4, 'blockSize', 30000, 'makeFigure', true, ...
        'meanTauQuantile', 0.6, 'maxMeanTauNs', [], 'outputDir', '', ...
        'outputRoot', '', ...
        'anchorModel', 'biexp', 'distOpts', struct());
    % 60 nodes over the range is ~2.5% spacing; parabolic refinement then takes
    % it well below the shot-noise width, so a finer scan only costs time.
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end

    analysisMat = immune_cell_MIET_biexp_resolve(source);
    fprintf('\nmeasure_slb_reference_lifetime\n  %s\n', analysisMat);
    loaded = load(analysisMat, 'result');
    result = loaded.result;
    dtNs = double(immune_cell_MIET_explorer_field(result, ...
        'bayesian.compact.dtNs'));
    periodNs = double(immune_cell_MIET_explorer_field(result, ...
        'bayesian.compact.pulsePeriodNs'));
    irf = double(immune_cell_MIET_explorer_field(result, 'irf.curve'));
    irf = max(irf(:), 0);
    if sum(irf) > 0; irf = irf / sum(irf); end

    cubeStruct = load(analysisMat, 'tcspc_pix');
    cube = cubeStruct.tcspc_pix;
    [nRow, nCol, nBin] = size(cube);
    basis = biexp_slb_basis(irf, dtNs, periodNs, nBin);

    slb = immune_cell_MIET_biexp_mask(result, opts.maskName, nRow, nCol);
    cell = immune_cell_MIET_biexp_mask(result, 'cellFootprint', nRow, nCol);
    % Stay clear of the cell edge: reconstruction spreads signal across a few
    % pixels, so bare-SLB pixels adjacent to the footprint are contaminated by
    % membrane photons and would drag the reference toward the cell value - the
    % very thing being tested.
    if opts.edgeMarginPix > 0
        k = 2 * opts.edgeMarginPix + 1;
        near = convn(double(cell), ones(k), 'same') > 0.5;
        slb = slb & ~near;
    end
    intensity = sum(double(cube), 3);
    slb = slb & intensity >= opts.minPhotons;
    nGeometric = nnz(slb);

    % Tighten to SLB-ONLY using the mean arrival-time map: the bilayer is the
    % shortest-lived thing in the field, so anything long-lived sitting outside
    % the cell is not SLB and should not set the reference.
    meanTau = immune_cell_MIET_explorer_field(result, ...
        'redMeanFlim.meanArrivalNs');
    if isempty(meanTau)
        meanTau = immune_cell_MIET_explorer_field(result, ...
            'blueMeanFlim.meanArrivalNs');
    end
    if ~isempty(meanTau) && isequal(size(meanTau), [nRow nCol])
        candidate = double(meanTau);
        vals = candidate(slb);
        vals = vals(isfinite(vals));
        cut = inf;
        if ~isempty(vals) && opts.meanTauQuantile < 1
            cut = quantileLocalBiexp(vals, opts.meanTauQuantile);
        end
        if ~isempty(opts.maxMeanTauNs)
            cut = min(cut, opts.maxMeanTauNs);
        end
        slb = slb & isfinite(candidate) & candidate <= cut;
        fprintf(['  mean-arrival cut at %.4f ns (quantile %.2f of the ' ...
            'candidates) kept %d of %d px\n'], cut, opts.meanTauQuantile, ...
            nnz(slb), nGeometric);
    else
        fprintf('  NOTE: no mean-arrival map found; using the geometric mask only\n');
    end

    idx = find(slb);
    fprintf(['  %d bare-SLB pixel(s) (%s + arrival cut, >= %d photons, %d px ' ...
        'clear of the cell), median %.0f photons\n'], numel(idx), ...
        opts.maskName, opts.minPhotons, opts.edgeMarginPix, ...
        median(intensity(idx)));
    if numel(idx) < 100
        error('measure_slb_reference_lifetime:TooFewPixels', ...
            'Only %d bare-SLB pixels survived the masks.', numel(idx));
    end

    tauGrid = logspace(log10(opts.tauRangeNs(1)), ...
        log10(opts.tauRangeNs(2)), opts.tauNodes);
    cubeFlat = reshape(cube, [], nBin);

    % ---- pooled: the anchor ---------------------------------------------
    pooled = sum(double(cubeFlat(idx, :)), 1)';
    fprintf('  POOLED bare-SLB decay, %.3g photons\n', sum(pooled));
    if strcmpi(opts.anchorModel, 'mono')
        pooledFit = fit_mono_poisson_batch(pooled, basis, tauGrid, struct());
        anchorNs = pooledFit.tauNs;
        fprintf('      mono tau = %.4f ns, redDev %.1f, B %.2f/bin\n', ...
            pooledFit.tauNs, pooledFit.reducedDeviance, pooledFit.background);
    else
        pooledFit = fit_biexp_distfluofit_seeded(pooled, irf, dtNs, ...
            periodNs, opts.distOpts);
        % The SLB is the SHORT component; the long one is the contamination
        % the mono fit used to smear into it.
        anchorNs = pooledFit.slbTauNs;
    end
    if pooledFit.reducedDeviance > 3
        fprintf(['      NOTE: even this model leaves reduced deviance %.1f, ' ...
            'so the anchor is the\n      best available summary rather than ' ...
            'a fully described decay.\n'], pooledFit.reducedDeviance);
    end

    % ---- per pixel at 1x1 -----------------------------------------------
    perPixel = runBlocks(cubeFlat, idx, basis, tauGrid, opts.blockSize, nBin);
    fprintf('  1x1     tau median %.4f ns, IQR [%.4f %.4f]\n', ...
        median(perPixel.tauNs), quantileLocalBiexp(perPixel.tauNs, 0.25), ...
        quantileLocalBiexp(perPixel.tauNs, 0.75));

    % ---- per pixel at the variation binning -----------------------------
    b = opts.binForVariation;
    binned = struct('tauNs', [], 'idx', []);
    if b > 1
        kern = ones(b);
        cubeB = convn(double(cube), kern, 'same');
        maskB = slb & (convn(double(slb), kern, 'same') >= b ^ 2 - 0.5);
        intenB = sum(cubeB, 3);
        maskB = maskB & intenB >= opts.minPhotons * b ^ 2;
        idxB = find(maskB);
        if numel(idxB) > 100
            flatB = reshape(cubeB, [], nBin);
            binned = runBlocks(flatB, idxB, basis, tauGrid, ...
                opts.blockSize, nBin);
            binned.idx = idxB;
            fprintf(['  %dx%d     tau median %.4f ns, IQR [%.4f %.4f] ' ...
                '(%d px, median %.0f photons)\n'], b, b, ...
                median(binned.tauNs), ...
                quantileLocalBiexp(binned.tauNs, 0.25), ...
                quantileLocalBiexp(binned.tauNs, 0.75), numel(idxB), ...
                median(intenB(idxB)));
        end
    end

    % ---- suggested prior -------------------------------------------------
    centre = anchorNs;
    if ~isempty(binned.tauNs)
        iqr = quantileLocalBiexp(binned.tauNs, 0.75) - ...
            quantileLocalBiexp(binned.tauNs, 0.25);
        sigma = max(iqr / 1.349, 0.005);
    else
        sigma = 0.02;
    end
    fprintf(['\n  SUGGESTED PRIOR for the cell fit: slbTauNs %.4f, ' ...
        'slbSigmaNs %.4f\n'], centre, sigma);
    fprintf(['      centre = SHORT component of the pooled bare-SLB ' ...
        'biexponential; width from\n      the %dx%d per-pixel spread, which ' ...
        'at %.0fx the photons is dominated by real\n      variation rather ' ...
        'than shot noise. A SMALL wiggle by design: the bilayer\n      ' ...
        'geometry does not change under the cell, so tau1 there should sit ' ...
        'at this\n      value.\n'], b, b, b ^ 2);

    out = struct('analysisMat', analysisMat, 'slbOnlyMask', slb, ...
        'geometricPixelCount', nGeometric, 'anchorNs', anchorNs, ...
        'anchorModel', opts.anchorModel, 'pooled', pooledFit, ...
        'perPixel', perPixel, 'binned', binned, 'pixelIndex', idx, ...
        'imageSize', [nRow nCol], 'maskName', opts.maskName, ...
        'suggestedPrior', struct('slbTauNs', centre, 'slbSigmaNs', sigma), ...
        'storedReferenceNs', double(immune_cell_MIET_explorer_field(result, ...
            'slbReference.fixedLifetimeNs')), 'opts', opts);

    stored = out.storedReferenceNs;
    if ~isempty(stored)
        fprintf(['  for comparison, the value stored in the analysis is ' ...
            '%.4f ns (%+.1f%%)\n'], stored(1), ...
            100 * (stored(1) - centre) / centre);
    end

    if isempty(opts.outputDir)
        if isempty(opts.outputRoot)
            opts.outputDir = fullfile(fileparts(analysisMat), 'slb_reference');
        else
            % Keep every result for the session under one root, named by
            % acquisition, rather than scattered through the data tree.
            [~, acqName] = fileparts(fileparts(fileparts(analysisMat)));
            opts.outputDir = fullfile(opts.outputRoot, ...
                sprintf('%s_slb_reference', acqName));
        end
    end
    if ~isfolder(opts.outputDir); mkdir(opts.outputDir); end
    out.dtNs = dtNs;
    out.periodNs = periodNs;
    % Kept so the figure can redraw the rejected mono fit and show WHICH pixels
    % were used against the real image, rather than as a bare 0/1 map.
    out.pooledIrf = irf;
    out.intensityMap = intensity;
    out.outputDir = opts.outputDir;
    save(fullfile(opts.outputDir, 'slb_reference_fit.mat'), 'out', '-v7.3');
    fprintf('  wrote %s\n', fullfile(opts.outputDir, 'slb_reference_fit.mat'));
    % makeFigure was previously declared and never acted on - the option existed
    % and produced nothing. It now writes the figure it always promised.
    if opts.makeFigure
        try
            out.figure = slb_reference_figure(out, opts.outputDir);
            fprintf('  wrote %s\n', out.figure);
        catch figError
            fprintf(['  WARNING: SLB reference figure failed (%s); the fit ' ...
                'is saved.\n'], figError.message);
        end
    end
end
