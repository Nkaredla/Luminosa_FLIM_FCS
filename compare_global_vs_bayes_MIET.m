function out = compare_global_vs_bayes_MIET(analysisMatFile, opts)
%COMPARE_GLOBAL_VS_BAYES_MIET Global pattern match vs per-pixel Bayesian maps.
%
% out = compare_global_vs_bayes_MIET(analysisMatFile)
% out = compare_global_vs_bayes_MIET(analysisMatFile, opts)
%
% Runs GlobalMultiExpPatternMatchFromTCSPC on the tcspc_pix cube already saved
% by immune_cell_MIET, and renders the per-pixel non-negative amplitude maps in
% the same three-panel layout as the Bayesian component figures so the two can
% be read side by side.
%
% WHY THIS COMPARISON
%
% The per-pixel fixed-SLB Bayesian route has three measured weaknesses on this
% data. The SLB lifetime is refitted per acquisition from an outside-cell
% region and came out between 0.108 and 0.399 ns across six nominally
% identical acquisitions, twice below the 0.16 ns bin width. Component
% counting is an M-open decision whose false-positive rate rose from 9% to 54%
% with photon count. And the SLB amplitude is a hard constraint carrying 17%
% calibration dispersion.
%
% Global analysis attacks all three at once:
%   - lifetimes come from the whole-image decay, of order 1e8 photons, so they
%     are far better determined than any per-pixel or sub-region fit
%   - each pixel estimates only NON-NEGATIVE AMPLITUDES on a fixed basis, so
%     "is a component present" becomes a continuous amplitude near zero rather
%     than a thresholded model-selection call
%   - the SLB becomes one basis pattern with a fitted amplitude, so nothing is
%     assumed, clipped, or calibrated
%
% WHAT IT COSTS
%
% Global analysis assumes the lifetimes are the same everywhere and only their
% amplitudes vary. For MIET that is not strictly true - lifetime varies
% continuously with emitter height, which is the signal. With a small basis
% that assumption is wrong; with a DENSE basis the per-pixel amplitudes become
% a lifetime distribution, which does represent continuous variation, but
% unmixing on a dense exponential basis is ill-conditioned (the Laplace
% inversion problem) and only recovers a smoothed distribution. Use
% opts.tau0 to choose where on that trade-off to sit, and read
% out.conditioning before trusting a dense result.
%
% opts fields (all optional)
%   tau0            initial lifetimes for the global fit.
%                   Default [0.35 1.7 3.3], i.e. the scale of the SLB and the
%                   two Bayesian membrane components on this data.
%   globalRegion    'all' | 'cell' | 'slb'  which pixels form the global decay.
%                   Default 'all'. 'cell' determines membrane lifetimes without
%                   the SLB dominating (the SLB is ~82% of the frame); 'slb'
%                   isolates the SLB reference itself.
%   mode            PatternMatchIm mode, default 'PIRLS'
%   useGPU          default false (the T1000 measured slower for this work)
%   displayMask     'cellMembrane' (default) | 'cellFootprint' | 'none'
%   outputDir       default alongside the analysis MAT
%   savePdf         also write a vector PDF, default true

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    defaults = struct( ...
        'tau0', [0.35 1.7 3.3], ...
        'globalRegion', 'all', ...
        'mode', 'PIRLS', ...
        'useGPU', false, ...
        'displayMask', 'cellMembrane', ...
        'outputDir', '', ...
        'savePdf', true);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if ~isfile(analysisMatFile)
        error('compare_global_vs_bayes_MIET:MissingFile', ...
            'Analysis MAT not found: %s', analysisMatFile);
    end
    for required = {'GlobalMultiExpPatternMatchFromTCSPC', 'Fluofit', ...
            'PatternMatchIm_matlab'}
        if exist(required{1}, 'file') ~= 2
            error('compare_global_vs_bayes_MIET:MissingDependency', ...
                'Required function %s.m is not on the MATLAB path.', ...
                required{1});
        end
    end
    if isempty(opts.outputDir)
        opts.outputDir = fileparts(analysisMatFile);
    end
    if ~isfolder(opts.outputDir)
        mkdir(opts.outputDir);
    end

    fprintf('\ncompare_global_vs_bayes_MIET\n  %s\n', analysisMatFile);
    store = matfile(analysisMatFile);
    variables = who(store);
    if ~ismember('tcspc_pix', variables)
        error('compare_global_vs_bayes_MIET:NoCube', ...
            ['This MAT has no tcspc_pix. Re-run with cfg.saveTcspcPix = ' ...
             'true.']);
    end
    result = store.result;
    cube = store.tcspc_pix;

    irf = double(result.irf.curve(:));
    dtNs = double(result.reassigned.dtNs);
    if ~isfinite(dtNs) || dtNs <= 0
        dtNs = double(result.config.tcspcBinNs);
    end
    pulsePeriodNs = double(result.bayesian.compact.pulsePeriodNs);
    masks = result.segmentation.masks;
    tauSlbBayes = double(result.slbReference.fixedLifetimeNs);

    fprintf(['  cube %s, dt %.4g ns, period %.4g ns, Bayesian tauSLB ' ...
        '%.4g ns\n'], mat2str(size(cube)), dtNs, pulsePeriodNs, tauSlbBayes);

    % The global decay is formed by summing the cube over pixels, so a region
    % is selected by zeroing everything outside it. The SLB covers most of the
    % frame, so 'all' determines the short component very precisely and the
    % membrane components only weakly; 'cell' reverses that.
    switch lower(opts.globalRegion)
        case 'all'
            regionMask = true(size(cube, 1), size(cube, 2));
        case 'cell'
            regionMask = logical(masks.cellMembrane);
        case 'slb'
            regionMask = logical(masks.slbReference);
        otherwise
            error('compare_global_vs_bayes_MIET:Region', ...
                'globalRegion must be all, cell or slb.');
    end
    fitCube = double(cube);
    fitCube(~repmat(regionMask, 1, 1, size(cube, 3))) = 0;
    fprintf('  global decay over %s: %d pixels, %.3g photons\n', ...
        opts.globalRegion, nnz(regionMask), sum(fitCube(:)));

    globalOpts = struct('useGPU', opts.useGPU, 'mode', opts.mode, ...
        'includeBackground', true, 'normalizePatterns', true, ...
        'sortLifetimes', true);
    started = tic;
    global_ = GlobalMultiExpPatternMatchFromTCSPC(fitCube, irf, ...
        pulsePeriodNs, dtNs, opts.tau0, globalOpts);
    elapsed = toc(started);

    fprintf('  global fit (%.1f s): tau = %s ns', elapsed, ...
        mat2str(round(global_.taufit(:).', 3)));
    if isfield(global_, 'chi') && ~isempty(global_.chi)
        fprintf(', chi2 = %.4g', global_.chi(1));
    end
    fprintf('\n');
    if isfield(global_, 'dtau') && ~isempty(global_.dtau)
        fprintf('  lifetime uncertainties: %s ns\n', ...
            mat2str(round(global_.dtau(:).', 4)));
    end
    fprintf('  Bayesian tauSLB was %.4g ns; global shortest is %.4g ns\n', ...
        tauSlbBayes, min(global_.taufit));

    % Amplitudes are non-negative and per pixel. Report them on the same
    % pixels the Bayesian component figures use, so the two are comparable.
    switch lower(opts.displayMask)
        case 'cellmembrane'
            showMask = logical(masks.cellMembrane);
        case 'cellfootprint'
            showMask = logical(masks.cellFootprint);
        otherwise
            showMask = true(size(cube, 1), size(cube, 2));
    end

    amplitude = double(global_.Amp);
    basisNames = global_.basisNames;
    % Drop the background column from the display panels; it is not a species.
    speciesColumns = 1:size(amplitude, 3);
    isBackground = false(1, numel(speciesColumns));
    for k = 1:numel(basisNames)
        if contains(lower(char(string(basisNames{k}))), 'back')
            isBackground(k) = true;
        end
    end
    speciesColumns = speciesColumns(~isBackground);

    photons = amplitude;
    figureFile = plotGlobalAmplitudes(photons, speciesColumns, basisNames, ...
        global_.taufit, showMask, opts, analysisMatFile);

    out = struct('global', global_, 'opts', opts, ...
        'tauSlbBayes', tauSlbBayes, 'dtNs', dtNs, ...
        'pulsePeriodNs', pulsePeriodNs, 'figure', figureFile, ...
        'seconds', elapsed);

    % Conditioning of the basis: exponentials are strongly non-orthogonal, so
    % report it rather than let the caller assume the unmixing is well posed.
    M = double(global_.patterns);
    out.conditioning = cond(M' * M);
    fprintf('  basis Gram condition number: %.3g\n', out.conditioning);
    if out.conditioning > 1e6
        fprintf(['  WARNING: the basis is close to degenerate. Per-pixel ' ...
            'amplitudes are\n  regularisation-dependent; treat them as a ' ...
            'smoothed distribution, not\n  as separated species.\n']);
    end
    fprintf('  wrote %s\n', figureFile);
end

function name = plotGlobalAmplitudes(photons, columns, basisNames, taufit, ...
        showMask, opts, analysisMatFile)
    panelCount = numel(columns);
    h = figure('Color', 'w', 'Visible', 'off', ...
        'Position', [60 60 460 * max(panelCount, 1), 620]);
    layout = tiledlayout(h, 1, max(panelCount, 1), 'Padding', 'compact', ...
        'TileSpacing', 'compact');

    for index = 1:panelCount
        column = columns(index);
        data = photons(:, :, column);
        mask = showMask & isfinite(data) & data > 0;
        ax = nexttile(layout);
        object = imagesc(ax, data.');
        axis(ax, 'image', 'off');
        set(ax, 'YDir', 'normal', 'Color', [0.015 0.015 0.025]);
        set(object, 'AlphaData', mask.');
        colormap(ax, turbo(256));
        values = data(mask);
        if isempty(values)
            clim(ax, [0 1]);
            detail = 'no pixels';
        else
            lower_ = prctile(values, 1);
            upper = prctile(values, 99);
            if ~isfinite(lower_) || ~isfinite(upper) || upper <= lower_
                lower_ = min(values); upper = max(values);
            end
            if upper <= lower_
                upper = lower_ + max(abs(lower_) * 1e-3, eps);
            end
            clim(ax, [lower_ upper]);
            detail = sprintf('%d px | median %.3g', nnz(mask), median(values));
        end
        bar = colorbar(ax);
        bar.Label.String = 'amplitude (photons)';
        label = sprintf('basis %d', column);
        if column <= numel(basisNames)
            label = char(string(basisNames{column}));
        end
        tauText = '';
        if column <= numel(taufit)
            tauText = sprintf(' | %.3g ns', taufit(column));
        end
        title(ax, {[label tauText], detail}, 'FontSize', 9);
    end

    [~, stem] = fileparts(analysisMatFile);
    title(layout, sprintf(['Global multiexponential amplitudes | region %s, ' ...
        'mode %s'], opts.globalRegion, opts.mode), 'FontWeight', 'bold');
    subtitle(layout, ['non-negative per-pixel amplitudes on a globally ' ...
        'fitted basis; no per-pixel model selection'], 'FontSize', 8);

    name = fullfile(opts.outputDir, sprintf('%s_global_amplitudes_%s.png', ...
        stem, lower(opts.globalRegion)));
    exportgraphics(h, name, 'Resolution', 200);
    if opts.savePdf
        exportgraphics(h, strrep(name, '.png', '.pdf'), ...
            'ContentType', 'vector');
    end
    close(h);
end
