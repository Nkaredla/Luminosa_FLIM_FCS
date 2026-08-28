function results = benchmark_biexp_speed(analysisMat, opts)
%BENCHMARK_BIEXP_SPEED Where the per-pixel fit's time actually goes.
%
% results = benchmark_biexp_speed(analysisMat)
% results = benchmark_biexp_speed(analysisMat, opts)
%
% Measures the three things that could make the fit faster, separately, so the
% decision is made on numbers rather than on which one sounds most promising:
%
%   1. BLOCK SIZE. The run measured 4.8 ms/pixel on 19000-pixel blocks but
%      2.0 ms/pixel on 4000-pixel blocks - same code. That points at memory
%      traffic rather than arithmetic, and it is free to fix.
%   2. SINGLE vs DOUBLE. The decays are photon counts in 156 bins; double buys
%      nothing here and costs twice the bandwidth.
%   3. GPU. This machine has a Quadro T1000 Max-Q. Its FP64 rate is 1/32 of
%      FP32, so in DOUBLE it is slower than the CPU - which is very likely why
%      an earlier attempt at GPU in this project measured 0.41-0.64x. The
%      question worth asking is whether SINGLE on the GPU wins, and that
%      depends on whether the 156-by-P arrays are large enough to cover the
%      transfer, so it has to be measured at several block sizes.
%
% Deliberately bounded: opts.maxPixels caps the sample so a benchmark cannot
% itself exhaust memory. This machine crashed during a full run with the page
% file at zero peak usage and C: 97% full, so peak allocation is a real
% constraint and not a theoretical one.

    if nargin < 2 || isempty(opts); opts = struct(); end
    defaults = struct('blockSizes', [1000 2000 4000 8000 16000], ...
        'maxPixels', 16000, 'tryGpu', true, 'slbTauNs', 0.0769, ...
        'slbSigmaNs', 0.0082);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end

    L = load(analysisMat, 'result');
    result = L.result;
    dtNs = double(immune_cell_MIET_explorer_field(result, ...
        'bayesian.compact.dtNs'));
    periodNs = double(immune_cell_MIET_explorer_field(result, ...
        'bayesian.compact.pulsePeriodNs'));
    irf = double(immune_cell_MIET_explorer_field(result, 'irf.curve'));
    irf = max(irf(:), 0); irf = irf / sum(irf);

    % Read only the pixels needed, through matfile, so the 452 MB double cube
    % is never materialised just to benchmark.
    % Load the cube ONCE as uint16 (113 MB) and index it, rather than pulling
    % pixels one at a time through matfile. The per-pixel route is what made an
    % earlier version of this benchmark time out before it fitted anything: the
    % measurement harness was slower than the thing being measured.
    S = load(analysisMat, 'tcspc_pix');
    cube = S.tcspc_pix;
    clear S;
    [nRow, nCol, nBin] = size(cube);
    mask = immune_cell_MIET_biexp_mask(result, 'cellFootprint', nRow, nCol);
    idx = find(mask);
    pick = idx(round(linspace(1, numel(idx), min(opts.maxPixels, numel(idx)))));
    flat = reshape(cube, [], nBin);
    Y = double(flat(pick, :))';
    clear cube flat;
    fprintf(['\nbenchmark_biexp_speed\n  %d pixel(s), %d bins, median %.0f ' ...
        'photons\n'], numel(pick), nBin, median(sum(Y, 1)));

    % One basis per precision. Casting only Y would be pointless: MATLAB
    % promotes to the wider operand, so a single Y against a double basis runs
    % entirely in double and the "single" row would measure nothing.
    basisDouble = biexp_slb_basis(irf, dtNs, periodNs, nBin, 'double');
    basisSingle = biexp_slb_basis(irf, dtNs, periodNs, nBin, 'single');
    base = struct('slbTauNs', opts.slbTauNs, 'slbSigmaNs', opts.slbSigmaNs, ...
        'gtol', 1e-3, 'fixSlbTau', true, 'innerSolver', 'irls');

    results = struct('label', {}, 'blockSize', {}, 'msPerPixel', {}, ...
        'tau2Median', {}, 'note', {});

    % ---- 1 + 2: block size, double then single --------------------------
    for useSingle = [false true]
        for b = opts.blockSizes
            if b > numel(pick); continue; end
            cols = 1:b;
            Yb = Y(:, cols);
            thisBasis = basisDouble;
            if useSingle; Yb = single(Yb); thisBasis = basisSingle; end
            t0 = tic;
            vp = biexp_slb_bfgs_batch(Yb, thisBasis, base);
            el = toc(t0);
            lbl = 'double';
            if useSingle; lbl = 'single'; end
            results(end + 1) = struct('label', lbl, 'blockSize', b, ...
                'msPerPixel', 1000 * el / b, ...
                'tau2Median', median(double(vp.tau2Ns)), 'note', ''); %#ok<AGROW>
            fprintf('  %-7s block %6d : %6.3f ms/px   tau2 %.4f ns\n', ...
                lbl, b, 1000 * el / b, median(double(vp.tau2Ns)));
        end
    end

    % ---- 3: GPU, single and double --------------------------------------
    if opts.tryGpu
        try
            g = gpuDevice;
            fprintf(['  GPU: %s, %.1f GB free, FP64 is 1/32 of FP32 on this ' ...
                'card\n'], g.Name, g.AvailableMemory / 1e9);
            for useSingle = [true false]
                b = min(max(opts.blockSizes), numel(pick));
                Yb = Y(:, 1:b);
                if useSingle; Yb = single(Yb); end
                try
                    Yg = gpuArray(Yb);
                    if useSingle; bg = basisSingle; else; bg = basisDouble; end
                    bg.C = gpuArray(bg.C);
                    bg.timeNs = gpuArray(bg.timeNs);
                    t0 = tic;
                    vp = biexp_slb_bfgs_batch(Yg, bg, base);
                    wait(g);
                    el = toc(t0);
                    lbl = 'gpu-double';
                    if useSingle; lbl = 'gpu-single'; end
                    results(end + 1) = struct('label', lbl, 'blockSize', b, ...
                        'msPerPixel', 1000 * el / b, ...
                        'tau2Median', median(double(gather(vp.tau2Ns))), ...
                        'note', ''); %#ok<AGROW>
                    fprintf('  %-10s block %6d : %6.3f ms/px   tau2 %.4f ns\n', ...
                        lbl, b, 1000 * el / b, ...
                        median(double(gather(vp.tau2Ns))));
                catch gpuError
                    fprintf('  GPU run failed (%s): %s\n', ...
                        class(Yb), gpuError.message);
                end
                clear Yg bg;
            end
        catch noGpu
            fprintf('  no usable GPU: %s\n', noGpu.message);
        end
    end

    % ---- summary --------------------------------------------------------
    if ~isempty(results)
        [~, best] = min([results.msPerPixel]);
        fprintf(['\n  FASTEST: %s at block %d, %.3f ms/px\n' ...
            '  the full run measured 4.8 ms/px at block 19000, so the ' ...
            'speedup available is %.1fx\n'], results(best).label, ...
            results(best).blockSize, results(best).msPerPixel, ...
            4.8 / results(best).msPerPixel);
        taus = [results.tau2Median];
        fprintf(['  tau2 across every configuration: %.4f to %.4f ns ' ...
            '(spread %.2f%%) - any\n  configuration that changes the ANSWER ' ...
            'is not an optimisation\n'], min(taus), max(taus), ...
            100 * (max(taus) - min(taus)) / median(taus));
    end
end
