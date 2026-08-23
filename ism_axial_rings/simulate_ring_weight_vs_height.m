function out = simulate_ring_weight_vs_height(opts)
%SIMULATE_RING_WEIGHT_VS_HEIGHT Detector-ring signature of emitter height.
%
% out = simulate_ring_weight_vs_height()
% out = simulate_ring_weight_vs_height(opts)
%
% Answers one question: for a molecule at height z above the focal plane,
% what fraction of its detected photons lands on each ring of the 23-channel
% SPAD array? And is that signature strong enough to separate membrane dye
% (in focus) from internalised dye (above focus) WITHOUT using lifetime?
%
% WHY THIS IS WORTH KNOWING
%
% The lifetime route cannot settle this. Beyond roughly 200 nm from the metal
% the MIET quenching curve flattens, so a long lifetime says "somewhere above
% 200 nm" and nothing more. The detector array carries an independent axial
% signal over that exact range: out-of-focus light spreads across more
% channels than in-focus light, so the DISTRIBUTION across channels encodes
% height regardless of lifetime. That makes the two reporters complementary
% rather than redundant.
%
% METHOD
%
% spadEffectivePSFArray gives hEff(y,x,z,k), the effective
% excitation-times-detection PSF for detector channel k, from the vectorial
% Richards-Wolf model in the Aberration codebase. Summing over the lateral
% scan gives the total signal a molecule at height z delivers to channel k:
%
%     W(k,z) = sum_{x,y} hEff(y,x,z,k)
%
% Channels are grouped into rings by their radial distance in the array, and
% two views are reported:
%
%   ABSOLUTE   W_ring(z)  - how much signal is lost with defocus
%   NORMALISED f_ring(z)  - the shape of the distribution, which is
%                           brightness-independent and is therefore what can
%                           actually separate heights from concentrations
%
% Separability is quantified with the Bhattacharyya coefficient between the
% ring distributions at two heights. For multinomial photon counts the
% classification error falls as BC^N, so the photons needed to distinguish
% two heights at error eps is about log(1/eps) / (-log BC). That is reported
% for both the ring-binned and the full 23-channel distributions, because
% binning into rings throws information away and the cost should be visible.
%
% WHAT THIS DOES NOT GIVE
%
% The PSF is symmetric about focus, so the ring signature measures |z|, not
% signed z. That is acceptable when the focus sits at the membrane and all
% the dye of interest is above it, but it cannot separate above from below.
%
% opts fields (all optional)
%   heightsUm        emitter heights to evaluate, default 0:0.05:1.0
%   lamExcUm         excitation wavelength, default 0.640
%   lamEmUm          emission wavelength, default 0.690
%   NA               objective NA, default 1.2
%   nMedium          sample refractive index, default 1.33
%   detPitchUm       detector pitch, default 0.18 (Luminosa honeycomb23)
%   fovUm            lateral field of view, default 3.0 (wide enough that
%                    defocused light is not truncated)
%   nx               lateral samples across the field, default 61
%   ringEdgesUm      radial ring boundaries; default derived from the layout
%   targetError      classification error for the photon-budget figure,
%                    default 0.05
%   mixtureHeightUm  the "internalised" height for the unmixing demo,
%                    default 0.4
%   outputDir        default alongside this file
%   makeFigure       default true

    thisDir = fileparts(mfilename('fullpath'));
    vendorDir = fullfile(thisDir, 'vendor', 'SPADarray_AberrationPSF');
    if isfolder(vendorDir) && ~contains([path pathsep], [vendorDir pathsep])
        addpath(vendorDir);
    end
    for required = {'defaultParams', 'psfBesselAirInterface', ...
            'spadEffectivePSFArrayInterface', 'detectorLayout', 'coeffStruct'}
        if exist(required{1}, 'file') ~= 2
            error('simulate_ring_weight_vs_height:MissingDependency', ...
                ['%s.m not found. It should be in ' ...
                 'ism_axial_rings/vendor/SPADarray_AberrationPSF.'], ...
                required{1});
        end
    end

    if nargin < 1 || isempty(opts); opts = struct(); end
    defaults = struct( ...
        'heightsUm', 0:0.05:1.0, ...
        'lamExcUm', 0.640, 'lamEmUm', 0.690, ...
        'NA', 1.45, 'nSample', 1.33, 'nImmersion', 1.52, ...
        'nGlass', 1.518, 'stageZUm', 0, ...
        'detPitchUm', 0.18, ...
        'fovUm', 3.0, 'nx', 61, ...
        'ringEdgesUm', [], ...
        'targetError', 0.05, ...
        'mixtureHeightUm', 0.4, ...
        'outputDir', thisDir, 'makeFigure', true);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if ~isfolder(opts.outputDir); mkdir(opts.outputDir); end

    heights = opts.heightsUm(:).';
    nz = numel(heights);

    % vectorialPSFBessel evaluates the axial phase as
    %   exp(1i*k0*nMedium*cosTheta*sim.z(iz))
    % so sim.z IS the emitter defocus and the whole sweep is one call.
    % NOTE beadBottomHeightUm does NOT do this - only airSurfaceBeadSlices
    % reads it, for the air-interface bead geometry, so driving the sweep
    % with it would return an identical PSF at every height.
    sim = defaultParams();
    sim.lamExc = opts.lamExcUm;
    sim.lamEm = opts.lamEmUm;
    sim.lamRef = opts.lamEmUm;
    % Full aperture, stratified sample. The homogeneous model would need
    % NA <= nMedium and so cannot express NA 1.45 in water, but the interface
    % model parameterises the pupil as q = NA*rho with
    % cosSample = sqrt(1-(q/nSample)^2), so supercritical components (q >
    % nSample) are retained with complex cosSample. Those are the
    % supercritical-angle fluorescence terms: near-field content of a dipole
    % close to the interface couples into propagating glass-side waves, the
    % reciprocal of TIRF. They are collected, and because they decay over
    % ~100-200 nm they are the steepest axial reporter available in exactly
    % the range the membrane occupies. Capping the NA would discard them.
    sim.NA = opts.NA;
    sim.nMedium = opts.nSample;
    sim.nSample = opts.nSample;
    sim.nGlass = opts.nGlass;
    sim.nDesignGlass = opts.nGlass;
    sim.nImmersion = opts.nImmersion;
    sim.sampleGeometry = 'airOnGlass';
    sim.detPitch = opts.detPitchUm;
    sim.detSize = opts.detPitchUm;
    sim.detectorLayout = 'honeycomb23';
    [sim.detXY, sim.detectorIndexGrid, sim.detectorLayoutInfo] = ...
        detectorLayout(sim.detectorLayout, sim.detPitch);
    sim.nDet = size(sim.detXY, 1);
    sim.fovXY = opts.fovUm;
    sim.nx = opts.nx;
    % defaultParams builds x, y and z at construction, so the grids must
    % be rebuilt after overriding the field of view - otherwise the
    % defaults silently persist and the sweep runs on the wrong grid.
    sim.x = linspace(-sim.fovXY / 2, sim.fovXY / 2, sim.nx);
    sim.y = sim.x;
    sim.dx = abs(sim.x(2) - sim.x(1));
    % In the stratified model the axial variable is the emitterHeight
    % ARGUMENT of psfBesselAirInterface, not sim.z. sim.z is left as a
    % single plane so the lateral grids stay consistent.
    sim.z = 0;
    sim.nz = 1;
    sim.nzRange = 0;

    nDet = sim.nDet;
    radius = sqrt(sum(sim.detXY .^ 2, 2));
    [ringIndex, ringEdges, ringRadius] = assignRings(radius, opts.ringEdgesUm);
    nRing = max(ringIndex);

    fprintf('\nsimulate_ring_weight_vs_height\n');
    fprintf(['  optics: NA %.3g (full aperture), n_imm %.4g, n_glass ' ...
        '%.4g, n_sample %.4g\n'], sim.NA, sim.nImmersion, ...
        sim.nGlass, sim.nSample);
    criticalNA = sim.nSample;
    fprintf(['    supercritical pupil fraction (q > n_sample): ' ...
        '%.1f%% of the radius\n'], ...
        100 * max(0, (sim.NA - criticalNA)) / sim.NA);
    fprintf('  stage focus: %.3g um relative to the interface\n', ...
        opts.stageZUm);
    fprintf('  wavelengths: exc %.3g um, em %.3g um; pitch %.3g um\n', ...
        sim.lamExc, sim.lamEm, sim.detPitch);
    fprintf('  field: %.3g um over %d samples; %d detectors\n', ...
        sim.fovXY, sim.nx, nDet);
    fprintf('  rings: %d, edges %s um, members %s\n', nRing, ...
        mat2str(round(ringEdges, 3)), mat2str(accumarray(ringIndex, 1)'));
    fprintf('  heights: %d from %.3g to %.3g um\n\n', nz, ...
        heights(1), heights(end));

    coeffs = coeffStruct(sim, zeros(numel(sim.modeOrder), 1));

    started = tic;
    hEff = spadEffectivePSFArrayInterface(sim, coeffs, ...
        opts.stageZUm, heights);
    fprintf('    forward model done in %.0f s, hEff %s\n', ...
        toc(started), mat2str(size(hEff)));
    % hEff is [ny nx nz nDet]; the lateral sum is the total signal a
    % molecule at that defocus delivers to each channel over the scan.
    channelWeight = reshape(sum(sum(double(hEff), 1), 2), nz, nDet).';
    if size(channelWeight, 1) ~= nDet || size(channelWeight, 2) ~= nz
        error('simulate_ring_weight_vs_height:Shape', ...
            'hEff collapsed to [%d %d]; expected [%d %d].', ...
            size(channelWeight, 1), size(channelWeight, 2), nDet, nz);
    end

    ringWeight = zeros(nRing, nz);
    for r = 1:nRing
        ringWeight(r, :) = sum(channelWeight(ringIndex == r, :), 1);
    end
    totalWeight = sum(ringWeight, 1);
    ringFraction = ringWeight ./ max(totalWeight, realmin);
    channelFraction = channelWeight ./ max(sum(channelWeight, 1), realmin);

    fprintf('\n  height   total signal   ring fractions\n');
    fprintf(  '  ------   ------------   --------------\n');
    for iz = 1:nz
        if mod(iz, 2) == 1 || iz == nz
            fprintf('  %6.3f   %12.4g   %s\n', heights(iz), ...
                totalWeight(iz) / max(totalWeight(1), realmin), ...
                mat2str(round(ringFraction(:, iz)', 4)));
        end
    end

    % Separability against the in-focus reference, for the ring-binned and
    % the full-channel distributions. Photons needed scales as
    % log(1/eps)/(-log BC).
    reference = 1;
    bcRing = ones(1, nz);
    bcChannel = ones(1, nz);
    for iz = 1:nz
        bcRing(iz) = bhattacharyya(ringFraction(:, reference), ...
            ringFraction(:, iz));
        bcChannel(iz) = bhattacharyya(channelFraction(:, reference), ...
            channelFraction(:, iz));
    end
    photonsRing = photonBudget(bcRing, opts.targetError);
    photonsChannel = photonBudget(bcChannel, opts.targetError);

    fprintf(['\n  photons needed to separate a molecule at height z from ' ...
        'one at focus (error %.3g)\n'], opts.targetError);
    fprintf('  height   using %d rings   using %d channels\n', nRing, nDet);
    for iz = 1:nz
        if heights(iz) == 0; continue; end
        if mod(iz, 2) == 1 || iz == nz
            fprintf('  %6.3f   %13s   %14s\n', heights(iz), ...
                budgetText(photonsRing(iz)), budgetText(photonsChannel(iz)));
        end
    end

    % Two-component unmixing: how precisely can the out-of-focus FRACTION be
    % recovered from the ring distribution alone? For a multinomial with
    % mixture p(phi) = (1-phi)*p0 + phi*p1, the Fisher information per photon
    % is sum_i (p1_i - p0_i)^2 / p_i, so the standard error on phi at N
    % photons is 1/sqrt(N*I).
    [~, mixtureIndex] = min(abs(heights - opts.mixtureHeightUm));
    p0 = ringFraction(:, reference);
    p1 = ringFraction(:, mixtureIndex);
    q0 = channelFraction(:, reference);
    q1 = channelFraction(:, mixtureIndex);
    phiGrid = [0.05 0.1 0.25 0.5];
    fisherRing = mixtureFisher(p0, p1, phiGrid);
    fisherChannel = mixtureFisher(q0, q1, phiGrid);

    fprintf(['\n  precision on the out-of-focus fraction, membrane at ' ...
        'focus vs dye at %.3g um\n'], heights(mixtureIndex));
    fprintf('  true phi   se at 1e3 photons (rings / channels)   at 1e4\n');
    for k = 1:numel(phiGrid)
        seRing3 = 1 / sqrt(1e3 * max(fisherRing(k), realmin));
        seChan3 = 1 / sqrt(1e3 * max(fisherChannel(k), realmin));
        seRing4 = 1 / sqrt(1e4 * max(fisherRing(k), realmin));
        fprintf('  %8.2f   %10.4f / %-10.4f            %10.4f\n', ...
            phiGrid(k), seRing3, seChan3, seRing4);
    end

    out = struct('heightsUm', heights, 'channelWeight', channelWeight, ...
        'ringWeight', ringWeight, 'ringFraction', ringFraction, ...
        'channelFraction', channelFraction, 'totalWeight', totalWeight, ...
        'ringIndex', ringIndex, 'ringEdgesUm', ringEdges, ...
        'ringRadiusUm', ringRadius, 'detXY', sim.detXY, ...
        'bhattacharyyaRing', bcRing, 'bhattacharyyaChannel', bcChannel, ...
        'photonsRing', photonsRing, 'photonsChannel', photonsChannel, ...
        'mixtureHeightUm', heights(mixtureIndex), ...
        'mixturePhiGrid', phiGrid, 'mixtureFisherRing', fisherRing, ...
        'mixtureFisherChannel', fisherChannel, 'sim', sim, 'opts', opts, ...
        'seconds', toc(started));

    if opts.makeFigure
        out.figure = plotRingWeights(out);
        fprintf('\n  wrote %s\n', out.figure);
    end

    fprintf(['\n  NOTE the PSF is symmetric about focus, so this measures ' ...
        '|z|, not signed z.\n  Optics default to NA %.3g in n %.4g - set ' ...
        'opts.NA and opts.nMedium to the real\n  objective before using ' ...
        'these numbers quantitatively.\n'], sim.NA, sim.nMedium);
end

% ------------------------------------------------------------------ helpers

function [ringIndex, edges, ringRadius] = assignRings(radius, requestedEdges)
    if ~isempty(requestedEdges)
        edges = reshape(sort(requestedEdges(:)), 1, []);
    else
        % Group by distinct radii in the honeycomb layout, tolerating
        % floating-point spread within a shell.
        uniqueRadii = uniquetol(radius, 1e-6);
        if numel(uniqueRadii) <= 1
            edges = [-inf inf];
        else
            midpoints = (uniqueRadii(1:end-1) + uniqueRadii(2:end)) / 2;
            edges = [-inf, midpoints(:)', inf];
        end
    end
    ringIndex = discretize(radius, edges);
    ringIndex(~isfinite(ringIndex)) = numel(edges) - 1;
    ringRadius = accumarray(ringIndex, radius, [], @mean);
    edges = edges(isfinite(edges));
end

function bc = bhattacharyya(p, q)
    p = max(double(p(:)), 0); q = max(double(q(:)), 0);
    p = p / max(sum(p), realmin); q = q / max(sum(q), realmin);
    bc = sum(sqrt(p .* q));
    bc = min(max(bc, realmin), 1);
end

function n = photonBudget(bc, targetError)
    % Chernoff-style bound: error <= BC^N, so N >= log(eps)/log(BC).
    n = log(targetError) ./ log(min(bc, 1 - eps));
    n(bc >= 1 - 1e-12) = Inf;
end

function text = budgetText(value)
    if ~isfinite(value)
        text = 'not separable';
    elseif value >= 1e6
        text = sprintf('%.1e', value);
    else
        text = sprintf('%.0f', value);
    end
end

function info = mixtureFisher(p0, p1, phiGrid)
    p0 = double(p0(:)); p1 = double(p1(:));
    info = zeros(1, numel(phiGrid));
    for k = 1:numel(phiGrid)
        phi = phiGrid(k);
        mix = (1 - phi) * p0 + phi * p1;
        good = mix > realmin;
        info(k) = sum((p1(good) - p0(good)) .^ 2 ./ mix(good));
    end
end

% ------------------------------------------------------------------ figures

function name = plotRingWeights(out)
    heights = out.heightsUm;
    h = figure('Color', 'w', 'Visible', 'off', 'Position', [60 60 1250 760]);
    layout = tiledlayout(h, 2, 3, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    colours = lines(size(out.ringWeight, 1));

    % (a) absolute signal per ring
    ax = nexttile(layout);
    hold(ax, 'on');
    for r = 1:size(out.ringWeight, 1)
        plot(ax, heights, out.ringWeight(r, :) / max(out.totalWeight(1), eps), ...
            '-', 'LineWidth', 1.5, 'Color', colours(r, :), ...
            'DisplayName', sprintf('ring %d (r=%.2f um)', r - 1, ...
                out.ringRadiusUm(r)));
    end
    plot(ax, heights, out.totalWeight / max(out.totalWeight(1), eps), 'k--', ...
        'LineWidth', 1.2, 'DisplayName', 'all channels');
    hold(ax, 'off'); grid(ax, 'on');
    xlabel(ax, 'emitter height [um]'); ylabel(ax, 'signal (focus = 1)');
    title(ax, 'Absolute collection vs height', 'FontSize', 9);
    legend(ax, 'Location', 'northeast', 'FontSize', 6, 'Box', 'off');

    % (b) normalised ring signature
    ax = nexttile(layout);
    hold(ax, 'on');
    for r = 1:size(out.ringFraction, 1)
        plot(ax, heights, out.ringFraction(r, :), '-', 'LineWidth', 1.5, ...
            'Color', colours(r, :));
    end
    hold(ax, 'off'); grid(ax, 'on'); ylim(ax, [0 1]);
    xlabel(ax, 'emitter height [um]'); ylabel(ax, 'fraction of photons');
    title(ax, 'Ring signature (brightness-independent)', 'FontSize', 9);

    % (c) outer-to-centre ratio
    ax = nexttile(layout);
    outerFraction = sum(out.ringFraction(2:end, :), 1);
    plot(ax, heights, outerFraction ./ max(out.ringFraction(1, :), realmin), ...
        'o-', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
    grid(ax, 'on');
    xlabel(ax, 'emitter height [um]'); ylabel(ax, 'outer / centre');
    title(ax, 'Single-number axial reporter', 'FontSize', 9);

    % (d) photon budget to separate from focus
    ax = nexttile(layout);
    keep = heights > 0;
    semilogy(ax, heights(keep), out.photonsRing(keep), 'o-', ...
        'LineWidth', 1.5, 'DisplayName', 'ring-binned');
    hold(ax, 'on');
    semilogy(ax, heights(keep), out.photonsChannel(keep), 's--', ...
        'LineWidth', 1.2, 'DisplayName', 'all 23 channels');
    hold(ax, 'off'); grid(ax, 'on');
    xlabel(ax, 'height difference from focus [um]');
    ylabel(ax, 'photons needed');
    title(ax, 'Cost of an axial call', 'FontSize', 9);
    legend(ax, 'Location', 'northeast', 'FontSize', 7, 'Box', 'off');

    % (e) detector map coloured by ring
    ax = nexttile(layout);
    scatter(ax, out.detXY(:, 1), out.detXY(:, 2), 220, out.ringIndex, ...
        'filled', 'MarkerEdgeColor', 'k');
    axis(ax, 'equal'); grid(ax, 'on');
    colormap(ax, colours);
    xlabel(ax, 'x [um]'); ylabel(ax, 'y [um]');
    title(ax, sprintf('honeycomb23 ring assignment (%d rings)', ...
        max(out.ringIndex)), 'FontSize', 9);

    % (f) precision on the out-of-focus fraction
    ax = nexttile(layout);
    photons = logspace(2, 5, 40);
    hold(ax, 'on');
    styles = {'-', '--', ':', '-.'};
    for k = 1:numel(out.mixturePhiGrid)
        se = 1 ./ sqrt(photons * max(out.mixtureFisherRing(k), realmin));
        plot(ax, photons, se, styles{min(k, numel(styles))}, ...
            'LineWidth', 1.4, 'DisplayName', ...
            sprintf('phi = %.2f', out.mixturePhiGrid(k)));
    end
    hold(ax, 'off');
    set(ax, 'XScale', 'log', 'YScale', 'log'); grid(ax, 'on');
    xlabel(ax, 'photons in the window');
    ylabel(ax, 'std error on out-of-focus fraction');
    title(ax, sprintf('Unmixing focus vs %.2f um', out.mixtureHeightUm), ...
        'FontSize', 9);
    legend(ax, 'Location', 'southwest', 'FontSize', 7, 'Box', 'off');

    title(layout, ['Detector-ring weight of a molecule versus height ' ...
        '(honeycomb23)'], 'FontWeight', 'bold');
    subtitle(layout, sprintf(['NA %.3g, n %.4g, exc %.3g um, em %.3g um, ' ...
        'pitch %.3g um - symmetric in z, so |z| only'], out.sim.NA, ...
        out.sim.nMedium, out.sim.lamExc, out.sim.lamEm, out.sim.detPitch), ...
        'FontSize', 8);

    name = fullfile(out.opts.outputDir, 'ism_ring_weight_vs_height.png');
    exportgraphics(h, name, 'Resolution', 200);
    exportgraphics(h, strrep(name, '.png', '.pdf'), 'ContentType', 'vector');
    close(h);
end
