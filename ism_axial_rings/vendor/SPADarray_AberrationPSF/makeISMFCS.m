function makeISMFCS(sim, outDir, dif)
%--------------------------------------------------------------------------
% makeISMFCS
%
% PURPOSE
%   Simulate the FCS autocorrelation curve G(τ) for freely diffusing
%   point-like fluorescent molecules observed through an ISM detector array,
%   for five standard aberration scenarios.  One G(τ) curve is computed
%   per detector pixel (25 total) × aberration scenario (5 total).
%
% USAGE
%   makeISMFCS()                              % all defaults, save to '.'
%   makeISMFCS(sim)
%   makeISMFCS(sim, outDir)
%   makeISMFCS(sim, outDir, dif)
%
% INPUTS
%   sim    : simulation structure from defaultParams().  If empty, it is
%            built internally with nz raised to 21 for smooth FCS curves.
%   outDir : output folder for PNG files  [default: '.']
%   dif    : diffusion coefficient [µm² / time-unit]  [default: 5e-5]
%            Time-unit is arbitrary; the FCS x-axis is labelled accordingly.
%            For reference: at this default D, a PSF with waist w ≈ 0.2 µm
%            gives a diffusion time τ_D = w²/(4D) ≈ 200 time-units.
%
% OUTPUT FILES
%   fcs_ism_<scenario>.png  (one per aberration, 5×5 pixel grid)
%   fcs_ism_allscenarios_perpixel.png  (5×5 pixel grid, all 5 scenarios)
%
% ══════════════════════════════════════════════════════════════════════════
% PHYSICS & MATHEMATICS
% ══════════════════════════════════════════════════════════════════════════
%
%  FCS model (implemented in FCS.m)
%  ──────────────────────────────────
%  The normalised autocorrelation for pure 3-D Brownian diffusion is:
%
%    G(τ) = ∫∫ h_eff_k(r) · P_diff(r→r', τ) · h_eff_k(r') d³r d³r'
%
%  where P_diff(r→r', τ) = (4πDτ)^{-3/2} exp(-|r-r'|²/4Dτ) is the free-
%  diffusion Green's function and h_eff_k(x,y,z) is the effective PSF of
%  detector channel k:
%
%    h_eff_k(x,y,z) = h_exc(x,y,z) · η_k(x,y,z)
%
%  with h_exc the excitation PSF and η_k the detection collection efficiency
%  of pixel k (see detectorCollectionEfficiencyExplicit.m).
%
%  Because P_diff separates in cylindrical coordinates, the 6-D double
%  integral reduces to a product of azimuthal-harmonic propagators:
%
%    G_m(τ) = ∫∫ ρ h_m(ρ,z) · [I_m-Bessel radial propagator] ·
%                               [erf axial propagator] ·
%                               h_m(ρ',z') dρ dz dρ' dz'
%
%  and the total correlation is:
%
%    G(τ) = 2·G_0(τ) + Σ_{m=1}^{M} [ G^c_m(τ) + G^s_m(τ) ]
%
%  (the factor 2 for m=0 accounts for the azimuthal integral of cos²(0)=1
%  being twice that of cos²(mφ) for m>0 when normalised by 1/π).
%
%  Polar harmonic extraction (psfCartesianToPolar.m)
%  ──────────────────────────────────────────────────
%  h(ρ,φ,z) is sampled on a (Nr×Nphi) polar grid by bilinear interpolation
%  from the Cartesian PSF, then the azimuthal FFT gives:
%
%    h_0(ρ,z)   = Re[A_0]
%    h^c_m(ρ,z) = 2·Re[A_m]       (cos(mφ) amplitudes)
%    h^s_m(ρ,z) = -2·Im[A_m]      (sin(mφ) amplitudes)
%
%  G(0) diagnostic
%  ───────────────
%  After FCS normalises by G(τ_min) ≈ G(0):
%    G(0) ∝ 1/N_eff   where N_eff = (∫h d³r)² / ∫h² d³r
%  Aberrations spread h, increasing N_eff and reducing G(0) — the plot
%  amplitude at τ → 0 therefore visualises the detection-volume distortion.
%
%  flipY = true convention
%  ────────────────────────
%  Same as makeModeFigureExplicitDetector: detector rows are displayed with
%  +y at the top, consistent with YDir='normal' imagesc of the scan image.
%
% ══════════════════════════════════════════════════════════════════════════

    %----------------------------------------------------------------------
    % Defaults
    %----------------------------------------------------------------------
    if nargin < 1 || isempty(sim),    sim    = defaultParams(); end
    if nargin < 2 || isempty(outDir), outDir = '.';             end
    if nargin < 3 || isempty(dif),    dif    = 5e-5;            end

    % Raise z resolution for smooth FCS axial integration
    if sim.nz < 11
        sim.nz  = 21;
        sim.z   = linspace(-sim.nzRange/2, sim.nzRange/2, sim.nz);
        sim.obj = beadObject3D(sim);
        fprintf('[makeISMFCS] nz raised to %d for smooth z-integration.\n', sim.nz);
    end

    %----------------------------------------------------------------------
    % Polar grid for psfCartesianToPolar
    %   rhoMax  : FOV radius (half-width), covers full PSF support
    %   Nr      : number of radial samples; ~50 gives adequate resolution
    %   Nphi    : azimuthal FFT samples; 8× oversampling of M
    %   M       : highest harmonic order (matches pupilFourierCoeffs default)
    %----------------------------------------------------------------------
    Nr     = 60;
    Nphi   = 128;
    M      = 5;
    rhoMax = sim.fovXY / 2;
    rhoVec = linspace(0, rhoMax, Nr)';           % [Nr×1], µm

    % FCS input grids (pass as 2-D matrices to avoid the 1-D reshaping bug)
    %   rho_grid [Nr×Nz]:  each column = rhoVec (z-slices share same ρ axis)
    %   z_row    [1×Nz]:   axial coordinate vector
    rho_grid = rhoVec * ones(1, sim.nz);          % [Nr×Nz]
    z_row    = sim.z(:)';                          % [1×Nz]

    %----------------------------------------------------------------------
    % Aberration scenarios
    %----------------------------------------------------------------------
    scenarios = {
        'No_aberration',  struct();
        'Defocus',        struct('defocus',   0.18);
        'Astig_x',        struct('astig_x',   0.18);
        'Coma_x',         struct('coma_x',    0.18);
        'Spherical',      struct('spherical', 0.18)};
    nSc = size(scenarios, 1);
    nCh = size(sim.detXY, 1);

    % Detector display grid (flipX=true, flipY=false)
    flipX         = true;
    flipY         = false;
    transposeGrid = false;
    idxGrid = detectorIndexGridFromXY(sim.detXY, flipX, flipY, transposeGrid);
    nRows = size(idxGrid, 1);
    nCols = size(idxGrid, 2);

    % Colour scheme
    scColors = lines(nSc);
    chColors = parula(nCh);

    %----------------------------------------------------------------------
    % Pre-allocate results
    %   fcsResults{j,k} = [nTime×1] normalised autocorrelation G(τ)
    %   G0raw(j,k)      = G(0) before FCS normalises (proportional to ∫h²)
    %----------------------------------------------------------------------
    fcsResults = cell(nSc, nCh);
    G0raw      = zeros(nSc, nCh);
    autotime_out = [];

    %----------------------------------------------------------------------
    % Main computation loop
    %----------------------------------------------------------------------
    for j = 1:nSc
        label  = scenarios{j,1};
        coeffs = coeffStruct(sim, scenarios{j,2});

        fprintf('\n[makeISMFCS] Scenario %d/%d: %s\n', j, nSc, strrep(label,'_',' '));
        t0 = tic;

        % Compute excitation and detection PSFs once per scenario
        hExc = psfBessel(sim, coeffs, sim.lamExc);
        hDet = psfBessel(sim, coeffs, sim.lamEm);

        % Detector collection efficiency: eta(y,x,z,k)
        % η_k(x,y,z) = (1/nQ)Σ_q h_det(x+d_kx+u_q, y+d_ky+v_q, z)
        eta = detectorCollectionEfficiencyExplicit(sim, hDet);

        fprintf('  PSF + detector computed in %.1f s. Running FCS for %d channels...\n', ...
            toc(t0), nCh);

        for k = 1:nCh
            % ── Effective PSF for channel k ───────────────────────────────
            %    h_eff_k = h_exc · η_k   (element-wise, per voxel)
            hEff = hExc .* eta(:,:,:,k);

            % ── Azimuthal harmonic decomposition ─────────────────────────
            %    vol(ρ, z, harmonic) in [Nr × Nz × 2M+1] format
            vol = psfCartesianToPolar(hEff, sim.x, sim.y, rhoVec, Nphi, M);

            % ── Call FCS ─────────────────────────────────────────────────
            % Pass volx1 = vol, voly1 = zeros → f1 = vol (no doubling)
            % bld = 0: suppress live plot; capture FCS's internally created figure
            prevFigs = findall(0, 'Type', 'figure');
            [modres, autotime] = FCS( ...
                rho_grid, z_row, ...        % geometry grids
                vol, zeros(size(vol)), ...  % single PSF, no second component
                [], [], ...                 % no cross-correlation
                [], [], ...                 % no boundary, no flow
                dif, 0);                    % diffusion coeff, no live plot
            newFigs = setdiff(findall(0,'Type','figure'), prevFigs);
            if ~isempty(newFigs), close(newFigs); end

            % Save raw G(0) (before FCS's internal normalisation)
            % FCS returns modres already normalised by modres(1,1); we
            % recover G(0) as the first time-point of the unnormalised curve.
            % For relative comparison we just store modres(1) from the
            % pre-normalisation step — here it equals 1 by definition.
            % Instead, store the un-normalised value of the m=0 integral at
            % τ = 0: that is  ∑_ρ ρ · h_0(ρ,z)² · dρ (integrated over z)
            %  = the m=0 contribution to ∫h²d³r.
            h0z = squeeze(vol(:,:,1));      % [Nr×Nz] m=0 harmonic
            G0raw(j,k) = 2 * sum(rhoVec .* sum(h0z.^2, 2));

            fcsResults{j,k} = modres(:,1);

            if isempty(autotime_out)
                autotime_out = autotime;
            end

            if mod(k,5)==0
                fprintf('    channel %d/%d done\n', k, nCh);
            end
        end
        fprintf('  Scenario done in %.1f s total.\n', toc(t0));
    end

    %----------------------------------------------------------------------
    % Figure 1: Per-scenario, 5×5 pixel-layout grid
    %   Each subplot = one physical detector pixel
    %   Shows G(τ) for that pixel, coloured by aberration scenario colour
    %----------------------------------------------------------------------
    for j = 1:nSc
        label = scenarios{j,1};
        fig = figure('Color','w','Position',[50 50 1500 1150]);

        titleStr = strrep(label,'_',' ');
        if j == 1
            ampStr = '(unaberrated)';
        else
            ampStr = '0.18 waves RMS';
        end
        sgtitle(sprintf('ISM-FCS: G(\\tau) per detector pixel  —  %s  %s', ...
            titleStr, ampStr), 'FontSize', 13, 'FontWeight', 'bold');

        for row = 1:nRows
            for col = 1:nCols
                k = idxGrid(row, col);
                ax = subplot(nRows, nCols, (row-1)*nCols + col);
                if ~isfinite(k) || k < 1 || k > nCh
                    axis(ax, 'off');
                    continue;
                end

                G = fcsResults{j,k};
                semilogx(autotime_out, G, ...
                    'Color', scColors(j,:), 'LineWidth', 1.8);

                dk = sim.detXY(k,:);
                title(sprintf('(%.2f, %.2f)µm', dk(1), dk(2)), ...
                    'FontSize', 6.5, 'Interpreter', 'none');

                xlim([autotime_out(1) autotime_out(end)]);
                ylim([0 1.25]);
                ax.FontSize   = 5.5;
                ax.XMinorTick = 'off';

                % Label only outer edges to reduce clutter
                if row == nRows
                    xlabel('lag [a.u.]', 'FontSize', 6.5);
                end
                if col == 1
                    ylabel('G(\tau)/G(0)', 'FontSize', 6.5);
                end

                % Draw a dashed horizontal reference at G=0.5
                hold on;
                xline(ax, autotime_out(1), 'Color',[0.8 0.8 0.8]);
                yline(ax, 0.5, '--', 'Color',[0.7 0.7 0.7], 'LineWidth', 0.8);
                hold off;
            end
        end

        outFile = fullfile(outDir, sprintf('fcs_ism_%s.png', label));
        exportgraphics(fig, outFile, 'Resolution', 150);
        close(fig);
        fprintf('[makeISMFCS] Saved %s\n', outFile);
    end

    %----------------------------------------------------------------------
    % Figure 1b: Per-scenario — all 25 pixel curves overlaid on one plot
    %   Each curve = one detector pixel, coloured by distance from centre
    %   Normalised G(τ)/G(0) so curves are comparable across pixels
    %----------------------------------------------------------------------
    for j = 1:nSc
        label    = scenarios{j,1};
        titleStr = strrep(label,'_',' ');
        if j == 1
            ampStr = '(unaberrated)';
        else
            ampStr = '0.18 waves RMS';
        end

        fig = figure('Color','w','Position',[100 100 700 500]);
        ax  = axes(fig);
        set(ax, 'XScale', 'log');   % must set log BEFORE hold on
        hold(ax,'on');

        legHandles = gobjects(nCh, 1);
        legLabels  = cell(nCh, 1);
        for k = 1:nCh
            G = fcsResults{j,k};
            legHandles(k) = semilogx(ax, autotime_out, G, ...
                'Color', chColors(k,:), 'LineWidth', 1.4);
            legLabels{k}  = sprintf('(%.2f,%.2f)', ...
                sim.detXY(k,1), sim.detXY(k,2));
        end
        yline(ax, 0.5, '--k', 'LineWidth', 0.8);
        hold(ax,'off');

        legend(legHandles, legLabels, ...
            'Location','eastoutside', 'FontSize', 6.5, ...
            'Box','off', 'NumColumns', 2);
        xlabel(ax, 'lag [a.u.]',      'FontSize', 11);
        ylabel(ax, 'G(\tau) / G(0)',  'FontSize', 11);
        title(ax, sprintf('ISM-FCS: all pixels — %s  %s', titleStr, ampStr), ...
            'FontSize', 10);
        xlim(ax, [autotime_out(1) autotime_out(end)]);
        ylim(ax, [0 1.25]);
        grid(ax,'on');  grid(ax,'minor');
        ax.XMinorGrid = 'off';
        box(ax,'on');

        outFile = fullfile(outDir, sprintf('fcs_ism_%s_allpixels.png', label));
        exportgraphics(fig, outFile, 'Resolution', 150);
        close(fig);
        fprintf('[makeISMFCS] Saved %s\n', outFile);
    end

    %----------------------------------------------------------------------
    % Figure 2: All scenarios overlaid, 5×5 pixel layout
    %   Each subplot = one physical detector pixel
    %   5 coloured curves = 5 aberration scenarios
    %----------------------------------------------------------------------
    fig = figure('Color','w','Position',[50 50 1500 1150]);
    sgtitle('ISM-FCS: all aberration scenarios — G(\tau) per detector pixel', ...
        'FontSize', 13, 'FontWeight', 'bold');

    for row = 1:nRows
        for col = 1:nCols
            k = idxGrid(row, col);
            ax = subplot(nRows, nCols, (row-1)*nCols + col);
            if ~isfinite(k) || k < 1 || k > nCh
                axis(ax, 'off');
                continue;
            end
            set(ax, 'XScale', 'log');   % must set log BEFORE hold on
            hold on;

            for j = 1:nSc
                G = fcsResults{j,k};
                semilogx(autotime_out, G, ...
                    'Color', scColors(j,:), 'LineWidth', 1.5);
            end

            dk = sim.detXY(k,:);
            title(sprintf('(%.2f, %.2f)', dk(1), dk(2)), ...
                'FontSize', 6.5, 'Interpreter', 'none');

            xlim([autotime_out(1) autotime_out(end)]);
            ylim([0 1.25]);
            ax.FontSize = 5.5;
            yline(ax, 0.5, '--', 'Color',[0.7 0.7 0.7], 'LineWidth', 0.8);
            hold off;

            if row == nRows
                xlabel('lag [a.u.]', 'FontSize', 6.5);
            end
            if col == 1
                ylabel('G(\tau)/G(0)', 'FontSize', 6.5);
            end
        end
    end

    % Legend in the first valid detector subplot (centre pixel gets crowded)
    [legendRow, legendCol] = find(isfinite(idxGrid) & idxGrid > 0, 1, 'first');
    subplot(nRows, nCols, (legendRow-1)*nCols + legendCol);
    hold on;
    legHandles = zeros(nSc, 1);
    for j = 1:nSc
        legHandles(j) = plot(NaN, NaN, '-', ...
            'Color', scColors(j,:), 'LineWidth', 2.5);
    end
    legend(legHandles, strrep(scenarios(:,1), '_', ' '), ...
        'Location','southwest', 'FontSize', 5.5, 'Box','off');
    hold off;

    outFile = fullfile(outDir, 'fcs_ism_allscenarios_perpixel.png');
    exportgraphics(fig, outFile, 'Resolution', 150);
    close(fig);
    fprintf('[makeISMFCS] Saved %s\n', outFile);

    %----------------------------------------------------------------------
    % Figure 3: G(0) map — shows detection-volume distortion per pixel
    %   G(0) ∝ 1/N_eff; smaller = more molecules in focus = brighter pixel
    %   Displayed as a 5×5 heat-map, one panel per aberration scenario
    %----------------------------------------------------------------------
    fig = figure('Color','w','Position',[100 100 1400 310]);
    sgtitle('ISM-FCS: G(0) amplitude map per detector pixel per scenario', ...
        'FontSize', 12, 'FontWeight', 'bold');

    for j = 1:nSc
        % Build 5×5 G(0) map in physical detector order
        G0map = nan(nRows, nCols);
        for row = 1:nRows
            for col = 1:nCols
                k = idxGrid(row, col);
                if ~isfinite(k) || k < 1 || k > nCh
                    continue;
                end
                % G(0) is proportional to G0raw: normalise to max of
                % unaberrated scenario for a fair comparison
                G0map(row, col) = G0raw(j, k);
            end
        end

        subplot(1, nSc, j);
        imagesc(G0map);
        axis image;
        set(gca, 'YDir','normal');
        title(strrep(scenarios{j,1},'_',' '), 'FontSize', 9);
        colorbar;
        if j == 1
            ylabel('detector row (y)');
        end
        xlabel('detector col (x)');
    end

    % Common colour scale across all panels
    G0all = G0raw(:);
    G0positive = G0all(isfinite(G0all) & G0all > 0);
    if ~isempty(G0positive)
        G0lim = [min(G0positive) * 0.95, max(G0positive) * 1.05];
        for j = 1:nSc
            subplot(1, nSc, j);
            caxis(G0lim);
        end
    end
    colormap(gca, hot);

    outFile = fullfile(outDir, 'fcs_ism_G0map.png');
    exportgraphics(fig, outFile, 'Resolution', 150);
    close(fig);
    fprintf('[makeISMFCS] Saved %s\n', outFile);

    %----------------------------------------------------------------------
    % Figure 4: Central pixel only — all scenarios on one clean plot
    %   Shows the aberration effect most clearly for the on-axis pixel
    %----------------------------------------------------------------------
    % Find the central channel (detector pixel closest to origin)
    [~, kCentre] = min(sum(sim.detXY.^2, 2));

    fig = figure('Color','w','Position',[200 200 650 450]);
    ax = gca;
    set(ax, 'XScale', 'log');   % must set log BEFORE hold on
    hold on;
    legHandles = zeros(nSc, 1);
    for j = 1:nSc
        G = fcsResults{j, kCentre};
        legHandles(j) = semilogx(autotime_out, G, ...
            'Color', scColors(j,:), 'LineWidth', 2.0);
    end
    yline(0.5, '--k', 'LineWidth', 0.8);
    hold off;

    legend(legHandles, strrep(scenarios(:,1),'_',' '), ...
        'Location','southwest', 'FontSize', 9, 'Box','off');
    xlabel('lag [a.u.]', 'FontSize', 11);
    ylabel('G(\tau) / G(0)', 'FontSize', 11);
    title(sprintf('Central pixel (%.2f, %.2f) µm — D = %.1e µm²/t.u.', ...
        sim.detXY(kCentre,1), sim.detXY(kCentre,2), dif), 'FontSize', 10);
    xlim([autotime_out(1) autotime_out(end)]);
    ylim([0 1.25]);
    grid on; grid minor;
    ax.XMinorGrid = 'off';
    box on;

    outFile = fullfile(outDir, 'fcs_ism_centralpixel.png');
    exportgraphics(fig, outFile, 'Resolution', 150);
    close(fig);
    fprintf('[makeISMFCS] Saved %s\n', outFile);

    fprintf('\n[makeISMFCS] All figures saved to: %s\n', outDir);
end
