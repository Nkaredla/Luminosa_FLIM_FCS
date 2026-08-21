function makePSF3DAllPixels(sim, outDir, tsh)
%--------------------------------------------------------------------------
% makePSF3DAllPixels
%
% PURPOSE
%   Generate 3-D isosurface visualizations of the effective PSF for every
%   detector pixel in the 5×5 ISM array, for each of the five standard
%   aberration scenarios.  One figure is produced per scenario; within each
%   figure the 25 sub-axes are laid out in the physical detector order.
%   Rendering is done via FocusImage3D (isosurface + Gouraud shading).
%
% USAGE
%   makePSF3DAllPixels()              % all defaults, save to current dir
%   makePSF3DAllPixels(sim)
%   makePSF3DAllPixels(sim, outDir)
%   makePSF3DAllPixels(sim, outDir, tsh)
%
% INPUTS
%   sim    : simulation structure from defaultParams().  If empty or absent,
%            defaultParams() is called internally.  nz is bumped to 21
%            automatically when < 11 to give smooth 3-D isosurfaces.
%   outDir : output directory for PNG files  [default: '.']
%   tsh    : 1×3 vector of isosurface levels as fractions of the channel
%            peak  [default: [0.5 0.2 0.05]  → 50 %, 20 %, 5 % of peak]
%
% OUTPUT FILES (one per scenario, in outDir)
%   psf3d_No_aberration.png
%   psf3d_Defocus.png
%   psf3d_Astig_x.png
%   psf3d_Coma_x.png
%   psf3d_Spherical.png
%
% ══════════════════════════════════════════════════════════════════════════
% PHYSICS / MATH SUMMARY
% ══════════════════════════════════════════════════════════════════════════
%
%  PSF computation  (scalarPSFBessel)
%  ──────────────────────────────────
%  Scalar Richards–Wolf diffraction with Abbe sine-condition apodization:
%
%    U(r,ψ,z) = Σ_m  2π · i^m · e^{imψ} ·
%               ∫₀¹ J_m(k₀·NA·r·ρ) · p_m(ρ) · e^{ik₀n·z·cosθ} · w(ρ) dρ
%
%  where  w(ρ) = (NA/n)² · ρ / cosθ  is the solid-angle Jacobian and
%  p_m(ρ) are the azimuthal Fourier coefficients of the aberrated pupil:
%
%    P(ρ,φ) = √cosθ · exp( i·2π·(λ_ref/λ)·W(ρ,φ) )
%    W      = Σ_k  a_k · Z_k(ρ,φ)         (Zernike expansion, OSA/ANSI)
%
%  PSF intensity:  h(x,y,z) = |U|²   (normalized to unit total power)
%
%  Detector collection efficiency  (detectorCollectionEfficiencyExplicit)
%  ────────────────────────────────────────────────────────────────────────
%  For detector pixel k centered at (d_kx, d_ky) in sample-equivalent
%  coordinates, the collection efficiency at sample point (x,y,z) is:
%
%    η_k(x,y,z) = (1/nQ) Σ_q  h_det( x + d_kx + uq,  y + d_ky + vq,  z )
%
%  where (uq,vq) are uniform quadrature points across the pixel area.
%  The (+x, +y) sign correctly accounts for microscope image inversion:
%  moving the emitter by +x shifts its image on the detector by −x, so
%  the detector pixel at +d_kx collects from emitters displaced by +d_kx
%  relative to the detection focus.
%
%  Effective PSF per detector channel
%  ───────────────────────────────────
%    h_eff_k(x,y,z) = h_exc(x,y,z) · η_k(x,y,z)
%
%  This is the quantity visualised here.  The scanned ISM image is the
%  correlation of the fluorescent object with h_eff_k (see
%  normalizedStackExplicitDetector for the full forward model).
%
%  Zernike modes (OSA/ANSI, normalised)
%  ─────────────────────────────────────
%    tilt_x/y   = 2ρ cos/sin(φ)
%    defocus    = √3 (2ρ²−1)
%    astig_x/y  = √6 ρ² cos/sin(2φ)
%    coma_x/y   = √8 (3ρ³−2ρ) cos/sin(φ)
%    spherical  = √5 (6ρ⁴−6ρ²+1)
%
% ══════════════════════════════════════════════════════════════════════════
% FLIP CONVENTION
% ══════════════════════════════════════════════════════════════════════════
%  flipX=false, flipY=true  (identical to makeModeFigureExplicitDetector)
%
%  With YDir='normal' in imagesc:  row 1 of the 5×5 grid → most-negative y
%  (bottom of image), row 5 → most-positive y (top).  This aligns the
%  detector-array orientation with the scan-image orientation.
%
%  Why flipY=true is necessary
%  ───────────────────────────
%  detectorIndexGridFromXY default (no flip): top row = most-positive y.
%  When displayed with imagesc + YDir='normal', row 1 appears at bottom,
%  so the default puts +y at bottom — OPPOSITE to the scan image.
%  Setting flipY=true inverts the row assignment so +y is again at top,
%  consistent with the scan image.
%--------------------------------------------------------------------------

    %----------------------------------------------------------------------
    % Defaults
    %----------------------------------------------------------------------
    if nargin < 1 || isempty(sim)
        sim = defaultParams();
    end
    if nargin < 2 || isempty(outDir)
        outDir = '.';
    end
    if nargin < 3 || isempty(tsh)
        tsh = [0.5, 0.2, 0.05];   % 50 %, 20 %, 5 % of channel peak
    end

    % Increase z resolution for smooth isosurfaces (default nz=5 is too coarse)
    if sim.nz < 11
        sim.nz  = 21;
        sim.z   = linspace(-sim.nzRange/2, sim.nzRange/2, sim.nz);
        sim.obj = beadObject3D(sim);   % recompute bead for new z grid
        fprintf('[makePSF3DAllPixels] nz set to %d for smooth 3-D display\n', sim.nz);
    end

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

    %----------------------------------------------------------------------
    % Detector display grid (flipY=true, flipX=false)
    %   idxGrid(row, col) = physical channel index k in sim.detXY
    %   row 1 → most-negative y (bottom)
    %   col 1 → most-negative x (left)
    %----------------------------------------------------------------------
    flipX        = false;
    flipY        = true;
    transposeGrid = false;
    idxGrid = detectorIndexGridFromXY(sim.detXY, flipX, flipY, transposeGrid);

    %----------------------------------------------------------------------
    % Pre-build 3-D coordinate meshgrids for FocusImage3D (Mode B)
    %
    %   MATLAB convention:  rows = y,  cols = x,  pages = z
    %   meshgrid(x,y,z) → X3D changes along dim-2, Y3D along dim-1
    %
    %   FocusImage3D rendering (ndims(rho)>2 branch):
    %     isosurface( R3D.*cos(Phi3D),  R3D.*sin(Phi3D),  Z3D,  PSF, lvl )
    %   = isosurface( X3D,             Y3D,               Z3D,  PSF, lvl )  ✓
    %----------------------------------------------------------------------
    [X3D, Y3D, Z3D] = meshgrid(sim.x, sim.y, sim.z);
    R3D   = sqrt(X3D.^2 + Y3D.^2);      % radial distance
    Phi3D = atan2(Y3D, X3D);             % azimuthal angle

    nRows = size(idxGrid, 1);
    nCols = size(idxGrid, 2);

    %----------------------------------------------------------------------
    % Main loop – one figure per aberration scenario
    %----------------------------------------------------------------------
    for j = 1:nSc
        label  = scenarios{j,1};
        coeffs = coeffStruct(sim, scenarios{j,2});

        fprintf('[makePSF3DAllPixels] Computing: %s ...\n', strrep(label,'_',' '));
        t0 = tic;

        % ── PSF at excitation wavelength ──────────────────────────────────
        hExc = psfBessel(sim, coeffs, sim.lamExc);

        % ── PSF at emission wavelength (for detection) ────────────────────
        hDet = psfBessel(sim, coeffs, sim.lamEm);

        % ── Explicit detector-pixel collection efficiency ─────────────────
        %    eta(y,x,z,k): fraction of emitted photons collected by pixel k
        %    when the emitter is at sample position (x,y,z)
        eta = detectorCollectionEfficiencyExplicit(sim, hDet);

        fprintf('  PSF + detector computed in %.1f s\n', toc(t0));

        % ── Figure layout ─────────────────────────────────────────────────
        fig = figure('Color', 'w', 'Position', [50 50 1500 1200]);

        titleStr = strrep(label, '_', ' ');
        if j == 1
            ampStr = '(unaberrated)';
        else
            ampStr = '(0.18 waves RMS)';
        end
        sgtitle(sprintf('Effective 3-D PSF per detector pixel  —  %s  %s', ...
            titleStr, ampStr), 'FontSize', 13, 'FontWeight', 'bold');

        % ── 25 subplots in physical detector order ────────────────────────
        for row = 1:nRows
            for col = 1:nCols
                ax = subplot(nRows, nCols, (row-1)*nCols + col);
                k = idxGrid(row, col);     % physical channel index
                if ~isfinite(k) || k < 1 || k > size(sim.detXY,1)
                    axis(ax, 'off');
                    continue;
                end

                % Effective PSF: excitation × detection collection efficiency
                %   h_eff_k(x,y,z) = h_exc(x,y,z) · η_k(x,y,z)
                hEff = hExc .* eta(:,:,:,k);

                % Normalize to peak = 1 for display
                pk = max(hEff(:));
                if pk > 0
                    hEffN = hEff / pk;
                else
                    hEffN = zeros(size(hEff));
                end

                %----------------------------------------------------------
                % FocusImage3D  –  Mode B (ndims(R3D) = 3)
                %
                %   Argument mapping:
                %     rho      → R3D    (r = √(x²+y²) grid, µm)
                %     z        → Z3D    (axial coordinate grid, µm)
                %     vol      → hEffN  (PSF values on same Cartesian grid)
                %     flag     → Phi3D  (φ = atan2(y,x) grid, rad)
                %     tsh      → isosurface fractions
                %     maxangle → 2π     (full azimuthal range)
                %     clfflag  → 1      (suppress lights inside FocusImage3D;
                %                        lights are added below per subplot)
                %
                %   Isosurface rendering:
                %     isosurface( R3D·cos(Phi3D),  R3D·sin(Phi3D),  Z3D, hEffN )
                %   = isosurface( X3D,             Y3D,             Z3D, hEffN )
                %   → correct Cartesian x/y/z in µm  ✓
                %----------------------------------------------------------
                FocusImage3D(R3D, Z3D, hEffN, Phi3D, tsh, 2*pi, 1);

                % Add lighting for this subplot
                light('Style', 'infinite', 'Position', [ 0  1  1]);
                light('Style', 'infinite', 'Position', [ 1 -1 -0.5]);

                % Title: physical detector coordinates (µm)
                dk = sim.detXY(k, :);
                title(sprintf('(%.2f, %.2f)', dk(1), dk(2)), ...
                    'FontSize', 6.5, 'Interpreter', 'none');

                % Compact axes (labels only on outer edges)
                ax = gca;
                ax.FontSize = 5.5;
                if row == nRows && col == ceil(nCols/2)
                    xlabel('x (\mum)', 'FontSize', 6);
                else
                    xlabel('');
                end
                if col == 1 && row == ceil(nRows/2)
                    ylabel('y (\mum)', 'FontSize', 6);
                else
                    ylabel('');
                end
                zlabel('');
            end
        end

        % ── Add a global colormap legend as a text annotation ─────────────
        annotation('textbox', [0.01 0.01 0.35 0.04], ...
            'String', sprintf('Isosurface levels: %.0f%% / %.0f%% / %.0f%% of peak', ...
                max(tsh)*100, median(tsh)*100, min(tsh)*100), ...
            'EdgeColor', 'none', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);

        % ── Save ──────────────────────────────────────────────────────────
        outFile = fullfile(outDir, sprintf('psf3d_%s.png', label));
        exportgraphics(fig, outFile, 'Resolution', 150);
        fprintf('  Saved: %s\n', outFile);
        close(fig);
    end

    fprintf('[makePSF3DAllPixels] All %d figures saved to: %s\n', nSc, outDir);
end
