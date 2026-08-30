function calib = make_rt_miet_calibration(cfg)
%MAKE_RT_MIET_CALIBRATION Calibration curve for the MemGlow-on-gold sessions.
%
% calib = make_rt_miet_calibration()
% calib = make_rt_miet_calibration(cfg)
%
% Builds the lifetime-to-height calibration used by the immune-cell MIET
% analyses and writes it beside the code as a two-column calibrationCurve MAT
% plus a diagnostic figure.
%
% PHOTOPHYSICS
%
%   tau_free = 2.0 ns, measured: the axial cross-sections in
%   free_dye_above_surface_summary.csv put the far-field plateau at 1.92-2.01 ns.
%
% The quantum yield is NOT assumed. It is derived from the bare-SLB lifetime,
% which is a genuine calibration point: the bilayer is a uniformly labelled
% single layer at a known small height. Set cfg.anchorSlbLifetimeNs and the
% quantum yield that reproduces it is solved for by fit_miet_slb_anchor.
%
% WHICH COMPONENT IS THE BILAYER: TWO READINGS
%
% A three-exponential tail fit of the bare-SLB pooled decay returns
% 0.329 +- 0.005 / 0.91 +- 0.22 / 2.39 +- 0.23 ns. Which of the first two is
% the supported bilayer is not settled by the fit, and the two readings give
% different photophysics:
%
%   H1  0.33 ns is the bilayer. The anchor then demands qy ~ 0.86 at isotropic
%       orientation and EXCLUDES an in-plane dipole outright.
%   H2  0.33 ns is substrate emission - gold photoluminescence under 640 nm
%       excitation, or dye collapsed against the silica - and 0.91 ns is the
%       bilayer. The anchor then returns qy = 0.34 with an IN-PLANE dipole,
%       which is what a membrane probe is expected to have.
%
% H2 is not a rescue of a failed fit: the independent literature parameters
% (qy 0.30, in-plane, 20 nm spacer, tau_0 2.0-2.085 ns) PREDICT 0.936-0.976 ns
% at z = 4 nm with no fitting, against a measured 0.910 +- 0.217 ns.
%
% An IRF artefact was ruled out quantitatively: the IRF's own tail decays with
% a 0.073-0.101 ns time constant, and the 0.33 ns component is unchanged when
% the fit window starts 1.60 ns after the peak, which is ~17 IRF time
% constants. Whatever the 0.33 ns component is, it is real emission.
%
% Only the QUANTUM YIELD and the free-space lifetime differ between H1 and H2.
% The Bayesian pipeline FIXES its short component and fits the longer one
% freely, so the fitted membrane lifetimes are identical under both readings
% and no refitting is needed - just a different calibration curve. Under H2 the
% bilayer's own emission beneath the cell (about 3% of the photons in the bare
% region) is left unmodelled, which biases the membrane lifetime slightly
% short; that is small next to the difference between the two readings.
%
% Both are settled by one control: a bare chip with spacer, no bilayer and no
% dye. If a 0.33 ns component carrying most of the counts still appears, H2 is
% right.
%
% ORIENTATION IS A CONSTRAINT, NOT A FREE CHOICE
%
% A strictly in-plane dipole cannot be quenched below tau_0/R_parallel, which
% for this stack at 4 nm is 0.418 ns. The measured SLB lifetimes are
% 0.313-0.390 ns, so the parallel assumption is excluded by the data - and so
% is a quantum yield of 0.30 at ANY orientation, since that would need
% <cos^2 theta> above 1. What the anchor leaves free is the pair (orientation,
% quantum yield) along a one-parameter locus; ONE measured lifetime cannot fix
% two unknowns.
%
% The anchor removes the degeneracy only at the SURFACE end. Anchored curves
% for <cos^2 theta> = 1/3 and 1 agree to 0.007 ns in their floor, because both
% are forced through the same measured SLB point. They do NOT agree higher up:
% the uniquely invertible range ends at 161.6 nm for isotropic and 121.7 nm for
% vertical, and on the RT session the median reported height moves by about
% 40 nm between the two. Orientation is therefore a leading systematic on
% absolute height, not a detail, and the two variants are meant to be run and
% compared rather than one of them believed. run_batch_immune_cell_MIET_height_maps
% tags every output with the <cos^2 theta> it used for exactly that reason.
%
% LAYER STACK (bottom to top; z = 0 at the top of the silica spacer, which is
% where the supported bilayer sits)
%
%   glass coverslip  n = 1.52   semi-infinite
%   titanium          2 nm      adhesion layer, index from metals.mat
%   gold             15 nm      index from metals.mat
%   titanium          2 nm      adhesion layer
%   silica spacer    20 nm      n = 1.46
%   water            n = 1.33   the molecule's half-space
%
% This is the stack recorded in Hela_MIET.mat, the only measured stack in the
% repository; it was confirmed by reproducing the existing calibrationCurve.mat
% from those settings.
%
% cfg fields (all optional)
%   tauFreeNs             2.0
%   quantumYield          used only when no SLB anchor is given
%   anchorSlbLifetimeNs   measured bare-SLB lifetime [ns]; when set, the
%                         quantum yield is solved for instead of assumed
%   anchorSlbHeightNm     4, the labelled bilayer treated as one layer
%   meanCosSquaredTheta   1/3 (isotropic wobble). 0 is in-plane, 1 is normal
%   lambdaNm              690, the ACO-ISM emission wavelength of the pipeline
%   goldThicknessNm       15
%   spacerThicknessNm     20
%   titaniumThicknessNm   2
%   zNm                   0.1:0.1:400
%   outputMat             auto-named from the resolved photophysics
%   outputFigure          same base name, .png
%   showFigure            false
%   writeFigure           true

    if nargin < 1 || isempty(cfg); cfg = struct(); end
    cfg = withDefaults(cfg, struct( ...
        'tauFreeNs', 2.0, ...
        'quantumYield', [], ...
        'anchorSlbLifetimeNs', [], ...
        'anchorSlbHeightNm', 4, ...
        'meanCosSquaredTheta', 1/3, ...
        'lambdaNm', 690, ...
        'goldThicknessNm', 15, ...
        'spacerThicknessNm', 20, ...
        'titaniumThicknessNm', 2, ...
        'glassIndex', 1.52, ...
        'spacerIndex', 1.46, ...
        'waterIndex', 1.33, ...
        'zNm', 0.1:0.1:400, ...
        'curveType', 'minimum', ...
        'outputMat', '', ...
        'outputFigure', '', ...
        'showFigure', false, ...
        'writeFigure', true));

    luminosa_miet_vendor_root();
    metals = load(fullfile(luminosa_miet_vendor_root(), 'subroutines', 'metals.mat'));
    goldIndex = metalIndexAt(metals.gold, metals.wavelength, cfg.lambdaNm, 'gold');
    titaniumIndex = metalIndexAt(metals.titan, metals.wavelength, cfg.lambdaNm, 'titanium');

    % Bottommost first, as the kernel wants: substrate, then upwards.
    params = struct();
    params.lambdaNm = cfg.lambdaNm;
    params.n0 = [cfg.glassIndex, titaniumIndex, goldIndex, titaniumIndex, cfg.spacerIndex];
    params.d0 = [cfg.titaniumThicknessNm, cfg.goldThicknessNm, ...
                 cfg.titaniumThicknessNm, cfg.spacerThicknessNm];
    params.n = cfg.waterIndex;
    params.d = 1000;      % irrelevant: n1 == n, so there is no upper interface
    params.n1 = cfg.waterIndex;
    params.d1 = [];
    params.tauFreeNs = cfg.tauFreeNs;
    params.dipoleOrientation = 'order';
    params.meanCosSquaredTheta = cfg.meanCosSquaredTheta;
    params.zNm = cfg.zNm;
    params.curveType = cfg.curveType;

    anchor = [];
    if ~isempty(cfg.anchorSlbLifetimeNs)
        anchorParams = params;
        anchorParams.slbHeightNm = cfg.anchorSlbHeightNm;
        anchorParams.slbLifetimeNs = cfg.anchorSlbLifetimeNs;
        anchorParams.quantumYield = [];
        anchorParams.verbose = false;
        anchor = fit_miet_slb_anchor(anchorParams);
        if ~all(anchor.feasible)
            error('make_rt_miet_calibration:AnchorInfeasible', ...
                ['The SLB anchor demands a quantum yield of %.3f, which is not ' ...
                 'physical. At <cos^2 theta> = %.3f nothing below %.4f ns is ' ...
                 'reachable, but tau_SLB = %.4f ns was measured. Raise ' ...
                 'meanCosSquaredTheta (>= %.3f keeps it physical), or revisit ' ...
                 'the spacer thickness, the SLB height or tau_free.'], ...
                max(anchor.quantumYield), cfg.meanCosSquaredTheta, ...
                cfg.tauFreeNs / anchor.rAtRequestedOrder, ...
                min(cfg.anchorSlbLifetimeNs), anchor.minimumFeasibleOrder);
        end
        cfg.quantumYield = anchor.quantumYield(1);
    elseif isempty(cfg.quantumYield)
        error('make_rt_miet_calibration:QuantumYield', ...
            ['Give either cfg.anchorSlbLifetimeNs, so the quantum yield can be ' ...
             'derived from the SLB, or cfg.quantumYield to assume one.']);
    end
    params.quantumYield = cfg.quantumYield;

    baseDir = fileparts(mfilename('fullpath'));
    if isempty(cfg.outputMat)
        cfg.outputMat = fullfile(baseDir, sprintf( ...
            'calibrationCurve_%dnm_f%03d_qy%03d_tau%gns.mat', ...
            round(cfg.lambdaNm), round(cfg.meanCosSquaredTheta * 100), ...
            round(cfg.quantumYield * 100), cfg.tauFreeNs));
    end
    if isempty(cfg.outputFigure)
        [figDir, figBase] = fileparts(cfg.outputMat);
        cfg.outputFigure = fullfile(figDir, [figBase '.png']);
    end

    orientationText = sprintf('<cos^2 th>=%.3f (theta_eff=%.1f deg)', ...
        cfg.meanCosSquaredTheta, acosd(sqrt(cfg.meanCosSquaredTheta)));
    if isempty(anchor)
        yieldText = sprintf('qy=%.3f assumed', cfg.quantumYield);
    else
        yieldText = sprintf('qy=%.3f from tau_SLB=%.4f ns at z=%g nm', ...
            cfg.quantumYield, cfg.anchorSlbLifetimeNs(1), cfg.anchorSlbHeightNm);
    end
    params.label = sprintf(['MemGlow on %g nm Au / %g nm SiO2, %s, %s, ' ...
        'tau_free=%.2f ns, %g nm emission'], cfg.goldThicknessNm, ...
        cfg.spacerThicknessNm, orientationText, yieldText, ...
        cfg.tauFreeNs, cfg.lambdaNm);

    calib = miet_calibration_curve(params);
    calib.sourceFunction = mfilename();
    calib.config = cfg;
    calib.anchor = anchor;

    % Saved under BOTH names so the existing loaders keep working: the
    % two-column `calibrationCurve` is what load_lifetime_height_calibration
    % looks for, `calibration` carries the provenance.
    calibrationCurve = calib.curve; %#ok<NASGU>
    calibration = calib; %#ok<NASGU>
    save(cfg.outputMat, 'calibrationCurve', 'calibration', '-v7');

    fprintf('make_rt_miet_calibration: %s\n', calib.label);
    fprintf('  usable height range   %.2f - %.2f nm (%d points)\n', ...
        calib.heightNm(1), calib.heightNm(end), numel(calib.heightNm));
    fprintf('  usable lifetime range %.4f - %.4f ns\n', ...
        calib.minLifetimeNs, calib.maxLifetimeNs);
    fprintf('  unambiguous up to     %.4f ns\n', calib.unambiguousLifetimeNs);
    if ~isempty(anchor)
        modelledSlb = interp1(calib.fullHeightNm, calib.fullLifetimeNs, ...
            cfg.anchorSlbHeightNm);
        fprintf(['  SLB anchor            qy %.4f reproduces tau_SLB %.4f ns ' ...
            'at z = %g nm (model gives %.4f ns)\n'], cfg.quantumYield, ...
            cfg.anchorSlbLifetimeNs(1), cfg.anchorSlbHeightNm, modelledSlb);
        fprintf(['  orientation limits    in-plane floor %.4f ns, vertical ' ...
            'floor %.4f ns; physical for <cos^2 th> >= %.3f\n'], ...
            anchor.parallelLimitLifetimeNs, anchor.minAttainableLifetimeNs, ...
            anchor.minimumFeasibleOrder);
    end
    fprintf('  saved                 %s\n', cfg.outputMat);

    if cfg.writeFigure
        writeCalibrationFigure(calib, cfg);
        fprintf('  figure                %s\n', cfg.outputFigure);
    end
end

function value = metalIndexAt(column, wavelengths, lambdaNm, name)
    exact = find(wavelengths == lambdaNm, 1);
    if ~isempty(exact)
        value = column(exact);
        return;
    end
    if lambdaNm < min(wavelengths) || lambdaNm > max(wavelengths)
        error('make_rt_miet_calibration:Wavelength', ...
            'metals.mat covers %g-%g nm; %g nm is outside it (%s).', ...
            min(wavelengths), max(wavelengths), lambdaNm, name);
    end
    value = interp1(double(wavelengths), double(column), lambdaNm, 'linear');
end

function writeCalibrationFigure(calib, cfg)
    visibility = 'off';
    if cfg.showFigure; visibility = 'on'; end
    h = figure('Color', 'w', 'Visible', visibility, 'Position', [80 80 1150 460]);
    layout = tiledlayout(h, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax = nexttile(layout);
    plot(ax, calib.fullHeightNm, calib.fullLifetimeNs, '-', ...
        'Color', [0.75 0.75 0.78], 'LineWidth', 1.2); hold(ax, 'on');
    plot(ax, calib.heightNm, calib.lifetimeNs, '-', ...
        'Color', [0.10 0.35 0.70], 'LineWidth', 2.2);
    yline(ax, calib.params.tauFreeNs, '--', 'free-space \tau_0', ...
        'Color', [0.4 0.4 0.4], 'LabelHorizontalAlignment', 'left');
    xlabel(ax, 'height above the spacer [nm]');
    ylabel(ax, 'lifetime [ns]');
    title(ax, 'MIET calibration curve');
    legend(ax, {'full simulation', 'uniquely invertible'}, 'Location', 'southeast');
    grid(ax, 'on'); box(ax, 'off');
    xlim(ax, [0 max(calib.fullHeightNm)]);

    ax = nexttile(layout);
    plot(ax, calib.lifetimeNs, calib.heightNm, '-', ...
        'Color', [0.10 0.35 0.70], 'LineWidth', 2.2);
    xlabel(ax, 'lifetime [ns]');
    ylabel(ax, 'height [nm]');
    title(ax, 'inverse used for the height maps');
    grid(ax, 'on'); box(ax, 'off');

    % The label is long enough to be clipped by a single-line tiledlayout
    % title, so it is wrapped at commas rather than shrunk to illegibility.
    wrapped = strtrim(split(string(calib.label), ', '));
    lines = strings(0, 1);
    current = "";
    for index = 1:numel(wrapped)
        candidate = current;
        if strlength(candidate) > 0; candidate = candidate + ", "; end
        candidate = candidate + wrapped(index);
        if strlength(candidate) > 72 && strlength(current) > 0
            lines(end + 1, 1) = current; %#ok<AGROW>
            current = wrapped(index);
        else
            current = candidate;
        end
    end
    lines(end + 1, 1) = current;
    title(layout, cellstr(lines), 'FontWeight', 'bold', ...
        'Interpreter', 'none', 'FontSize', 10);
    subtitle(layout, sprintf(['unambiguous up to %.3f ns (%.1f nm); measured ' ...
        'lifetimes above that cannot be inverted and are left blank'], ...
        calib.unambiguousLifetimeNs, calib.maxHeightNm));
    exportgraphics(h, cfg.outputFigure, 'Resolution', 250);
    if ~cfg.showFigure; close(h); end
end

function cfg = withDefaults(cfg, defaults)
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(cfg, names{index}) || isempty(cfg.(names{index}))
            cfg.(names{index}) = defaults.(names{index});
        end
    end
end
