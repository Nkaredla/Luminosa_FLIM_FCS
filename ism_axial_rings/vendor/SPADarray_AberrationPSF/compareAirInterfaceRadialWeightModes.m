function report = compareAirInterfaceRadialWeightModes(varargin)
%COMPAREAIRINTERFACERADIALWEIGHTMODES Compare corrected and legacy PSF weights.
%
%   report = compareAirInterfaceRadialWeightModes()
%
%   Computes water/glass/air scalar PSFs with:
%       sampleSolidAngle  - current default, sample-side angular weighting
%       immersionDebye    - legacy immersion-side Debye weighting
%
%   and writes XY/XZ/YZ comparison figures. This is a fast core-model
%   diagnostic; it does not load experimental bead data.

    p = inputParser;
    p.FunctionName = 'compareAirInterfaceRadialWeightModes';
    addParameter(p, 'outputDir', '');
    addParameter(p, 'figureVisible', 'off');
    addParameter(p, 'stageZUm', -0.8:0.05:0.8);
    addParameter(p, 'emitterHeightUm', 0);
    addParameter(p, 'objectiveNA', 1.2);
    addParameter(p, 'immersionRefractiveIndex', 1.33);
    addParameter(p, 'glassRefractiveIndex', 1.518);
    addParameter(p, 'sampleRefractiveIndex', 1.0003);
    parse(p, varargin{:});
    opts = p.Results;
    opts.figureVisible = char(opts.figureVisible);

    if isempty(opts.outputDir)
        root = fileparts(fileparts(mfilename('fullpath')));
        opts.outputDir = fullfile(root, 'output_matlab', ...
            'air_interface_radial_weight_comparison');
    else
        opts.outputDir = char(opts.outputDir);
    end
    if exist(opts.outputDir, 'dir') ~= 7
        mkdir(opts.outputDir);
    end

    sim = defaultParams();
    sim.sampleGeometry = 'airOnGlass';
    sim.NA = opts.objectiveNA;
    sim.nImmersion = opts.immersionRefractiveIndex;
    sim.nMedium = sim.nImmersion;
    sim.nGlass = opts.glassRefractiveIndex;
    sim.nSample = opts.sampleRefractiveIndex;
    sim.nDesignGlass = opts.glassRefractiveIndex;
    sim.coverslipThicknessUm = 190;
    sim.designCoverslipThicknessUm = 190;
    sim.nx = 41;
    sim.ny = 41;
    sim.x = linspace(-1.0, 1.0, sim.nx);
    sim.y = sim.x;
    sim.dx = sim.x(2) - sim.x(1);
    sim.Nr = 72;
    sim.Nphi = 144;

    stageZ = double(opts.stageZUm(:)).';
    emitterHeight = double(opts.emitterHeightUm);
    if isscalar(emitterHeight)
        emitterHeight = repmat(emitterHeight, size(stageZ));
    end

    sampleSim = sim;
    sampleSim.interfaceRadialWeightMode = 'sampleSolidAngle';
    legacySim = sim;
    legacySim.interfaceRadialWeightMode = 'immersionDebye';

    samplePSF = scalarPSFBesselAirInterface( ...
        sampleSim, struct(), sampleSim.lamEm, stageZ, emitterHeight);
    legacyPSF = scalarPSFBesselAirInterface( ...
        legacySim, struct(), legacySim.lamEm, stageZ, emitterHeight);
    samplePSF = samplePSF / sum(samplePSF(:));
    legacyPSF = legacyPSF / sum(legacyPSF(:));
    diffPSF = samplePSF - legacyPSF;

    rho = linspace(0,1,sim.Nr);
    quadrature = ones(size(rho));
    quadrature([1 end]) = 0.5;
    quadrature = quadrature * (rho(2)-rho(1));
    [~, weightDiag] = airInterfaceRadialWeights(sampleSim, rho, quadrature);

    report = struct();
    report.sampleSolidAnglePSF = samplePSF;
    report.legacyImmersionDebyePSF = legacyPSF;
    report.differencePSF = diffPSF;
    report.relativeL2Difference = norm(diffPSF(:)) / max(norm(legacyPSF(:)), eps);
    report.maxAbsDifference = max(abs(diffPSF(:)));
    report.radialWeightDiagnostics = weightDiag;
    report.radialWeightProfile = analyzeRadialWeightProfile(sim);
    report.outputDir = opts.outputDir;

    projections = {'xy','xz','yz'};
    for k = 1:numel(projections)
        writeComparisonFigure(samplePSF, legacyPSF, diffPSF, projections{k}, ...
            opts, fullfile(opts.outputDir, sprintf( ...
            'air_interface_%s_sampleSolidAngle_vs_immersionDebye.png', ...
            upper(projections{k}))));
    end
    writeRadialProfileFigure(report.radialWeightProfile, opts, ...
        fullfile(opts.outputDir, 'air_interface_radial_weight_saf_profile.png'));
    save(fullfile(opts.outputDir, 'air_interface_radial_weight_comparison.mat'), ...
        'report', '-v7.3');

    prof = report.radialWeightProfile;
    fprintf('\nAir-interface radial weight comparison\n');
    fprintf('  sample-vs-legacy PSF relative L2 difference: %.4g\n', ...
        report.relativeL2Difference);
    fprintf('  critical angle (air->glass): %.2f deg in glass; critical pupil radius rho_c = %.3f\n', ...
        prof.criticalAngleGlassDeg, prof.criticalPupilRadius);
    fprintf('  SAF pupil area fraction: %.3f\n', ...
        weightDiag.supercriticalAreaFraction);
    fprintf('  SAF radial-weight fraction: sample %.3f, legacy %.3f\n', ...
        prof.safWeightFractionSample, prof.safWeightFractionLegacy);
    fprintf('  sample-vs-legacy radial-weight abs difference: %.4g\n', ...
        weightDiag.sampleVsLegacyRelativeAbsDifference);
    fprintf('  outputs: %s\n\n', opts.outputDir);
end

function prof = analyzeRadialWeightProfile(sim)
%   Rigorous, theory-anchored radial diagnostics for the air-interface pupil.
%   Everything here is exact given the scalar layered model: the critical
%   angle/radius, the Fresnel s/p power transmittances, and the model's own
%   radial weighting recomputed on the pupil grid. These validate that the
%   SAF (supercritical) annulus sits where theory says and that the model's
%   weighting behaves correctly there. The one thing they do NOT assert is
%   the absolute SAF magnitude of a real dipole -- that needs a vectorial
%   reference, exposed as dipoleReference below (a labelled hook, not a
%   fabricated curve), because near the critical angle the emission is
%   strongly p-polarized and the scalar 1/2(s+p) average mis-weights it.
    rho = linspace(0, 1, sim.Nr);
    quad = ones(size(rho));
    quad([1 end]) = 0.5;
    quad = quad * (rho(2) - rho(1));

    [~, dSample] = airInterfaceRadialWeights( ...
        setRadialMode(sim, 'sampleSolidAngle'), rho, quad);
    [~, dLegacy] = airInterfaceRadialWeights( ...
        setRadialMode(sim, 'immersionDebye'), rho, quad);

    q = sim.NA * rho;                          % k_parallel / k0 (invariant)
    nImm = sim.nImmersion; nGlass = sim.nGlass; nSample = sim.nSample;
    cosImm = rootLocal(1 - (q/nImm).^2);
    cosGlass = rootLocal(1 - (q/nGlass).^2);
    cosSample = rootLocal(1 - (q/nSample).^2);
    [tpIG, tsIG] = fresnelTLocal(nImm, cosImm, nGlass, cosGlass);
    [tpGS, tsGS] = fresnelTLocal(nGlass, cosGlass, nSample, cosSample);
    Tp = abs(tpIG .* tpGS).^2;                 % |t_p|^2 through the stack
    Ts = abs(tsIG .* tsGS).^2;                 % |t_s|^2 through the stack

    superMask = q > nSample;                   % supercritical (SAF) annulus
    wSample = abs(dSample.selectedWeight);
    wLegacy = abs(dLegacy.selectedWeight);

    prof = struct();
    prof.rho = rho;
    prof.qOverK0 = q;
    prof.criticalPupilRadius = nSample / sim.NA;
    prof.criticalAngleGlassDeg = real(asind(nSample / nGlass));
    prof.sampleWeight = wSample;
    prof.legacyWeight = wLegacy;
    prof.fresnelTp2 = Tp;
    prof.fresnelTs2 = Ts;
    prof.supercriticalMask = superMask;
    prof.supercriticalAreaFraction = dSample.supercriticalAreaFraction;
    prof.safWeightFractionSample = safFraction(wSample, superMask);
    prof.safWeightFractionLegacy = safFraction(wLegacy, superMask);
    prof.dipoleReference = dipoleBfpReferenceStub();
end

function s = dipoleBfpReferenceStub()
%   Placeholder for a vectorial dipole-at-interface back-focal-plane radiance
%   (Axelrod 2012 / Enderlein), orientation-averaged. Left unfilled on
%   purpose: shipping an unverified prefactor as "ground truth" would be
%   worse than none. Fill .profile (same rho grid) with your trusted
%   convention to overlay it; the figure picks it up automatically.
    s = struct('available', false, 'profile', [], ...
        'note', ['Supply an orientation-averaged Axelrod/Enderlein dipole ' ...
        'BFP radiance on the same rho grid to overlay the vectorial SAF reference.']);
end

function sim = setRadialMode(sim, mode)
    sim.interfaceRadialWeightMode = mode;
end

function [tp, ts] = fresnelTLocal(n1, c1, n2, c2)
    tp = 2*n1*c1 ./ (n2*c1 + n1*c2);
    ts = 2*n1*c1 ./ (n1*c1 + n2*c2);
end

function root = rootLocal(value)
    root = sqrt(complex(value, 0));
    root(imag(root) < 0) = conj(root(imag(root) < 0));
end

function frac = safFraction(weight, superMask)
    total = sum(abs(weight));
    if total <= 0
        frac = 0;
    else
        frac = sum(abs(weight(superMask))) / total;
    end
end

function writeRadialProfileFigure(prof, opts, outFile)
    fig = figure('Visible', opts.figureVisible, 'Color', 'w', ...
        'Position', [80 80 900 520]);
    ax = axes(fig); hold(ax, 'on');
    nrm = @(v) v / max(max(abs(v)), eps);
    plot(ax, prof.rho, nrm(prof.sampleWeight), 'LineWidth', 1.8, ...
        'DisplayName', 'sample-side weight (model default)');
    plot(ax, prof.rho, nrm(prof.legacyWeight), '--', 'LineWidth', 1.4, ...
        'DisplayName', 'legacy immersion weight');
    plot(ax, prof.rho, nrm(prof.fresnelTp2), ':', 'LineWidth', 1.6, ...
        'DisplayName', '|t_p|^2 (p-polarized)');
    plot(ax, prof.rho, nrm(prof.fresnelTs2), ':', 'LineWidth', 1.6, ...
        'DisplayName', '|t_s|^2 (s-polarized)');
    if isstruct(prof.dipoleReference) && prof.dipoleReference.available && ...
            ~isempty(prof.dipoleReference.profile)
        plot(ax, prof.rho, nrm(prof.dipoleReference.profile), 'k-', ...
            'LineWidth', 2, 'DisplayName', 'vectorial dipole reference');
    end
    yl = ylim(ax);
    xc = prof.criticalPupilRadius;
    plot(ax, [xc xc], yl, 'k-.', 'LineWidth', 1, ...
        'DisplayName', sprintf('critical radius \\rho_c=%.2f', xc));
    if xc < 1
        patch(ax, [xc 1 1 xc], [yl(1) yl(1) yl(2) yl(2)], [1 0.9 0.6], ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'SAF annulus');
    end
    hold(ax, 'off');
    xlabel(ax, 'normalized pupil radius \rho'); ylabel(ax, 'normalized radial weight');
    legend(ax, 'Location', 'northwest'); grid(ax, 'on');
    title(ax, sprintf(['Air-interface pupil weighting and SAF annulus ' ...
        '(\\theta_c=%.1f deg, SAF weight frac sample %.2f)'], ...
        prof.criticalAngleGlassDeg, prof.safWeightFractionSample), ...
        'Interpreter', 'tex');
    exportgraphics(fig, outFile, 'Resolution', 180);
    close(fig);
end

function writeComparisonFigure(samplePSF, legacyPSF, diffPSF, planeName, opts, outFile)
    fig = figure('Visible', opts.figureVisible, 'Color', 'w', ...
        'Position', [80 80 1200 900]);
    tl = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    data = {projectOne(samplePSF, planeName), ...
        projectOne(legacyPSF, planeName), projectOne(diffPSF, planeName)};
    titles = {'sampleSolidAngle weighting', 'legacy immersionDebye weighting', ...
        'sampleSolidAngle - legacy'};
    for k = 1:3
        ax = nexttile(tl);
        imagesc(ax, data{k});
        axis(ax, 'image');
        axis(ax, 'off');
        title(ax, titles{k}, 'Interpreter', 'none');
        if k < 3
            colormap(ax, hot);
        else
            lim = max(abs(data{k}(:)), [], 'omitnan');
            if ~isfinite(lim) || lim <= 0, lim = 1; end
            caxis(ax, [-lim lim]);
            colormap(ax, redBlueMap(256));
        end
        colorbar(ax);
    end
    title(tl, sprintf('Air-interface radial weighting: %s projection', ...
        upper(planeName)), 'Interpreter', 'none');
    exportgraphics(fig, outFile, 'Resolution', 180);
    close(fig);
end

function image = projectOne(vol, planeName)
    switch lower(planeName)
        case 'xy'
            image = squeeze(max(vol, [], 3));
        case 'xz'
            image = squeeze(max(vol, [], 1)).';
        case 'yz'
            image = squeeze(max(vol, [], 2)).';
        otherwise
            error('compareAirInterfaceRadialWeightModes:BadProjection', ...
                'Projection must be xy, xz, or yz.');
    end
end

function cmap = redBlueMap(n)
    x = linspace(-1,1,n).';
    cmap = [max(0,-x), 1-abs(x), max(0,x)];
    cmap = max(min(cmap,1),0);
end
