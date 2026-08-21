function report = validateAirInterfaceForwardModel(varargin)
%VALIDATEAIRINTERFACEFORWARDMODEL Numerical checks for the layered PSF core.

    p = inputParser;
    addParameter(p,'homogeneousTolerance',1e-6);
    addParameter(p,'minimumContrastWeightDifference',0.05);
    addParameter(p,'interfaceRadialWeightMode','sampleSolidAngle');
    addParameter(p,'airInterfaceStageMedium','immersion');
    addParameter(p,'writeOutput','');
    parse(p,varargin{:});
    opts = p.Results;

    sim = defaultParams();
    sim.sampleGeometry = 'airOnGlass';
    sim.nx = 19;
    sim.ny = 19;
    sim.x = linspace(-0.8,0.8,sim.nx);
    sim.y = sim.x;
    sim.dx = sim.x(2)-sim.x(1);
    sim.z = [-0.2 0 0.2];
    sim.nz = numel(sim.z);
    sim.Nr = 36;
    sim.Nphi = 72;
    sim.interfaceRadialWeightMode = char(opts.interfaceRadialWeightMode);
    sim.airInterfaceStageMedium = char(opts.airInterfaceStageMedium);

    homogeneous = sim;
    homogeneous.nImmersion = 1.33;
    homogeneous.nGlass = 1.33;
    homogeneous.nSample = 1.33;
    homogeneous.nDesignGlass = 1.33;
    homogeneous.nMedium = 1.33;
    homogeneous.coverslipThicknessUm = 190;
    homogeneous.designCoverslipThicknessUm = 190;
    validateScalarOptics(homogeneous);

    reference = scalarPSFBessel(homogeneous,struct(),homogeneous.lamEm);
    layered = scalarPSFBesselAirInterface(homogeneous,struct(), ...
        homogeneous.lamEm,homogeneous.z,0);
    reference = reference/sum(reference(:));
    layered = layered/sum(layered(:));
    homogeneousRelativeError = norm(reference(:)-layered(:))/ ...
        max(norm(reference(:)),eps);

    air = sim;
    air.nImmersion = 1.33;
    air.nGlass = 1.518;
    air.nSample = 1.0003;
    air.nDesignGlass = 1.518;
    air.nMedium = air.nImmersion;
    air.coverslipThicknessUm = 190;
    air.designCoverslipThicknessUm = 190;
    validateScalarOptics(air);
    interfaceField = scalarPSFBesselAirInterface( ...
        air,struct(),air.lamEm,[0 0],[0 0.16]);
    rho = linspace(0,1,air.Nr);
    quadrature = ones(size(rho));
    if numel(rho) > 1
        quadrature([1 end]) = 0.5;
        quadrature = quadrature*(rho(2)-rho(1));
    end
    [~, radialWeightDiagnostics] = airInterfaceRadialWeights( ...
        air, rho, quadrature);

    finiteGate = all(isfinite(interfaceField(:)));
    nonnegativeGate = all(interfaceField(:)>=0);
    surfacePeak = max(interfaceField(:,:,1),[],'all');
    elevatedPeak = max(interfaceField(:,:,2),[],'all');
    evanescentDecayGate = elevatedPeak < surfacePeak;
    contrastWeightModeGate = strcmpi(radialWeightDiagnostics.mode, ...
        'sampleSolidAngle');
    finiteWeightGate = radialWeightDiagnostics.finiteSelectedWeight;
    supercriticalGate = radialWeightDiagnostics.supercriticalAreaFraction > 0;
    contrastWeightDifferenceGate = ...
        radialWeightDiagnostics.sampleVsLegacyRelativeAbsDifference >= ...
        opts.minimumContrastWeightDifference;

    report = struct();
    report.homogeneousRelativeError = homogeneousRelativeError;
    report.homogeneousLimitGate = ...
        homogeneousRelativeError <= opts.homogeneousTolerance;
    report.finiteGate = finiteGate;
    report.nonnegativeGate = nonnegativeGate;
    report.surfacePeak = surfacePeak;
    report.elevatedPeak = elevatedPeak;
    report.evanescentDecayGate = evanescentDecayGate;
    report.radialWeightDiagnostics = radialWeightDiagnostics;
    report.contrastWeightModeGate = contrastWeightModeGate;
    report.finiteWeightGate = finiteWeightGate;
    report.supercriticalGate = supercriticalGate;
    report.contrastWeightDifferenceGate = contrastWeightDifferenceGate;
    report.allPassed = report.homogeneousLimitGate && finiteGate && ...
        nonnegativeGate && evanescentDecayGate && contrastWeightModeGate && ...
        finiteWeightGate && supercriticalGate && contrastWeightDifferenceGate;

    if ~isempty(opts.writeOutput)
        save(char(opts.writeOutput),'report');
    end
    fprintf('\nAir-interface forward-model validation\n');
    fprintf('  homogeneous relative error: %.4g (pass %d)\n', ...
        homogeneousRelativeError,report.homogeneousLimitGate);
    fprintf('  finite/nonnegative: %d/%d\n',finiteGate,nonnegativeGate);
    fprintf('  surface/elevated peak: %.4g / %.4g (decay pass %d)\n\n', ...
        surfacePeak,elevatedPeak,evanescentDecayGate);
    fprintf(['  radial weighting: %s; SAF area fraction %.3f; ' ...
        'sample-vs-legacy diff %.3f\n'], ...
        radialWeightDiagnostics.mode, ...
        radialWeightDiagnostics.supercriticalAreaFraction, ...
        radialWeightDiagnostics.sampleVsLegacyRelativeAbsDifference);
    fprintf(['  contrast gates: mode %d, finite weights %d, ' ...
        'supercritical %d, contrast difference %d; allPassed %d\n\n'], ...
        contrastWeightModeGate, finiteWeightGate, supercriticalGate, ...
        contrastWeightDifferenceGate, report.allPassed);
end
