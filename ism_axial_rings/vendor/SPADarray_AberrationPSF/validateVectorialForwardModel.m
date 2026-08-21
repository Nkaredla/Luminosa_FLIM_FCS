function report = validateVectorialForwardModel(varargin)
%VALIDATEVECTORIALFORWARDMODEL Numerical checks for the vectorial PSF core.
%
%   report = validateVectorialForwardModel()
%
%   Gates:
%     homogeneousLimitGate     vectorial air-interface with no index contrast
%                              reduces to the homogeneous vectorial model.
%     lowNAScalarGate          at low NA the vectorial PSF approaches the
%                              scalar PSF (longitudinal field negligible).
%     highNAVectorialGate      at high NA the vectorial PSF differs from the
%                              scalar one by a meaningful amount (the
%                              longitudinal field is actually contributing).
%     evanescentDecayGate      air-side SAF/evanescent field decays with
%                              emitter height above the interface.
%     excitationCollectionDirectionGate
%                              excitation and collection use distinct
%                              interface transfer/weighting paths.
%     finiteGate/nonnegativeGate  outputs are finite and nonnegative.

    p = inputParser;
    addParameter(p, 'homogeneousTolerance', 1e-6);
    addParameter(p, 'lowNATolerance', 0.05);
    addParameter(p, 'minHighNADifference', 0.02);
    addParameter(p, 'minimumDirectionDifference', 1e-3);
    addParameter(p, 'interfaceRadialWeightMode', 'sampleSolidAngle');
    addParameter(p, 'airInterfaceStageMedium', 'immersion');
    addParameter(p, 'vectorialPolarizationMode', 'xyAverage');
    addParameter(p, 'excitationPolarizationMode', 'x');
    addParameter(p, 'collectionPolarizationMode', 'xyAverage');
    addParameter(p, 'writeOutput', '');
    parse(p, varargin{:});
    opts = p.Results;

    base = defaultParams();
    base.sampleGeometry = 'airOnGlass';
    base.nx = 19; base.ny = 19;
    base.x = linspace(-0.8, 0.8, base.nx); base.y = base.x;
    base.dx = base.x(2) - base.x(1);
    base.z = [-0.2 0 0.2]; base.nz = numel(base.z);
    base.Nr = 36; base.Nphi = 72;
    base.includesVectorialPolarization = true;
    base.interfaceRadialWeightMode = char(opts.interfaceRadialWeightMode);
    base.airInterfaceStageMedium = char(opts.airInterfaceStageMedium);
    base.vectorialPolarizationMode = char(opts.vectorialPolarizationMode);
    base.excitationPolarizationMode = char(opts.excitationPolarizationMode);
    base.collectionPolarizationMode = char(opts.collectionPolarizationMode);

    % ---- homogeneous-limit: vectorial air-interface -> vectorial homogeneous
    homog = base;
    homog.nImmersion = 1.33; homog.nGlass = 1.33; homog.nSample = 1.33;
    homog.nDesignGlass = 1.33; homog.nMedium = 1.33;
    homog.coverslipThicknessUm = 190; homog.designCoverslipThicknessUm = 190;
    validateScalarOptics(homog);
    vecHomog = vectorialPSFBessel(homog, struct(), homog.lamEm);
    vecLayered = vectorialPSFBesselAirInterface(homog, struct(), ...
        homog.lamEm, homog.z, 0);
    vecHomog = vecHomog / sum(vecHomog(:));
    vecLayered = vecLayered / sum(vecLayered(:));
    homogeneousRelativeError = norm(vecHomog(:) - vecLayered(:)) / ...
        max(norm(vecHomog(:)), eps);

    % ---- low-NA: vectorial approaches scalar -----------------------------
    lowNA = homog; lowNA.NA = 0.4;
    validateScalarOptics(lowNA);
    vecLow = normVol(vectorialPSFBessel(lowNA, struct(), lowNA.lamEm));
    scaLow = normVol(scalarPSFBessel(lowNA, struct(), lowNA.lamEm));
    lowNADifference = norm(vecLow(:) - scaLow(:)) / max(norm(scaLow(:)), eps);

    % ---- high-NA: vectorial differs from scalar (longitudinal field) -----
    highNA = homog; highNA.NA = 1.2;
    validateScalarOptics(highNA);
    vecHigh = normVol(vectorialPSFBessel(highNA, struct(), highNA.lamEm));
    scaHigh = normVol(scalarPSFBessel(highNA, struct(), highNA.lamEm));
    highNADifference = norm(vecHigh(:) - scaHigh(:)) / max(norm(scaHigh(:)), eps);

    % ---- air interface: SAF/evanescent decay with emitter height ---------
    air = base;
    air.nImmersion = 1.33; air.nGlass = 1.518; air.nSample = 1.0003;
    air.nDesignGlass = 1.518; air.nMedium = air.nImmersion;
    air.coverslipThicknessUm = 190; air.designCoverslipThicknessUm = 190;
    validateScalarOptics(air);
    airCollection = air;
    airCollection.airInterfaceDirection = 'collection';
    airCollection = setVectorialPolarizationForPath(airCollection, 'collection');
    interfaceField = vectorialPSFBesselAirInterface( ...
        airCollection, struct(), air.lamEm, [0 0], [0 0.16]);
    surfacePeak = max(interfaceField(:,:,1), [], 'all');
    elevatedPeak = max(interfaceField(:,:,2), [], 'all');

    airExcitation = air;
    airExcitation.airInterfaceDirection = 'excitation';
    airExcitation = setVectorialPolarizationForPath(airExcitation, 'excitation');
    excitationField = vectorialPSFBesselAirInterface( ...
        airExcitation, struct(), air.lamExc, 0, 0.02);
    collectionField = vectorialPSFBesselAirInterface( ...
        airCollection, struct(), air.lamEm, 0, 0.02);
    excitationField = normVol(excitationField);
    collectionField = normVol(collectionField);
    directionDifference = norm(excitationField(:) - collectionField(:)) / ...
        max(norm(collectionField(:)), eps);

    % ---- gates -----------------------------------------------------------
    finiteGate = all(isfinite(interfaceField(:))) && ...
        all(isfinite(excitationField(:))) && all(isfinite(vecHigh(:)));
    nonnegativeGate = all(interfaceField(:) >= 0) && ...
        all(excitationField(:) >= 0) && all(vecHigh(:) >= 0);

    report = struct();
    report.homogeneousRelativeError = homogeneousRelativeError;
    report.homogeneousLimitGate = homogeneousRelativeError <= opts.homogeneousTolerance;
    report.lowNADifference = lowNADifference;
    report.lowNAScalarGate = lowNADifference <= opts.lowNATolerance;
    report.highNADifference = highNADifference;
    report.highNAVectorialGate = highNADifference >= opts.minHighNADifference;
    report.surfacePeak = surfacePeak;
    report.elevatedPeak = elevatedPeak;
    report.evanescentDecayGate = elevatedPeak < surfacePeak;
    report.directionDifference = directionDifference;
    report.excitationCollectionDirectionGate = ...
        directionDifference >= opts.minimumDirectionDifference;
    report.finiteGate = finiteGate;
    report.nonnegativeGate = nonnegativeGate;
    report.allPassed = report.homogeneousLimitGate && report.lowNAScalarGate && ...
        report.highNAVectorialGate && report.evanescentDecayGate && ...
        report.excitationCollectionDirectionGate && ...
        finiteGate && nonnegativeGate;

    if ~isempty(opts.writeOutput)
        save(char(opts.writeOutput), 'report');
    end
    fprintf('\nVectorial forward-model validation\n');
    fprintf('  homogeneous-limit error: %.4g (pass %d)\n', ...
        homogeneousRelativeError, report.homogeneousLimitGate);
    fprintf('  low-NA vs scalar: %.4g (pass %d); high-NA vs scalar: %.4g (pass %d)\n', ...
        lowNADifference, report.lowNAScalarGate, highNADifference, ...
        report.highNAVectorialGate);
    fprintf('  SAF surface/elevated peak: %.4g / %.4g (decay pass %d)\n', ...
        surfacePeak, elevatedPeak, report.evanescentDecayGate);
    fprintf('  excitation/collection normalized-field difference: %.4g (pass %d)\n', ...
        directionDifference, report.excitationCollectionDirectionGate);
    fprintf('  finite/nonnegative: %d/%d; allPassed %d\n\n', ...
        finiteGate, nonnegativeGate, report.allPassed);
end

function v = normVol(v)
    s = sum(v(:));
    if s > 0, v = v / s; end
end
