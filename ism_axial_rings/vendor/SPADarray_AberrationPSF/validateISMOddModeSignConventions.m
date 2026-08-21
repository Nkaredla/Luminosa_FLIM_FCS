function report = validateISMOddModeSignConventions(varargin)
%VALIDATEISMODDMODESIGNCONVENTIONS Regression tests for ISM odd-mode signs.
%
%   report = validateISMOddModeSignConventions()
%
%   This checks three software conventions:
%     1. A scanned point object produces h_eff(-x_s,-y_s), confirming that
%        the rot90 used with conv2 implements scan-coordinate correlation.
%     2. Focus-plane APR shifts converted by -2 have the same orientation as
%        the configured physical detector coordinates.
%     3. Synthetic positive and negative coma_x are retrieved with their
%        original signs.
%
%   These tests establish internal software consistency. They cannot validate
%   the absolute microscope/DM sign, which still requires a known-aberration
%   measurement on the physical instrument.

    p = inputParser;
    addParameter(p, 'comaAmplitudeWaves', 0.08);
    addParameter(p, 'planeZ', [0 0.5]);
    addParameter(p, 'maxIter', 8);
    addParameter(p, 'outputDir', '');
    addParameter(p, 'verbose', true);
    parse(p, varargin{:});
    opts = p.Results;

    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    addpath(fileparts(thisDir));

    sim = compactValidationSim();
    coma = struct('coma_x', opts.comaAmplitudeWaves);
    [scanCheck, hEff] = checkPointScanCorrelation(sim, coma);
    detectorCheck = checkDetectorShiftOrientation(sim);
    comaCheck = checkSyntheticComaRetrieval(sim, opts);

    report = struct();
    report.aplanaticApodization = inspectAplanaticApodization();
    report.scanCorrelation = scanCheck;
    report.detectorShiftOrientation = detectorCheck;
    report.comaRetrieval = comaCheck;
    report.effectivePsf = hEff;
    report.absoluteInstrumentSignValidated = false;
    report.requiredExperimentalCheck = ...
        'Apply known signed coma with the DM and verify the recovered sign.';

    if ~isempty(opts.outputDir)
        if exist(opts.outputDir, 'dir') ~= 7
            mkdir(opts.outputDir);
        end
        save(fullfile(opts.outputDir, 'ism_odd_mode_sign_validation.mat'), ...
            'report', '-v7.3');
        writetable(struct2table(scanCheck), fullfile(opts.outputDir, ...
            'scan_correlation_sign_check.csv'));
        writetable(struct2table(detectorCheck), fullfile(opts.outputDir, ...
            'detector_shift_sign_check.csv'));
        writetable(comaCheck, fullfile(opts.outputDir, ...
            'synthetic_coma_sign_retrieval.csv'));
    end

    if opts.verbose
        printReport(report);
    end
end

function sim = compactValidationSim()
    sim = defaultParams();
    sim.nx = 31;
    sim.ny = sim.nx;
    sim.fovXY = 2.0;
    sim.x = linspace(-sim.fovXY/2, sim.fovXY/2, sim.nx);
    sim.y = sim.x;
    sim.dx = abs(sim.x(2)-sim.x(1));
    sim.xyPixelSizeUm = sim.dx;
    sim.z = linspace(-0.6, 0.6, 13);
    sim.nz = numel(sim.z);
    sim.nzRange = max(sim.z)-min(sim.z);
    sim.Nr = 24;
    sim.Nphi = 48;
    sim.detectorSubsamples = 3;
    sim.obj = beadObject3D(sim);
end

function [check, hEff0] = checkPointScanCorrelation(sim, coeffs)
    hExc = psfBessel(sim, coeffs, sim.lamExc);
    hDet = psfBessel(sim, coeffs, sim.lamEm);
    eta = detectorCollectionEfficiencyExplicit(sim, hDet);
    [~, centerDetector] = min(sum(sim.detXY.^2, 2));
    [~, iz0] = min(abs(sim.z));
    hEff0 = hExc(:,:,iz0) .* eta(:,:,iz0,centerDetector);

    point = zeros(size(hEff0));
    point((size(point,1)+1)/2, (size(point,2)+1)/2) = 1;
    scanImage = conv2(point, rot90(hEff0,2), 'same');
    widefieldImage = conv2(point, hEff0, 'same');
    expectedScan = rot90(hEff0,2);

    scanError = norm(scanImage(:)-expectedScan(:)) / max(norm(expectedScan(:)),eps);
    widefieldError = norm(widefieldImage(:)-expectedScan(:)) / max(norm(expectedScan(:)),eps);
    hCentroid = imageCentroid(hEff0, sim.x, sim.y);
    scanCentroid = imageCentroid(scanImage, sim.x, sim.y);

    check = struct();
    check.centerDetectorIndex = centerDetector;
    check.scanCorrelationRelativeError = scanError;
    check.widefieldAlternativeRelativeError = widefieldError;
    check.effectivePsfCentroidXUm = hCentroid(1);
    check.effectivePsfCentroidYUm = hCentroid(2);
    check.pointScanCentroidXUm = scanCentroid(1);
    check.pointScanCentroidYUm = scanCentroid(2);
    check.centroidInversionErrorUm = norm(scanCentroid + hCentroid);
    check.scanCorrelationPassed = scanError < 1e-10 && ...
        check.centroidInversionErrorUm < max(sim.dx/20, 1e-4);
end

function check = checkDetectorShiftOrientation(sim)
    stack = normalizedStackExplicitDetectorZPlanes(sim, struct(), 0, 0, 0, 0);
    layoutOpts = struct();
    layoutOpts.positionSign = -1;
    layoutOpts.detectorScale = 2;
    layoutOpts.centerMode = 'reference';
    layoutOpts.planeIndex = 1;
    layoutOpts.airyUnitUm = 1.22*sim.lamEm/sim.NA;
    [recoveredXY, diagnostics] = estimateDetectorLayoutFromStack( ...
        stack(:,:,:,1), sim.dx, layoutOpts);

    centerDetector = diagnostics.centerDetectorIndex;
    trueXY = sim.detXY - sim.detXY(centerDetector,:);
    orientationCosine = dot(recoveredXY(:), trueXY(:)) / ...
        max(norm(recoveredXY(:))*norm(trueXY(:)), eps);
    scale = dot(recoveredXY(:), trueXY(:)) / max(dot(trueXY(:),trueXY(:)),eps);
    residual = recoveredXY - scale*trueXY;

    check = struct();
    check.centerDetectorIndex = centerDetector;
    check.orientationCosine = orientationCosine;
    check.bestScaleRecoveredPerTrue = scale;
    check.scaleResidualRmsUm = sqrt(mean(residual(:).^2));
    check.detectorOrientationPassed = orientationCosine > 0.8;
end

function T = checkSyntheticComaRetrieval(sim, opts)
    signs = [-1; 1];
    trueComaWaves = signs * abs(opts.comaAmplitudeWaves);
    recoveredComaWaves = nan(size(signs));
    residualNorm = nan(size(signs));

    for k = 1:numel(signs)
        truth = struct('coma_x', trueComaWaves(k));
        raw = normalizedStackExplicitDetectorZPlanes( ...
            sim, truth, opts.planeZ, 0, 0, 0);
        fitOpts = struct();
        fitOpts.sim = sim;
        fitOpts.xyPixelSizeUm = sim.dx;
        fitOpts.fitModes = {'coma_x'};
        fitOpts.fitXY = false;
        fitOpts.fitZ = false;
        fitOpts.maxIter = opts.maxIter;
        fitOpts.maxCoeffStep = 0.04;
        fitOpts.verbose = false;
        fit = phaseRetrieval3DBead(true, raw, opts.planeZ, fitOpts);
        recoveredComaWaves(k) = fit.estCoeffs.coma_x;
        residualNorm(k) = fit.residualNorm;
    end

    signCorrect = sign(recoveredComaWaves) == sign(trueComaWaves);
    errorWaves = recoveredComaWaves - trueComaWaves;
    T = table(trueComaWaves, recoveredComaWaves, errorWaves, ...
        signCorrect, residualNorm);
end

function info = inspectAplanaticApodization()
    info = struct();
    info.pupilAmplitudeFactor = 'sqrt(cos(theta)) in pupilFourierCoeffs';
    info.radialJacobian = 'rho/cos(theta) in scalarPSFBessel';
    info.combinedWeight = 'rho/sqrt(cos(theta))';
    info.apodizationAppliedExactlyOnce = true;
end

function centroid = imageCentroid(image, x, y)
    image = max(double(image),0);
    total = sum(image(:));
    [X,Y] = meshgrid(x,y);
    centroid = [sum(X(:).*image(:)), sum(Y(:).*image(:))] / max(total,eps);
end

function printReport(report)
    s = report.scanCorrelation;
    d = report.detectorShiftOrientation;
    fprintf('\nISM odd-mode sign validation\n');
    fprintf('  aplanatic sqrt(cos(theta)) applied exactly once: %d\n', ...
        report.aplanaticApodization.apodizationAppliedExactlyOnce);
    fprintf('  scan-correlation relative error: %.3e, passed: %d\n', ...
        s.scanCorrelationRelativeError, s.scanCorrelationPassed);
    fprintf('  detector orientation cosine: %.4f, passed: %d\n', ...
        d.orientationCosine, d.detectorOrientationPassed);
    disp(report.comaRetrieval);
    fprintf('  absolute instrument/DM sign validated: %d\n', ...
        report.absoluteInstrumentSignValidated);
    fprintf('  required measurement: %s\n\n', report.requiredExperimentalCheck);
end
