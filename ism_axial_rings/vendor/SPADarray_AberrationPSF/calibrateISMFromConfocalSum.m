function calibration = calibrateISMFromConfocalSum(stackInput, varargin)
%CALIBRATEISMFROMCONFOCALSUM Two-stage, decoupled ISM calibration from one bead.
%
%   calibration = calibrateISMFromConfocalSum(stackInput)
%
%   Stage A (wavefront, confocal):
%       Sum all 23 detector pixels into a single confocal 3-D bead PSF and
%       retrieve the objective wavefront from that volume alone -- the
%       classic Gustafsson-style 3-D phase-retrieval problem. No detector
%       geometry, no SPAD crosstalk (summing recovers the confocal-equivalent
%       PSF), no scan-registration coupling. The retrieved coma sign is
%       physical (PSF asymmetry), not tied to the detector convention chain.
%
%   Stage B (detector collection efficiency + magnification, at fixed wavefront):
%       For the SAME aberrations, forward-generate the per-detector intensity
%       distribution at the CENTRE point of the bead (the detection PSF on the
%       array, which carries the full aberration) as a function of detector
%       pixel coordinates (set by the sample-space pitch) and match it to the
%       measured centre-point microimage across recorded planes. Solves 23
%       relative collection-efficiency values + sample-space detector pitch
%       (+ per-plane photon scale), with the measured background fixed, and
%       checks whether the focal-plane collection efficiency also describes the
%       defocused planes (collectionEfficiencyPlaneStability).
%
%   The saved calibration {wavefront, collection efficiency, detector pitch,
%   background} applies to any other bead on the same system.
%
%   Forward model: water/glass/air layered vectorial PSF (air-interface).
%
%   Key options (name-value):
%       'fitModes'        wavefront modes (default no defocus; z0 covers it)
%       'darkFile','channelIDs','backgroundMode'   (loader pass-through)
%       'objectiveNA','immersionRefractiveIndex','glassRefractiveIndex',
%       'sampleRefractiveIndex','coverslipThicknessUm'  (optics)
%       'beadDiameterUm','airBeadAxialSamples'           (bead model)
%       'detectorLayout','detectorPitchProfileUm'        (Stage B pitch grid)
%       'hardwarePitchUm' physical SPAD pitch; if given, magnification is
%                         reported as hardwarePitchUm / sampleSpacePitch
%       'outputDir','writeOutputs','verbose'

    if nargin < 1 || isempty(stackInput)
%         stackInput = ['D:\Luminosa\Data\PSF_batch_outputs\ISM_Aberation_72\' ...
%             'xz_yz_plots\test_20260515_144001_frame_alignment.csv'];

%         stackInput = "D:\Luminosa\Data\PSF_batch_outputs\ISM_Aberation2_73\xz_yz_plots\x0_4uW_0_19collar_80mmlens_20260515_155248_frame_alignment.csv";
        stackInput = "D:\Luminosa\Data\PSF_batch_outputs\xz_yz_plots\x0_4uW_0_19collar_80mmlens_20260515_155744_frame_alignment.csv";

    end
    addpath(fileparts(mfilename('fullpath')));
    opts = parseOptions(varargin{:});
    % Accept either a path/numeric stack or an already-loaded data struct. The
    % latter lets an outer loop (e.g. the effective-NA profile) load the frames
    % ONCE and reuse them, instead of re-reading every PTU plane each iteration.
    preloaded = isstruct(stackInput) && isfield(stackInput, 'rawCounts');
    stackRef = stackInputRef(stackInput);
    opts.outputDir = resolveOutputDir(stackRef, opts.outputDir);
    if opts.writeOutputs && exist(opts.outputDir, 'dir') ~= 7
        mkdir(opts.outputDir);
    end

    % ---------------------------------------------------------------- load
    if preloaded
        data = stackInput;                       % reuse in-memory stack; no reload
        if opts.verbose
            fprintf(['[confocal-sum calibration] Using pre-loaded data struct ' ...
                '(%d planes); skipping frame load.\n'], size(data.rawCounts, 4));
        end
    else
        data = loadFullStackISMData(stackInput, ...
            'channelIDs', opts.channelIDs, 'darkFile', opts.darkFile, ...
            'darkScale', opts.darkScale, 'backgroundMode', opts.backgroundMode, ...
            'xyPixelSizeUm', opts.xyPixelSizeUm, 'reuseCache', opts.reuseCache, ...
            'verbose', opts.verbose);
    end
    [ny, nx, nCh, nPlane] = size(data.rawCounts);
    stageZ = data.stageZUm(:).';
    sim = configureSim(data, opts);

    % ============================ STAGE A: confocal wavefront ============
    bgPerPix = reshape(data.backgroundPerPixel, nCh, nPlane);     % kept for Stage B
    alignShifts = resolveAlignmentShifts(stackRef, data);         % previously-computed
    [confocal, goodPlaneMask, frameInfo] = buildConfocalVolume(data, alignShifts, opts);
    [z0Guess, focusIdx] = focusPlaneZ(confocal, stageZ, goodPlaneMask);
    if opts.verbose
        fprintf('\n[confocal-sum calibration] Stage A: wavefront from summed 3-D PSF\n');
        fprintf('  inter-plane shifts applied: %d (%s); max |shift| %.2f px; using %d of %d planes (%dx%d)\n', ...
            frameInfo.alignmentApplied, frameInfo.shiftSource, frameInfo.maxAbsShiftPx, ...
            nnz(goodPlaneMask), nPlane, ny, nx);
        if any(~goodPlaneMask)
            fprintf('  excluded shifted/bad planes: %s\n', num2str(find(~goodPlaneMask)));
        end
        fprintf('  focus guess z0 = %.3f um\n', z0Guess);
    end
    switch lower(opts.wavefrontMethod)
        case 'hanser'
            wavefront = runHanserStageA(confocal, stageZ, goodPlaneMask, z0Guess, sim, opts);
        case 'poissonmle'
            wavefront = runPoissonStageA(data, confocal, stageZ, frameInfo, bgPerPix, ...
                goodPlaneMask, z0Guess, sim, opts);
        otherwise
            wavefront = fitConfocalWavefront(confocal, stageZ, sim, z0Guess, ...
                goodPlaneMask, opts);
    end
    wavefront.excludedPlanes = find(~goodPlaneMask);
    wavefront.frameAlignment = frameInfo;

    if opts.verbose
        printCoeffs('  wavefront', wavefront.coeffs, opts.fitModes, sim.lamRef);
        fprintf('  fitted z0 = %.3f um; lateral (%.3f, %.3f) um; rel. residual %.4f\n', ...
            wavefront.z0Um, wavefront.xyUm(1), wavefront.xyUm(2), ...
            wavefront.relResidual);
        fprintf('  empirical model blur disabled: XY %.0f nm, Z %.0f nm\n', ...
            1000*wavefront.modelBlurXYUm, 1000*wavefront.modelBlurZUm);
    end

    % ===================== STAGE B: detector collection efficiency + pitch =====
    % Per-detector intensity distribution at the fitted bead centre in the
    % detector-summed focal image. This microimage is the detection PSF on the
    % array, so it carries the full aberration (coma included); the scan sum is
    % the exc-det correlation and is centro-symmetric (coma cancels).
    [measPattern, microimageSelection] = centrePointMicroimage( ...
        data.rawCounts, bgPerPix, focusIdx, opts);   % [23 x nPlane]
    if opts.verbose
        fprintf('[confocal-sum calibration] Stage B: per-detector collection efficiency + sample-space pitch\n');
        fprintf(['  centre microimage scan pixel row=%d, col=%d; fitted bead centre ' ...
            'x=%.3f, y=%.3f px (%s)\n'], ...
            microimageSelection.centerPixelYX(1), microimageSelection.centerPixelYX(2), ...
            microimageSelection.centerXY(1), microimageSelection.centerXY(2), ...
            microimageSelection.method);
    end
    % well-registered planes (small alignment shift) for the CE-vs-defocus check
    if isempty(alignShifts)
        smallShiftPlanes = true(1, nPlane);
    else
        smallShiftPlanes = hypot(alignShifts(:,1), alignShifts(:,2)).' <= opts.maxAlignShiftPx;
    end
    detector = fitDetectorCollectionEfficiencyAndPitch(measPattern, data, focusIdx, stageZ, sim, ...
        wavefront.coeffs, wavefront.z0Um, smallShiftPlanes, opts);
    detector.microimageSelection = microimageSelection;

    if opts.verbose
        fprintf('  sample-space detector pitch = %.4f um (%s; residual %.4f)\n', ...
            detector.detectorPitchSampleUm, detector.pitchSource, detector.relResidual);
        if ~isempty(detector.crossCorrPitch)
            fprintf(['    cross-correlation: %.4f um (reassignment scale %.2f); ' ...
                'scan-ruler implies scale %.2f\n'], ...
                detector.crossCorrPitch.detectorPitchSampleUm, ...
                detector.crossCorrPitch.reassignmentScale, ...
                detector.impliedReassignmentScale);
        end
        if ~isempty(detector.scanPitch) && isfinite(detector.scanPitch.detectorPitchSampleUm)
            fprintf(['    scan-ruler: %.4f um from %d scan points; ' ...
                'rotation %.1f deg, anisotropy %.3f, linearity residual %.3f\n'], ...
                detector.scanPitch.detectorPitchSampleUm, ...
                detector.scanPitch.nScanPointsUsed, ...
                detector.scanPitch.rotationDeg, detector.scanPitch.anisotropy, ...
                detector.scanPitch.fitRelResidual);
        end
        if ~isempty(detector.patternMatchPitch)
            fprintf('    pattern-match cross-check: %.4f um\n', ...
                detector.patternMatchPitch.detectorPitchSampleUm);
        end
        if isfinite(detector.magnification)
            fprintf('  net magnification = %.4g x (hardware pitch %.4g um)\n', ...
                detector.magnification, opts.hardwarePitchUm);
        end
        fprintf('  relative collection efficiency range %.3f to %.3f\n', ...
            min(detector.relativeCollectionEfficiency), max(detector.relativeCollectionEfficiency));
        if isfield(detector, 'collectionEfficiencyPlaneStability') && ...
                ~isempty(detector.collectionEfficiencyPlaneStability)
            st = detector.collectionEfficiencyPlaneStability.summary;
            fprintf(['  focal-plane collection-efficiency transfer across z: median corr %.3f, ' ...
                'median RMS diff %.3f, worst RMS diff %.3f\n'], ...
                st.medianCorrelationWithFocalCE, ...
                st.medianRmsDiffFromFocalCE, st.maxRmsDiffFromFocalCE);
        end
    end

    % ============================ assemble calibration ===================
    calibration = struct();
    calibration.version = 1;
    calibration.createdAt = datestr(now, 30);
    calibration.method = ['confocal-sum wavefront + fixed-wavefront ' ...
        'detector collection-efficiency/pitch'];
    calibration.sourceStack = char(stackInput);
    calibration.lambdaRefUm = sim.lamRef;
    calibration.optics = struct('objectiveNA', sim.NA, ...
        'immersionRefractiveIndex', sim.nImmersion, ...
        'glassRefractiveIndex', sim.nGlass, ...
        'sampleRefractiveIndex', sim.nSample, ...
        'coverslipThicknessUm', sim.coverslipThicknessUm, ...
        'sampleGeometry', sim.sampleGeometry, ...
        'interfaceRadialWeightMode', sim.interfaceRadialWeightMode, ...
        'diffractionModel', sim.diffractionModel, ...
        'includesVectorialPolarization', sim.includesVectorialPolarization, ...
        'vectorialPolarizationMode', sim.vectorialPolarizationMode, ...
        'excitationPolarizationMode', sim.excitationPolarizationMode, ...
        'collectionPolarizationMode', sim.collectionPolarizationMode, ...
        'airInterfaceStageMedium', sim.airInterfaceStageMedium);
    calibration.wavefront = wavefront;
    calibration.detector = detector;
    calibration.background = struct('perChannelPerPixel', bgPerPix, ...
        'channelIDs', data.channelIDs(:).', 'source', data.backgroundSource);
    calibration.centerCoordinateFitDefaults = struct( ...
        'initialCoeffs', wavefront.coeffs, ...
        'fitModes', {opts.fitModes}, ...
        'detectorCollectionEfficiency', detector.relativeCollectionEfficiency(:).', ...
        'detectorPitchSampleUm', detector.detectorPitchSampleUm, ...
        'modelBlurSigmaXYUm', 0, ...
        'modelBlurSigmaZUm', 0, ...
        'diffractionModel', sim.diffractionModel, ...
        'interfaceRadialWeightMode', sim.interfaceRadialWeightMode, ...
        'vectorialPolarizationMode', sim.vectorialPolarizationMode, ...
        'excitationPolarizationMode', sim.excitationPolarizationMode, ...
        'collectionPolarizationMode', sim.collectionPolarizationMode, ...
        'airInterfaceStageMedium', sim.airInterfaceStageMedium, ...
        'estimateDetectorLayout', false, ...
        'fitDetectorPitchScale', false);
    calibration.options = opts;
    calibration.outputDir = opts.outputDir;

    if opts.writeOutputs
        writeOutputs(calibration, confocal, wavefront, measPattern, detector, ...
            stageZ, sim, opts);
    end
end

% ========================================================================
function opts = parseOptions(varargin)
    p = inputParser; p.FunctionName = 'calibrateISMFromConfocalSum';
    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'darkFile', 'D:\Luminosa\Data\ISMdark_counts.ptu');
    addParameter(p, 'darkScale', 1);
    addParameter(p, 'backgroundMode', 'auto');
    addParameter(p, 'xyPixelSizeUm', []);
    addParameter(p, 'reuseCache', true);

    addParameter(p, 'objectiveNA', 1.1);
    addParameter(p, 'excitationWavelengthUm', 0.640);
    addParameter(p, 'emissionWavelengthUm', 0.690);
    addParameter(p, 'immersionRefractiveIndex', 1.33);
    addParameter(p, 'glassRefractiveIndex', 1.518);
    addParameter(p, 'sampleRefractiveIndex', 1.0003);
    addParameter(p, 'designGlassRefractiveIndex', 1.518);
    addParameter(p, 'coverslipThicknessUm', 190);
    addParameter(p, 'designCoverslipThicknessUm', 190);
    addParameter(p, 'beadDiameterUm', 0.10);
    addParameter(p, 'airBeadAxialSamples', 5);
    addParameter(p, 'beadSubsamples', [3 3 3]);
    addParameter(p, 'interfaceRadialWeightMode', 'sampleSolidAngle');
    addParameter(p, 'diffractionModel', 'vectorial');
    addParameter(p, 'vectorialPolarizationMode', 'xyAverage');
    addParameter(p, 'excitationPolarizationMode', 'x');
    addParameter(p, 'collectionPolarizationMode', 'xyAverage');
    addParameter(p, 'airInterfaceStageMedium', 'immersion');

    addParameter(p, 'detectorLayout', 'honeycomb23');
    addParameter(p, 'detectorPixelShape', 'hex');
    addParameter(p, 'detectorSubsamples', 5);
    addParameter(p, 'detFillRatio', 1.0);
    addParameter(p, 'detectorPitchProfileUm', 0.10:0.02:0.30);
    addParameter(p, 'hardwarePitchUm', []);
    addParameter(p, 'microimageSelectionSignalFraction', 0.20);
    addParameter(p, 'pitchMethod', 'crosscorr');  % 'crosscorr' | 'scan' | 'patternMatch'
    addParameter(p, 'scanPitchSignalFraction', 0.1);
    addParameter(p, 'hotPixelFactor', 8);
    % the focal-CE-describes-defocus check uses only well-registered planes
    addParameter(p, 'maxAlignShiftPx', 2.0);
    % cross-correlation reassignment layout (reused estimateDetectorLayoutFromStack)
    addParameter(p, 'detectorReassignmentScale', 2);   % position = scale x shift (~2-2.1)
    addParameter(p, 'detectorLayoutPositionSign', -1);
    addParameter(p, 'detectorShiftUpsample', 20);

    % Stage-A wavefront method: 'parametric' uses the configured vectorial
    % forward model. 'hanser' is available as scalar 3-D-PSF -> pupil phase
    % retrieval and Zernike decomposition.
    addParameter(p, 'wavefrontMethod', 'hanser'); %'parametric'
    addParameter(p, 'hanserPadSize', 96);
    addParameter(p, 'hanserIterations', 40);
    addParameter(p, 'hanserPropagatorMode', 'airInterface');  % 'airInterface' | 'homogeneous'
    addParameter(p, 'fitModes', {'astig_x','astig_y','coma_x','coma_y', ...
        'trefoil_x','trefoil_y','spherical', ...
        'secondary_astig_x','secondary_astig_y', ...
        'quadrafoil_x','quadrafoil_y'});
    addParameter(p, 'planeHalfRangeUm', []);      % [] -> use all (good) planes
    % Frame-shift handling for the confocal volume: apply the saved alignment
    % shifts, then exclude planes whose centroid still jumps (failed frames).
    addParameter(p, 'applyAlignmentShifts', true);
    addParameter(p, 'excludeShiftedFrames', false);   % keep all planes; shifts align them
    addParameter(p, 'frameTrajWindow', 7);        % moving-median window [planes]
    addParameter(p, 'frameOutlierFactor', 4);     % robust-sigma multiple
    addParameter(p, 'minShiftOutlierPx', 1.5);    % floor on the exclusion threshold
    addParameter(p, 'maxIter', 12);
    addParameter(p, 'fdCoeff', 0.02);
    addParameter(p, 'fdZ0Um', 0.03);
    addParameter(p, 'fdXYUm', 0.02);
    addParameter(p, 'coeffBoundsWaves', [-0.4 0.4]);
    % Poisson maximum-likelihood refinement (wavefrontMethod='poissonMLE').
    % Fits the air-interface forward model to the RAW summed photon counts
    % by Fisher scoring on the Poisson deviance; initialised by the
    % effective-NA Hanser retrieval. Photon scale, background scale and
    % (optionally) the effective NA are fitted jointly; CRLB error bars
    % come from the Fisher info.
    addParameter(p, 'poissonMaxIter', 30);
    addParameter(p, 'poissonFitEffectiveNA', true);
    addParameter(p, 'fdNA', 0.01);            % finite-difference step for effective NA
    addParameter(p, 'poissonNAMin', 0.6);     % lower bound on the fitted effective NA
    % Empirical blur is intentionally disabled by default. A fitted axial
    % blur can hide an incorrect interface/vector model and produce meaningless
    % values such as a 200+ nm "model blur".
    addParameter(p, 'fitModelBlur', false);
    addParameter(p, 'initialBlurXYUm', 0);
    addParameter(p, 'initialBlurZUm', 0);
    addParameter(p, 'blurBoundsUm', [0 0]);
    addParameter(p, 'fdBlurUm', 0);
    addParameter(p, 'alsIterations', 25);

    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'verbose', true);
    parse(p, varargin{:});
    opts = p.Results;
    if ischar(opts.fitModes) || isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    end
    opts.fitModes = opts.fitModes(:).';
    opts.diffractionModel = char(opts.diffractionModel);
    opts.vectorialPolarizationMode = char(opts.vectorialPolarizationMode);
    opts.excitationPolarizationMode = char(opts.excitationPolarizationMode);
    opts.collectionPolarizationMode = char(opts.collectionPolarizationMode);
    opts.airInterfaceStageMedium = char(opts.airInterfaceStageMedium);
    opts.interfaceRadialWeightMode = char(opts.interfaceRadialWeightMode);
    opts.fitModelBlur = false;
    opts.initialBlurXYUm = 0;
    opts.initialBlurZUm = 0;
    opts.blurBoundsUm = [0 0];
    opts.fdBlurUm = 0;
end

% ------------------------------------------------------------------------
function sim = configureSim(data, opts)
    sim = defaultParams();
    sim.sampleGeometry = 'airOnGlass';
    sim.NA = opts.objectiveNA;
    sim.lamExc = opts.excitationWavelengthUm;
    sim.lamEm = opts.emissionWavelengthUm;
    sim.lamRef = opts.emissionWavelengthUm;
    sim.nImmersion = opts.immersionRefractiveIndex;
    sim.nMedium = sim.nImmersion;
    sim.nGlass = opts.glassRefractiveIndex;
    sim.nSample = opts.sampleRefractiveIndex;
    sim.nDesignGlass = opts.designGlassRefractiveIndex;
    sim.coverslipThicknessUm = opts.coverslipThicknessUm;
    sim.designCoverslipThicknessUm = opts.designCoverslipThicknessUm;
    sim.beadRadius = opts.beadDiameterUm/2;
    sim.beadSubsamples = opts.beadSubsamples;
    sim.airBeadAxialSamples = opts.airBeadAxialSamples;
    sim.interfaceRadialWeightMode = opts.interfaceRadialWeightMode;
    sim.airInterfaceStageMedium = opts.airInterfaceStageMedium;
    sim.diffractionModel = opts.diffractionModel;
    sim.includesVectorialPolarization = usesVectorialOption(opts.diffractionModel);
    sim.vectorialPolarizationMode = opts.vectorialPolarizationMode;
    sim.excitationPolarizationMode = opts.excitationPolarizationMode;
    sim.collectionPolarizationMode = opts.collectionPolarizationMode;
    sim.detectorSubsamples = opts.detectorSubsamples;
    sim.detectorPixelShape = opts.detectorPixelShape;
    sim.detectorImageInverted = true;
    sim.x = data.xUm(:).'; sim.y = data.yUm(:).';
    sim.nx = numel(sim.x); sim.ny = numel(sim.y);
    sim.dx = data.xyPixelSizeUm;
end

function tf = usesVectorialOption(model)
    model = lower(char(model));
    tf = contains(model, 'vector') && ~contains(model, 'scalar');
end

% ------------------------------------------------------------------------
function [z0, idx] = focusPlaneZ(confocal, stageZ, mask)
    axialSignal = squeeze(sum(sum(confocal,1),2)).';
    if nargin >= 3 && ~isempty(mask)
        axialSignal(~mask) = -inf;       % focus must be a kept plane
    end
    [~, idx] = max(axialSignal);
    z0 = stageZ(idx);
end

% ------------------------------------------------------------------------
function [confocal, goodMask, info] = buildConfocalVolume(data, shifts, opts)
    % Detector-summed confocal volume, background-subtracted, with the
    % previously-computed per-plane inter-plane shifts applied (same
    % fourierShift2D convention as the batch alignment). Planes whose centroid
    % still jumps off the smooth trajectory (failed registration) are removed.
    [ny, nx, nCh, nPlane] = size(data.rawCounts);
    bgPerPix = reshape(data.backgroundPerPixel, nCh, nPlane);
    bgSum = sum(bgPerPix, 1);
    hasShift = opts.applyAlignmentShifts && ~isempty(shifts) && ...
        size(shifts,1) == nPlane;
    maxAbs = 0;
    confocal = zeros(ny, nx, nPlane);
    for iz = 1:nPlane
        plane = max(sum(data.rawCounts(:,:,:,iz), 3) - bgSum(iz), 0);
        if hasShift
            sx = shifts(iz,1); sy = shifts(iz,2);
            if isfinite(sx) && isfinite(sy) && (abs(sx) > 0 || abs(sy) > 0)
                plane = max(fourierShift2D(plane, sx, sy), 0);
                maxAbs = max(maxAbs, max(abs(sx), abs(sy)));
            end
        end
        confocal(:,:,iz) = plane;
    end
    info = struct('alignmentApplied', hasShift, 'maxAbsShiftPx', maxAbs, ...
        'shiftSource', shiftSourceLabel(shifts, hasShift), 'shifts', shifts);
    if opts.excludeShiftedFrames
        [goodMask, dinfo] = detectGoodPlanes(confocal, opts);
    else
        goodMask = true(1, nPlane); dinfo = struct();
    end
    info.detect = dinfo;
    info.goodMask = goodMask;
    info.excludedPlanes = find(~goodMask);
end

% ------------------------------------------------------------------------
function [good, info] = detectGoodPlanes(confocal, opts)
    [ny, nx, nPlane] = size(confocal);
    [X, Y] = meshgrid(1:nx, 1:ny);
    cx = nan(1,nPlane); cy = nan(1,nPlane); tot = zeros(1,nPlane);
    for iz = 1:nPlane
        w = confocal(:,:,iz);
        w = max(w - 0.1*max(w(:)), 0);        % suppress noise floor for the centroid
        s = sum(w(:)); tot(iz) = s;
        if s > 0
            cx(iz) = sum(sum(w.*X))/s;
            cy(iz) = sum(sum(w.*Y))/s;
        end
    end
    win = max(3, round(opts.frameTrajWindow));
    cxs = movmedian(cx, win, 'omitnan');
    cys = movmedian(cy, win, 'omitnan');
    r = hypot(cx - cxs, cy - cys);          % centroid jump off the trajectory [px]
    rr = r(isfinite(r));
    sigma = 1.4826 * median(abs(rr - median(rr)));
    thr = max(opts.frameOutlierFactor * sigma, opts.minShiftOutlierPx);
    good = isfinite(r) & (r <= thr) & (tot > 0);
    info = struct('centroidResidualPx', r, 'thresholdPx', thr, ...
        'centroidX', cx, 'centroidY', cy, 'totalSignal', tot);
end

% ------------------------------------------------------------------------
function shifts = resolveAlignmentShifts(stackInput, data)
    % Previously-computed inter-plane shifts: prefer the batch frame_alignment
    % CSV's total_shift_x/y (authoritative, cache-independent), else the shifts
    % already attached to the loaded data.
    shifts = [];
    nPlane = size(data.rawCounts, 4);
    s = char(stackInput);
    csv = '';
    if exist(s, 'file') == 2 && endsWith(lower(s), '.csv')
        csv = s;
    elseif exist(s, 'dir') == 7
        f = dir(fullfile(s, '**', '*frame_alignment.csv'));
        if ~isempty(f)
            [~, i] = max([f.datenum]);
            csv = fullfile(f(i).folder, f(i).name);
        end
    end
    if ~isempty(csv) && exist(csv, 'file') == 2
        try
            T = readtable(csv);
            if all(ismember({'total_shift_x','total_shift_y'}, ...
                    T.Properties.VariableNames)) && height(T) == nPlane
                shifts = [double(T.total_shift_x(:)), double(T.total_shift_y(:))];
            end
        catch
        end
    end
    if isempty(shifts) && isfield(data, 'alignmentShiftXPx') && ...
            numel(data.alignmentShiftXPx) == nPlane && ...
            isfield(data, 'alignmentShiftYPx') && numel(data.alignmentShiftYPx) == nPlane
        shifts = [data.alignmentShiftXPx(:), data.alignmentShiftYPx(:)];
    end
end

function s = shiftSourceLabel(shifts, hasShift)
    if hasShift && ~isempty(shifts)
        s = 'previously-computed shifts';
    else
        s = 'none / zero';
    end
end

% ============================ Stage A ===================================
function wavefront = runHanserStageA(confocal, stageZ, goodMask, z0, sim, opts)
    % Hanser 3-D-PSF -> pupil phase retrieval on the registered, good planes.
    [sel, planeZsel] = selectPlanes(stageZ, z0, goodMask, opts);
    params = struct('NA', sim.NA, 'nMedium', sim.nImmersion, ...
        'wavelengthUm', sim.lamExc, 'dxUm', sim.dx, 'lamRefUm', sim.lamRef, ...
        'padSize', opts.hanserPadSize, 'nIter', opts.hanserIterations, ...
        'propagatorMode', opts.hanserPropagatorMode, ...
        'nImmersion', sim.nImmersion, 'nGlass', sim.nGlass, ...
        'nSample', sim.nSample, 'nDesignGlass', sim.nDesignGlass, ...
        'coverslipUm', sim.coverslipThicknessUm, ...
        'designCoverslipUm', sim.designCoverslipThicknessUm, ...
        'emitterHeightUm', sim.beadBottomHeightUm, ...
        'fitModes', {opts.fitModes}, 'verbose', opts.verbose);
    res = hanserPhaseRetrieval(confocal(:,:,sel), planeZsel - z0, params);

    wavefront = struct();
    wavefront.method = 'hanser';
    wavefront.coeffs = res.coeffs;
    wavefront.fitModes = opts.fitModes;
    wavefront.coeffVectorWaves = res.coeffVector;
    wavefront.coeffNmRms = res.coeffNmRms;
    wavefront.z0Um = z0;
    wavefront.xyUm = [0 0];
    wavefront.modelBlurXYUm = 0;
    wavefront.modelBlurZUm = 0;
    wavefront.globalPhotonScale = 1;
    wavefront.relResidual = res.relResidual;
    wavefront.selectedPlaneIndices = sel;
    wavefront.modelVolume = res.modelStack;       % already in the input frame
    wavefront.lambdaExcUm = sim.lamExc;
    wavefront.lambdaRefUm = sim.lamRef;
    wavefront.hanser = res;
end

function wavefront = fitConfocalWavefront(confocal, stageZ, sim, z0Guess, goodMask, opts)
    % Damped Gauss-Newton least-squares fit of the summed (confocal) volume.
    % Parameters: [coeffs(fitModes), z0, x0, y0, blurXY, blurZ]; global photon
    % scale is profiled out analytically at each evaluation. Only good (non-
    % shifted) planes within the near-focus window are fitted.
    [sel, planeZsel] = selectPlanes(stageZ, z0Guess, goodMask, opts);
    y = confocal(:,:,sel);
    y = y(:);

    nMode = numel(opts.fitModes);
    % p = [coeffs(1:nMode), z0, x0, y0, blurXYum, blurZum]
    if opts.fitModelBlur
        initBlurXY = opts.initialBlurXYUm;
        initBlurZ = opts.initialBlurZUm;
    else
        initBlurXY = 0;
        initBlurZ = 0;
    end
    p = [zeros(1,nMode), z0Guess, 0, 0, initBlurXY, initBlurZ];
    blurStep = opts.fdBlurUm * double(opts.fitModelBlur);   % 0 freezes the blur
    steps = [opts.fdCoeff*ones(1,nMode), opts.fdZ0Um, opts.fdXYUm, opts.fdXYUm, ...
        blurStep, blurStep];
    nParam = numel(p);

    [r, dev, scale] = evalConfocal(p, y, planeZsel, sim, opts);
    reg = 1e-3;
    for it = 1:opts.maxIter
        J = zeros(numel(r), nParam);
        for q = 1:nParam
            pp = p; pp(q) = pp(q) + steps(q);
            pp = clampConfocal(pp, nMode, opts);
            rp = evalConfocal(pp, y, planeZsel, sim, opts);
            J(:,q) = (rp - r) / max(steps(q), eps);
        end
        H = J.'*J; g = J.'*r;
        delta = -(H + reg*diag(diag(H)+eps)) \ g;
        delta = delta(:).';
        accepted = false; ls = 1;
        for s = 1:10
            trial = clampConfocal(p + ls*delta, nMode, opts);
            [rt, dt, st] = evalConfocal(trial, y, planeZsel, sim, opts);
            if dt < dev
                p = trial; r = rt; dev = dt; scale = st; accepted = true; break;
            end
            ls = ls/2;
        end
        if opts.verbose
            fprintf('    [confocal GN] iter %2d: rel.residual %.5f\n', ...
                it, sqrt(dev)/max(norm(y),eps));
        end
        if ~accepted || norm(ls*delta) < 1e-4
            break;
        end
    end

    coeffs = makeCoeffs(opts.fitModes, p(1:nMode));
    wavefront = struct();
    wavefront.coeffs = coeffs;
    wavefront.fitModes = opts.fitModes;
    wavefront.coeffVectorWaves = p(1:nMode);
    wavefront.coeffNmRms = p(1:nMode) * sim.lamRef * 1000;
    wavefront.z0Um = p(nMode+1);
    wavefront.xyUm = p(nMode+2:nMode+3);
    wavefront.modelBlurXYUm = p(nMode+4);
    wavefront.modelBlurZUm = p(nMode+5);
    wavefront.globalPhotonScale = scale;
    wavefront.relResidual = sqrt(dev)/max(norm(y),eps);
    wavefront.selectedPlaneIndices = sel;
    wavefront.lambdaExcUm = sim.lamExc;
    wavefront.lambdaRefUm = sim.lamRef;
end

function [r, dev, scale] = evalConfocal(p, y, planeZsel, sim, opts)
    nMode = numel(opts.fitModes);
    coeffs = makeCoeffs(opts.fitModes, p(1:nMode));
    z0 = p(nMode+1); x0 = p(nMode+2); y0 = p(nMode+3);
    sigXY = p(nMode+4); sigZ = p(nMode+5);
    vol = buildConfocalModel(sim, coeffs, planeZsel - z0, x0, y0, sigXY, sigZ);
    m = vol(:);
    scale = profileScale(y, m);
    r = y - scale*m;
    dev = sum(r.^2);
end

function vol = buildConfocalModel(sim, coeffs, planeZ, x0, y0, sigXYum, sigZum)
    % Air-interface confocal volume, then an empirical effective blur so the
    % model can reach the measured scanned-bead width.
    vol = confocalAirInterfaceVolume(sim, coeffs, planeZ, x0, y0);
    sigXpix = max(sigXYum, 0) / sim.dx;
    if numel(planeZ) > 1
        dz = median(abs(diff(planeZ)));
    else
        dz = 1;
    end
    sigZpl = 0;
    if dz > 0, sigZpl = max(sigZum, 0) / dz; end
    vol = gaussianBlur3(vol, sigXpix, sigXpix, sigZpl);
    total = sum(vol(:));
    if total > 0, vol = vol / total; end
end

function vol = gaussianBlur3(vol, sigXpix, sigYpix, sigZpix)
    if sigXpix > 0
        kx = gaussKernel(sigXpix); vol = convn(vol, reshape(kx,1,[]), 'same');
    end
    if sigYpix > 0
        ky = gaussKernel(sigYpix); vol = convn(vol, reshape(ky,[],1), 'same');
    end
    if sigZpix > 0 && size(vol,3) > 1
        kz = gaussKernel(sigZpix); vol = convn(vol, reshape(kz,1,1,[]), 'same');
    end
end

function k = gaussKernel(sigPix)
    half = max(1, ceil(3*sigPix));
    x = -half:half;
    k = exp(-(x.^2)/(2*sigPix^2));
    k = k / sum(k);
end

function scale = profileScale(y, m)
    den = sum(m.^2);
    if den <= 0, scale = 0; else, scale = max(sum(y.*m)/den, 0); end
end

function vol = confocalAirInterfaceVolume(sim, coeffs, planeZ, x0, y0)
    % Confocal (open-pinhole / excitation) air-interface bead volume: the
    % detector-summed model with the array aperture treated as fully open.
    stageZ = double(planeZ(:).');
    [objectSlices, heights] = airSurfaceBeadSlices(sim);
    nPlane = numel(stageZ); nHeight = numel(heights);
    stageConditions = repelem(stageZ(:), nHeight);
    heightConditions = repmat(heights(:), nPlane, 1);
    simExc = sim;
    simExc.airInterfaceDirection = 'excitation';
    simExc = setVectorialPolarizationForPath(simExc, 'excitation');
    hExc = psfBesselAirInterface(simExc, coeffs, sim.lamExc, ...
        stageConditions, heightConditions);
    vol = zeros(numel(sim.y), numel(sim.x), nPlane);
    for ip = 1:nPlane
        for ih = 1:nHeight
            cond = (ip-1)*nHeight + ih;
            vol(:,:,ip) = vol(:,:,ip) + ...
                conv2(objectSlices(:,:,ih), rot90(hExc(:,:,cond),2), 'same');
        end
    end
    if (x0 ~= 0 || y0 ~= 0)
        for ip = 1:nPlane
            vol(:,:,ip) = max(fourierShift2D(vol(:,:,ip), x0/sim.dx, y0/sim.dx), 0);
        end
    end
    total = sum(vol(:));
    if total > 0, vol = vol/total; end
end

function p = clampConfocal(p, nMode, opts)
    b = opts.coeffBoundsWaves;
    p(1:nMode) = min(max(p(1:nMode), b(1)), b(2));
    if opts.fitModelBlur
        bb = opts.blurBoundsUm;
        p(nMode+4) = min(max(p(nMode+4), bb(1)), bb(2));   % blur XY
        p(nMode+5) = min(max(p(nMode+5), bb(1)), bb(2));   % blur Z
    else
        p(nMode+4:nMode+5) = 0;
    end
end

% ===================== Stage A: Poisson maximum likelihood ===============
function wavefront = runPoissonStageA(data, confocal, stageZ, frameInfo, bgPerPix, ...
        goodMask, z0Guess, sim, opts)
    % Poisson-MLE wavefront. Initialise with the effective-NA Hanser retrieval,
    % then fit the air-interface forward model to the RAW summed photon counts
    % by Fisher scoring on the Poisson deviance. Photon scale, background scale
    % and (optionally) the effective NA are fitted jointly; the Fisher
    % information at the optimum gives CRLB error bars on the coefficients.
    init = runHanserStageA(confocal, stageZ, goodMask, z0Guess, sim, opts);

    [sel, planeZsel] = selectPlanes(stageZ, z0Guess, goodMask, opts);
    [counts, bgPlane] = buildRawConfocalSummed(data, bgPerPix, frameInfo, sel);

    init0 = struct('coeffs', init.coeffVectorWaves, 'z0', z0Guess, 'xy', [0 0]);
    if isfield(init,'effectiveNA') && ~isempty(init.effectiveNA) && isfinite(init.effectiveNA)
        init0.effectiveNA = init.effectiveNA;
    else
        init0.effectiveNA = sim.NA;
    end

    fit = fitConfocalWavefrontPoisson(counts, bgPlane, planeZsel, sim, init0, opts);

    wavefront = struct();
    wavefront.method = 'poissonMLE';
    wavefront.coeffs = makeCoeffs(opts.fitModes, fit.coeffs);
    wavefront.fitModes = opts.fitModes;
    wavefront.coeffVectorWaves = fit.coeffs;
    wavefront.coeffNmRms = fit.coeffs * sim.lamRef * 1000;
    wavefront.coeffStdErrWaves = fit.coeffStdErr;
    wavefront.coeffStdErrNm = fit.coeffStdErr * sim.lamRef * 1000;
    wavefront.z0Um = fit.z0;
    wavefront.xyUm = fit.xy;
    wavefront.modelBlurXYUm = 0;
    wavefront.modelBlurZUm = 0;
    wavefront.effectiveNA = fit.effectiveNA;
    wavefront.photonScale = fit.photonScale;
    wavefront.globalPhotonScale = fit.photonScale;   % figure/back-compat
    wavefront.backgroundScale = fit.backgroundScale;
    wavefront.relResidual = fit.relResidual;
    wavefront.poissonDeviance = fit.deviance;
    wavefront.reducedDeviance = fit.reducedDeviance;
    wavefront.logLikelihood = fit.logLik;
    wavefront.selectedPlaneIndices = sel;
    wavefront.modelVolume = fit.sourceVolume;        % s*m (source only) [counts]
    wavefront.lambdaExcUm = sim.lamExc;
    wavefront.lambdaRefUm = sim.lamRef;
    wavefront.hanserInit = init;
    wavefront.naScan = init.naScan;
    if opts.verbose
        fprintf('  [poisson MLE] effective NA %.3f, photon scale %.3g, bg scale %.3f\n', ...
            fit.effectiveNA, fit.photonScale, fit.backgroundScale);
        fprintf('  [poisson MLE] deviance %.4g (reduced %.3f), rel.residual %.4f\n', ...
            fit.deviance, fit.reducedDeviance, fit.relResidual);
    end
end

function [counts, bgPlane] = buildRawConfocalSummed(data, bgPerPix, frameInfo, sel)
    % Raw detector-summed counts on the kept planes, registered with the same
    % inter-plane shifts as the confocal volume but WITHOUT background
    % subtraction or clamping to zero, so the data keep their Poisson statistics.
    % The flat per-plane background level is returned separately for the model.
    [ny, nx, ~, ~] = size(data.rawCounts);
    bgSum = sum(bgPerPix, 1);                       % [1 nPlane] per-pixel bg level
    shifts = frameInfo.shifts;
    hasShift = frameInfo.alignmentApplied && ~isempty(shifts);
    counts = zeros(ny, nx, numel(sel));
    bgPlane = zeros(1, numel(sel));
    for k = 1:numel(sel)
        iz = sel(k);
        plane = sum(data.rawCounts(:,:,:,iz), 3);   % sum of independent Poissons
        if hasShift
            sx = shifts(iz,1); sy = shifts(iz,2);
            if isfinite(sx) && isfinite(sy) && (sx ~= 0 || sy ~= 0)
                plane = fourierShift2D(plane, sx, sy);   % registration shift
            end
        end
        counts(:,:,k) = max(plane, 0);              % counts are non-negative
        bgPlane(k) = bgSum(iz);
    end
end

function fit = fitConfocalWavefrontPoisson(counts, bgPlane, planeZsel, sim, init0, opts)
    % Fisher-scoring (IRLS) Poisson maximum-likelihood fit of the air-interface
    % forward model. mu_i = s*m_i(theta) + beta*bg_i ; minimise the Poisson
    % deviance. Parameter vector:
    %   [coeffs(1:nMode), z0, x0, y0, (effNA), logS, beta]
    nMode = numel(opts.fitModes);
    [ny, nx, nPl] = size(counts);
    n  = counts(:);
    bg = reshape(repmat(bgPlane(:).', ny*nx, 1), ny, nx, nPl);
    bg = bg(:);

    fitNA = opts.poissonFitEffectiveNA;
    iZ0 = nMode+1; iX = nMode+2; iY = nMode+3;
    if fitNA, iNA = nMode+4; iLogS = nMode+5; iBeta = nMode+6;
    else,     iNA = 0;        iLogS = nMode+4; iBeta = nMode+5;
    end
    nParam = iBeta;

    NAmin = opts.poissonNAMin; NAmax = sim.NA;
    effNA = min(max(init0.effectiveNA, NAmin), NAmax);

    p = zeros(1, nParam);
    p(1:nMode) = init0.coeffs(:).';
    p(iZ0) = init0.z0; p(iX) = init0.xy(1); p(iY) = init0.xy(2);
    if fitNA, p(iNA) = effNA; end

    % seed the scale parameters from the initial model shape
    m0 = modelShape(p, planeZsel, sim, opts, nMode, iZ0, iX, iY, iNA, effNA);
    s0 = max(sum(n) - sum(bg), eps) / max(sum(m0(:)), eps);
    p(iLogS) = log(max(s0, eps));
    p(iBeta) = 1;

    stepPSF = [opts.fdCoeff*ones(1,nMode), opts.fdZ0Um, opts.fdXYUm, opts.fdXYUm];
    if fitNA, stepPSF(iNA) = opts.fdNA; end
    psfIdx = 1:iY; if fitNA, psfIdx = [1:iY, iNA]; end

    [mu, mv, dev] = evalPoisson(p, planeZsel, sim, opts, nMode, ...
        iZ0, iX, iY, iNA, iLogS, iBeta, n, bg, effNA);
    reg = 1e-3;
    for it = 1:opts.poissonMaxIter
        [J, ~] = poissonJacobian(p, mv, planeZsel, sim, opts, nMode, ...
            iZ0, iX, iY, iNA, iLogS, iBeta, bg, effNA, stepPSF, psfIdx, ...
            NAmin, NAmax, nParam);
        w  = 1 ./ sqrt(max(mu, eps));
        rP = (n - mu) .* w;
        Jw = J .* w;
        H = Jw.'*Jw; g = Jw.'*rP;
        delta = (H + reg*diag(diag(H)+eps)) \ g;
        delta = delta(:).';
        accepted = false; ls = 1;
        for ssearch = 1:12
            trial = clampPoisson(p + ls*delta, nMode, iNA, iBeta, NAmin, NAmax, opts);
            [mut, mvt, devt] = evalPoisson(trial, planeZsel, sim, opts, nMode, ...
                iZ0, iX, iY, iNA, iLogS, iBeta, n, bg, effNA);
            if devt < dev
                p = trial; mu = mut; mv = mvt; dev = devt; accepted = true; break;
            end
            ls = ls/2;
        end
        if opts.verbose
            fprintf('    [poisson IRLS] iter %2d: deviance %.6g\n', it, dev);
        end
        if ~accepted || norm(ls*delta) < 1e-5
            break;
        end
    end

    % Fisher information at the optimum -> CRLB on the coefficients (waves)
    Jw = poissonJacobian(p, mv, planeZsel, sim, opts, nMode, ...
        iZ0, iX, iY, iNA, iLogS, iBeta, bg, effNA, stepPSF, psfIdx, ...
        NAmin, NAmax, nParam) ./ sqrt(max(mu, eps));
    FIM = Jw.'*Jw;
    C = pinv(FIM);
    dC = max(diag(C), 0);
    stdErr = sqrt(dC(1:nMode)).';

    nObs = numel(n); dof = max(nObs - nParam, 1);
    fit = struct();
    fit.coeffs = p(1:nMode);
    fit.coeffStdErr = stdErr;
    fit.z0 = p(iZ0); fit.xy = [p(iX) p(iY)];
    if fitNA, fit.effectiveNA = p(iNA); else, fit.effectiveNA = effNA; end
    fit.photonScale = exp(p(iLogS));
    fit.backgroundScale = p(iBeta);
    fit.deviance = dev;
    fit.reducedDeviance = dev / dof;
    fit.logLik = poissonLogLik(n, mu);
    src = fit.photonScale * mv;
    fit.sourceVolume = reshape(src, ny, nx, nPl);
    meas = max(n - p(iBeta)*bg, 0);
    fit.relResidual = sqrt(sum((meas - src).^2)) / max(sqrt(sum(meas.^2)), eps);
end

function [J, mv] = poissonJacobian(p, mv, planeZsel, sim, opts, nMode, ...
        iZ0, iX, iY, iNA, iLogS, iBeta, bg, effNA, stepPSF, psfIdx, NAmin, NAmax, nParam)
    % d(mu)/d(theta): finite difference on the model shape for the PSF
    % parameters, analytic for the photon scale (logS) and background (beta).
    s = exp(p(iLogS));
    J = zeros(numel(mv), nParam);
    for q = psfIdx
        h = stepPSF(q); if h == 0, h = 1e-3; end
        pp = p; pp(q) = pp(q) + h;
        pp = clampPoisson(pp, nMode, iNA, iBeta, NAmin, NAmax, opts);
        dq = pp(q) - p(q);
        if abs(dq) < 1e-12, continue; end           % step killed by a bound
        mq = modelShape(pp, planeZsel, sim, opts, nMode, iZ0, iX, iY, iNA, effNA);
        J(:,q) = s * (mq(:) - mv) / dq;
    end
    J(:, iLogS) = s * mv;                            % mu = exp(logS)*m + beta*bg
    J(:, iBeta) = bg;
end

function [mu, mv, dev] = evalPoisson(p, planeZsel, sim, opts, nMode, ...
        iZ0, iX, iY, iNA, iLogS, iBeta, n, bg, effNA)
    m  = modelShape(p, planeZsel, sim, opts, nMode, iZ0, iX, iY, iNA, effNA);
    mv = m(:);
    s = exp(p(iLogS)); beta = p(iBeta);
    mu = max(s*mv + beta*bg, 1e-9);
    dev = poissonDeviance(n, mu);
end

function m = modelShape(p, planeZsel, sim, opts, nMode, iZ0, iX, iY, iNA, effNAfixed)
    coeffs = makeCoeffs(opts.fitModes, p(1:nMode));
    z0 = p(iZ0); x0 = p(iX); y0 = p(iY);
    if iNA > 0, NA = p(iNA); else, NA = effNAfixed; end
    s2 = sim; s2.NA = NA;
    m = confocalAirInterfaceVolume(s2, coeffs, planeZsel - z0, x0, y0);  % sums to 1
end

function p = clampPoisson(p, nMode, iNA, iBeta, NAmin, NAmax, opts)
    b = opts.coeffBoundsWaves;
    p(1:nMode) = min(max(p(1:nMode), b(1)), b(2));
    if iNA > 0, p(iNA) = min(max(p(iNA), NAmin), NAmax); end
    p(iBeta) = max(p(iBeta), 0);                     % background scale >= 0
end

function D = poissonDeviance(n, mu)
    nn = max(n, 0);
    term = -(nn - mu);
    pos = nn > 0;
    term(pos) = term(pos) + nn(pos) .* log(nn(pos) ./ mu(pos));
    D = 2 * sum(term);
end

function L = poissonLogLik(n, mu)
    nn = max(n, 0);
    L = sum(nn .* log(max(mu, 1e-12)) - mu);         % up to the constant -log(n!)
end

% ============================ Stage B ===================================
function detector = fitDetectorCollectionEfficiencyAndPitch(measPattern, data, focusIdx, stageZ, ...
        sim, coeffs, z0, smallShiftPlanes, opts)
    % Detector pitch first (model-free scan ruler by default), then per-channel
    % collection efficiency + per-plane photon scale by ALS at that fixed pitch.
    relZ = stageZ - z0;
    % Cross-correlation reassignment is the primary estimate; the factor-free
    % scan-ruler is always computed as a cross-check on the x2 assumption.
    scanPitch = estimateDetectorPitchFromScan(data, focusIdx, opts);
    crossPitch = []; gridPitch = [];
    if any(strcmpi(opts.pitchMethod, {'crosscorr','crosscorrelation'}))
        crossPitch = crossCorrPitch(data, focusIdx, sim, opts);
    end
    if any(strcmpi(opts.pitchMethod, {'patternmatch'}))
        gridPitch = patternMatchPitch(measPattern, sim, coeffs, relZ, opts);
    end
    switch lower(opts.pitchMethod)
        case {'crosscorr','crosscorrelation'}
            usedPitch = crossPitch.detectorPitchSampleUm;
            src = 'cross-correlation reassignment';
        case 'scan'
            usedPitch = scanPitch.detectorPitchSampleUm;
            src = 'scan-ruler';
        case 'patternmatch'
            usedPitch = gridPitch.detectorPitchSampleUm;
            src = 'pattern-match';
        otherwise
            error('calibrateISMFromConfocalSum:BadPitchMethod', ...
                'pitchMethod must be ''crosscorr'', ''scan'', or ''patternMatch''.');
    end
    if ~isfinite(usedPitch) || usedPitch <= 0
        error('calibrateISMFromConfocalSum:BadPitch', ...
            'Detector pitch estimate failed (%s).', src);
    end

    M = perDetectorModel(sim, coeffs, relZ, usedPitch, opts);
    [ceRaw, s, sse] = alsRankOne(measPattern, M, opts.alsIterations);
    meanCE = mean(ceRaw(ceRaw > 0));
    if ~isfinite(meanCE) || meanCE <= 0, meanCE = 1; end

    mag = NaN;
    if ~isempty(opts.hardwarePitchUm) && isfinite(opts.hardwarePitchUm) && ...
            opts.hardwarePitchUm > 0
        mag = opts.hardwarePitchUm / usedPitch;
    end

    detector = struct();
    detector.detectorPitchSampleUm = usedPitch;
    detector.pitchSource = src;
    detector.crossCorrPitch = crossPitch;
    detector.scanPitch = scanPitch;
    detector.patternMatchPitch = gridPitch;
    if ~isempty(crossPitch) && isfinite(scanPitch.detectorPitchSampleUm) && ...
            crossPitch.detectorPitchSampleUm > 0
        % effective reassignment scale implied by the factor-free cross-check
        detector.impliedReassignmentScale = opts.detectorReassignmentScale * ...
            scanPitch.detectorPitchSampleUm / crossPitch.detectorPitchSampleUm;
    else
        detector.impliedReassignmentScale = NaN;
    end
    detector.relativeCollectionEfficiency = (ceRaw(:) / meanCE).';        % normalized to mean 1
    detector.channelIDs = data.channelIDs(:).';
    detector.planePhotonScale = s(:).' * meanCE;
    detector.magnification = mag;
    detector.hardwarePitchUm = opts.hardwarePitchUm;
    detector.detectorLayout = opts.detectorLayout;
    detector.relResidual = sqrt(sse) / max(norm(measPattern(:)), eps);
    detector.modelPattern = M;
    detector.measuredPattern = measPattern;
    detector.collectionEfficiencyPlaneStability = ...
        collectionEfficiencyPlaneStability(measPattern, M, detector.relativeCollectionEfficiency(:), ...
        focusIdx, stageZ, smallShiftPlanes);
end

function info = crossCorrPitch(data, focusIdx, sim, opts)
    % Detector pitch from the average reassignment shift between detector scan
    % images, via FFT phase-correlation (reuses estimateDetectorLayoutFromStack).
    % shift_k = d_k * sigma_exc^2/(sigma_exc^2+sigma_det^2) ~ d_k/2, so the
    % position scale (~2) is the one model assumption; the scan-ruler validates it.
    layoutOpts = struct( ...
        'positionSign', opts.detectorLayoutPositionSign, ...
        'detectorScale', opts.detectorReassignmentScale, ...
        'centerMode', 'reference', ...
        'upsampleReg', opts.detectorShiftUpsample, ...
        'airyUnitUm', 1.22*sim.lamEm / max(min(sim.NA, sim.nSample), eps));
    [detXY, diagnostics] = estimateDetectorLayoutFromStack( ...
        data.rawCounts(:,:,:,focusIdx), data.xyPixelSizeUm, layoutOpts);
    info = struct();
    info.method = 'cross-correlation reassignment shifts';
    info.detXY = detXY;
    info.detectorPitchSampleUm = medianNearestPitch(detXY);
    info.reassignmentScale = opts.detectorReassignmentScale;
    info.planeIndex = focusIdx;
    info.diagnostics = diagnostics;
end

function pitch = medianNearestPitch(xy)
    n = size(xy,1);
    nearest = nan(n,1);
    for k = 1:n
        d = hypot(xy(:,1)-xy(k,1), xy(:,2)-xy(k,2));
        d(k) = inf;
        nearest(k) = min(d);
    end
    nearest = nearest(isfinite(nearest) & nearest > 0);
    if isempty(nearest), pitch = NaN; else, pitch = median(nearest); end
end

function info = estimateDetectorPitchFromScan(data, focusIdx, opts)
    % Model-free detector pitch from the calibrated scan ruler.
    % In descanned ISM the descanned image of a fixed bead moves on the array
    % by exactly the (negative) scan displacement in sample-equivalent units.
    % So the per-scan-position detector-microimage centroid, expressed in
    % honeycomb-layout units (nearest-neighbour = 1), shifts linearly with scan
    % position: d(centroid)/d(scan_um) = A = -inv(D), where D maps layout units
    % to sample um. The pitch (and rotation/anisotropy vs the scan axes) reads
    % off D = -inv(A). The slope is aberration-independent -- coma/astig only
    % add a constant centroid offset, which differentiation removes.
    raw = double(data.rawCounts(:,:,:,focusIdx));            % [ny nx 23]
    [nyy, nxx, nCh] = size(raw);
    bg = reshape(data.backgroundPerPixel(:,:,:,focusIdx), 1, 1, nCh);
    raw = max(raw - bg, 0);

    bgVec = reshape(data.backgroundPerPixel(:,:,:,focusIdx), nCh, 1);
    medBg = median(bgVec);
    includeCh = bgVec <= opts.hotPixelFactor * max(medBg, eps);   % drop hot pixels
    if nnz(includeCh) < 6, includeCh = true(nCh,1); end

    L = detectorLayout(opts.detectorLayout, 1);              % [23 x 2], nn = 1
    Lx = reshape(L(includeCh,1), 1, 1, []);
    Ly = reshape(L(includeCh,2), 1, 1, []);
    rawInc = raw(:,:,includeCh);
    total = sum(rawInc, 3);                                   % [ny nx]
    Cx = sum(rawInc .* Lx, 3) ./ max(total, eps);             % centroid, layout units
    Cy = sum(rawInc .* Ly, 3) ./ max(total, eps);

    thr = opts.scanPitchSignalFraction * max(total(:));
    idx = find(total > max(thr, 0));
    if numel(idx) < 4
        [~, order] = sort(total(:), 'descend');
        idx = order(1:min(12, numel(order)));
    end

    dx = data.xyPixelSizeUm;
    [SX, SY] = meshgrid((0:nxx-1)*dx, (0:nyy-1)*dx);          % scan position [um]
    w = total(idx);
    G = [SX(idx), SY(idx), ones(numel(idx),1)];
    WG = w(:) .* G;
    coefX = (G.' * WG) \ (G.' * (w(:) .* Cx(idx)));           % [a11 a12 bx]
    coefY = (G.' * WG) \ (G.' * (w(:) .* Cy(idx)));           % [a21 a22 by]
    A = [coefX(1) coefX(2); coefY(1) coefY(2)];

    % Linearity of centroid-vs-scan is the built-in check that this is a
    % descanned geometry with the assumed channel->layout mapping. A small
    % residual = trustworthy pitch; a large one flags a wrong assumption.
    resid = [Cx(idx) - G*coefX; Cy(idx) - G*coefY];
    spread = [Cx(idx) - mean(Cx(idx)); Cy(idx) - mean(Cy(idx))];
    fitRelResidual = norm(resid) / max(norm(spread), eps);

    info = struct();
    info.method = 'scan-ruler centroid';
    info.fitRelResidual = fitRelResidual;
    info.nScanPointsUsed = numel(idx);
    info.includedChannels = find(includeCh).';
    info.centroidPerScanMatrix = A;
    info.scanPixelSizeUm = dx;
    info.scanCentroidLayoutX = Cx;       % [ny nx] for the shift-map figure
    info.scanCentroidLayoutY = Cy;
    info.scanX = SX;
    info.scanY = SY;
    info.scanMaskIdx = idx;
    info.scanTotalSignal = total;
    info.affineCoefX = coefX;
    info.affineCoefY = coefY;
    if abs(det(A)) < eps
        info.detectorPitchSampleUm = NaN;
        info.latticeVectorsUm = nan(2);
        info.anisotropy = NaN; info.rotationDeg = NaN;
        warning('calibrateISMFromConfocalSum:DegenerateScanPitch', ...
            'Scan-ruler centroid fit is degenerate; check signal coverage.');
        return;
    end
    D = -inv(A);                                             % layout unit -> sample um
    sv = svd(D);
    info.latticeVectorsUm = D;
    info.detectorPitchSampleUm = sqrt(sv(1) * sv(2));        % geometric-mean pitch
    info.anisotropy = sv(1) / max(sv(2), eps);
    [U, ~, Vt] = svd(D); Rm = U * Vt.';
    info.rotationDeg = atan2d(Rm(2,1), Rm(1,1));
end

function info = patternMatchPitch(measPattern, sim, coeffs, relZ, opts)
    % Optional model-dependent cross-check: pitch that best matches the
    % per-detector pattern shape (degenerate with PSF width -- see notes).
    candidates = opts.detectorPitchProfileUm(:).';
    candidates = candidates(isfinite(candidates) & candidates > 0);
    bestSSE = inf; bestPitch = NaN;
    profile = nan(numel(candidates), 2);
    for c = 1:numel(candidates)
        M = perDetectorModel(sim, coeffs, relZ, candidates(c), opts);
        [~, ~, sse] = alsRankOne(measPattern, M, opts.alsIterations);
        profile(c,:) = [candidates(c), sse];
        if sse < bestSSE, bestSSE = sse; bestPitch = candidates(c); end
    end
    info = struct('detectorPitchSampleUm', bestPitch, ...
        'pitchProfile', profile, 'sse', bestSSE);
end

function M = perDetectorModel(sim, coeffs, relZ, pitch, opts)
    % Per-detector intensity distribution at the CENTRE point of the bead
    % (scanning focus on the bead) = the detection PSF eta_k sampled by the
    % array. Carries the full aberration; the scan-summed signal would be the
    % exc-det correlation (centro-symmetric, coma cancels) instead.
    sim.detXY = detectorLayout(opts.detectorLayout, pitch);
    sim.nDet = size(sim.detXY, 1);
    sim.detPitch = pitch;
    sim.detSize = opts.detFillRatio * pitch;
    sim.detectorHexRadius = sim.detSize / sqrt(3);
    stack = normalizedStackAirInterfaceZPlanes(sim, coeffs, relZ, 0, 0, 0);
    [~, xc] = min(abs(sim.x));            % bead-centre scan pixel (x = y = 0)
    [~, yc] = min(abs(sim.y));
    M = reshape(stack(yc, xc, :, :), size(sim.detXY,1), numel(relZ));
end

function [P, info] = centrePointMicroimage(rawCounts, bgPerPix, focusIdx, opts)
    % The 23 detector values at the bead-centre scan pixel from a Gaussian fit
    % seeded by the centre of mass in the detector-summed focal image.
    [~, ~, nCh, nPlane] = size(rawCounts);
    focusIdx = min(max(round(focusIdx), 1), nPlane);
    bgFocus = reshape(bgPerPix(:, focusIdx), 1, 1, nCh);
    focusPlane = max(double(rawCounts(:,:,:,focusIdx)) - bgFocus, 0);
    sumImage = sum(focusPlane, 3);
    [centerXY, fitInfo] = gaussianCenterFromCom(sumImage, ...
        opts.microimageSelectionSignalFraction);
    centerYX = [min(max(round(centerXY(2)), 1), size(sumImage, 1)), ...
        min(max(round(centerXY(1)), 1), size(sumImage, 2))];
    yc = centerYX(1);
    xc = centerYX(2);

    P = zeros(nCh, nPlane);
    for z = 1:nPlane
        vals = reshape(rawCounts(yc, xc, :, z), nCh, 1);
        P(:, z) = max(vals - bgPerPix(:, z), 0);
    end

    info = struct();
    info.method = fitInfo.method;
    info.focusIndex = focusIdx;
    info.centerXY = centerXY;
    info.centerPixelYX = centerYX;
    info.fit = fitInfo;
end

function [centerXY, info] = centerOfMassXY(img, thresholdFraction)
    positive = max(double(img) - min(img(:)), 0);
    threshold = thresholdFraction * max(positive(:));
    weights = positive;
    weights(weights < threshold) = 0;
    if sum(weights(:)) <= 0
        [~, idx] = max(img(:));
        [cy, cx] = ind2sub(size(img), idx);
        centerXY = [cx, cy];
        info = struct('method', 'summed focus image peak fallback');
    else
        [yy, xx] = ndgrid(1:size(img,1), 1:size(img,2));
        mass = sum(weights(:));
        centerXY = [sum(xx(:).*weights(:)), sum(yy(:).*weights(:))] / mass;
        info = struct('method', 'summed focus image center of mass', ...
            'mass', mass, 'threshold', threshold);
    end
end

function [centerXY, info] = gaussianCenterFromCom(img, thresholdFraction)
    [seedXY, comInfo] = centerOfMassXY(img, thresholdFraction);
    img = double(img);
    [ny, nx] = size(img);
    cx0 = min(max(seedXY(1), 1), nx);
    cy0 = min(max(seedXY(2), 1), ny);
    radius = max(4, ceil(min([ny nx]) / 6));
    xlo = max(1, floor(cx0 - radius));
    xhi = min(nx, ceil(cx0 + radius));
    ylo = max(1, floor(cy0 - radius));
    yhi = min(ny, ceil(cy0 + radius));
    patch = img(ylo:yhi, xlo:xhi);
    [yy, xx] = ndgrid(ylo:yhi, xlo:xhi);

    bg0 = median(patch(:));
    amp0 = max(patch(:)) - bg0;
    if ~isfinite(amp0) || amp0 <= 0
        centerXY = seedXY;
        info = comInfo;
        info.method = 'summed focus image center of mass fallback';
        return;
    end

    sigma0 = max(1.5, radius / 2);
    p0 = [bg0, amp0, cx0, cy0, log(sigma0), log(sigma0)];
    objective = @(p) gaussianSse(p, xx, yy, patch);
    options = optimset('Display', 'off', 'MaxIter', 200, 'MaxFunEvals', 800);
    p = fminsearch(objective, p0, options);
    centerXY = [min(max(p(3), 1), nx), min(max(p(4), 1), ny)];
    info = struct('method', 'summed focus image Gaussian fit from COM seed', ...
        'centerOfMassXY', seedXY, 'fitSigmaXY', exp(p(5:6)), ...
        'fitAmplitude', p(2), 'fitBackground', p(1), ...
        'fitWindowYX', [ylo yhi xlo xhi]);
end

function sse = gaussianSse(p, xx, yy, patch)
    sx = max(exp(p(5)), 0.5);
    sy = max(exp(p(6)), 0.5);
    model = p(1) + p(2) * exp(-0.5 * (((xx - p(3)) / sx).^2 + ((yy - p(4)) / sy).^2));
    residual = model - patch;
    sse = sum(residual(:).^2);
    if ~isfinite(sse)
        sse = realmax;
    end
end

function [CE, s, sse] = alsRankOne(P, M, nIter)
    % Solve P_k(z) ~ CE_k * s(z) * M_k(z) by alternating least squares, where
    % CE_k is the per-detector collection efficiency and s(z) the per-plane
    % photon scale. Weighting by M is implicit (the model carries the shape).
    [nCh, nZ] = size(P);
    s = ones(1, nZ);
    CE = ones(nCh, 1);
    for it = 1:nIter
        % update CE_k = sum_z P_k M_k s_z / sum_z (M_k s_z)^2
        Ms = M .* s;                          % [nCh x nZ]
        num = sum(P .* Ms, 2);
        den = sum(Ms.^2, 2);
        CE = num ./ max(den, eps);
        CE = max(CE, 0);
        % update s_z = sum_k P_k M_k CE_k / sum_k (M_k CE_k)^2
        MQ = M .* CE;                         % [nCh x nZ]
        nums = sum(P .* MQ, 1);
        dens = sum(MQ.^2, 1);
        s = nums ./ max(dens, eps);
        s = max(s, 0);
    end
    fit = (CE * s) .* M;
    sse = sum((P(:) - fit(:)).^2);
end

function stability = collectionEfficiencyPlaneStability(P, M, allCE, focusIdx, stageZ, goodPlanes)
    % Estimate relative detector collection efficiency independently at each z
    % plane and compare it with the focal-plane estimate. If the efficiency is
    % truly detector-only, the normalized per-plane estimates are constant with
    % z; structured changes indicate optical-model mismatch or axial dependence.
    % Only well-registered planes (goodPlanes, small alignment shift) are used.
    [nCh, nZ] = size(P);
    if nargin < 4 || isempty(focusIdx) || focusIdx < 1 || focusIdx > nZ
        [~, focusIdx] = max(sum(P, 1));
    end
    if nargin < 5 || isempty(stageZ)
        stageZ = 1:nZ;
    end
    stageZ = stageZ(:).';
    if nargin < 6 || isempty(goodPlanes)
        goodPlanes = true(1, nZ);
    end
    goodPlanes = logical(goodPlanes(:).');

    perPlaneCollectionEfficiency = nan(nCh, nZ);
    corrFocal = nan(1, nZ);
    rmsFocal = nan(1, nZ);
    corrAll = nan(1, nZ);
    rmsAll = nan(1, nZ);
    relResidualAllCE = nan(1, nZ);
    nUsed = zeros(1, nZ);
    allCE = normalizePositiveVector(allCE(:));

    for iz = 1:nZ
        if ~goodPlanes(iz) && iz ~= focusIdx   % keep focus as the reference
            continue;
        end
        pz = P(:,iz);
        mz = M(:,iz);
        include = isfinite(pz) & isfinite(mz) & mz > 0.01*max(mz);
        if nnz(include) < max(6, round(0.4*nCh))
            include = isfinite(pz) & isfinite(mz) & mz > 0;
        end
        nUsed(iz) = nnz(include);
        if nnz(include) < 3
            continue;
        end

        q = nan(nCh,1);
        q(include) = pz(include) ./ max(mz(include), eps);
        q = normalizePositiveVector(q);
        perPlaneCollectionEfficiency(:,iz) = q;

        pred = allCE .* mz;
        alpha = profileScale(pz(include), pred(include));
        relResidualAllCE(iz) = norm(pz(include) - alpha*pred(include)) / ...
            max(norm(pz(include)), eps);
    end

    focalCE = perPlaneCollectionEfficiency(:,focusIdx);
    if all(~isfinite(focalCE))
        focalCE = allCE;
    end
    focalCE = normalizePositiveVector(focalCE);
    for iz = 1:nZ
        q = perPlaneCollectionEfficiency(:,iz);
        include = isfinite(q) & isfinite(focalCE);
        if nnz(include) >= 3
            corrFocal(iz) = corrNoToolbox(q(include), focalCE(include));
            rmsFocal(iz) = sqrt(mean((q(include) - focalCE(include)).^2));
        end
        include = isfinite(q) & isfinite(allCE);
        if nnz(include) >= 3
            corrAll(iz) = corrNoToolbox(q(include), allCE(include));
            rmsAll(iz) = sqrt(mean((q(include) - allCE(include)).^2));
        end
    end

    T = table((1:nZ).', stageZ(:), nUsed(:), corrFocal(:), rmsFocal(:), ...
        corrAll(:), rmsAll(:), relResidualAllCE(:), ...
        'VariableNames', {'planeIndex','stageZUm','nDetectorsUsed', ...
        'corrWithFocalCE','rmsDiffFromFocalCE', ...
        'corrWithAllPlaneCE','rmsDiffFromAllPlaneCE', ...
        'relResidualUsingAllPlaneCE'});

    finiteRms = rmsFocal(isfinite(rmsFocal));
    finiteCorr = corrFocal(isfinite(corrFocal));
    stability = struct();
    stability.method = ['Independent per-plane relative collection efficiency ' ...
        '= measured/model, mean-normalized; compared to the focal-plane and ' ...
        'all-plane collection efficiency.'];
    stability.focusPlaneIndex = focusIdx;
    stability.stageZUm = stageZ;
    stability.goodPlanes = goodPlanes;
    stability.nGoodPlanesUsed = nnz(goodPlanes);
    stability.perPlaneCollectionEfficiency = perPlaneCollectionEfficiency;
    stability.focalPlaneCollectionEfficiency = focalCE;
    stability.allPlaneCollectionEfficiency = allCE;
    stability.table = T;
    stability.summary = struct( ...
        'medianCorrelationWithFocalCE', medianOrNaN(finiteCorr), ...
        'medianRmsDiffFromFocalCE', medianOrNaN(finiteRms), ...
        'maxRmsDiffFromFocalCE', maxOrNaN(finiteRms), ...
        'medianRelResidualUsingAllPlaneCE', ...
            medianOrNaN(relResidualAllCE(isfinite(relResidualAllCE))));
end

function q = normalizePositiveVector(q)
    q = double(q(:));
    good = isfinite(q) & q > 0;
    m = mean(q(good));
    if ~isfinite(m) || m <= 0
        m = 1;
    end
    q = q / m;
end

function c = corrNoToolbox(a, b)
    a = a(:) - mean(a(:));
    b = b(:) - mean(b(:));
    d = sqrt(sum(a.^2) * sum(b.^2));
    if d <= 0
        c = NaN;
    else
        c = sum(a.*b) / d;
    end
end

function v = medianOrNaN(x)
    if isempty(x), v = NaN; else, v = median(x); end
end

function v = maxOrNaN(x)
    if isempty(x), v = NaN; else, v = max(x); end
end

% ============================ shared ====================================
function [sel, planeZsel] = selectPlanes(stageZ, z0Guess, goodMask, opts)
    goodMask = logical(goodMask(:).');
    inRange = true(size(stageZ));
    if ~isempty(opts.planeHalfRangeUm)
        inRange = abs(stageZ - z0Guess) <= opts.planeHalfRangeUm;
    end
    sel = find(goodMask & inRange);
    if numel(sel) < 5                          % don't starve the fit
        cand = find(goodMask);
        if numel(cand) < 5, cand = 1:numel(stageZ); end
        sel = cand;
    end
    planeZsel = stageZ(sel);
end

function c = makeCoeffs(modes, vals)
    c = struct();
    for i = 1:numel(modes), c.(modes{i}) = vals(i); end
end

function printCoeffs(label, coeffs, modes, lamRef)
    fprintf('%s:\n', label);
    for k = 1:numel(modes)
        v = 0; if isfield(coeffs, modes{k}), v = coeffs.(modes{k}); end
        fprintf('     %-10s %+.4f waves (%+.1f nm RMS)\n', modes{k}, v, v*lamRef*1000);
    end
end

function outDir = resolveOutputDir(stackInput, outDir)
    if ~isempty(outDir), return; end
    root = fileparts(fileparts(mfilename('fullpath')));
    [~, stem] = fileparts(char(stackInput));
    if isempty(stem), stem = 'preloaded_stack'; end
    outDir = fullfile(root, 'output_matlab', 'confocal_sum_calibration', stem);
end

function ref = stackInputRef(stackInput)
    % String reference for naming/alignment lookup. For an already-loaded data
    % struct, fall back to its source CSV (the inter-plane shifts are already
    % attached to the struct, so an empty reference is fine).
    if isstruct(stackInput)
        ref = '';
        if isfield(stackInput, 'alignmentCsv') && ~isempty(stackInput.alignmentCsv)
            ref = char(stackInput.alignmentCsv);
        elseif isfield(stackInput, 'sourceMat') && ~isempty(stackInput.sourceMat)
            ref = char(stackInput.sourceMat);
        end
    elseif isnumeric(stackInput)
        ref = '';
    else
        ref = char(stackInput);
    end
end

function writeOutputs(calibration, confocal, wavefront, measPattern, detector, ...
        stageZ, sim, opts)
    outDir = calibration.outputDir;
    save(fullfile(outDir, 'confocal_sum_calibration.mat'), 'calibration', '-v7.3');

    modes = opts.fitModes(:);
    coeffWaves = wavefront.coeffVectorWaves(:);
    writetable(table(modes, coeffWaves, coeffWaves*sim.lamRef*1000, ...
        'VariableNames', {'mode','waves','nmRms'}), ...
        fullfile(outDir, 'wavefront_coefficients.csv'));
    writetable(table(detector.channelIDs(:), detector.relativeCollectionEfficiency(:), ...
        'VariableNames', {'detectorIndex','relativeCollectionEfficiency'}), ...
        fullfile(outDir, 'detector_relative_collection_efficiency.csv'));
    if isfield(detector, 'collectionEfficiencyPlaneStability') && ...
            ~isempty(detector.collectionEfficiencyPlaneStability)
        writetable(detector.collectionEfficiencyPlaneStability.table, ...
            fullfile(outDir, 'detector_collection_efficiency_plane_stability.csv'));
        ceNames = cellstr(strcat('detector_', string(detector.channelIDs(:).')));
        cePlane = array2table(detector.collectionEfficiencyPlaneStability.perPlaneCollectionEfficiency.', ...
            'VariableNames', matlab.lang.makeValidName(ceNames));
        cePlane = addvars(cePlane, ...
            detector.collectionEfficiencyPlaneStability.stageZUm(:), ...
            'Before', 1, 'NewVariableNames', 'stageZUm');
        writetable(cePlane, fullfile(outDir, ...
            'detector_relative_collection_efficiency_by_plane.csv'));
    end
    if ~isempty(detector.patternMatchPitch) && ...
            ~isempty(detector.patternMatchPitch.pitchProfile)
        writetable(array2table(detector.patternMatchPitch.pitchProfile, ...
            'VariableNames', {'detectorPitchSampleUm','sse'}), ...
            fullfile(outDir, 'detector_pitch_profile.csv'));
    end

    figs = { ...
        @() writeSummaryFigure(confocal, wavefront, measPattern, detector, ...
            stageZ, sim, fullfile(outDir, 'confocal_sum_calibration_summary.png')), ...
        @() writePsfProjectionFigure(confocal, wavefront, sim, stageZ, ...
            fullfile(outDir, 'stageA_psf_fit_projections.png')), ...
        @() writeCollectionEfficiencyMapFigure(detector, opts, ...
            fullfile(outDir, 'stageB_detector_collection_efficiency_map.png')), ...
        @() writeCollectionEfficiencyDiagnosticFigure(detector, opts, ...
            fullfile(outDir, 'stageB_expected_measured_collection_efficiency.png')), ...
        @() writeShiftMapFigure(detector, ...
            fullfile(outDir, 'stageB_scan_shift_map.png'))};
    if isfield(wavefront, 'hanser') && ~isempty(wavefront.hanser)
        figs{end+1} = @() writeHanserFigure(wavefront.hanser, ...
            fullfile(outDir, 'stageA_hanser_pupil.png'));
    end
    for f = 1:numel(figs)
        try
            figs{f}();
        catch err
            warning('calibrateISMFromConfocalSum:Figure', ...
                'A diagnostic figure failed: %s', err.message);
        end
    end
end

function writeHanserFigure(res, outFile)
    fig = figure('Visible','off','Color','w','Position',[60 60 980 420]);
    tl = tiledlayout(fig, 1, 2, 'Padding','compact','TileSpacing','compact');
    ax = nexttile(tl);
    ph = res.pupilPhase; ph(~res.aperture) = NaN;
    imagesc(ax, ph, 'AlphaData', res.aperture); axis(ax,'image'); set(ax,'YDir','normal');
    colormap(ax, parula); colorbar(ax); title(ax,'retrieved pupil phase [rad]');
    ax = nexttile(tl);
    plot(ax, res.errorHistory, 'o-'); grid(ax,'on');
    xlabel(ax,'iteration'); ylabel(ax,'amplitude error'); title(ax,'convergence');
    parts = cell(1, numel(res.fitModes));
    for k = 1:numel(res.fitModes)
        parts{k} = sprintf('%s %+.3f', res.fitModes{k}, res.coeffVector(k));
    end
    title(tl, sprintf('Hanser pupil retrieval: %s [waves]', strjoin(parts, '  ')), ...
        'FontWeight','bold');
    exportgraphics(fig, outFile, 'Resolution', 160); close(fig);
end

function writePsfProjectionFigure(confocal, wavefront, sim, stageZ, outFile)
    % Measured vs fitted confocal 3-D PSF, three orthogonal max-projections.
    sel = wavefront.selectedPlaneIndices;
    Vmeas = confocal(:,:,sel);
    z0 = wavefront.z0Um; xy = wavefront.xyUm;
    if isfield(wavefront, 'modelVolume') && ~isempty(wavefront.modelVolume)
        Vmodel = wavefront.modelVolume;            % Hanser retrieved-pupil model
        for p = 1:size(Vmodel,3)                   % scale per plane for display
            sm = sum(sum(Vmeas(:,:,p))); sf = sum(sum(Vmodel(:,:,p)));
            if sf > 0, Vmodel(:,:,p) = Vmodel(:,:,p) * sm / sf; end
        end
    else
        Vmodel = wavefront.globalPhotonScale * buildConfocalModel( ...
            sim, wavefront.coeffs, stageZ(sel) - z0, xy(1), xy(2), ...
            wavefront.modelBlurXYUm, wavefront.modelBlurZUm);
    end
    names = {'XY','XZ','YZ'};
    fig = figure('Visible','off','Color','w','Position',[60 60 1080 980]);
    tl = tiledlayout(fig, 3, 3, 'Padding','compact','TileSpacing','compact');
    for r = 1:3
        Pm = projectVol(Vmeas, names{r});
        Pf = projectVol(Vmodel, names{r});
        cmax = max([Pm(:); Pf(:)]); if ~isfinite(cmax) || cmax <= 0, cmax = 1; end
        ax = nexttile(tl); imagesc(ax, Pm, [0 cmax]); axis(ax,'image'); set(ax,'YDir','normal');
        colormap(ax, hot); colorbar(ax); ylabel(ax, names{r}, 'FontWeight','bold');
        if r==1, title(ax,'measured'); end
        ax = nexttile(tl); imagesc(ax, Pf, [0 cmax]); axis(ax,'image'); set(ax,'YDir','normal');
        colormap(ax, hot); colorbar(ax);
        if r==1, title(ax,'fitted (air-interface)'); end
        ax = nexttile(tl); Dd = Pm - Pf; lim = max(abs(Dd(:)));
        if ~isfinite(lim) || lim <= 0, lim = 1; end
        imagesc(ax, Dd, [-lim lim]); axis(ax,'image'); set(ax,'YDir','normal');
        colormap(ax, redBlue(256)); colorbar(ax);
        if r==1, title(ax,'measured - fitted'); end
    end
    title(tl, sprintf(['Stage A: confocal 3-D PSF fit  ' ...
        '(z0=%.3f um, rel.residual %.3f)\n%s'], ...
        z0, wavefront.relResidual, coefficientSummary(wavefront)), ...
        'FontWeight','bold');
    exportgraphics(fig, outFile, 'Resolution', 160); close(fig);
end

function txt = coefficientSummary(wavefront)
    parts = cell(1, numel(wavefront.fitModes));
    for k = 1:numel(wavefront.fitModes)
        m = wavefront.fitModes{k};
        parts{k} = sprintf('%s %+.3f', m, getc(wavefront.coeffs, m));
    end
    txt = strjoin(parts, '  ');
    if numel(txt) > 180
        txt = [txt(1:180) ' ...'];
    end
    txt = [txt ' [waves]'];
end

function writeCollectionEfficiencyMapFigure(detector, opts, outFile)
    detXY = detectorLayout(opts.detectorLayout, detector.detectorPitchSampleUm);
    fig = figure('Visible','off','Color','w','Position',[80 80 980 430]);
    tl = tiledlayout(fig, 1, 2, 'Padding','compact','TileSpacing','compact');
    ax = nexttile(tl);
    plotDetectorHexMap(detXY, detector.relativeCollectionEfficiency, 'Parent', ax);
    colormap(ax, parula); colorbar(ax);
    title(ax,'relative collection efficiency map (mean 1)');
    ax = nexttile(tl);
    bar(ax, detector.relativeCollectionEfficiency); grid(ax,'on'); yline(ax, 1, 'k--');
    xlabel(ax,'detector index'); ylabel(ax,'relative collection efficiency');
    title(ax,'per-channel collection efficiency');
    title(tl, sprintf('Detector collection efficiency   (pitch %.4f um, %s)', ...
        detector.detectorPitchSampleUm, detector.pitchSource), 'FontWeight','bold');
    exportgraphics(fig, outFile, 'Resolution', 160); close(fig);
end

function writeShiftMapFigure(detector, outFile)
    sp = detector.scanPitch;
    if isempty(sp) || ~isfield(sp,'scanCentroidLayoutX'), return; end
    idx = sp.scanMaskIdx;
    sx = sp.scanX(idx); sy = sp.scanY(idx);
    cx = sp.scanCentroidLayoutX(idx); cy = sp.scanCentroidLayoutY(idx);
    fig = figure('Visible','off','Color','w','Position',[60 60 1240 420]);
    tl = tiledlayout(fig, 1, 3, 'Padding','compact','TileSpacing','compact');

    ax = nexttile(tl);
    quiver(ax, sx(:), sy(:), cx(:)-mean(cx), cy(:)-mean(cy), 0, 'LineWidth', 0.8);
    axis(ax,'image'); grid(ax,'on'); set(ax,'YDir','reverse');
    xlabel(ax,'scan x [um]'); ylabel(ax,'scan y [um]');
    title(ax,'detector-image centroid shift vs scan');

    ax = nexttile(tl);
    plot(ax, sx(:), cx(:), '.'); hold(ax,'on');
    xf = [min(sx) max(sx)];
    bx = median(cx(:) - sp.affineCoefX(1)*sx(:));
    plot(ax, xf, sp.affineCoefX(1)*xf + bx, 'r-', 'LineWidth', 1.4); grid(ax,'on');
    xlabel(ax,'scan x [um]'); ylabel(ax,'centroid x [layout units]');
    title(ax, sprintf('slope %.3f/um  ->  %.3f um', sp.affineCoefX(1), ...
        1/max(abs(sp.affineCoefX(1)),eps)));

    ax = nexttile(tl);
    plot(ax, sy(:), cy(:), '.'); hold(ax,'on');
    yf = [min(sy) max(sy)];
    by = median(cy(:) - sp.affineCoefY(2)*sy(:));
    plot(ax, yf, sp.affineCoefY(2)*yf + by, 'r-', 'LineWidth', 1.4); grid(ax,'on');
    xlabel(ax,'scan y [um]'); ylabel(ax,'centroid y [layout units]');
    title(ax, sprintf('slope %.3f/um  ->  %.3f um', sp.affineCoefY(2), ...
        1/max(abs(sp.affineCoefY(2)),eps)));

    title(tl, sprintf(['Scan-ruler pitch: %.4f um  (linearity residual %.3f, ' ...
        'rotation %.1f deg, anisotropy %.3f)'], sp.detectorPitchSampleUm, ...
        sp.fitRelResidual, sp.rotationDeg, sp.anisotropy), 'FontWeight','bold');
    exportgraphics(fig, outFile, 'Resolution', 160); close(fig);
end

function writeCollectionEfficiencyDiagnosticFigure(detector, opts, outFile)
    % Expected (model) vs measured per-detector centre-point microimage at
    % focus, and the collection efficiency CE = measured/model. If CE
    % anti-correlates with the measured fingerprint, CE is absorbing an
    % optical-model mismatch (e.g. SPAD crosstalk), not pure detector response.
    detXY = detectorLayout(opts.detectorLayout, detector.detectorPitchSampleUm);
    P = detector.measuredPattern; M = detector.modelPattern;
    CE = detector.relativeCollectionEfficiency(:);
    [~, fp] = max(sum(P, 1));                 % brightest (focal) plane
    pm = P(:,fp) / max(sum(P(:,fp)), eps);
    mm = M(:,fp) / max(sum(M(:,fp)), eps);

    fig = figure('Visible','off','Color','w','Position',[60 60 1180 380]);
    tl = tiledlayout(fig, 1, 3, 'Padding','compact','TileSpacing','compact');
    ax = nexttile(tl); plotDetectorHexMap(detXY, mm, 'Parent', ax);
    colormap(ax, hot); colorbar(ax); title(ax,'expected on detector (model)');
    ax = nexttile(tl); plotDetectorHexMap(detXY, pm, 'Parent', ax);
    colormap(ax, hot); colorbar(ax); title(ax,'measured');
    ax = nexttile(tl); plotDetectorHexMap(detXY, CE, 'Parent', ax);
    colormap(ax, parula); colorbar(ax);
    title(ax,'relative collection efficiency = measured / model');

    c = corrVec(CE - mean(CE), pm - mean(pm));
    title(tl, sprintf(['Per-detector: expected vs measured vs collection efficiency   ' ...
        '(corr(CE, measured) = %+.2f)\nnear -1  this CE is absorbing an optical ' ...
        'mismatch (e.g. crosstalk), not pure sensitivity'], c), 'FontWeight','bold');
    exportgraphics(fig, outFile, 'Resolution', 160); close(fig);
end

function c = corrVec(a, b)
    a = a(:); b = b(:);
    d = sqrt(sum(a.^2) * sum(b.^2));
    if d <= 0, c = 0; else, c = sum(a.*b) / d; end
end

function P = projectVol(V, name)
    switch upper(name)
        case 'XY', P = max(V, [], 3);
        case 'XZ', P = squeeze(max(V, [], 1)).';
        case 'YZ', P = squeeze(max(V, [], 2)).';
        otherwise, error('projectVol:bad','unknown projection');
    end
end

function v = getc(coeffs, name)
    if isfield(coeffs, name), v = coeffs.(name); else, v = 0; end
end

function cmap = redBlue(n)
    x = linspace(-1,1,n).';
    cmap = max(min([max(0,x), 1-abs(x), max(0,-x)], 1), 0);
end

function writeSummaryFigure(confocal, wavefront, measPattern, detector, ...
        stageZ, sim, outFile)
    fig = figure('Visible','off','Color','w','Position',[60 60 1300 760]);
    tl = tiledlayout(fig, 2, 3, 'Padding','compact','TileSpacing','compact');

    ax = nexttile(tl);
    plot(ax, stageZ, squeeze(sum(sum(confocal,1),2)), 'ko-'); grid(ax,'on');
    xline(ax, wavefront.z0Um, '--', 'z0');
    xlabel(ax,'stage z [um]'); ylabel(ax,'confocal photons');
    title(ax,'Stage A: confocal axial trace');

    ax = nexttile(tl);
    [~,fp] = min(abs(stageZ - wavefront.z0Um));
    imagesc(ax, sim.x, sim.y, confocal(:,:,fp)); axis(ax,'image'); set(ax,'YDir','normal');
    colormap(ax, hot); colorbar(ax); title(ax,'confocal focal plane (measured)');

    ax = nexttile(tl);
    bar(ax, wavefront.coeffVectorWaves); grid(ax,'on');
    set(ax,'XTick',1:numel(wavefront.fitModes),'XTickLabel',wavefront.fitModes, ...
        'XTickLabelRotation',40);
    ylabel(ax,'waves RMS'); title(ax,'retrieved wavefront');

    ax = nexttile(tl);
    if ~isempty(detector.patternMatchPitch) && ...
            ~isempty(detector.patternMatchPitch.pitchProfile)
        pp = detector.patternMatchPitch.pitchProfile;
        plot(ax, pp(:,1), pp(:,2), 'o-'); grid(ax,'on');
        xline(ax, detector.detectorPitchSampleUm, '--', 'used');
        xlabel(ax,'sample-space pitch [um]'); ylabel(ax,'SSE');
        title(ax,'Stage B: pitch profile (pattern-match)');
    else
        text(ax, 0.5, 0.6, sprintf('pitch = %.4f um', detector.detectorPitchSampleUm), ...
            'HorizontalAlignment','center','Units','normalized','FontWeight','bold');
        if ~isempty(detector.scanPitch)
            text(ax, 0.5, 0.4, sprintf(['scan ruler: %d points\\newlinerotation %.1f' ...
                ' deg, anisotropy %.3f'], detector.scanPitch.nScanPointsUsed, ...
                detector.scanPitch.rotationDeg, detector.scanPitch.anisotropy), ...
                'HorizontalAlignment','center','Units','normalized');
        end
        axis(ax,'off'); title(ax,'Stage B: pitch (scan ruler)');
    end

    ax = nexttile(tl);
    bar(ax, detector.relativeCollectionEfficiency); grid(ax,'on');
    xlabel(ax,'detector index'); ylabel(ax,'relative collection efficiency');
    title(ax,'detector collection efficiency (mean 1)');

    ax = nexttile(tl);
    fit = (detector.relativeCollectionEfficiency(:)*detector.planePhotonScale) .* detector.modelPattern;
    plot(ax, measPattern(:), fit(:), '.'); grid(ax,'on'); axis(ax,'equal');
    hold(ax,'on'); lim=[0 max(measPattern(:))]; plot(ax,lim,lim,'r-'); hold(ax,'off');
    xlabel(ax,'measured per-detector'); ylabel(ax,'fitted'); title(ax,'Stage B fit');

    title(tl, sprintf(['Confocal-sum calibration: z0=%.3f um, pitch=%.4f um, ' ...
        'collection efficiency %.2f-%.2f'], wavefront.z0Um, detector.detectorPitchSampleUm, ...
        min(detector.relativeCollectionEfficiency), max(detector.relativeCollectionEfficiency)), 'FontWeight','bold');
    exportgraphics(fig, outFile, 'Resolution', 160);
    close(fig);
end
