function result = estimateISMAstigFromDetectorPattern(stackInput, varargin)
%--------------------------------------------------------------------------
% estimateISMAstigFromDetectorPattern
%
% PURPOSE
%   Easy, scale-free estimate of low-order lateral aberrations (astigmatism,
%   optionally coma) from an ISM SPAD-array z-stack, by matching the SHAPE of
%   the per-plane detector pattern instead of absolute Poisson counts.
%
% WHY THIS IS ROBUST
%   The detector-integrated signal of channel k (summing the scan image over
%   x,y) is proportional to the overlap integral of the excitation PSF with
%   that channel's detection efficiency:
%
%        P_k(z)  =  sum_xy [ h_exc(z) .* eta_k(z) ]   ~   integral h_exc * eta_k
%
%   This is independent of the bead's lateral position and of the absolute
%   per-plane brightness. Normalising each plane's 23-channel pattern to
%   sum 1 therefore removes exactly the axial-trace / scale mismatch that
%   makes the full absolute-count fit (estimateFullStackISMWavefront) fail on
%   the air-interface data, while keeping the astigmatism signature (the
%   pattern elongates along one axis below focus and the orthogonal axis
%   above focus).
%
%   z0 is fitted (cheap, 1 parameter) so a wrong focus guess cannot bias the
%   result. Spherical/defocus are NOT fitted (defocus is supplied by each
%   plane's own z). A homogeneous scalar PSF is used for simplicity.
%
% USAGE
%   result = estimateISMAstigFromDetectorPattern();      % default dataset
%   result = estimateISMAstigFromDetectorPattern(alignmentCsv);
%   result = estimateISMAstigFromDetectorPattern([], 'fitModes', ...
%       {'astig_x','astig_y','coma_x','coma_y'});
%
% KEY OPTIONS (name-value)
%   'fitModes'          modes to fit  (default {astig_x,astig_y,coma_x,coma_y})
%   'planeHalfRangeUm'  +/- window about focus used         (default 0.3)
%   'maxPlanes'         cap on number of planes used         (default 9)
%   'focusZUm'          force focus z [um]    (default [] -> brightness peak)
%   'fitZ0'             also fit the focus position          (default true)
%   'fitDetScale'       fit an isotropic detector-scale that (default true)
%                       absorbs the PSF-width / pitch mismatch
%                       between the simple homogeneous model
%                       and the real (broader) air PSF
%   'detScaleBounds'    bounds for that scale                (default [0.3 2.5])
%   'NA','lamExc','lamEm','nMedium'           optics         (1.2/.488/.520/1.33)
%   'detectorLayout'    'fixed' = nominal PDA23 honeycomb    (default 'fixed')
%                       with no data-derived geometry;
%                       'regular' = honeycomb aligned to the
%                       data; 'measured' = data-derived
%   'coeffBoundsWaves'  coefficient bounds                   (default [-0.3 0.3])
%   'z0HalfRangeUm'     z0 search half-range about peak      (default 0.3)
%   'hotPixelFactor'    exclude channels with background >   (default 8)
%                       factor x median (auto-handles the hot pixel)
%   'excludeChannelIDs' channel IDs to drop from the fit     (default [])
%   ... loader passthrough: 'channelIDs','darkFile','darkScale',
%       'backgroundMode','xyPixelSizeUm','reuseCache','cacheFile'
%   'outputDir','writeOutputs','verbose'
%
% OUTPUT (result struct)
%   estCoeffs / estCoeffVector / estCoeffNm   fitted coefficients
%   estZ0Um                                    fitted focus
%   selectedPlaneIndices / relativeZUm         planes used
%   measuredPattern / modelPattern             normalised 23 x nPlane patterns
%   patternResidual / perPlaneCorrelation      fit quality
%--------------------------------------------------------------------------

    addpath(fileparts(mfilename('fullpath')));
    if nargin < 1 || isempty(stackInput)
        stackInput = ['D:\Luminosa\Data\PSF_batch_outputs\ISM_Aberation2_73\' ...
            'xz_yz_plots\x0_4uW_0_19collar_80mmlens_20260515_155744_frame_alignment.csv'];
    end
    opts = parseInputs(varargin{:});

    % ----- load raw detector stack (reuse the existing loader) ------------
    data = loadFullStackISMData(stackInput, ...
        'channelIDs', opts.channelIDs, 'darkFile', opts.darkFile, ...
        'darkScale', opts.darkScale, 'backgroundMode', opts.backgroundMode, ...
        'xyPixelSizeUm', opts.xyPixelSizeUm, 'reuseCache', opts.reuseCache, ...
        'cacheFile', opts.cacheFile, 'verbose', opts.verbose);

    [ny, nx, nCh, nPlane] = size(data.rawCounts);
    stageZ = data.stageZUm(:).';

    % ----- measured detector pattern: integrate scan image over x,y -------
    rawSumXY  = reshape(sum(sum(data.rawCounts,1),2), nCh, nPlane);
    bgPerPix  = reshape(data.backgroundPerPixel, nCh, nPlane);
    measPattern = max(rawSumXY - bgPerPix*nx*ny, 0);     % [nCh x nPlane]
    axialSignal = sum(measPattern, 1);

    % ----- focus plane + symmetric near-focus window ----------------------
    if isempty(opts.focusZUm)
        [~, focusIdx] = max(axialSignal);
    else
        [~, focusIdx] = min(abs(stageZ - opts.focusZUm));
    end
    z0Guess = stageZ(focusIdx);
    sel = find(abs(stageZ - z0Guess) <= opts.planeHalfRangeUm);
    if numel(sel) > opts.maxPlanes
        keep = round(linspace(1, numel(sel), opts.maxPlanes));
        sel = sel(keep);
    end
    if numel(sel) < 3
        [~, order] = sort(abs(stageZ - z0Guess));
        sel = sort(order(1:min(5,nPlane)));
    end
    relZ = stageZ(sel) - z0Guess;

    % ----- channels to use (drop hot / requested pixels) ------------------
    medBg = median(bgPerPix(:, focusIdx));
    hot = bgPerPix(:, focusIdx) > opts.hotPixelFactor * max(medBg, eps);
    excludeByID = ismember(data.channelIDs(:), opts.excludeChannelIDs(:));
    includeCh = ~(hot | excludeByID);
    if nnz(includeCh) < 6
        includeCh = true(nCh,1);     % safety: never drop everything
    end

    % ----- minimal simulation config (homogeneous scalar PSF) -------------
    sim = configureSim(data, focusIdx, opts);

    % measured patterns normalised over included channels
    measNorm = normalizePattern(measPattern(:, sel), includeCh);

    % ----- coarse search to seed the optimiser ----------------------------
    % First the isotropic detector scale (absorbs PSF-width / pitch mismatch
    % between the simple homogeneous model and the real air-interface PSF),
    % then astigmatism at that scale.
    c0 = zeros(1, numel(opts.fitModes));
    detScale0 = opts.detScaleInit;
    if opts.fitDetScale
        detScale0 = coarseDetScaleSearch(c0, opts.fitModes, sim, relZ, ...
            measNorm, includeCh, opts);
    end
    c0 = coarseAstigSearch(c0, opts.fitModes, sim, relZ, measNorm, ...
        includeCh, opts, detScale0);

    % ----- refine with Nelder-Mead (base MATLAB, no toolbox) --------------
    theta0 = c0;
    if opts.fitZ0, theta0 = [theta0, z0Guess]; end
    if opts.fitDetScale, theta0 = [theta0, detScale0]; end
    objFun = @(t) patternObjective(t, opts.fitModes, opts.fitZ0, z0Guess, ...
        opts.fitDetScale, opts.detScaleInit, sim, stageZ(sel), measNorm, ...
        includeCh, opts);
    fmOpts = optimset('Display','off', 'MaxFunEvals', opts.maxFunEvals, ...
        'MaxIter', opts.maxIter, 'TolX', 1e-4, 'TolFun', 1e-6);
    [thetaHat, fval] = fminsearch(objFun, theta0, fmOpts);

    % ----- unpack result --------------------------------------------------
    [coeffVec, z0Hat, detScaleHat] = unpackTheta(thetaHat, opts.fitModes, ...
        opts.fitZ0, z0Guess, opts.fitDetScale, opts.detScaleInit);
    coeffVec = min(max(coeffVec, opts.coeffBoundsWaves(1)), opts.coeffBoundsWaves(2));
    detScaleHat = min(max(detScaleHat, opts.detScaleBounds(1)), opts.detScaleBounds(2));
    estCoeffs = makeCoeffs(opts.fitModes, coeffVec);
    relZHat = stageZ(sel) - z0Hat;
    modelNorm = normalizePattern( ...
        forwardPattern(sim, estCoeffs, relZHat, detScaleHat), includeCh);

    perPlaneCorr = zeros(1, numel(sel));
    for p = 1:numel(sel)
        a = measNorm(includeCh, p); b = modelNorm(includeCh, p);
        perPlaneCorr(p) = corr2safe(a, b);
    end

    result = struct();
    result.estCoeffs        = estCoeffs;
    result.estCoeffVector   = coeffVec;
    result.estCoeffNm       = coeffVec * sim.lamRef * 1000;
    result.fitModes         = opts.fitModes;
    result.estZ0Um          = z0Hat;
    result.estDetScale      = detScaleHat;
    result.focusGuessZUm    = z0Guess;
    result.selectedPlaneIndices = sel;
    result.relativeZUm      = relZHat;
    result.measuredPattern  = measNorm;
    result.modelPattern     = modelNorm;
    result.patternResidual  = measNorm - modelNorm;
    result.perPlaneCorrelation = perPlaneCorr;
    result.meanCorrelation  = mean(perPlaneCorr);
    result.objectiveValue   = fval;
    result.includedChannels = find(includeCh).';
    result.excludedChannelIDs = data.channelIDs(~includeCh);
    result.detXY            = sim.detXY;
    result.sim              = sim;
    result.options          = opts;
    result.outputDir        = resolveOutputDir(opts);

    if opts.verbose, printSummary(result, data); end
    if opts.writeOutputs, writeOutputs(result, data); end
end

% ========================================================================
function opts = parseInputs(varargin)
    p = inputParser; p.FunctionName = 'estimateISMAstigFromDetectorPattern';
    addParameter(p, 'fitModes', {'astig_x','astig_y','coma_x','coma_y'});
    addParameter(p, 'planeHalfRangeUm', 0.15);
    addParameter(p, 'maxPlanes', 9);
    addParameter(p, 'focusZUm', []);
    addParameter(p, 'fitZ0', true);
    addParameter(p, 'NA', 1.2);
    addParameter(p, 'lamExc', 0.488);
    addParameter(p, 'lamEm', 0.520);
    addParameter(p, 'nMedium', 1.33);
    addParameter(p, 'beadRadiusUm', 0.08);
    addParameter(p, 'detectorSubsamples', 5);
    addParameter(p, 'detFillRatio', 1.0);
    addParameter(p, 'detectorLayout', 'fixed');  % 'fixed' | 'regular' | 'measured'
    addParameter(p, 'detectorLayoutScale', 2);
    addParameter(p, 'detectorLayoutPositionSign', -1);
    addParameter(p, 'detectorShiftUpsample', 20);
    addParameter(p, 'detectorPitchSampleUm', []);
    addParameter(p, 'detectorHardwarePitchUm', []);
    addParameter(p, 'detectorTotalMagnification', []);
    addParameter(p, 'coeffBoundsWaves', [-0.3 0.3]);
    addParameter(p, 'z0HalfRangeUm', 0.3);
    addParameter(p, 'fitDetScale', true);
    addParameter(p, 'detScaleInit', 1);
    addParameter(p, 'detScaleBounds', [0.3 2.5]);
    addParameter(p, 'detScaleCoarse', [0.4 0.55 0.7 0.85 1.0 1.2 1.5 2.0]);
    addParameter(p, 'coarseGridWaves', -0.15:0.075:0.15);
    addParameter(p, 'maxFunEvals', 900);
    addParameter(p, 'maxIter', 700);
    addParameter(p, 'hotPixelFactor', 8);
    addParameter(p, 'excludeChannelIDs', []);
    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'darkFile', 'D:\Luminosa\Data\ISMdark_counts.ptu');
    addParameter(p, 'darkScale', 1);
    addParameter(p, 'backgroundMode', 'auto');
    addParameter(p, 'xyPixelSizeUm', []);
    addParameter(p, 'reuseCache', false);
    addParameter(p, 'cacheFile', '');
    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'verbose', true);
    parse(p, varargin{:});
    opts = p.Results;
    if ischar(opts.fitModes) || isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    end
    opts.detectorLayout = lower(char(opts.detectorLayout));
    if ~ismember(opts.detectorLayout, {'fixed','regular','measured'})
        error('estimateISMAstigFromDetectorPattern:BadDetectorLayout', ...
            'detectorLayout must be ''fixed'', ''regular'', or ''measured''.');
    end
end

% ------------------------------------------------------------------------
function sim = configureSim(data, focusIdx, opts)
    sim = defaultParams();
    sim.NA = opts.NA; sim.nMedium = opts.nMedium;
    sim.lamExc = opts.lamExc; sim.lamEm = opts.lamEm; sim.lamRef = opts.lamEm;
    sim.sampleGeometry = 'homogeneous';
    sim.beadRadius = opts.beadRadiusUm;
    sim.x = data.xUm(:).'; sim.y = data.yUm(:).';
    sim.nx = numel(sim.x); sim.ny = numel(sim.y);
    sim.dx = data.xyPixelSizeUm;
    sim.detectorSubsamples = opts.detectorSubsamples;
    sim.detectorImageInverted = true;
    sim.detectorPixelShape = 'hex';
    detectorPitch = resolveDetectorPitchSampleUmAstig(opts, sim.detPitch);

    switch opts.detectorLayout
        case 'fixed'
            detXY = detectorLayout('honeycomb23', detectorPitch);
            detMeasured = detXY;
        case {'regular','measured'}
            layoutOpts = struct('positionSign', opts.detectorLayoutPositionSign, ...
                'detectorScale', opts.detectorLayoutScale, 'centerMode', 'reference', ...
                'upsampleReg', opts.detectorShiftUpsample, ...
                'airyUnitUm', 1.22*sim.lamEm/sim.NA);
            detMeasured = estimateDetectorLayoutFromStack( ...
                data.rawCounts(:,:,:,focusIdx), data.xyPixelSizeUm, layoutOpts);
            if strcmpi(opts.detectorLayout, 'regular')
                % Keep the layout a clean honeycomb (no per-detector shifts), but
                % aligned to the data's own pitch/orientation so channel identity and
                % astigmatism orientation are preserved.
                detXY = regularizeDetectorLayout(detMeasured);
            else
                detXY = detMeasured;
            end
    end
    sim.detXY = detXY;
    sim.detXYMeasured = detMeasured;
    sim.detectorLayoutMode = opts.detectorLayout;
    sim.detectorPitchSource = detectorPitchSourceAstig(opts);
    sim.nDet = size(detXY,1);
    pitch = medianNearestPitch(detXY);
    sim.detPitch = pitch;
    sim.detSize = opts.detFillRatio * pitch;
    sim.detectorHexRadius = sim.detSize / sqrt(3);
end

% ------------------------------------------------------------------------
function detReg = regularizeDetectorLayout(detMeasured)
    % Snap the measured detector cloud onto an ideal honeycomb of the same
    % pitch/orientation, preserving each channel's index. Robust to arbitrary
    % rotation and reflection (the SPAD frame may be flipped vs scan axes).
    template = detectorLayout('honeycomb23', 1);
    template = template - mean(template, 1);
    centroid = mean(detMeasured, 1);
    M = detMeasured - centroid;
    pitch = medianNearestPitch(detMeasured);
    if ~isfinite(pitch) || pitch <= 0, detReg = detMeasured; return; end
    theta0 = estimateLatticeAngleDeg(detMeasured);

    bestCost = inf; detReg = detMeasured;
    for f = [1 -1]                                  % no-reflection / reflection
        Tf = template; Tf(:,1) = f * Tf(:,1);
        for extra = 0:60:300                        % hex-symmetry orientations
            th = deg2rad(theta0 + extra);
            R = [cos(th) -sin(th); sin(th) cos(th)];
            T = (pitch * Tf) * R.';
            [assign, cost] = greedyAssign(M, T);
            if cost < bestCost
                bestCost = cost;
                detReg = T(assign, :) + centroid;   % regular site per channel
            end
        end
    end
end

% ------------------------------------------------------------------------
function [assign, cost] = greedyAssign(M, T)
    % Unique nearest assignment of each row of M to a row of T (n small).
    n = size(M,1);
    D = zeros(n,n);
    for i = 1:n
        D(i,:) = hypot(T(:,1)-M(i,1), T(:,2)-M(i,2)).';
    end
    [~, order] = sort(D(:));
    assign = zeros(n,1);
    usedT = false(n,1); doneM = false(n,1); cost = 0; filled = 0;
    for q = 1:numel(order)
        [i, j] = ind2sub([n n], order(q));
        if ~doneM(i) && ~usedT(j)
            assign(i) = j; doneM(i) = true; usedT(j) = true;
            cost = cost + D(i,j); filled = filled + 1;
            if filled == n, break; end
        end
    end
end

% ------------------------------------------------------------------------
function thDeg = estimateLatticeAngleDeg(xy)
    % Dominant nearest-neighbour bond angle, folded to [0,60) (hex lattice).
    n = size(xy,1); ang = zeros(n,1);
    for k = 1:n
        d = xy - xy(k,:); dist = hypot(d(:,1), d(:,2)); dist(k) = inf;
        [~, j] = min(dist);
        ang(k) = atan2(d(j,2), d(j,1));
    end
    z = mean(exp(1i*6*ang));            % 60-degree periodicity
    thDeg = mod(rad2deg(angle(z))/6, 60);
end

% ------------------------------------------------------------------------
function P = forwardPattern(sim, coeffs, relZ, detScale)
    % Detector-integrated, per-channel ISM pattern at each relative z.
    % detScale isotropically rescales the detector geometry to absorb any
    % PSF-width / pitch-calibration mismatch (model PSF vs real air PSF).
    if nargin < 4 || isempty(detScale), detScale = 1; end
    sim.detXY = sim.detXY * detScale;
    sim.detSize = sim.detSize * detScale;
    sim.detectorHexRadius = sim.detectorHexRadius * detScale;
    sim.z = relZ(:).';
    hExc = psfBessel(sim, coeffs, sim.lamExc);                % [ny nx nP]
    hDet = psfBessel(sim, coeffs, sim.lamEm);
    eta  = detectorCollectionEfficiencyExplicit(sim, hDet);  % [ny nx nP nCh]
    nCh = size(sim.detXY,1); nP = numel(relZ);
    P = zeros(nCh, nP);
    for ip = 1:nP
        he = hExc(:,:,ip);
        P(:,ip) = reshape(sum(sum(eta(:,:,ip,:) .* he, 1), 2), nCh, 1);
    end
end

% ------------------------------------------------------------------------
function Pn = normalizePattern(P, includeCh)
    Pn = zeros(size(P));
    for p = 1:size(P,2)
        s = sum(P(includeCh, p));
        if s > 0, Pn(:,p) = P(:,p) / s; end
    end
end

% ------------------------------------------------------------------------
function f = patternObjective(theta, modes, fitZ0, z0Guess, fitDetScale, ...
        detScaleInit, sim, planeZ, measNorm, includeCh, opts)
    [coeffVec, z0, detScale] = unpackTheta(theta, modes, fitZ0, z0Guess, ...
        fitDetScale, detScaleInit);
    penalty = boundsPenalty(coeffVec, z0, z0Guess, detScale, opts);
    coeffVec = min(max(coeffVec, opts.coeffBoundsWaves(1)), opts.coeffBoundsWaves(2));
    detScale = min(max(detScale, opts.detScaleBounds(1)), opts.detScaleBounds(2));
    coeffs = makeCoeffs(modes, coeffVec);
    modelNorm = normalizePattern( ...
        forwardPattern(sim, coeffs, planeZ - z0, detScale), includeCh);
    d = measNorm(includeCh, :) - modelNorm(includeCh, :);
    f = sum(d(:).^2) + penalty;
end

% ------------------------------------------------------------------------
function c0 = coarseAstigSearch(c0, modes, sim, relZ, measNorm, includeCh, ...
        opts, detScale)
    ix = find(strcmp(modes, 'astig_x'), 1);
    iy = find(strcmp(modes, 'astig_y'), 1);
    if isempty(ix) && isempty(iy), return; end
    grid = opts.coarseGridWaves;
    best = inf; bestC = c0;
    for ax = grid
        for ay = grid
            c = c0;
            if ~isempty(ix), c(ix) = ax; end
            if ~isempty(iy), c(iy) = ay; end
            modelNorm = normalizePattern( ...
                forwardPattern(sim, makeCoeffs(modes, c), relZ, detScale), includeCh);
            d = measNorm(includeCh, :) - modelNorm(includeCh, :);
            val = sum(d(:).^2);
            if val < best, best = val; bestC = c; end
        end
    end
    c0 = bestC;
end

% ------------------------------------------------------------------------
function s0 = coarseDetScaleSearch(c0, modes, sim, relZ, measNorm, includeCh, opts)
    best = inf; s0 = opts.detScaleInit;
    for s = opts.detScaleCoarse
        modelNorm = normalizePattern( ...
            forwardPattern(sim, makeCoeffs(modes, c0), relZ, s), includeCh);
        d = measNorm(includeCh, :) - modelNorm(includeCh, :);
        val = sum(d(:).^2);
        if val < best, best = val; s0 = s; end
    end
end

% ------------------------------------------------------------------------
function [coeffVec, z0, detScale] = unpackTheta(theta, modes, fitZ0, ...
        z0Guess, fitDetScale, detScaleInit)
    nM = numel(modes);
    coeffVec = theta(1:nM);
    idx = nM;
    if fitZ0, idx = idx + 1; z0 = theta(idx); else, z0 = z0Guess; end
    if fitDetScale, idx = idx + 1; detScale = theta(idx); else, detScale = detScaleInit; end
end

% ------------------------------------------------------------------------
function pen = boundsPenalty(coeffVec, z0, z0Guess, detScale, opts)
    over = max(abs(coeffVec) - opts.coeffBoundsWaves(2), 0);
    zover = max(abs(z0 - z0Guess) - opts.z0HalfRangeUm, 0);
    sover = max(detScale - opts.detScaleBounds(2), 0) + ...
        max(opts.detScaleBounds(1) - detScale, 0);
    pen = 50 * (sum(over.^2) + zover.^2 + sover.^2);
end

% ------------------------------------------------------------------------
function c = makeCoeffs(modes, vals)
    c = struct();
    for i = 1:numel(modes), c.(modes{i}) = vals(i); end
end

% ------------------------------------------------------------------------
function r = corr2safe(a, b)
    a = a - mean(a); b = b - mean(b);
    d = sqrt(sum(a.^2) * sum(b.^2));
    if d <= 0, r = NaN; else, r = sum(a.*b) / d; end
end

% ------------------------------------------------------------------------
function pitch = medianNearestPitch(xy)
    n = size(xy,1); nearest = nan(n,1);
    for k = 1:n
        d = hypot(xy(:,1)-xy(k,1), xy(:,2)-xy(k,2)); d(k) = inf;
        nearest(k) = min(d);
    end
    nearest = nearest(isfinite(nearest) & nearest > 0);
    if isempty(nearest), pitch = 1; else, pitch = median(nearest); end
end

% ------------------------------------------------------------------------
function pitch = resolveDetectorPitchSampleUmAstig(opts, defaultPitch)
    pitch = defaultPitch;
    if isnumeric(opts.detectorPitchSampleUm) && ...
            isscalar(opts.detectorPitchSampleUm) && ...
            isfinite(opts.detectorPitchSampleUm) && opts.detectorPitchSampleUm > 0
        pitch = double(opts.detectorPitchSampleUm);
        return;
    end
    physical = opts.detectorHardwarePitchUm;
    magnification = opts.detectorTotalMagnification;
    if isnumeric(physical) && isscalar(physical) && isfinite(physical) && ...
            physical > 0 && isnumeric(magnification) && isscalar(magnification) && ...
            isfinite(magnification) && magnification > 0
        pitch = double(physical) / double(magnification);
    end
end

% ------------------------------------------------------------------------
function source = detectorPitchSourceAstig(opts)
    if isnumeric(opts.detectorPitchSampleUm) && ...
            isscalar(opts.detectorPitchSampleUm) && ...
            isfinite(opts.detectorPitchSampleUm) && opts.detectorPitchSampleUm > 0
        source = 'detectorPitchSampleUm';
    elseif isnumeric(opts.detectorHardwarePitchUm) && ...
            isscalar(opts.detectorHardwarePitchUm) && ...
            isfinite(opts.detectorHardwarePitchUm) && ...
            opts.detectorHardwarePitchUm > 0 && ...
            isnumeric(opts.detectorTotalMagnification) && ...
            isscalar(opts.detectorTotalMagnification) && ...
            isfinite(opts.detectorTotalMagnification) && ...
            opts.detectorTotalMagnification > 0
        source = 'detectorHardwarePitchUm / detectorTotalMagnification';
    else
        source = 'defaultParams.detPitch';
    end
end

% ------------------------------------------------------------------------
function outDir = resolveOutputDir(opts)
    if ~isempty(opts.outputDir)
        outDir = opts.outputDir;
    else
        root = fileparts(fileparts(mfilename('fullpath')));
        outDir = fullfile(root, 'output_matlab', 'ism_astig_detector_pattern');
    end
end

% ------------------------------------------------------------------------
function printSummary(result, data) %#ok<INUSD>
    fprintf('\n[estimateISMAstigFromDetectorPattern] scale-free pattern fit\n');
    fprintf('  focus guess z0 = %.3f um; fitted z0 = %.3f um\n', ...
        result.focusGuessZUm, result.estZ0Um);
    fprintf('  detector layout: %s; pitch source: %s\n', ...
        result.sim.detectorLayoutMode, result.sim.detectorPitchSource);
    fprintf('  fitted detector scale: %.3f (1 = nominal; absorbs PSF/pitch width)\n', ...
        result.estDetScale);
    fprintf('  planes used: %d (relative z %.2f..%.2f um)\n', ...
        numel(result.selectedPlaneIndices), min(result.relativeZUm), ...
        max(result.relativeZUm));
    if ~isempty(result.excludedChannelIDs)
        fprintf('  excluded channel IDs: %s\n', num2str(result.excludedChannelIDs(:).'));
    end
    fprintf('  mean per-plane pattern correlation: %.4f\n', result.meanCorrelation);
    for k = 1:numel(result.fitModes)
        fprintf('     %-10s %+.4f waves (%+.1f nm RMS)\n', ...
            result.fitModes{k}, result.estCoeffVector(k), result.estCoeffNm(k));
    end
    fprintf('\n');
end

% ------------------------------------------------------------------------
function writeOutputs(result, data) %#ok<INUSD>
    outDir = result.outputDir;
    if exist(outDir,'dir') ~= 7, mkdir(outDir); end
    T = table(result.fitModes(:), result.estCoeffVector(:), result.estCoeffNm(:), ...
        'VariableNames', {'mode','estimateWaves','estimateNmRms'});
    writetable(T, fullfile(outDir, 'astig_pattern_coefficients.csv'));
    save(fullfile(outDir, 'ism_astig_pattern_fit.mat'), 'result', '-v7.3');
    writeFigure(result, fullfile(outDir, 'astig_pattern_fit.png'));
end

% ------------------------------------------------------------------------
function writeFigure(result, outFile)
    sel = result.selectedPlaneIndices; nSel = numel(sel);
    detXY = result.detXY; inc = false(size(detXY,1),1);
    inc(result.includedChannels) = true;
    meas = result.measuredPattern; model = result.modelPattern;
    resid = result.patternResidual;
    meas(~inc,:) = NaN; model(~inc,:) = NaN; resid(~inc,:) = NaN;
    rmax = max(abs(resid(:)), [], 'omitnan'); if isempty(rmax)||rmax==0, rmax=1; end
    divmap = bluewhitered();

    fig = figure('Visible','off','Color','w', ...
        'Position',[80 80 720 max(180*nSel,360)]);
    tl = tiledlayout(fig, nSel, 3, 'Padding','compact','TileSpacing','compact');
    for p = 1:nSel
        cmax = max([meas(:,p);model(:,p)], [], 'omitnan');
        if isempty(cmax)||cmax<=0, cmax = 1; end
        ax = nexttile(tl); plotDetectorHexMap(detXY, meas(:,p), 'Parent', ax);
        caxis(ax,[0 cmax]); colormap(ax, parula);
        if p==1, title(ax,'measured'); end
        ylabel(ax, sprintf('z=%+.2f', result.relativeZUm(p)), ...
            'Rotation',0,'HorizontalAlignment','right','Visible','on');
        ax = nexttile(tl); plotDetectorHexMap(detXY, model(:,p), 'Parent', ax);
        caxis(ax,[0 cmax]); colormap(ax, parula);
        if p==1, title(ax,'fitted'); end
        ax = nexttile(tl); plotDetectorHexMap(detXY, resid(:,p), 'Parent', ax);
        caxis(ax,[-rmax rmax]); colormap(ax, divmap);
        if p==1, title(ax,'residual'); end
    end
    parts = cell(1, numel(result.fitModes));
    for k = 1:numel(result.fitModes)
        parts{k} = sprintf('%s %+.3f\\lambda', result.fitModes{k}, ...
            result.estCoeffVector(k));
    end
    coeffStr = strjoin(parts, ',  ');
    title(tl, {sprintf('Scale-free detector-pattern fit: %s', coeffStr), ...
        sprintf('z0 = %.3f um, det scale %.2f, mean corr %.3f', ...
        result.estZ0Um, result.estDetScale, result.meanCorrelation)}, ...
        'FontWeight','bold');
    exportgraphics(fig, outFile, 'Resolution', 160);
    close(fig);
end

% ------------------------------------------------------------------------
function cmap = bluewhitered()
    n = 128;
    top = [linspace(0,1,n).' linspace(0,1,n).' ones(n,1)];
    bot = [ones(n,1) linspace(1,0,n).' linspace(1,0,n).'];
    cmap = [top; bot];
end
