function res = phaseRetrieval3DBead(useRealData, rawData, planeZ, opts)
%--------------------------------------------------------------------------
% phaseRetrieval3DBead
%
% PURPOSE
%   Retrieve a multi-mode pupil phase from a detector-resolved, multi-plane
%   bead measurement using the configured SPAD-array forward model.
%
% USAGE
%   res = phaseRetrieval3DBead(false);                    % simulate 3 planes
%   res = phaseRetrieval3DBead(false, [], [-0.2 0 0.2]);  % simulate custom z
%   res = phaseRetrieval3DBead(true, rawData, planeZ);     % fit real data
%
% INPUTS
%   useRealData : false = simulate noisy data, true = fit supplied rawData
%   rawData     : measured data. Accepted layouts:
%                   [ny nx nDetector nPlane]
%                   [ny nx nPlane nDetector]
%                   [ny nx nDetector] for one plane
%                 or an nDetector-element cell array, each cell [ny nx nPlane].
%   planeZ      : recorded bead-centre z positions [um]. Default:
%                   [-0.1 0 0.1] for +/-100 nm plus focus.
%   opts        : optional struct. Common fields:
%                   sim, fitModes, maxIter, photons, bg, trueCoeffs,
%                   trueXY, trueZOffset, initialCoeffs, initialXY,
%                   initialZOffset, fitXY, fitZ, autoZGrid, modelDz,
%                   modelZPadding, xyPixelSizeUm, verbose.
%
% OUTPUT
%   res : structure containing normalized data, fitted coefficients, fitted
%         phase map, fitted model stack, residual history, and simulation
%         truth fields when useRealData is false.
%
% NOTES
%   Coefficients are fitted in waves RMS at sim.lamRef. The default fit modes
%   exclude tilt because x/y bead shift already absorbs most tilt-like image
%   displacement in intensity data. z-position fitting is off by default
%   because axial offset and defocus are strongly correlated.
%--------------------------------------------------------------------------

    addParentFolderToPath();

    if nargin < 1 || isempty(useRealData), useRealData = false; end
    if ischar(useRealData) || isstring(useRealData)
        useRealData = any(strcmpi(char(useRealData), {'real','data','measured','true'}));
    end
    if nargin < 2, rawData = []; end
    if nargin < 3 || isempty(planeZ), planeZ = [-0.1 0 0.1]; end
    if nargin < 4 || isempty(opts), opts = struct(); end

    planeZ = planeZ(:).';
    sim = getOption(opts, 'sim', defaultParams());
    sim = prepareMultiplaneSim(sim, planeZ, opts);
    fdXYWasProvided = isfield(opts, 'fdXY') && ~isempty(opts.fdXY);

    fitModes = getOption(opts, 'fitModes', ...
        {'defocus','astig_x','astig_y','coma_x','coma_y','spherical'});
    if ischar(fitModes)
        fitModes = {fitModes};
    end
    if isstring(fitModes)
        fitModes = cellstr(fitModes);
    end
    fitModes = fitModes(:).';

    opts.fitXY = getOption(opts, 'fitXY', true);
    opts.fitZ = getOption(opts, 'fitZ', false);
    opts.maxIter = getOption(opts, 'maxIter', 8);
    opts.fdCoeff = getOption(opts, 'fdCoeff', 0.01);
    opts.fdXY = getOption(opts, 'fdXY', sim.dx/4);
    opts.fdZ = getOption(opts, 'fdZ', 0.02);
    opts.regCoeff = getOption(opts, 'regCoeff', 1e-5);
    opts.regXY = getOption(opts, 'regXY', 1e-6);
    opts.regZ = getOption(opts, 'regZ', 1e-6);
    opts.maxCoeffStep = getOption(opts, 'maxCoeffStep', 0.05);
    opts.maxXYStep = getOption(opts, 'maxXYStep', 0.05);
    opts.maxZStep = getOption(opts, 'maxZStep', 0.05);
    opts.verbose = getOption(opts, 'verbose', true);

    if useRealData
        if isempty(rawData)
            error('phaseRetrieval3DBead:MissingData', ...
                'rawData is required when useRealData is true.');
        end
        raw4 = standardizeSpadData(rawData, sim, numel(planeZ));
        sim = matchSimLateralGridToData(sim, raw4, opts);
        if ~fdXYWasProvided
            opts.fdXY = sim.dx/4;
        end
    else
        trueCoeffs = getOption(opts, 'trueCoeffs', defaultTrueCoeffs());
        trueXY = getOption(opts, 'trueXY', [0.02 -0.015]);
        trueZOffset = getOption(opts, 'trueZOffset', 0);
        photons = getOption(opts, 'photons', 3e5);
        bg = getOption(opts, 'bg', 0.1);

        trueStack = normalizedStackExplicitDetectorZPlanes( ...
            sim, trueCoeffs, planeZ, trueXY(1), trueXY(2), trueZOffset);
        raw4 = poissonSample(photons*trueStack + bg);
    end

    [dataN, beta] = normalizeSpadVolume(raw4);

    paramNames = buildParamNames(fitModes, opts.fitXY, opts.fitZ);
    p = initialParameterVector(sim, fitModes, opts);
    step = finiteDifferenceSteps(numel(fitModes), opts.fitXY, opts.fitZ, opts);
    reg = regularizationVector(numel(fitModes), opts.fitXY, opts.fitZ, opts);
    maxStep = maxUpdateVector(numel(fitModes), opts.fitXY, opts.fitZ, opts);

    history = zeros(opts.maxIter, 2);

    for it = 1:opts.maxIter
        m0 = modelFromParams(sim, fitModes, planeZ, p, opts.fitXY, opts.fitZ);
        r = dataN(:) - m0(:);
        J = zeros(numel(m0), numel(p));

        for q = 1:numel(p)
            pp = p;
            pm = p;
            pp(q) = pp(q) + step(q);
            pm(q) = pm(q) - step(q);

            mp = modelFromParams(sim, fitModes, planeZ, pp, opts.fitXY, opts.fitZ);
            mm = modelFromParams(sim, fitModes, planeZ, pm, opts.fitXY, opts.fitZ);
            J(:,q) = (mp(:) - mm(:)) / (2*step(q));
        end

        H = J.'*J + diag(reg);
        g = J.'*r;

        if any(~isfinite(H(:))) || any(~isfinite(g(:)))
            error('phaseRetrieval3DBead:NonFiniteSystem', ...
                'Non-finite values encountered in Gauss-Newton system.');
        end

        if rcond(H) < 1e-12
            delta = pinv(H) * g;
        else
            delta = H \ g;
        end

        delta = max(min(delta(:).', maxStep), -maxStep);
        p = p + delta;

        history(it,1) = norm(r);
        history(it,2) = norm(delta);

        if opts.verbose
            fprintf('[phaseRetrieval3DBead] iter %02d  residual %.4e  step %.4e\n', ...
                it, history(it,1), history(it,2));
        end

        if norm(delta) < getOption(opts, 'tolStep', 1e-4)
            history = history(1:it,:);
            break;
        end
    end

    fitStack = modelFromParams(sim, fitModes, planeZ, p, opts.fitXY, opts.fitZ);
    [estCoeffs, estXY, estZOffset] = unpackParams(sim, fitModes, p, opts.fitXY, opts.fitZ);

    res.useRealData = logical(useRealData);
    res.rawData = raw4;
    res.dataN = dataN;
    res.background = beta;
    res.planeZ = planeZ;
    res.fitModes = fitModes;
    res.paramNames = paramNames;
    res.paramVector = p;
    res.estCoeffs = estCoeffs;
    res.estCoeffVector = coeffStructToVector(sim, estCoeffs);
    res.estXY = estXY;
    res.estZOffset = estZOffset;
    res.fitStack = fitStack;
    res.residual = dataN - fitStack;
    res.residualNorm = norm(res.residual(:));
    res.history = history;
    res.phase = zernikePhaseMap(sim, estCoeffs, sim.lamEm);
    res.sim = sim;
    res.options = opts;

    if ~useRealData
        res.trueCoeffs = trueCoeffs;
        res.trueCoeffVector = coeffStructToVector(sim, trueCoeffs);
        res.trueXY = trueXY;
        res.trueZOffset = trueZOffset;
        res.trueStack = trueStack;
        res.truePhase = zernikePhaseMap(sim, trueCoeffs, sim.lamEm);
    end
end

function addParentFolderToPath()
    thisDir = fileparts(mfilename('fullpath'));
    parentDir = fileparts(thisDir);

    if exist(fullfile(parentDir, 'defaultParams.m'), 'file') == 2
        addpath(parentDir);
    end
end

function v = getOption(opts, name, defaultValue)
    if isfield(opts, name) && ~isempty(opts.(name))
        v = opts.(name);
    else
        v = defaultValue;
    end
end

function sim = prepareMultiplaneSim(sim, planeZ, opts)
    if ~getOption(opts, 'autoZGrid', true)
        return;
    end

    dzTarget = getOption(opts, 'modelDz', 0.05);          % 50 nm
    zPadding = getOption(opts, 'modelZPadding', 0.50);    % +/-500 nm support

    zMin = min([planeZ(:); 0]) - zPadding;
    zMax = max([planeZ(:); 0]) + zPadding;
    nZ = ceil((zMax - zMin) / dzTarget) + 1;
    if mod(nZ, 2) == 0
        nZ = nZ + 1;
    end

    if numel(sim.z) > 1
        currentDz = max(abs(diff(sim.z)));
    else
        currentDz = Inf;
    end

    needsGridUpdate = currentDz > 1.25*dzTarget || ...
        min(sim.z) > zMin || max(sim.z) < zMax;

    if needsGridUpdate
        sim.z = linspace(zMin, zMax, nZ);
        sim.nz = numel(sim.z);
        sim.nzRange = zMax - zMin;
        sim.obj = beadObject3D(sim);

        if getOption(opts, 'verbose', true)
            fprintf('[phaseRetrieval3DBead] model z grid: %d planes, %.3f to %.3f um\n', ...
                sim.nz, sim.z(1), sim.z(end));
        end
    end
end

function sim = matchSimLateralGridToData(sim, raw4, opts)
    ny = size(raw4, 1);
    nx = size(raw4, 2);

    if ny ~= nx
        error('phaseRetrieval3DBead:NonSquareData', ...
            'The current forward model expects square y-x images.');
    end

    pixelSizeUm = getOption(opts, 'xyPixelSizeUm', []);
    if ~isfinitePositiveScalarLocal(pixelSizeUm) && ...
            isfield(sim, 'xyPixelSizeUm') && isfinitePositiveScalarLocal(sim.xyPixelSizeUm)
        pixelSizeUm = sim.xyPixelSizeUm;
    end
    if ~isfinitePositiveScalarLocal(pixelSizeUm)
        if ~getOption(opts, 'allowUncalibratedFOV', false)
            error('phaseRetrieval3DBead:MissingLateralScale', ...
                ['Real-data fitting requires opts.xyPixelSizeUm or ' ...
                'opts.sim.xyPixelSizeUm. The default 1.8 um simulation FOV ' ...
                'must not be silently applied to measured data.']);
        end
        pixelSizeUm = sim.fovXY / max(nx-1, 1);
    end

    gridAlreadyMatches = numel(sim.y) == ny && numel(sim.x) == nx && ...
        numel(sim.x) > 1 && abs(abs(sim.x(2)-sim.x(1)) - pixelSizeUm) <= ...
        max(1e-12, 1e-6*pixelSizeUm);
    if gridAlreadyMatches
        sim.xyPixelSizeUm = pixelSizeUm;
        return;
    end

    if ~getOption(opts, 'autoXYGrid', true)
        currentPixelSizeUm = NaN;
        if numel(sim.x) > 1
            currentPixelSizeUm = abs(sim.x(2)-sim.x(1));
        end
        error('phaseRetrieval3DBead:GridMismatch', ...
            ['rawData is %dx%d at %.6g um/pixel, but sim.x/sim.y are ' ...
            '%dx%d at %.6g um/pixel.'], ny, nx, pixelSizeUm, ...
            numel(sim.y), numel(sim.x), currentPixelSizeUm);
    end

    sim.nx = nx;
    sim.ny = ny;
    sim.xyPixelSizeUm = pixelSizeUm;
    sim.fovX = pixelSizeUm * max(nx-1, 0);
    sim.fovY = pixelSizeUm * max(ny-1, 0);
    sim.fovXY = max(sim.fovX, sim.fovY);
    sim.x = ((1:nx) - (nx+1)/2) * pixelSizeUm;
    sim.y = ((1:ny) - (ny+1)/2) * pixelSizeUm;
    sim.dx = pixelSizeUm;
    sim.obj = beadObject3D(sim);

    if getOption(opts, 'verbose', true)
        fprintf(['[phaseRetrieval3DBead] lateral grid matched to data: ' ...
            '%d x %d pixels at %.3f nm/pixel\n'], ny, nx, 1000*pixelSizeUm);
    end
end

function tf = isfinitePositiveScalarLocal(value)
    tf = isnumeric(value) && isscalar(value) && isfinite(value) && value > 0;
end

function coeffs = defaultTrueCoeffs()
    coeffs = struct( ...
        'defocus', 0.035, ...
        'astig_x', 0.025, ...
        'coma_y', -0.020, ...
        'spherical', 0.015);
end

function names = buildParamNames(fitModes, fitXY, fitZ)
    names = fitModes;
    if fitXY
        names = [names {'x_shift','y_shift'}];
    end
    if fitZ
        names = [names {'z_offset'}];
    end
end

function p = initialParameterVector(sim, fitModes, opts)
    coeffs0 = getOption(opts, 'initialCoeffs', struct());
    coeffs0 = coeffStruct(sim, coeffs0);

    p = zeros(1, numel(fitModes));
    for k = 1:numel(fitModes)
        if isfield(coeffs0, fitModes{k})
            p(k) = coeffs0.(fitModes{k});
        end
    end

    if opts.fitXY
        xy0 = getOption(opts, 'initialXY', [0 0]);
        p = [p xy0(1) xy0(2)];
    end

    if opts.fitZ
        z0 = getOption(opts, 'initialZOffset', 0);
        p = [p z0];
    end
end

function step = finiteDifferenceSteps(nModes, fitXY, fitZ, opts)
    step = opts.fdCoeff * ones(1, nModes);
    if fitXY
        step = [step opts.fdXY opts.fdXY];
    end
    if fitZ
        step = [step opts.fdZ];
    end
end

function reg = regularizationVector(nModes, fitXY, fitZ, opts)
    reg = opts.regCoeff * ones(1, nModes);
    if fitXY
        reg = [reg opts.regXY opts.regXY];
    end
    if fitZ
        reg = [reg opts.regZ];
    end
end

function maxStep = maxUpdateVector(nModes, fitXY, fitZ, opts)
    maxStep = opts.maxCoeffStep * ones(1, nModes);
    if fitXY
        maxStep = [maxStep opts.maxXYStep opts.maxXYStep];
    end
    if fitZ
        maxStep = [maxStep opts.maxZStep];
    end
end

function stack = modelFromParams(sim, fitModes, planeZ, p, fitXY, fitZ)
    [coeffs, xy, zOffset] = unpackParams(sim, fitModes, p, fitXY, fitZ);
    stack = normalizedStackExplicitDetectorZPlanes( ...
        sim, coeffs, planeZ, xy(1), xy(2), zOffset);
end

function [coeffs, xy, zOffset] = unpackParams(sim, fitModes, p, fitXY, fitZ)
    fullVec = zeros(1, numel(sim.modeOrder));

    for k = 1:numel(fitModes)
        idx = find(strcmp(sim.modeOrder, fitModes{k}), 1, 'first');
        if isempty(idx)
            error('phaseRetrieval3DBead:UnknownMode', ...
                'Unknown fit mode "%s".', fitModes{k});
        end
        fullVec(idx) = p(k);
    end

    coeffs = coeffStruct(sim, fullVec);
    next = numel(fitModes) + 1;

    xy = [0 0];
    if fitXY
        xy = [p(next) p(next+1)];
        next = next + 2;
    end

    zOffset = 0;
    if fitZ
        zOffset = p(next);
    end
end

function v = coeffStructToVector(sim, coeffs)
    coeffs = coeffStruct(sim, coeffs);
    v = zeros(1, numel(sim.modeOrder));

    for k = 1:numel(sim.modeOrder)
        name = sim.modeOrder{k};
        if isfield(coeffs, name)
            v(k) = coeffs.(name);
        end
    end
end

function raw4 = standardizeSpadData(rawData, sim, nPlane)
    nCh = size(sim.detXY, 1);

    if iscell(rawData)
        if numel(rawData) ~= nCh
            error('phaseRetrieval3DBead:BadCellData', ...
                'Cell-array rawData must contain one volume per detector channel.');
        end

        first = rawData{1};
        if ndims(first) == 2
            first = reshape(first, size(first,1), size(first,2), 1);
        end

        [ny, nx, nz] = size(first);
        if nz ~= nPlane
            error('phaseRetrieval3DBead:PlaneMismatch', ...
                'rawData has %d planes, but planeZ has %d values.', nz, nPlane);
        end

        raw4 = zeros(ny, nx, nCh, nPlane);
        for k = 1:nCh
            vol = rawData{k};
            if ndims(vol) == 2
                vol = reshape(vol, size(vol,1), size(vol,2), 1);
            end
            if size(vol,1) ~= ny || size(vol,2) ~= nx || size(vol,3) ~= nPlane
                error('phaseRetrieval3DBead:BadCellVolume', ...
                    'All channel volumes must have size [ny nx nPlane].');
            end
            raw4(:,:,k,:) = reshape(vol, ny, nx, 1, nPlane);
        end
        return;
    end

    dims = size(rawData);
    if ndims(rawData) == 3
        if dims(3) ~= nCh
            error('phaseRetrieval3DBead:BadDataShape', ...
                '3D rawData must have size [ny nx %d].', nCh);
        end
        if nPlane ~= 1
            error('phaseRetrieval3DBead:PlaneMismatch', ...
                '3D rawData has one plane, but planeZ has %d values.', nPlane);
        end
        raw4 = reshape(rawData, dims(1), dims(2), nCh, 1);
        return;
    end

    if ndims(rawData) ~= 4
        error('phaseRetrieval3DBead:BadDataShape', ...
            'rawData must be [ny nx %d nPlane], [ny nx nPlane %d], or a %d-cell array.', ...
            nCh, nCh, nCh);
    end

    if dims(3) == nCh && dims(4) == nPlane
        raw4 = rawData;
    elseif dims(4) == nCh && dims(3) == nPlane
        raw4 = permute(rawData, [1 2 4 3]);
    else
        error('phaseRetrieval3DBead:BadDataShape', ...
            'Could not match rawData dimensions to %d channels and %d z planes.', nCh, nPlane);
    end
end

function [dataN, beta] = normalizeSpadVolume(data)
    border = [reshape(data(1,:,:,:), [], 1); ...
              reshape(data(end,:,:,:), [], 1); ...
              reshape(data(:,1,:,:), [], 1); ...
              reshape(data(:,end,:,:), [], 1)];

    beta = max(median(border), 0);
    d = max(data - beta, 0);
    s = sum(d(:));

    if s > 0
        dataN = d / s;
    else
        dataN = d;
    end
end
