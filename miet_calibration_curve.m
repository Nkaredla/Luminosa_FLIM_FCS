function calib = miet_calibration_curve(params)
%MIET_CALIBRATION_CURVE Lifetime-versus-height calibration from a layer stack.
%
% calib = miet_calibration_curve(params)
%
% This is the calculation the MIET-GUI performs when you press "calculate
% calibration curve", called programmatically so that a curve can be
% regenerated and version-controlled instead of being clicked out of a window
% and saved to an anonymous MAT file. It uses the vendored MIET-GUI kernel
% (LifetimeLSimpsExp) unchanged, so the numbers are the GUI's numbers.
%
% LAYER-ORDER CONVENTION - THE ONE THING THAT IS EASY TO GET WRONG
%
% LifetimeL/LifetimeLSimpsExp take n0 BOTTOMMOST-FIRST: n0(1) is the
% semi-infinite substrate, n0(end) is the layer touching the molecule's layer,
% and d0(k) is the thickness of n0(k+1), so numel(d0) == numel(n0)-1. The
% parameter WINDOW of the GUI lists the same stack the other way round
% (topmost first) and pads its tables with a trailing zero row, which is why a
% saved GUI settings file cannot be fed to the kernel verbatim. Getting this
% backwards does not error - it silently returns the curve for a mirrored
% stack. It was checked here by reproducing the existing calibrationCurve.mat
% from the Hela_MIET.mat settings: bottommost-first gives a 3.3% relative
% residual, topmost-first gives 38%.
%
% REQUIRED FIELDS
%   lambdaNm            emission wavelength in vacuum [nm]
%   n0                  refractive indices below the molecule, BOTTOMMOST first
%   d0                  thicknesses [nm] of n0(2:end)
%   n                   refractive index of the molecule's layer (real)
%   d                   thickness [nm] of the molecule's layer
%   n1                  refractive indices above the molecule, topmost LAST
%   d1                  thicknesses [nm] of n1(2:end); [] for a half-space
%   quantumYield        free-space quantum yield, 0-1 (>1 is read as percent)
%   tauFreeNs           free-space lifetime [ns]
%
% OPTIONAL FIELDS
%   dipoleOrientation   'parallel' (in-plane, e.g. a membrane probe lying in
%                       the bilayer), 'vertical', 'fast_rotating' (default,
%                       isotropic average), 'random_fixed', or 'order' - see
%                       meanCosSquaredTheta
%   meanCosSquaredTheta <cos^2(theta)> of the dipole against the surface
%                       NORMAL, used when dipoleOrientation is 'order'. This
%                       is the one continuous knob that contains all three
%                       named cases: 0 is parallel, 1/3 is isotropic
%                       (identical to 'fast_rotating'), 1 is vertical. A real
%                       membrane probe is neither of the extremes - it wobbles
%                       in a cone - and the named labels force a choice
%                       between two ends of a continuum that the data can
%                       actually constrain. Relates to the usual second-rank
%                       order parameter as S = (3*<cos^2 theta> - 1)/2.
%   zNm                 height grid [nm], default 0.1:0.1:400
%   curveType           'minimum' (default) keeps only heights whose lifetime
%                       is unique, 'maximum' keeps up to the first maximum,
%                       'manual' keeps the whole grid
%   label               free-text description stored with the curve
%
% OUTPUT
%   calib.heightNm / calib.lifetimeNs      the usable, truncated curve
%   calib.fullHeightNm / calib.fullLifetimeNs   before truncation
%   calib.unambiguousLifetimeNs            largest uniquely invertible lifetime
%   calib.maxHeightNm                      height at that lifetime
%   calib.curve                            [height lifetime], the two-column
%                                          layout load_lifetime_height_calibration
%                                          expects
%   calib.params                           the resolved inputs
%
% See also LUMINOSA_MIET_CALC, LOAD_LIFETIME_HEIGHT_CALIBRATION.

    params = fillDefaults(params);
    ensureMietKernelOnPath();

    n0 = params.n0(:).';
    d0 = params.d0(:).';
    n1 = params.n1(:).';
    d1 = params.d1(:).';
    if numel(d0) ~= numel(n0) - 1
        error('miet_calibration_curve:LowerStack', ...
            ['The lower stack needs one more refractive index than thickness ' ...
             '(numel(n0)=%d, numel(d0)=%d). Remember the GUI pads its tables ' ...
             'with a trailing zero row.'], numel(n0), numel(d0));
    end
    if numel(d1) ~= numel(n1) - 1
        error('miet_calibration_curve:UpperStack', ...
            ['The upper stack needs one more refractive index than thickness ' ...
             '(numel(n1)=%d, numel(d1)=%d).'], numel(n1), numel(d1));
    end
    if ~isreal(params.n)
        error('miet_calibration_curve:MoleculeIndex', ...
            'The molecule layer index must be real; LifetimeLSimpsExp requires it.');
    end

    z = double(params.zNm(:)).';
    k = 2 * pi / double(params.lambdaNm);
    [~, ~, ~, ~, qvd, qvu, qpd, qpu] = LifetimeLSimpsExp( ...
        k * z, n0, params.n, n1, k * d0, k * params.d, k * d1);

    sv = (qvu + qvd).';     % vertical dipole, total emitted power
    sp = (qpu + qpd).';     % parallel dipole, total emitted power
    lifetimeNs = orientedLifetime(sv, sp, params.n, params.quantumYield, ...
        params.tauFreeNs, params.dipoleOrientation, params.meanCosSquaredTheta);

    [usableHeight, usableLifetime, limitLifetime] = ...
        truncateCurve(z(:), lifetimeNs, params.curveType);

    calib = struct();
    calib.label = params.label;
    calib.heightNm = usableHeight;
    calib.lifetimeNs = usableLifetime;
    calib.curve = [usableHeight usableLifetime];
    calib.fullHeightNm = z(:);
    calib.fullLifetimeNs = lifetimeNs;
    calib.unambiguousLifetimeNs = limitLifetime;
    calib.maxHeightNm = usableHeight(end);
    calib.minLifetimeNs = usableLifetime(1);
    calib.maxLifetimeNs = usableLifetime(end);
    calib.quenchedSurfaceLifetimeNs = usableLifetime(1);
    calib.params = params;
    calib.kernel = 'LifetimeLSimpsExp (vendored MIET-GUI)';
    calib.createdOn = datestr(now, 'yyyy-mm-ddTHH:MM:SS'); %#ok<TNOW1,DATST>
end

function params = fillDefaults(params)
    if ~isstruct(params)
        error('miet_calibration_curve:Params', 'params must be a struct.');
    end
    required = {'lambdaNm', 'n0', 'd0', 'n', 'd', 'n1', 'quantumYield', 'tauFreeNs'};
    for index = 1:numel(required)
        if ~isfield(params, required{index}) || isempty(params.(required{index}))
            error('miet_calibration_curve:MissingField', ...
                'params.%s is required.', required{index});
        end
    end
    if ~isfield(params, 'd1'); params.d1 = []; end
    if ~isfield(params, 'dipoleOrientation') || isempty(params.dipoleOrientation)
        params.dipoleOrientation = 'fast_rotating';
    end
    if ~isfield(params, 'meanCosSquaredTheta'); params.meanCosSquaredTheta = []; end
    % 'order' is local, not a MIET-GUI label, so it must be recognised before
    % luminosa_miet_normalize_orientation folds every unknown string into
    % 'fast_rotating' and silently discards the order parameter.
    if strcmpi(char(params.dipoleOrientation), 'order')
        params.dipoleOrientation = 'order';
        if isempty(params.meanCosSquaredTheta)
            error('miet_calibration_curve:MeanCosSquared', ...
                'dipoleOrientation ''order'' requires params.meanCosSquaredTheta.');
        end
        if params.meanCosSquaredTheta < 0 || params.meanCosSquaredTheta > 1
            error('miet_calibration_curve:MeanCosSquared', ...
                '<cos^2 theta> must lie in [0, 1]; got %g.', ...
                params.meanCosSquaredTheta);
        end
    else
        params.dipoleOrientation = ...
            luminosa_miet_normalize_orientation(params.dipoleOrientation);
        params.meanCosSquaredTheta = namedOrderParameter(params.dipoleOrientation);
    end
    if ~isfield(params, 'zNm') || isempty(params.zNm)
        params.zNm = 0.1:0.1:400;
    end
    if ~isfield(params, 'curveType') || isempty(params.curveType)
        params.curveType = 'minimum';
    end
    if ~isfield(params, 'label'); params.label = ''; end
    if params.quantumYield > 1
        params.quantumYield = params.quantumYield / 100;
    end
    if params.quantumYield <= 0 || params.quantumYield > 1
        error('miet_calibration_curve:QuantumYield', ...
            'quantumYield must lie in (0, 1]; got %g.', params.quantumYield);
    end
    if params.tauFreeNs <= 0
        error('miet_calibration_curve:TauFree', ...
            'tauFreeNs must be positive; got %g.', params.tauFreeNs);
    end
end

function value = namedOrderParameter(mode)
% <cos^2 theta> equivalent of each named mode, recorded so that every curve
% carries the same continuous descriptor regardless of how it was requested.
% 'random_fixed' has no exact equivalent - it averages LIFETIMES, not rates -
% so it is reported as NaN rather than as a misleading 1/3.
    switch mode
        case 'parallel';      value = 0;
        case 'vertical';      value = 1;
        case 'fast_rotating'; value = 1/3;
        otherwise;            value = NaN;
    end
end

function lifetimeNs = orientedLifetime(sv, sp, n, qy, tauFree, mode, meanCosSq)
% Normalisation follows the kernel: the free-dipole emission rate is 4/3*n,
% and only the radiative part of the decay is modified, hence the (1-qy) term.
    switch mode
        case 'order'
            % Fast wobble within a cone: the RATES average over orientation,
            % so the mean of cos^2(theta) is the only thing the decay sees.
            % meanCosSq = 0 reproduces 'parallel', 1/3 reproduces
            % 'fast_rotating' exactly, and 1 reproduces 'vertical'.
            sr = meanCosSq * sv + (1 - meanCosSq) * sp;
            lifetimeNs = tauFree ./ ((1 - qy) + qy * sr / (4/3 * n));
        case 'vertical'
            lifetimeNs = tauFree ./ ((1 - qy) + qy * sv / (4/3 * n));
        case 'parallel'
            lifetimeNs = tauFree ./ ((1 - qy) + qy * sp / (4/3 * n));
        case 'random_fixed'
            % Static isotropic ensemble: average the LIFETIME over theta with
            % the sin(theta) solid-angle weight. This is NOT the same as the
            % fast-rotating case, where a single molecule samples all angles
            % within its lifetime and the RATES average instead.
            theta = linspace(0, pi/2, 181);
            weights = sin(theta);
            sr = sv * (cos(theta).^2) + sp * (sin(theta).^2);
            lifeByTheta = tauFree ./ ((1 - qy) + qy * sr / (4/3 * n));
            lifetimeNs = trapz(theta, lifeByTheta .* weights, 2) ./ trapz(theta, weights);
        otherwise
            lifetimeNs = tauFree ./ ((1 - qy) + qy * (sv + 2 * sp) / (4 * n));
    end
    lifetimeNs = lifetimeNs(:);
end

function [heightNm, lifetimeNs, limitLifetime] = truncateCurve(z, life, curveType)
% The MIET curve rises with height, turns over, and then oscillates. Only the
% stretch below the first turning point is invertible, and even there a
% lifetime above the first post-peak minimum recurs higher up, so the honest
% limit is that minimum, not the peak.
    peak = find(diff(life) < 0, 1);
    if isempty(peak)
        limitLifetime = life(end);
        heightNm = z;
        lifetimeNs = life;
        return;
    end
    switch lower(char(curveType))
        case 'maximum'
            limitLifetime = life(peak);
            keep = 1:peak;
        case 'manual'
            limitLifetime = min(life(peak:end));
            keep = 1:numel(z);
        otherwise
            limitLifetime = min(life(peak:end));
            firstRepeat = find(life > limitLifetime, 1);
            if isempty(firstRepeat)
                keep = 1:numel(z);
            else
                keep = 1:(firstRepeat - 1);
            end
    end
    heightNm = z(keep);
    lifetimeNs = life(keep);
end

function ensureMietKernelOnPath()
    if exist('LifetimeLSimpsExp', 'file') == 2
        return;
    end
    root = luminosa_miet_vendor_root();
    if isempty(root) || ~isfolder(root)
        error('miet_calibration_curve:Kernel', ...
            'The vendored MIET-GUI tree was not found; LifetimeLSimpsExp is unavailable.');
    end
    addpath(genpath(root));
    if exist('LifetimeLSimpsExp', 'file') ~= 2
        error('miet_calibration_curve:Kernel', ...
            'LifetimeLSimpsExp was not found under %s.', root);
    end
end
