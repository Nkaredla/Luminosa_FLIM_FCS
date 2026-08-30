function anchor = fit_miet_slb_anchor(params)
%FIT_MIET_SLB_ANCHOR Use the bare-SLB lifetime to pin the MIET photophysics.
%
% anchor = fit_miet_slb_anchor(params)
%
% The supported bilayer is the one object in a MIET experiment whose height is
% known independently: it is a single labelled layer sitting a fixed, small
% distance above the spacer. Its measured lifetime is therefore a calibration
% POINT, not just a fitted number, and it can be used to determine the dye
% photophysics that the calibration curve otherwise has to assume.
%
% WHAT ONE ANCHOR CAN AND CANNOT DO
%
% Near a metal the decay rate of a dipole is
%
%     1/tau(z) = (1/tau_0) * [ (1 - q) + q * S(z, f) / (4/3 n) ]
%
% with q the free-space quantum yield and f = <cos^2 theta> the mean squared
% projection of the dipole on the surface normal. Writing R = S/(4/3 n),
%
%     q = (tau_0/tau_SLB - 1) / (R(f) - 1)
%
% so ONE measured lifetime fixes ONE combination of q and f, not both. This
% routine therefore does not return "the" quantum yield. It returns the whole
% feasible locus q(f), and solves for q at whichever f the caller commits to.
% Presenting it as a single fitted number would hide a real degeneracy: a dim
% dye standing up and a bright dye lying down produce the same SLB lifetime,
% and they produce DIFFERENT calibration curves higher up.
%
% The locus is bounded. R(f) is largest for a vertical dipole, so
%
%     tau_SLB >= tau_0 / R(1)
%
% is a hard floor: a lifetime shorter than that cannot be produced by this
% stack at ANY orientation or quantum yield, and it means the geometry (spacer
% thickness, SLB height) or tau_0 is wrong rather than the photophysics.
%
% Read the other way round, the same relation caps tau_0. Since q <= 1,
%
%     tau_0 <= tau_SLB * R(f)
%
% so the SLB measurement is an INDEPENDENT UPPER BOUND on the free-space
% lifetime, reported as tauFreeMaxNs. That is worth having because tau_0 sets
% the whole calibration's ceiling and is otherwise measured somewhere else
% entirely (the axial cross-sections). If a tau_0 taken from elsewhere exceeds
% this bound, the two measurements are inconsistent and one of them is wrong -
% which is a much more useful thing to learn than a silently rescaled height
% map. It also rules out reading tau_0 off the upper tail of the cell's own
% lifetime map, which would be circular anyway: defining tau_0 as the longest
% lifetime present guarantees that nothing exceeds the ceiling, so the
% "above ceiling" diagnostic would report success by construction.
%
% REQUIRED FIELDS
%   lambdaNm, n0, d0, n, d, n1, d1   the layer stack, as in miet_calibration_curve
%   tauFreeNs                        free-space lifetime [ns]
%   slbHeightNm                      height of the labelled SLB layer [nm]
%   slbLifetimeNs                    measured SLB lifetime [ns]; may be a
%                                    vector, one entry per acquisition
%
% OPTIONAL FIELDS
%   meanCosSquaredTheta   the f to solve q at (default 1/3, isotropic)
%   quantumYield          if given, f is solved for instead of q
%   orderGrid             f values for the reported locus, default 0:0.01:1
%   labels                cell array of names for the slbLifetimeNs entries
%   verbose               true; set false to suppress the printed report
%
% OUTPUT
%   anchor.quantumYield          q at the requested f, per lifetime
%   anchor.feasible              q lies in (0, 1]
%   anchor.minAttainableLifetimeNs   tau_0/R(1), the hard floor
%   anchor.parallelLimitLifetimeNs   tau_0/R(0), the floor if strictly in-plane
%   anchor.tauFreeMaxNs          tau_SLB*R(f), the largest tau_0 compatible
%                                with q <= 1 at the requested orientation
%   anchor.tauFreeMaxByOrientationNs  the same bound for f = 0, 1/3 and 1
%   anchor.locus                 table of f versus required q, per lifetime
%   anchor.impliedOrder          f at the requested quantumYield, if given
%
% See also MIET_CALIBRATION_CURVE, MAKE_RT_MIET_CALIBRATION.

    params = fillDefaults(params);
    if exist('LifetimeLSimpsExp', 'file') ~= 2
        luminosa_miet_vendor_root();
    end

    n0 = params.n0(:).';
    d0 = params.d0(:).';
    n1 = params.n1(:).';
    d1 = params.d1(:).';
    k = 2 * pi / double(params.lambdaNm);

    [~, ~, ~, ~, qvd, qvu, qpd, qpu] = LifetimeLSimpsExp( ...
        k * double(params.slbHeightNm), n0, params.n, n1, ...
        k * d0, k * params.d, k * d1);
    freeRate = 4/3 * params.n;
    rVertical = (qvu(1) + qvd(1)) / freeRate;
    rParallel = (qpu(1) + qpd(1)) / freeRate;

    tauFree = double(params.tauFreeNs);
    tauSlb = double(params.slbLifetimeNs(:));
    ratio = tauFree ./ tauSlb;              % required total rate enhancement

    f = double(params.meanCosSquaredTheta);
    rAtF = orderRate(rVertical, rParallel, f);
    quantumYield = (ratio - 1) ./ (rAtF - 1);

    anchor = struct();
    anchor.slbHeightNm = params.slbHeightNm;
    anchor.slbLifetimeNs = tauSlb;
    anchor.tauFreeNs = tauFree;
    anchor.labels = params.labels;
    anchor.rVertical = rVertical;
    anchor.rParallel = rParallel;
    anchor.meanCosSquaredTheta = f;
    anchor.rAtRequestedOrder = rAtF;
    anchor.quantumYield = quantumYield;
    anchor.feasible = quantumYield > 0 & quantumYield <= 1;
    % q <= 1 turns the anchor into an upper bound on tau_0, per acquisition.
    anchor.tauFreeMaxNs = tauSlb * rAtF;
    anchor.tauFreeMaxByOrientationNs = struct( ...
        'inPlane', tauSlb * rParallel, ...
        'isotropic', tauSlb * orderRate(rVertical, rParallel, 1/3), ...
        'vertical', tauSlb * rVertical);
    anchor.tauFreeExceedsAnchorBound = tauFree > anchor.tauFreeMaxNs;
    anchor.minAttainableLifetimeNs = tauFree / rVertical;
    anchor.parallelLimitLifetimeNs = tauFree / rParallel;
    anchor.isotropicLimitLifetimeNs = ...
        tauFree / orderRate(rVertical, rParallel, 1/3);

    % The locus: for every orientation, the quantum yield the anchor demands.
    grid = double(params.orderGrid(:));
    rGrid = orderRate(rVertical, rParallel, grid);
    anchor.orderGrid = grid;
    anchor.quantumYieldLocus = (ratio.' - 1) ./ (rGrid - 1);   % nF x nTau
    anchor.locusFeasible = anchor.quantumYieldLocus > 0 & ...
        anchor.quantumYieldLocus <= 1;
    % Smallest f that keeps every acquisition physical: below it the anchor
    % demands a quantum yield above one somewhere in the session.
    allFeasible = all(anchor.locusFeasible, 2);
    if any(allFeasible)
        anchor.minimumFeasibleOrder = grid(find(allFeasible, 1));
    else
        anchor.minimumFeasibleOrder = NaN;
    end

    if ~isempty(params.quantumYield)
        q = double(params.quantumYield);
        rNeeded = 1 + (ratio - 1) ./ q;
        impliedOrder = (rNeeded - rParallel) ./ (rVertical - rParallel);
        anchor.requestedQuantumYield = q;
        anchor.impliedOrder = impliedOrder;
        anchor.impliedOrderFeasible = impliedOrder >= 0 & impliedOrder <= 1;
    end

    if params.verbose
        printAnchor(anchor);
    end
end

function r = orderRate(rVertical, rParallel, f)
% Fast wobble averages RATES, so the orientation enters only through
% <cos^2 theta>. Linear in f by construction, which is why the locus can be
% inverted in closed form in both directions.
    r = f * rVertical + (1 - f) * rParallel;
end

function params = fillDefaults(params)
    required = {'lambdaNm', 'n0', 'd0', 'n', 'd', 'n1', 'tauFreeNs', ...
        'slbHeightNm', 'slbLifetimeNs'};
    for index = 1:numel(required)
        if ~isfield(params, required{index}) || isempty(params.(required{index}))
            error('fit_miet_slb_anchor:MissingField', ...
                'params.%s is required.', required{index});
        end
    end
    if ~isfield(params, 'd1'); params.d1 = []; end
    if ~isfield(params, 'meanCosSquaredTheta') || isempty(params.meanCosSquaredTheta)
        params.meanCosSquaredTheta = 1/3;
    end
    if ~isfield(params, 'quantumYield'); params.quantumYield = []; end
    if ~isfield(params, 'orderGrid') || isempty(params.orderGrid)
        params.orderGrid = 0:0.01:1;
    end
    if ~isfield(params, 'verbose') || isempty(params.verbose)
        params.verbose = true;
    end
    if ~isfield(params, 'labels') || isempty(params.labels)
        params.labels = arrayfun(@(i) sprintf('#%d', i), ...
            1:numel(params.slbLifetimeNs), 'UniformOutput', false);
    end
end

function printAnchor(anchor)
    fprintf('fit_miet_slb_anchor: SLB at %.1f nm, tau_0 = %.3f ns\n', ...
        anchor.slbHeightNm, anchor.tauFreeNs);
    fprintf(['  rate enhancement at the SLB: parallel %.3f, isotropic %.3f, ' ...
        'vertical %.3f\n'], anchor.rParallel, ...
        orderRate(anchor.rVertical, anchor.rParallel, 1/3), anchor.rVertical);
    fprintf(['  shortest SLB lifetime this stack can ever give: %.4f ns ' ...
        '(vertical, qy = 1)\n'], anchor.minAttainableLifetimeNs);
    fprintf('  if strictly in-plane, nothing below %.4f ns is reachable\n', ...
        anchor.parallelLimitLifetimeNs);
    fprintf(['  largest tau_0 compatible with qy <= 1 at this orientation: ' ...
        '%.4f ns\n'], min(anchor.tauFreeMaxNs));
    if any(anchor.tauFreeExceedsAnchorBound)
        fprintf(2, ['  tau_0 = %.4f ns EXCEEDS that bound for %d of %d ' ...
            'acquisitions: the independently measured tau_0 and the SLB ' ...
            'anchor disagree.\n'], anchor.tauFreeNs, ...
            nnz(anchor.tauFreeExceedsAnchorBound), ...
            numel(anchor.tauFreeExceedsAnchorBound));
    end
    fprintf('  solving at <cos^2 theta> = %.4f\n', anchor.meanCosSquaredTheta);
    for index = 1:numel(anchor.slbLifetimeNs)
        flag = '';
        if ~anchor.feasible(index); flag = '   <-- UNPHYSICAL'; end
        fprintf('    %-14s tau_SLB %.4f ns  ->  qy = %.4f%s\n', ...
            anchor.labels{index}, anchor.slbLifetimeNs(index), ...
            anchor.quantumYield(index), flag);
    end
    if isfinite(anchor.minimumFeasibleOrder)
        fprintf(['  every acquisition stays physical only for ' ...
            '<cos^2 theta> >= %.3f\n'], anchor.minimumFeasibleOrder);
    else
        fprintf(2, ['  NO orientation keeps every acquisition physical: the ' ...
            'geometry or tau_0 is wrong, not the photophysics.\n']);
    end
end
