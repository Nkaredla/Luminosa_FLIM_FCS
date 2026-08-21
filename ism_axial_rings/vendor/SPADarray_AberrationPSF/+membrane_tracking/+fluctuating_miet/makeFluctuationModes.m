function modes = makeFluctuationModes(opts)
    import membrane_tracking.fluctuating_miet.*

%MAKEFLUCTUATIONMODES Helfrich mode basis over the observed patch.
%
%   Half-plane cosine/sine basis. Var(a_m) = Var(b_m) = S(q_m) with
%       S(q) = 2 / ( A * ( kappaBend*q^4 + tension*q^2 ) )     [kT = 1]
%   so that Var(dh) = sum_m S(q_m). Only modes resolvable across the patch
%   are retained: the longest wavelength is the patch diameter (longer
%   modes are pure piston and are absorbed into h0), and the shortest is
%   set by nModeShells.

    modes = struct();
    patchRadius = opts.membraneRadiusUm;
    patchArea = pi * patchRadius^2;
    qFundamental = pi / patchRadius;

    if ~opts.enableFluctuations
        modes.qVectors = zeros(0, 2);
        modes.variance = zeros(0, 1);
        modes.relaxationTimeS = zeros(0, 1);
        modes.nModes = 0;
        modes.fieldRmsUm = 0;
        modes.patchAreaUm2 = patchArea;
        modes.qFundamentalPerUm = qFundamental;
        modes.confoundedWithCurvature = zeros(0, 1);
        return;
    end

    n = opts.nModeShells;
    [ix, iy] = meshgrid(-n:n, -n:n);
    ix = ix(:); iy = iy(:);
    % Retain one representative per +/- pair: the half plane.
    keep = (ix > 0) | (ix == 0 & iy > 0);
    ix = ix(keep); iy = iy(keep);
    qVectors = qFundamental * [ix, iy];
    qMagnitude = hypot(qVectors(:,1), qVectors(:,2));
    inBand = qMagnitude <= qFundamental * n + 1e-9;
    qVectors = qVectors(inBand, :);
    qMagnitude = qMagnitude(inBand);

    variance = 2 ./ (patchArea * (opts.bendingModulusKT * qMagnitude.^4 + ...
        opts.tensionKTPerUm2 * qMagnitude.^2));

    % Guard against an unphysically soft membrane producing a field larger
    % than the mean corrugation, which would invalidate the height-graph
    % (small-slope) representation.
    rms = sqrt(sum(variance));
    if rms > opts.maxFieldRmsUm
        variance = variance * (opts.maxFieldRmsUm / rms)^2;
        rms = opts.maxFieldRmsUm;
        warning('simulateFluctuatingMIETMembraneTracking:FieldRescaled', ...
            ['The Helfrich spectrum implied an rms height of more than ' ...
             'maxFieldRmsUm; mode variances were rescaled to keep the ' ...
             'small-slope surface model valid.']);
    end

    % tau(q) ~ 1/(kappaBend*q^3 + tension*q), normalised at the fundamental.
    rate = opts.bendingModulusKT * qMagnitude.^3 + ...
        opts.tensionKTPerUm2 * qMagnitude;
    rateFundamental = opts.bendingModulusKT * qFundamental^3 + ...
        opts.tensionKTPerUm2 * qFundamental;
    relaxationTime = opts.fluctuationRelaxationTimeS * ...
        rateFundamental ./ max(rate, realmin);

    modes.qVectors = qVectors;
    modes.qMagnitudePerUm = qMagnitude;
    modes.variance = variance;
    modes.relaxationTimeS = relaxationTime;
    modes.nModes = size(qVectors, 1);
    modes.fieldRmsUm = rms;
    modes.patchAreaUm2 = patchArea;
    modes.qFundamentalPerUm = qFundamental;
    % A mode is quadratic across the observed disk, and therefore
    % confounded with kappa, when q*R is of order one.
    modes.confoundedWithCurvature = ...
        qMagnitude * opts.maxLocalizationRadiusUm < pi;
    modes.confoundedVarianceFraction = ...
        sum(variance(modes.confoundedWithCurvature)) / max(sum(variance), realmin);
end
