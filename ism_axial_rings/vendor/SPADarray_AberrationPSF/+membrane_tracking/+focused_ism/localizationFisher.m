function fisher = localizationFisher(A, b, prob, dPdx, dPdy, opts)
    import membrane_tracking.focused_ism.*

    mu = max(b + A * prob(:), opts.minExpectedCount);
    J = [A * dPdx(:), A * dPdy(:), prob(:), ones(numel(prob), 1)];
    weightedJ = bsxfun(@rdivide, J, mu);
    I = J.' * weightedJ;
    I = 0.5 * (I + I.');

    Ixy = I(1:2, 1:2);
    Inuis = I(3:4, 3:4);
    Icross = I(1:2, 3:4);
    [covNuis, nuisValid] = invertPositiveDefinite(Inuis);
    if nuisValid
        Ieff = Ixy - Icross * covNuis * Icross.';
        Ieff = 0.5 * (Ieff + Ieff.');
        [covXY, xyValid] = invertPositiveDefinite(Ieff);
    else
        Ieff = nan(2);
        covXY = nan(2);
        xyValid = false;
    end

    fisher = struct();
    fisher.informationFull = I;
    fisher.informationXYMarginal = Ieff;
    fisher.covarianceXY = covXY;
    fisher.isValid = nuisValid && xyValid;
end
