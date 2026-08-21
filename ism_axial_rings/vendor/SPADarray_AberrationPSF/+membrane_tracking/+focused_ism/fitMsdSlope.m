function [D, intercept, beta] = fitMsdSlope(tau, y, nPairs)
    import membrane_tracking.focused_ism.*

    tau = tau(:);
    y = y(:);
    nPairs = nPairs(:);

    if numel(tau) < 2
        D = NaN;
        intercept = NaN;
        beta = [NaN; NaN];
        return;
    end

    X = [ones(numel(tau), 1), tau];
    w = sqrt(max(nPairs, 1) / max(sum(nPairs), 1));
    Xw = bsxfun(@times, X, w);
    yw = y .* w;
    beta = Xw \ yw;

    intercept = beta(1);
    D = beta(2) / 4;
end
