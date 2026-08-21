function [hessian, ok] = numericalHessian2(objective, point, bounds)
    import membrane_tracking.curved_miet.*

    point = double(point(:).');
    baseStep = [max(1e-5, 1e-3*max(abs(point(1)),0.1)), ...
        max(1e-4, 1e-3*max(abs(point(2)),0.1))];
    distanceToBounds = [min(point(1)-bounds(1,1), ...
        bounds(1,2)-point(1)), min(point(2)-bounds(2,1), ...
        bounds(2,2)-point(2))];
    step = min(baseStep, 0.2 * distanceToBounds);
    hessian = nan(2);
    ok = all(step > 10*eps) && all(isfinite(step));
    if ~ok
        return;
    end

    e1 = [step(1) 0];
    e2 = [0 step(2)];
    f0 = objective(point);
    f1p = objective(point+e1);
    f1m = objective(point-e1);
    f2p = objective(point+e2);
    f2m = objective(point-e2);
    fpp = objective(point+e1+e2);
    fpm = objective(point+e1-e2);
    fmp = objective(point-e1+e2);
    fmm = objective(point-e1-e2);

    hessian(1,1) = (f1p - 2*f0 + f1m) / step(1)^2;
    hessian(2,2) = (f2p - 2*f0 + f2m) / step(2)^2;
    hessian(1,2) = (fpp - fpm - fmp + fmm) / ...
        (4 * step(1) * step(2));
    hessian(2,1) = hessian(1,2);
    hessian = 0.5 * (hessian + hessian.');
    ok = all(isfinite(hessian(:)));
end
