function [vxx, vyy, vxy, determinant] = diffusionStepCovariance( ...
        D, noiseCov, stepDtS)
    import membrane_tracking.focused_ism.*

    brownianVariance = 2 * D * stepDtS;
    vxx = brownianVariance + reshape(noiseCov(1,1,:), [], 1);
    vyy = brownianVariance + reshape(noiseCov(2,2,:), [], 1);
    vxy = reshape(noiseCov(1,2,:), [], 1);
    determinant = max(vxx .* vyy - vxy.^2, realmin);
end
