function transform = tangentTransformFromGradient(gradient, metricScale)
    import membrane_tracking.fluctuating_miet.*

    gradientNorm = norm(gradient);
    if gradientNorm > 10*eps
        uphill = gradient/gradientNorm;
    else
        uphill = [1 0];
    end
    crossSlope = [-uphill(2), uphill(1)];
    transform = [sqrt(metricScale)*uphill; crossSlope];
end
