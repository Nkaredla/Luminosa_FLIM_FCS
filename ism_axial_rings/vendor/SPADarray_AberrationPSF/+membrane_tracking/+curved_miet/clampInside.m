function value = clampInside(value, bounds)
    import membrane_tracking.curved_miet.*

    margin = 1e-6 * diff(bounds);
    value = min(max(value, bounds(1)+margin), bounds(2)-margin);
end
