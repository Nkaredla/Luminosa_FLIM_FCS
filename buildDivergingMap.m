function map = buildDivergingMap()
%BUILDDIVERGINGMAP Blue-white-red map for signed quantities.
%
% Written out rather than using a toolbox colormap so the prior-pull panel reads
% correctly: a signed quantity on a sequential map hides the sign, which is the
% only thing that panel is for.
    n = 128;
    ramp = linspace(0, 1, n)';
    lower = [ramp * 0 + 0.23 + 0.77 * ramp, ...
             ramp * 0 + 0.30 + 0.70 * ramp, ...
             ramp * 0 + 0.75 + 0.25 * ramp];
    upper = [ones(n, 1), 1 - 0.70 * ramp, 1 - 0.75 * ramp];
    map = [lower; upper];
    map = min(max(map, 0), 1);
end
