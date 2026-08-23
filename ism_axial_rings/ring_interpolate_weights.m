function w = ring_interpolate_weights(weights, height)
%RING_INTERPOLATE_WEIGHTS Ring collection weights at an arbitrary height.
%
% w = ring_interpolate_weights(weights, height)
%
% weights.table    [nRing x nz] ring weight w(ring, z)
% weights.heights  [1 x nz] the heights the table is tabulated on, um
%
% The height is CLAMPED to the tabulated range. Clamping rather than
% extrapolating matters: linear extrapolation of a decaying weight curve goes
% negative, and a negative collection weight is not a physical model. Any
% caller that searches over height must be aware the objective is flat outside
% the table, so a fitted height sitting exactly on an endpoint means "at or
% beyond the edge of the model", not a measurement.

    if ~isstruct(weights) || ~isfield(weights, 'table') || ...
            ~isfield(weights, 'heights')
        error('ring_interpolate_weights:Weights', ...
            'weights must be a struct with fields table and heights.');
    end
    grid = weights.heights(:);
    height = min(max(height, grid(1)), grid(end));
    w = interp1(grid, weights.table', height, 'linear')';
    w = max(w(:), 0);
end
