function beta = poisson_nnls2_exact(G, b)
%POISSON_NNLS2_EXACT Exact non-negative least squares for TWO parameters.
%
% beta = poisson_nnls2_exact(G, b)
%
% G is 3-by-P holding each pixel's symmetric 2-by-2 normal matrix as
% [g11; g12; g22]; b is 2-by-P. Returns the exact minimiser of
%
%     q(beta) = beta' * G * beta - 2 * beta' * b     subject to  beta >= 0
%
% column by column, with no iteration.
%
% Same argument as poisson_nnls3_exact, one dimension smaller: the constrained
% minimiser of a strictly convex quadratic sits on some face of the non-negative
% quadrant and is the unconstrained minimiser there, and with two parameters
% there are only four faces. Used for a background-plus-single-exponential fit,
% where the two parameters are [B; a].

    P = size(G, 2);
    g11 = G(1, :); g12 = G(2, :); g22 = G(3, :);
    b1 = b(1, :); b2 = b(2, :);

    scale = max(max(abs(G), [], 1), eps);
    guard = (1e-12 * scale) .^ 2;

    beta = zeros(2, P);
    best = zeros(1, P);              % the origin, always feasible, q = 0

    % ---- face {1} : background only -------------------------------------
    v1 = b1 ./ max(g11, eps);
    v = [v1; zeros(1, P)];
    ok = g11 > 0 & all(v >= 0, 1) & all(isfinite(v), 1);
    q = v(1, :) .^ 2 .* g11 - 2 * v(1, :) .* b1;
    take = ok & q < best;
    best(take) = q(take); beta(:, take) = v(:, take);

    % ---- face {2} : exponential only ------------------------------------
    v2 = b2 ./ max(g22, eps);
    v = [zeros(1, P); v2];
    ok = g22 > 0 & all(v >= 0, 1) & all(isfinite(v), 1);
    q = v(2, :) .^ 2 .* g22 - 2 * v(2, :) .* b2;
    take = ok & q < best;
    best(take) = q(take); beta(:, take) = v(:, take);

    % ---- face {1,2} : both free -----------------------------------------
    det2 = g11 .* g22 - g12 .^ 2;
    ok = abs(det2) > guard;
    den = det2; den(~ok) = 1;
    v = [(g22 .* b1 - g12 .* b2) ./ den; (g11 .* b2 - g12 .* b1) ./ den];
    v(:, ~ok) = -1;
    feasible = all(v >= 0, 1) & all(isfinite(v), 1);
    q = v(1, :) .^ 2 .* g11 + v(2, :) .^ 2 .* g22 ...
        + 2 * v(1, :) .* v(2, :) .* g12 ...
        - 2 * (v(1, :) .* b1 + v(2, :) .* b2);
    take = feasible & q < best;
    beta(:, take) = v(:, take);
end
