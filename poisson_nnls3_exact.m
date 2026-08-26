function beta = poisson_nnls3_exact(G, b)
%POISSON_NNLS3_EXACT Exact non-negative least squares for THREE parameters.
%
% beta = poisson_nnls3_exact(G, b)
%
% G is 6-by-P holding each pixel's symmetric 3-by-3 normal matrix as
% [g11; g12; g13; g22; g23; g33]; b is 3-by-P. Returns the exact minimiser of
%
%     q(beta) = beta' * G * beta - 2 * beta' * b     subject to  beta >= 0
%
% for every column independently, with no iteration and no call to lsqnonneg.
%
% WHY ENUMERATION IS EXACT
%
% For a strictly convex quadratic the constrained minimiser lies on some face of
% the non-negative octant, and ON that face it is the UNCONSTRAINED minimiser
% over the free coordinates. With three parameters there are only 2^3 = 8 faces,
% so solving all eight in closed form and keeping the feasible one of lowest
% objective is exact.
%
% The point of doing it this way is that it VECTORISES. lsqnonneg cannot be
% batched over pixels, and that is the practical reason the original solver
% reached for the Gram-matrix shortcut - lsqnonneg(x'*w*x, x'*w*y) - which
% minimises the normal-equation residual instead of the weighted residual and is
% therefore wrong exactly when a bound is active. Here the constrained solve is
% eight closed-form 3-by-3 solves across all pixels at once, so no shortcut is
% needed.
%
% Each face is evaluated by masking the normal matrix to the identity on its
% inactive coordinates, which turns every face into the same 3-by-3 solve and
% forces the inactive components to zero. Near-singular faces are rejected by a
% determinant guard rather than being allowed to emit a huge spurious beta, and
% the origin is always feasible, so no column can come back empty.

    P = size(G, 2);
    g11 = G(1, :); g12 = G(2, :); g13 = G(3, :);
    g22 = G(4, :); g23 = G(5, :); g33 = G(6, :);
    b1 = b(1, :); b2 = b(2, :); b3 = b(3, :);

    scale = max(max(abs(G), [], 1), eps);
    guard = (1e-12 * scale) .^ 3;

    beta = zeros(3, P);
    best = zeros(1, P);        % the origin: q = 0, always feasible

    faces = logical([ ...
        1 0 0; 0 1 0; 0 0 1; ...
        1 1 0; 1 0 1; 0 1 1; ...
        1 1 1]);

    for k = 1:size(faces, 1)
        m = faces(k, :);

        % Mask to identity on the inactive coordinates.
        if m(1); h11 = g11; else; h11 = ones(1, P); end
        if m(2); h22 = g22; else; h22 = ones(1, P); end
        if m(3); h33 = g33; else; h33 = ones(1, P); end
        if m(1) && m(2); h12 = g12; else; h12 = zeros(1, P); end
        if m(1) && m(3); h13 = g13; else; h13 = zeros(1, P); end
        if m(2) && m(3); h23 = g23; else; h23 = zeros(1, P); end
        if m(1); c1 = b1; else; c1 = zeros(1, P); end
        if m(2); c2 = b2; else; c2 = zeros(1, P); end
        if m(3); c3 = b3; else; c3 = zeros(1, P); end

        % 3-by-3 symmetric solve by the adjugate.
        a11 = h22 .* h33 - h23 .^ 2;
        a12 = h13 .* h23 - h12 .* h33;
        a13 = h12 .* h23 - h13 .* h22;
        a22 = h11 .* h33 - h13 .^ 2;
        a23 = h13 .* h12 - h11 .* h23;
        a33 = h11 .* h22 - h12 .^ 2;
        det3 = h11 .* a11 + h12 .* a12 + h13 .* a13;

        ok = abs(det3) > guard;
        den = det3;
        den(~ok) = 1;
        v = [(a11 .* c1 + a12 .* c2 + a13 .* c3) ./ den; ...
             (a12 .* c1 + a22 .* c2 + a23 .* c3) ./ den; ...
             (a13 .* c1 + a23 .* c2 + a33 .* c3) ./ den];
        v(~m, :) = 0;                       % inactive coordinates are exactly 0
        v(:, ~ok) = -1;                     % rejected: fails the feasibility test

        feasible = all(v >= 0, 1) & all(isfinite(v), 1);
        if ~any(feasible); continue; end

        % Objective against the ORIGINAL G and b, not the masked ones.
        q = v(1, :) .^ 2 .* g11 + v(2, :) .^ 2 .* g22 + v(3, :) .^ 2 .* g33 ...
            + 2 * (v(1, :) .* v(2, :) .* g12 + v(1, :) .* v(3, :) .* g13 ...
            + v(2, :) .* v(3, :) .* g23) ...
            - 2 * (v(1, :) .* b1 + v(2, :) .* b2 + v(3, :) .* b3);

        take = feasible & (q < best);
        if any(take)
            best(take) = q(take);
            beta(:, take) = v(:, take);
        end
    end
end
