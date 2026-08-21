function res = singleModeFit(sim, modeName, amp, x0, y0, photons, bg, maxIter)
%--------------------------------------------------------------------------
% singleModeFit
%
% PURPOSE
%   Simulate noisy detector-stack data for a single aberration mode and
%   retrieve that mode amplitude using a local Gauss-Newton fit.
%
% INPUTS
%   sim      : simulation structure
%   modeName : one mode from sim.modeOrder
%   amp      : true aberration amplitude
%   x0,y0    : true bead shift [um]
%   photons  : photon budget
%   bg       : background level
%   maxIter  : number of Gauss-Newton iterations
%
% OUTPUT
%   res : structure containing true and estimated quantities
%
% DESCRIPTION
%   The fit estimates:
%       - one selected mode amplitude
%       - x bead shift
%       - y bead shift
%
%   using finite-difference derivatives of the normalized detector stack.
%--------------------------------------------------------------------------

    % Default true bead shift and noise parameters
    if nargin < 4, x0 = 0.02; end
    if nargin < 5, y0 = -0.015; end
    if nargin < 6, photons = 2e5; end
    if nargin < 7, bg = 0.1; end
    if nargin < 8, maxIter = 12; end

    % Find the requested mode index
    midx = find(strcmp(sim.modeOrder, modeName));

    % True aberration coefficient vector
    trueA = zeros(1, numel(sim.modeOrder));
    trueA(midx) = amp;

    % Simulate noiseless normalized stack for the true parameters
    mTrue = normalizedStack(sim, coeffStruct(sim, trueA), x0, y0);

    % Add Poisson noise and constant background
    data = poissonSample(photons*mTrue + bg);

    % Normalize noisy data for fitting
    [dataN, beta] = normalizeDataStack(data);

    % Parameter vector:
    %   p(1:numModes) = aberration coefficients
    %   p(end-1)      = x shift
    %   p(end)        = y shift
    p = zeros(1, numel(sim.modeOrder)+2);

    % Finite-difference step for aberration coefficient
    fd = 0.01;

    % Finite-difference step for x,y position
    stepXY = sim.dx/4;

    % Gauss-Newton iterations
    for it = 1:maxIter

        % Current model prediction
        m0 = normalizedStack(sim, coeffStruct(sim, p(1:numel(sim.modeOrder))), p(end-1), p(end));

        % Residual vector
        r = dataN(:) - m0(:);

        % Jacobian for:
        %   column 1 = selected aberration mode
        %   column 2 = x shift
        %   column 3 = y shift
        J = zeros(numel(m0), 3);

        % Finite difference for aberration mode
        pp = p; pm = p;
        pp(midx)=pp(midx)+fd;
        pm(midx)=pm(midx)-fd;
        mp = normalizedStack(sim, coeffStruct(sim, pp(1:numel(sim.modeOrder))), pp(end-1), pp(end));
        mm = normalizedStack(sim, coeffStruct(sim, pm(1:numel(sim.modeOrder))), pm(end-1), pm(end));
        J(:,1) = (mp(:)-mm(:))/(2*fd);

        % Finite differences for x and y position
        for q = 1:2
            pp = p; pm = p;
            pp(end-2+q)=pp(end-2+q)+stepXY;
            pm(end-2+q)=pm(end-2+q)-stepXY;
            mp = normalizedStack(sim, coeffStruct(sim, pp(1:numel(sim.modeOrder))), pp(end-1), pp(end));
            mm = normalizedStack(sim, coeffStruct(sim, pm(1:numel(sim.modeOrder))), pm(end-1), pm(end));
            J(:,1+q) = (mp(:)-mm(:))/(2*stepXY);
        end

        % Regularized normal equations
        H = J.'*J + diag([1e-4 1e-5 1e-5]);
        g = J.'*r;

        % Guard against NaN / Inf in the system
        if any(~isfinite(H(:))) || any(~isfinite(g(:)))
            error('Non-finite values encountered in Gauss-Newton system.');
        end

        % Solve using pseudoinverse if matrix is nearly singular
        if rcond(H) < 1e-12
            delta = pinv(H) * g;
        else
            delta = H \ g;
        end

        % Update only the chosen mode and x,y shifts
        p(midx) = p(midx) + delta(1);
        p(end-1) = p(end-1) + delta(2);
        p(end) = p(end) + delta(3);
    end

    % Pack outputs
    res.modeName = modeName;
    res.trueAmp = amp;
    res.estAmp = p(midx);
    res.trueXY = [x0 y0];
    res.estXY = [p(end-1) p(end)];
    res.beta = beta;
end