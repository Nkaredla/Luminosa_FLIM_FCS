function fits = immune_cell_MIET_explorer_fit(decay, irf, dtNs, periodNs)
%IMMUNE_CELL_MIET_EXPLORER_FIT Bi- and tri-exponential fits of a pooled decay.
%
% fits = immune_cell_MIET_explorer_fit(decay, irf, dtNs, periodNs)
%
% Fits the same two models the pipeline compares - two and three exponentials -
% to a pooled TCSPC decay, and returns both with their residuals so the GUI can
% draw them together. Returned struct array fields:
%
%   nExp        2 or 3
%   tauNs       fitted lifetimes, ascending
%   amplitude   species fractions matching tauNs, summing to 1
%   offset      fitted constant term
%   shiftNs     fitted IRF colour shift
%   model       the fitted curve, on the same bins as decay
%   residual    (decay - model) / sqrt(model), i.e. Poisson sigmas
%   chiSquare   reduced chi-square against Poisson counts
%   photons     total photons in the decay
%   ok          false if the fit failed, with the reason in .message
%
% THE MODEL COMES FROM FLUOFIT, IT IS NOT REBUILT HERE
%
% Fluofit's 8th output zz holds the fitted basis already scaled by the
% amplitudes, so sum(zz, 2) IS the fitted curve and the residual computed from
% it is the residual of the fit that actually happened.
%
% An earlier version rebuilt the curve from the returned lifetimes and
% amplitudes. That was wrong twice over and the errors were large: Fluofit
% shifts the IRF by a fitted COLOUR SHIFT (Fluofit.m:117) and NORMALISES each
% basis column before solving the amplitudes (Fluofit.m:143), so amplitudes
% belonging to a normalised, time-shifted basis were being applied to an
% un-normalised, unshifted one. Reduced chi-square came out at 47 to 1383 and
% residuals reached 125 sigma - not a bad fit, a mis-reconstructed curve.
%
% This is the second time in this project that drawing a model separately from
% the one that was fitted has produced a figure disagreeing with its own
% residuals, so the rule is now: never rebuild, always use the fitter's curve.
%
% WHY IRF RECONVOLUTION AND NOT A TAIL FIT
%
% The question the GUI exists to answer is whether a pixel needs two or three
% components, and the short component - the quenched SLB near 0.34 ns - is
% comparable to the IRF width. A tail fit discards exactly the early bins that
% carry it.
%
% BOTH MODELS ARE ALWAYS SHOWN, AND THE RESIDUALS ARE THE POINT
%
% Three exponentials can only reduce the residual, so a smaller chi-square is
% not evidence for three components. What distinguishes them is STRUCTURE in the
% two-exponential residual: if that is flat within noise, the third component is
% not needed for that pixel whatever the chi-squares say. This project has
% already measured how badly likelihood alone behaves here - the false
% three-component rate RISES with photon count, from 9% to 54% - so the residual
% panels are the honest display and the chi-square is reported without being
% interpreted.

    fits = struct('nExp', {}, 'tauNs', {}, 'amplitude', {}, 'offset', {}, ...
        'shiftNs', {}, 'model', {}, 'residual', {}, 'chiSquare', {}, ...
        'photons', {}, 'ok', {}, 'message', {});
    if exist('Fluofit', 'file') ~= 2
        error('immune_cell_MIET_explorer_fit:NoFluofit', ...
            'Fluofit.m must be on the MATLAB path.');
    end

    decay = double(decay(:));
    irf = double(irf(:));
    if numel(irf) ~= numel(decay)
        irf = interp1(linspace(0, 1, numel(irf))', irf, ...
            linspace(0, 1, numel(decay))', 'linear', 0);
    end
    irf = max(irf, 0);
    if sum(irf) > 0; irf = irf / sum(irf); end
    photons = sum(decay);

    for nExp = [2 3]
        entry = struct('nExp', nExp, 'tauNs', nan(1, nExp), ...
            'amplitude', nan(1, nExp), 'offset', NaN, 'shiftNs', NaN, ...
            'model', [], 'residual', [], 'chiSquare', NaN, ...
            'photons', photons, 'ok', false, 'message', '');
        if photons < 50
            entry.message = sprintf('only %.0f photons', photons);
            fits(end + 1) = entry; %#ok<AGROW>
            continue;
        end
        seed = logspace(log10(0.35), log10(3.5), nExp);
        try
            [shift, offset, amplitude, tau, ~, ~, ~, basis] = ...
                Fluofit(irf, decay, periodNs, dtNs, seed, [], 0, 'mle', false);
            tau = double(tau(:)');
            amplitude = double(amplitude(:)');

            % sum(basis, 2) is the fitted curve: Fluofit scales each basis
            % column - including the leading constant - by its amplitude.
            model = sum(double(basis), 2);
            if numel(model) ~= numel(decay)
                entry.message = sprintf(['Fluofit returned a %d-bin fit for ' ...
                    '%d bins of data'], numel(model), numel(decay));
                fits(end + 1) = entry; %#ok<AGROW>
                continue;
            end

            % The amplitude vector carries a leading constant term; the
            % remaining entries pair with the lifetimes.
            if numel(amplitude) > numel(tau)
                constant = amplitude(1);
                amplitude = amplitude(end - numel(tau) + 1:end);
            else
                constant = double(offset(1));
            end
            [tau, order] = sort(tau, 'ascend');
            amplitude = amplitude(order);

            positive = model > 0;
            dof = max(1, nnz(positive) - (2 * nExp + 2));
            entry.tauNs = tau;
            entry.amplitude = amplitude / max(sum(amplitude), eps);
            entry.offset = constant;
            entry.shiftNs = double(shift(1));
            entry.model = model;
            entry.residual = (decay - model) ./ sqrt(max(model, 1));
            entry.chiSquare = sum((decay(positive) - model(positive)) .^ 2 ./ ...
                model(positive)) / dof;
            entry.ok = all(isfinite(tau)) && all(tau > 0) && ...
                all(isfinite(model));
            if ~entry.ok
                entry.message = 'fit returned a non-physical lifetime';
            end
        catch fitError
            entry.message = fitError.message;
        end
        fits(end + 1) = entry; %#ok<AGROW>
    end
end
