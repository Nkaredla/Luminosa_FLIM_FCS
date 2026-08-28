function [P, total] = pie_pattern_columns(basis, tauNs, pulseIndex)
%PIE_PATTERN_COLUMNS Unit-sum patterns for lifetimes attached to a given pulse.
%
% [P, total] = pie_pattern_columns(basis, tauNs, pulseIndex)
%
% BASIS comes from pie_two_pulse_basis. TAUNS is 1-by-K. PULSEINDEX selects
% which excitation pulse these components are excited by (1 = the reference the
% IRF is aligned to). P is nBin-by-K, each column normalised to UNIT SUM so its
% fitted amplitude is that component's PHOTON COUNT. TOTAL is the
% pre-normalisation sum, which converts an amplitude back to a pre-exponential
% weight.
%
% Periodic re-excitation is handled by the circular convolution itself: over one
% full period the wrap IS the previous pulse's remaining decay, so no explicit
% multi-period sum is needed. That is a real simplification over the gated
% patterns, which had to sum several periods precisely because they only saw a
% fraction of one.

    tauNs = double(tauNs(:))';
    if any(tauNs <= 0)
        error('pie_pattern_columns:NonPositiveTau', ...
            'Lifetimes must be positive; min is %.4g.', min(tauNs));
    end
    if nargin < 3 || isempty(pulseIndex); pulseIndex = 1; end
    if pulseIndex < 1 || pulseIndex > numel(basis.shifts)
        error('pie_pattern_columns:BadPulse', ...
            'pulseIndex %d is outside the %d pulse(s) in the basis.', ...
            pulseIndex, numel(basis.shifts));
    end

    t = basis.timeNs;
    E = exp(-t ./ tauNs);
    U = basis.C * E;
    S = sum(U, 1);
    S(S <= 0) = eps;
    P = U ./ S;
    total = S;

    shift = basis.shifts(pulseIndex);
    if shift ~= 0
        P = circshift(P, shift, 1);
    end
end
