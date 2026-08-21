function [rho, ms, pm] = pupilFourierCoeffs(sim, coeffs, wavelength)
%--------------------------------------------------------------------------
% pupilFourierCoeffs
%
% PURPOSE
%   Compute the azimuthal Fourier coefficients of the aberrated pupil.
%
% INPUTS
%   sim        : simulation parameters
%   coeffs     : structure containing aberration amplitudes
%   wavelength : wavelength at which the pupil phase is evaluated
%
% OUTPUTS
%   rho : radial pupil samples
%   ms  : azimuthal harmonic indices kept, from -M to +M
%   pm  : Fourier coefficients pm(rho,m)
%
% DESCRIPTION
%   The pupil is:
%       pupil(rho,phi) = sqrt(cos(theta)) * exp(i * phase(rho,phi))
%
%   where the phase is constructed from Zernike modes. The azimuthal FFT
%   over phi is then used to extract coefficients for harmonics m.
%--------------------------------------------------------------------------

    validateScalarOptics(sim);
    if ~isnumeric(wavelength) || ~isscalar(wavelength) || ...
            ~isfinite(wavelength) || wavelength <= 0
        error('pupilFourierCoeffs:InvalidWavelength', ...
            'wavelength must be a finite positive scalar in micrometers.');
    end

    % Radial and angular pupil sampling
    rho = linspace(0, 1, sim.Nr);
    phi = linspace(0, 2*pi, sim.Nphi+1);
    phi(end) = [];  % remove duplicate 2*pi sample

    % Form 2D pupil coordinate grids
    [R, P] = ndgrid(rho, phi);

    % Initialize wavefront phase map
    W = zeros(size(R));

    % Evaluate all basis modes
    Z = zernikeModes(R, P);

    % Add weighted aberration modes into the total OPD W
    fn = fieldnames(coeffs);
    for k = 1:numel(fn)
        W = W + coeffs.(fn{k}) * Z.(fn{k});
    end

    % Sine-theta maximum: normalized radial extent in the medium
    smax = sim.NA / sim.nMedium;

    % Cos(theta) factor across the pupil
    cosT = sqrt(max(1 - (smax*R).^2, 0));

    % Aberrated pupil function
    % lamRef/wavelength rescales the OPD coefficient convention
    pupil = sqrt(cosT) .* exp(1i*2*pi*(sim.lamRef/wavelength).*W);

    % Fourier transform along azimuthal coordinate phi
    A = fft(pupil, [], 2) / sim.Nphi;

    % Retain harmonic orders m = -M ... +M
    ms = -sim.M:sim.M;

    % Allocate Fourier coefficient table
    pm = complex(zeros(sim.Nr, numel(ms)));

    % Extract the requested harmonics
    for ii = 1:numel(ms)
        idx = mod(ms(ii), sim.Nphi) + 1;
        pm(:,ii) = A(:,idx);
    end
end
