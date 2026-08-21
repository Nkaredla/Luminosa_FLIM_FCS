function psf = scalarPSFBessel(sim, coeffs, wavelength)
%--------------------------------------------------------------------------
% scalarPSFBessel
%
% PURPOSE
%   Compute a scalar 3D PSF using a Bessel / azimuthal harmonic expansion.
%
% INPUTS
%   sim        : simulation parameters
%   coeffs     : aberration coefficient structure
%   wavelength : wavelength at which the PSF is computed
%
% OUTPUT
%   psf(y, x, z) : normalized 3D PSF
%
% DESCRIPTION
%   The pupil is first decomposed into azimuthal harmonics m. Then the
%   image-plane field is evaluated using Bessel functions J_m and an axial
%   phase factor. The intensity PSF is |U|^2.
%
% NOTES
%   This is a scalar model, not a full vectorial dipole-emission model.
%--------------------------------------------------------------------------

    % Get radial pupil samples, harmonic indices, and pupil harmonics
    [rho, ms, pm] = pupilFourierCoeffs(sim, coeffs, wavelength);

    % Build lateral image-plane coordinate grid
    [X, Y] = meshgrid(sim.x, sim.y);
    r = sqrt(X.^2 + Y.^2);     % radial coordinate in image plane
    psi = atan2(Y, X);         % azimuthal angle in image plane

    % Flatten radius array for vectorized evaluation
    rv = r(:);

    % Reduce repeated radius evaluations by working only on unique radii
    ur = unique(round(rv, 12));
    [~, inv] = ismember(round(rv, 12), ur);

    % Wavenumber
    k0 = 2*pi / wavelength;

    % Maximum pupil sine angle
    smax = sim.NA / sim.nMedium;

    % cos(theta) across the pupil
    cosT = sqrt(max(1 - (smax*rho).^2, 1e-12));

    % Debye solid-angle Jacobian. This geometric measure is rho/cosT.
    % pupilFourierCoeffs separately includes the aplanatic sqrt(cosT)
    % amplitude exactly once, giving a combined rho/sqrt(cosT) factor.
    % Replacing this measure by rho/sqrt(cosT) would double-count apodization.
    weight = (smax^2) * rho ./ cosT;

    % Precompute Bessel arguments for all unique image radii and pupil rho
    arg = k0 * sim.NA * (ur * rho);

    % Precompute Bessel tables for each harmonic order m
    Jtab = cell(numel(ms),1);
    for ii = 1:numel(ms)
        m = ms(ii);
        J = besselj(abs(m), arg);

        % For negative m, use parity relation
        if m < 0
            J = ((-1)^(-m)) * J;
        end

        Jtab{ii} = J;
    end

    % Allocate output PSF
    psf = zeros(numel(sim.y), numel(sim.x), numel(sim.z));

    % Flatten image azimuths
    psiV = psi(:);

    % Loop over z planes
    for iz = 1:numel(sim.z)

        % Axial phase factor exp(i*k*n*z*cos(theta))
        axial = exp(1i * k0 * sim.nMedium * sim.z(iz) * cosT);

        % Radial part of the field for each harmonic m
        radial = complex(zeros(numel(ur), numel(ms)));

        % Perform radial integral for each harmonic
        for ii = 1:numel(ms)
            basis = (pm(:,ii).').*axial.*weight;
            radial(:,ii) = trapz(rho, Jtab{ii} .* basis, 2);
        end

        % Reconstruct complex field in the image plane
        U = complex(zeros(size(rv)));
        for ii = 1:numel(ms)
            m = ms(ii);
            U = U + (2*pi) * (1i^m) * exp(1i*m*psiV) .* radial(inv,ii);
        end

        % Convert field to intensity
        psf(:,:,iz) = abs(reshape(U, size(r))).^2;
    end

    % Normalize PSF to unit total power
    psf = psf / sum(psf(:));
end
