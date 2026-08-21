function vol = psfCartesianToPolar(h, x, y, rhoVec, Nphi, M)
%--------------------------------------------------------------------------
% psfCartesianToPolar
%
% PURPOSE
%   Convert a 3-D PSF given on a Cartesian grid h(y,x,z) into the
%   azimuthal-harmonic format [Nr × Nz × (2M+1)] expected by FCS.
%
% INPUTS
%   h      : [Ny × Nx × Nz]  PSF (or any 3-D function) on Cartesian grid
%   x      : [1 × Nx] or [Nx × 1]  x coordinates [µm]
%   y      : [1 × Ny] or [Ny × 1]  y coordinates [µm]
%   rhoVec : [Nr × 1]  radial sample points [µm], 0 ≤ rho ≤ rho_max
%   Nphi   : number of azimuthal sample points for the FFT  (default 64)
%   M      : highest harmonic order to retain  (default 5)
%
% OUTPUT
%   vol : [Nr × Nz × (2M+1)]  azimuthal Fourier coefficients
%
%   Harmonic encoding (matching FCS.m convention):
%     vol(:,:,1)       = h_0(rho,z)              DC (m=0)
%     vol(:,:,1+m)     = h_c_m(rho,z)            cos(m·phi) amplitude
%     vol(:,:,M+1+m)   = h_s_m(rho,z)            sin(m·phi) amplitude
%
%   so that the full reconstruction is:
%     h(rho,phi,z) = h_0 + sum_{m=1}^{M} [ h_c_m·cos(m·phi)
%                                          + h_s_m·sin(m·phi) ]
%
% MATHEMATICAL DERIVATION
%   1. Interpolate h onto a polar grid (rho_i, phi_j, z_k):
%        h_polar(i,j,k) = h( rho_i·cos(phi_j), rho_i·sin(phi_j), z_k )
%
%   2. Discrete azimuthal FFT:
%        A_m(rho_i, z_k) = (1/Nphi) sum_j h_polar(i,j,k) exp(-i·m·phi_j)
%
%   3. Extract real harmonics:
%        h_0    = real( A_0 )
%        h_c_m  = 2·real( A_m )    m = 1..M
%        h_s_m  = -2·imag( A_m )   m = 1..M
%
%   This is the standard cosine-sine decomposition: the factor of 2 in
%   h_c_m and h_s_m accounts for the negative-frequency mirror images in
%   the one-sided representation.
%
% NOTE ON rho = 0
%   At rho = 0 all azimuthal components collapse to the same value; the
%   interpolation naturally returns h(0,0,z) for all phi, so only the m=0
%   component is non-zero — which is physically correct.
%--------------------------------------------------------------------------

    if nargin < 5 || isempty(Nphi), Nphi = 64; end
    if nargin < 6 || isempty(M),    M    = 5;  end

    Nr = numel(rhoVec);
    Nz = size(h, 3);
    rhoVec = rhoVec(:);                       % ensure column [Nr×1]
    x = x(:)';                                % ensure row [1×Nx]
    y = y(:)';                                % ensure row [1×Ny]

    phiVec = (0:Nphi-1) / Nphi * 2*pi;       % [1×Nphi], does not repeat 2π

    vol = zeros(Nr, Nz, 2*M+1);

    for iz = 1:Nz
        slice = h(:,:,iz);                    % [Ny × Nx]

        % ── Sample PSF on polar grid ──────────────────────────────────────
        % xp, yp: [Nr × Nphi] query point grids
        xp = rhoVec * cos(phiVec);            % [Nr × Nphi]
        yp = rhoVec * sin(phiVec);            % [Nr × Nphi]

        % interp2 with 'linear' extrapolates to 0 outside grid
        h_polar = interp2(x, y, slice, xp, yp, 'linear', 0);  % [Nr × Nphi]

        % ── Azimuthal FFT (over columns = phi direction) ──────────────────
        A = fft(h_polar, [], 2) / Nphi;       % [Nr × Nphi], complex

        % ── Harmonic extraction ───────────────────────────────────────────
        vol(:, iz, 1) = real(A(:, 1));                  % m = 0

        for m = 1:M
            vol(:, iz, 1+m)   =  2 * real(A(:, m+1));  % cos(m·phi)
            vol(:, iz, M+1+m) = -2 * imag(A(:, m+1));  % sin(m·phi)
        end
    end
end
