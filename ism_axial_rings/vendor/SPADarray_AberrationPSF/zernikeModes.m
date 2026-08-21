function Z = zernikeModes(rho, phi)
%--------------------------------------------------------------------------
% zernikeModes
%
% PURPOSE
%   Return named Zernike-like mode values evaluated on the pupil.
%
% INPUTS
%   rho : radial pupil coordinate, normalized to [0,1]
%   phi : azimuthal pupil angle
%
% OUTPUT
%   Z : structure with named aberration mode arrays
%
% NOTES
%   Coefficients multiplying these modes are RMS waves. The original low
%   order repository basis is preserved and extended through radial order 7
%   for calibration use: trefoil, quadrafoil, secondary astigmatism/coma,
%   pentafoil, secondary spherical, hexafoil, and the radial-order-7 group
%   (tertiary coma/trefoil, secondary pentafoil, heptafoil).
%
%   The astigmatism suffixes are the legacy cosine/sine quadrature names:
%   astig_x = cos(2phi), with axes at 0/90 degrees; astig_y = sin(2phi),
%   the same astigmatism rotated by 45 degrees.
%--------------------------------------------------------------------------

    Z.tilt_x = 2*rho.*cos(phi);
    Z.tilt_y = 2*rho.*sin(phi);

    Z.astig_x = sqrt(6)*rho.^2.*cos(2*phi);
    Z.astig_y = sqrt(6)*rho.^2.*sin(2*phi);

    Z.defocus = sqrt(3)*(2*rho.^2 - 1);

    Z.coma_x = sqrt(8)*(3*rho.^3 - 2*rho).*cos(phi);
    Z.coma_y = sqrt(8)*(3*rho.^3 - 2*rho).*sin(phi);

    Z.trefoil_x = sqrt(8)*rho.^3.*cos(3*phi);
    Z.trefoil_y = sqrt(8)*rho.^3.*sin(3*phi);

    Z.spherical = sqrt(5)*(6*rho.^4 - 6*rho.^2 + 1);

    Z.secondary_astig_x = sqrt(10)*(4*rho.^4 - 3*rho.^2).*cos(2*phi);
    Z.secondary_astig_y = sqrt(10)*(4*rho.^4 - 3*rho.^2).*sin(2*phi);
    Z.quadrafoil_x = sqrt(10)*rho.^4.*cos(4*phi);
    Z.quadrafoil_y = sqrt(10)*rho.^4.*sin(4*phi);

    Z.secondary_coma_x = sqrt(12)*(10*rho.^5 - 12*rho.^3 + 3*rho).*cos(phi);
    Z.secondary_coma_y = sqrt(12)*(10*rho.^5 - 12*rho.^3 + 3*rho).*sin(phi);
    Z.secondary_trefoil_x = sqrt(12)*(5*rho.^5 - 4*rho.^3).*cos(3*phi);
    Z.secondary_trefoil_y = sqrt(12)*(5*rho.^5 - 4*rho.^3).*sin(3*phi);
    Z.pentafoil_x = sqrt(12)*rho.^5.*cos(5*phi);
    Z.pentafoil_y = sqrt(12)*rho.^5.*sin(5*phi);

    Z.secondary_spherical = sqrt(7)*(20*rho.^6 - 30*rho.^4 + 12*rho.^2 - 1);
    Z.tertiary_astig_x = sqrt(14)*(15*rho.^6 - 20*rho.^4 + 6*rho.^2).*cos(2*phi);
    Z.tertiary_astig_y = sqrt(14)*(15*rho.^6 - 20*rho.^4 + 6*rho.^2).*sin(2*phi);
    Z.secondary_quadrafoil_x = sqrt(14)*(6*rho.^6 - 5*rho.^4).*cos(4*phi);
    Z.secondary_quadrafoil_y = sqrt(14)*(6*rho.^6 - 5*rho.^4).*sin(4*phi);
    Z.hexafoil_x = sqrt(14)*rho.^6.*cos(6*phi);
    Z.hexafoil_y = sqrt(14)*rho.^6.*sin(6*phi);

    % Radial order 7 (normalization sqrt(2*(n+1)) = 4 for these m != 0 modes).
    Z.tertiary_coma_x = 4*(35*rho.^7 - 60*rho.^5 + 30*rho.^3 - 4*rho).*cos(phi);
    Z.tertiary_coma_y = 4*(35*rho.^7 - 60*rho.^5 + 30*rho.^3 - 4*rho).*sin(phi);
    Z.tertiary_trefoil_x = 4*(21*rho.^7 - 30*rho.^5 + 10*rho.^3).*cos(3*phi);
    Z.tertiary_trefoil_y = 4*(21*rho.^7 - 30*rho.^5 + 10*rho.^3).*sin(3*phi);
    Z.secondary_pentafoil_x = 4*(7*rho.^7 - 6*rho.^5).*cos(5*phi);
    Z.secondary_pentafoil_y = 4*(7*rho.^7 - 6*rho.^5).*sin(5*phi);
    Z.heptafoil_x = 4*rho.^7.*cos(7*phi);
    Z.heptafoil_y = 4*rho.^7.*sin(7*phi);
end
