function out = fourierShift2D(img, shiftXPx, shiftYPx)
%--------------------------------------------------------------------------
% fourierShift2D
%
% PURPOSE
%   Apply a subpixel shift to a 2D image using the Fourier shift theorem.
%
% INPUTS
%   img       : input 2D image
%   shiftXPx  : x shift in pixels
%   shiftYPx  : y shift in pixels
%
% OUTPUT
%   out : shifted image
%
% NOTES
%   This method preserves interpolation quality better than naive spatial
%   interpolation for many optics problems.
%--------------------------------------------------------------------------

    [ny, nx] = size(img);

    % Frequency coordinates
    ux = ifftshift((-floor(nx/2):ceil(nx/2)-1)/nx);
    uy = ifftshift((-floor(ny/2):ceil(ny/2)-1)/ny);
    [Ux, Uy] = meshgrid(ux, uy);

    % Fourier-domain phase ramp corresponding to the desired shift
    phase = exp(-2i*pi*(shiftXPx*Ux + shiftYPx*Uy));

    % Shifted image
    out = real(ifft2(fft2(img).*phase));
end