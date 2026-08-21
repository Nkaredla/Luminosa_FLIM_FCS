
function out = boxBlurStack(stack, detPix)
%--------------------------------------------------------------------------
% boxBlurStack
%
% PURPOSE
%   Approximate finite detector-pixel integration by convolving each z-slice
%   with a square box kernel.
%
% INPUTS
%   stack  : 3D stack (y,x,z)
%   detPix : detector size in image pixels
%
% OUTPUT
%   out : blurred 3D stack
%
% NOTES
%   The detector pixel width is forced to be odd so the blur kernel remains
%   centered.
%--------------------------------------------------------------------------

    % Enforce odd-sized kernel
    if mod(detPix,2)==0
        detPix = detPix+1;
    end

    % Uniform box kernel
    ker = ones(detPix, detPix);
    ker = ker / sum(ker(:));

    % Allocate output
    out = zeros(size(stack));

    % Convolve each z slice independently
    for iz = 1:size(stack,3)
        out(:,:,iz) = conv2(stack(:,:,iz), ker, 'same');
    end
end
