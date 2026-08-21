function stack = normalizedStackExplicitDetector(sim, coeffs, x0, y0)
%--------------------------------------------------------------------------
% normalizedStackExplicitDetector
%
% PURPOSE
%   Generate the full detector-resolved ISM stack using explicit
%   detector-plane pixel integration rather than the shifted-PSF shortcut.
%
% MODEL
%   For each detector pixel k:
%
%     eta_k(x,y,z) = \int_{\Omega_k} h_det(u + x, v + y, z) du dv
%
%   where h_det is the centered detection intensity PSF in sample-equivalent
%   detector coordinates, and the (+x,+y) argument accounts for the normal
%   microscope image inversion between sample and detector plane.
%
%   The effective channel PSF is then
%
%     h_eff,k(x,y,z) = h_exc(x,y,z) .* eta_k(x,y,z)
%
%   and the scanned image is the correlation of the object with h_eff,k:
%
%     I_k(x_s) = integral object(r) h_eff,k(r-x_s) dr
%
%   The rot90 below is therefore intentional. MATLAB conv2 performs a
%   convolution and flips its kernel internally; pre-rotating h_eff converts
%   that operation into the scan-coordinate cross-correlation above.
%
% NOTES
%   - This removes the per-channel renormalization that was distorting the
%     relative brightness of the 25 channels.
%   - It is still a scalar model, but the detector integration is now done
%     explicitly in the detector plane.
%--------------------------------------------------------------------------

    if nargin < 3 || isempty(x0), x0 = 0; end
    if nargin < 4 || isempty(y0), y0 = 0; end

    coeffs = coeffStruct(sim, coeffs);

    % Excitation and centered detection PSFs
    hExc = psfBessel(sim, coeffs, sim.lamExc);
    hDet = psfBessel(sim, coeffs, sim.lamEm);

    % Explicit detector collection-efficiency stacks
    eta = detectorCollectionEfficiencyExplicit(sim, hDet);

    nCh = size(sim.detXY,1);
    chans = zeros(sim.nx, sim.nx, nCh);

    for k = 1:nCh
        hEff = hExc .* eta(:,:,:,k);   % no per-channel renormalization

        img = zeros(sim.nx, sim.nx);
        for iz = 1:size(sim.obj,3)
            % Scanned image formation is a correlation with hEff
            img = img + conv2(sim.obj(:,:,iz), rot90(hEff(:,:,iz), 2), 'same');
        end

        chans(:,:,k) = img;
    end

    % Optional overall sample shift
    if abs(x0) > 0 || abs(y0) > 0
        shiftedStack = zeros(size(chans));
        for k = 1:size(chans,3)
            shiftedStack(:,:,k) = max(fourierShift2D(chans(:,:,k), x0/sim.dx, y0/sim.dx), 0);
        end
        chans = shiftedStack;
    end

    % Only one normalization at the very end
    s = sum(chans(:));
    if s > 0
        stack = chans / s;
    else
        stack = chans;
    end
end
