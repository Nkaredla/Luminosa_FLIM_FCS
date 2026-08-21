function stack = normalizedStackExplicitDetectorZPlanes(sim, coeffs, planeZ, x0, y0, z0)
%--------------------------------------------------------------------------
% normalizedStackExplicitDetectorZPlanes
%
% PURPOSE
%   Generate a detector-resolved bead scan for multiple recorded axial
%   planes using explicit detector-plane SPAD integration.
%
% INPUTS
%   sim    : simulation structure from defaultParams()
%   coeffs : aberration coefficient struct or numeric vector
%   planeZ : recorded bead-centre z positions [um], e.g. [-0.1 0 0.1]
%   x0,y0  : global lateral bead shift [um]
%   z0     : global axial bead-offset correction [um]
%
% OUTPUT
%   stack(y, x, detectorChannel, zPlane) normalized to sum 1 over all planes
%
% MODEL
%   For each detector channel k and recorded plane p:
%
%       h_eff,k = h_exc .* eta_k
%       image_k,p(x_s) = integral object(r) h_eff,k(r-x_s) dr
%
%   This is a cross-correlation because x_s is the scanning-focus position,
%   while h_eff is parameterized by emitter displacement from that focus.
%   MATLAB conv2 flips its second input, so conv2(obj,rot90(h_eff,2)) is the
%   required scan correlation. Removing rot90 would model a widefield camera
%   convolution and reverse odd-aberration signs in scan coordinates.
%   eta_k is computed by detectorCollectionEfficiencyExplicit().
%--------------------------------------------------------------------------

    if nargin < 3 || isempty(planeZ), planeZ = [-0.1 0 0.1]; end
    if nargin < 4 || isempty(x0),     x0 = 0;              end
    if nargin < 5 || isempty(y0),     y0 = 0;              end
    if nargin < 6 || isempty(z0),     z0 = 0;              end

    planeZ = planeZ(:).';
    coeffs = coeffStruct(sim, coeffs);

    hExc = psfBessel(sim, coeffs, sim.lamExc);
    hDet = psfBessel(sim, coeffs, sim.lamEm);
    eta = detectorCollectionEfficiencyExplicit(sim, hDet);

    ny = numel(sim.y);
    nx = numel(sim.x);
    nCh = size(sim.detXY, 1);
    nPlane = numel(planeZ);
    stack = zeros(ny, nx, nCh, nPlane);
    hEffAll = zeros(size(hExc,1), size(hExc,2), size(hExc,3), nCh);

    for k = 1:nCh
        hEffAll(:,:,:,k) = hExc .* eta(:,:,:,k);
    end

    for ip = 1:nPlane
        obj = beadObject3DAtZ(sim, planeZ(ip) + z0);

        for k = 1:nCh
            img = zeros(ny, nx);
            for iz = 1:size(obj, 3)
                img = img + conv2(obj(:,:,iz), rot90(hEffAll(:,:,iz,k), 2), 'same');
            end

            stack(:,:,k,ip) = img;
        end
    end

    if abs(x0) > 0 || abs(y0) > 0
        shiftedStack = zeros(size(stack));
        for ip = 1:nPlane
            for k = 1:nCh
                shiftedStack(:,:,k,ip) = max( ...
                    fourierShift2D(stack(:,:,k,ip), x0/sim.dx, y0/sim.dx), 0);
            end
        end
        stack = shiftedStack;
    end

    s = sum(stack(:));
    if s > 0
        stack = stack / s;
    end
end
