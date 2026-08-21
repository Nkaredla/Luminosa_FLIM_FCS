function eta = detectorCollectionEfficiencyExplicit(sim, hDet)
%--------------------------------------------------------------------------
% detectorCollectionEfficiencyExplicit
%
% PURPOSE
%   Compute detector-pixel collection efficiency explicitly by integrating
%   the centered detection PSF over each finite detector pixel.
%
% INPUTS
%   sim  : simulation structure
%   hDet : centered 3D detection intensity PSF, size [ny nx nz]
%
% OUTPUT
%   eta(y,x,z,k) : collection-efficiency stack for detector channel k
%
% REQUIRED / OPTIONAL sim FIELDS
%   sim.detXY                : [nCh x 2] detector centers in sample-
%                              equivalent coordinates
%   sim.detSize              : detector pixel width in same units
%   sim.detectorPixelShape   : 'hex' or 'square', default from defaultParams
%   sim.detectorHexRadius    : center-to-vertex radius for hex pixels
%   sim.detectorSubsamples   : optional, default 5
%   sim.detectorImageInverted: optional, default true
%
% DESCRIPTION
%   For an emitter displaced by (x,y), the detector-plane image is inverted.
%   In sample-equivalent detector coordinates, the detector pixel centered
%   at d_k therefore samples the centered detection PSF at:
%
%       (u, v) = (d_kx + x, d_ky + y)
%
%   This explicit detector-plane integration is more faithful than shifting
%   a sample-space PSF and blurring it there.
%--------------------------------------------------------------------------

    if isfield(sim, 'detectorSubsamples') && ~isempty(sim.detectorSubsamples)
        nSub = sim.detectorSubsamples;
    else
        nSub = 5;
    end

    if isfield(sim, 'detectorImageInverted') && ~isempty(sim.detectorImageInverted)
        detectorImageInverted = logical(sim.detectorImageInverted);
    else
        detectorImageInverted = true;
    end

    [X, Y] = meshgrid(sim.x, sim.y);
    nCh = size(sim.detXY,1);
    nz  = size(hDet,3);

    eta = zeros(numel(sim.y), numel(sim.x), nz, nCh);

    % Uniform quadrature points across the configured detector pixel shape.
    [DU, DV] = detectorPixelOffsets(sim, nSub);
    nQ = numel(DU);

    for k = 1:nCh
        dkx = sim.detXY(k,1);
        dky = sim.detXY(k,2);

        for iz = 1:nz
            acc = zeros(size(X));

            for q = 1:nQ
                if detectorImageInverted
                    % Sample -> detector image inversion
                    Uq = X + dkx + DU(q);
                    Vq = Y + dky + DV(q);
                else
                    Uq = dkx + DU(q) - X;
                    Vq = dky + DV(q) - Y;
                end

                acc = acc + interp2(sim.x, sim.y, hDet(:,:,iz), Uq, Vq, 'linear', 0);
            end

            eta(:,:,iz,k) = acc / nQ;
        end
    end
end
