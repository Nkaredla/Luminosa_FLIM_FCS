function obj = beadObject3DAtZ(sim, zCenter)
%--------------------------------------------------------------------------
% beadObject3DAtZ
%
% PURPOSE
%   Build a normalized fractional-voxel fluorescent bead whose centre is
%   axially displaced relative to the focal plane.
%
% INPUTS
%   sim     : simulation structure from defaultParams()
%   zCenter : bead centre position [um] in the same coordinate system as
%             sim.z. Use zCenter = 0 for an in-focus bead.
%
% OUTPUT
%   obj(y, x, z) : normalized 3D fluorescence object
%--------------------------------------------------------------------------

    if nargin < 2 || isempty(zCenter)
        zCenter = 0;
    end

    obj = softBeadObject3D(sim,zCenter);
end
