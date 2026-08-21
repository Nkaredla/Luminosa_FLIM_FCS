
function obj = beadObject3D(sim)
%--------------------------------------------------------------------------
% beadObject3D
%
% PURPOSE
%   Build a fractional-voxel sphere representing a fluorescent bead.
%
% OUTPUT
%   obj(y, x, z) : normalized 3D fluorescence object
%
% DESCRIPTION
%   The sphere is subvoxel-integrated to avoid plane-dependent changes in
%   apparent bead volume as its centre moves relative to the axial grid.
%--------------------------------------------------------------------------

    obj = softBeadObject3D(sim,0);
end

