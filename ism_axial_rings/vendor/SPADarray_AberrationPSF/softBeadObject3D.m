function obj = softBeadObject3D(sim, zCenter)
%SOFTBEADOBJECT3D Fractional-voxel fluorescent sphere on the simulation grid.
%
%   obj = softBeadObject3D(sim, zCenter)
%
%   sim.beadSubsamples controls subvoxel integration as either one integer
%   or [nx ny nz]. The default is [3 3 3].

    if nargin < 2 || isempty(zCenter)
        zCenter = 0;
    end
    if ~isfield(sim,'beadRadius') || ~isscalar(sim.beadRadius) || ...
            ~isfinite(sim.beadRadius) || sim.beadRadius <= 0
        error('softBeadObject3D:InvalidRadius', ...
            'sim.beadRadius must be a finite positive scalar.');
    end

    nSub = [3 3 3];
    if isfield(sim,'beadSubsamples') && ~isempty(sim.beadSubsamples)
        value = double(sim.beadSubsamples(:)).';
        if isscalar(value)
            nSub = repmat(value,1,3);
        elseif numel(value) == 3
            nSub = value;
        else
            error('softBeadObject3D:InvalidSubsamples', ...
                'sim.beadSubsamples must be a scalar or [nx ny nz].');
        end
    end
    if any(~isfinite(nSub)) || any(nSub < 1) || any(nSub ~= round(nSub))
        error('softBeadObject3D:InvalidSubsamples', ...
            'Bead subvoxel counts must be positive integers.');
    end

    x = double(sim.x(:)).';
    y = double(sim.y(:)).';
    z = double(sim.z(:)).';
    dx = gridSpacing(x,fieldOr(sim,'dx',NaN));
    dy = gridSpacing(y,dx);
    dz = gridSpacing(z,min(dx,dy));
    ox = centeredSubpixelOffsets(dx,nSub(1));
    oy = centeredSubpixelOffsets(dy,nSub(2));
    oz = centeredSubpixelOffsets(dz,nSub(3));

    [X,Y] = meshgrid(x,y);
    obj = zeros(numel(y),numel(x),numel(z));
    radius2 = sim.beadRadius^2;
    nSamples = prod(nSub);

    for iz = 1:numel(z)
        axialCenter = z(iz)-zCenter;
        if abs(axialCenter) > sim.beadRadius+dz/2
            continue;
        end
        occupancy = zeros(size(X));
        for kz = 1:numel(oz)
            axial = axialCenter+oz(kz);
            lateralRadius2 = radius2-axial^2;
            if lateralRadius2 < 0
                continue;
            end
            for ky = 1:numel(oy)
                y2 = (Y+oy(ky)).^2;
                for kx = 1:numel(ox)
                    occupancy = occupancy + ...
                        ((X+ox(kx)).^2+y2 <= lateralRadius2);
                end
            end
        end
        obj(:,:,iz) = occupancy/nSamples;
    end

    total = sum(obj(:));
    if total <= 0
        error('softBeadObject3D:BeadOutsideGrid', ...
            'The bead does not overlap the configured simulation grid.');
    end
    obj = obj/total;
end

function spacing = gridSpacing(axisValues,fallback)
    if numel(axisValues) > 1
        spacing = median(abs(diff(axisValues)));
    else
        spacing = fallback;
    end
    if ~isfinite(spacing) || spacing <= 0
        error('softBeadObject3D:InvalidGrid', ...
            'Simulation axes must have finite positive spacing.');
    end
end

function offsets = centeredSubpixelOffsets(spacing,n)
    offsets = (((1:n)-0.5)/n-0.5)*spacing;
end

function value = fieldOr(S,name,fallback)
    value = fallback;
    if isfield(S,name) && ~isempty(S.(name))
        value = S.(name);
    end
end
