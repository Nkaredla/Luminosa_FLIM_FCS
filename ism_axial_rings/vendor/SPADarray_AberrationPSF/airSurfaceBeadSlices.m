function [objectSlices,heightsUm,weights] = airSurfaceBeadSlices(sim)
%AIRSURFACEBEADSLICES Quadrature model of a sphere resting on glass in air.

    nAxial = 3;
    if isfield(sim,'airBeadAxialSamples') && ~isempty(sim.airBeadAxialSamples)
        nAxial = double(sim.airBeadAxialSamples);
    end
    if ~isscalar(nAxial) || ~isfinite(nAxial) || ...
            nAxial < 1 || nAxial ~= round(nAxial)
        error('airSurfaceBeadSlices:BadAxialSamples', ...
            'sim.airBeadAxialSamples must be a positive integer.');
    end
    bottom = 0;
    if isfield(sim,'beadBottomHeightUm') && ~isempty(sim.beadBottomHeightUm)
        bottom = double(sim.beadBottomHeightUm);
    end
    if ~isscalar(bottom) || ~isfinite(bottom) || bottom < 0
        error('airSurfaceBeadSlices:BadBottomHeight', ...
            'sim.beadBottomHeightUm must be finite and nonnegative.');
    end

    [nodes,quadratureWeights] = gaussLegendre(nAxial);
    radius = sim.beadRadius;
    heightsUm = bottom+radius*(nodes(:).'+1);
    axialWeights = radius*quadratureWeights(:).';
    centerHeight = bottom+radius;

    lateralSubsamples = 3;
    if isfield(sim,'beadSubsamples') && ~isempty(sim.beadSubsamples)
        value = double(sim.beadSubsamples(:)).';
        if isscalar(value)
            lateralSubsamples = value;
        else
            lateralSubsamples = max(value(1:2));
        end
    end
    lateralSubsamples = max(1,round(lateralSubsamples));
    dx = median(abs(diff(sim.x)));
    dy = median(abs(diff(sim.y)));
    ox = (((1:lateralSubsamples)-0.5)/lateralSubsamples-0.5)*dx;
    oy = (((1:lateralSubsamples)-0.5)/lateralSubsamples-0.5)*dy;
    [X,Y] = meshgrid(sim.x,sim.y);

    objectSlices = zeros(numel(sim.y),numel(sim.x),nAxial);
    weights = zeros(1,nAxial);
    for iz = 1:nAxial
        crossRadius2 = max(radius^2-(heightsUm(iz)-centerHeight)^2,0);
        occupancy = zeros(size(X));
        for ky = 1:numel(oy)
            for kx = 1:numel(ox)
                occupancy = occupancy+ ...
                    ((X+ox(kx)).^2+(Y+oy(ky)).^2 <= crossRadius2);
            end
        end
        occupancy = occupancy/(numel(ox)*numel(oy));
        objectSlices(:,:,iz) = occupancy*axialWeights(iz);
        weights(iz) = sum(objectSlices(:,:,iz),'all');
    end
    total = sum(objectSlices,'all');
    if total <= 0
        error('airSurfaceBeadSlices:EmptyBead', ...
            'The surface bead does not overlap the lateral simulation grid.');
    end
    objectSlices = objectSlices/total;
    weights = weights/sum(weights);
end

function [nodes,weights] = gaussLegendre(n)
    if n == 1
        nodes = 0;
        weights = 2;
        return;
    end
    beta = (1:n-1)./sqrt(4*(1:n-1).^2-1);
    J = diag(beta,1)+diag(beta,-1);
    [V,D] = eig(J);
    [nodes,order] = sort(diag(D));
    V = V(:,order);
    weights = 2*(V(1,:).^2).';
end
