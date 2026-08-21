function mesh = makePeriodicVoronoiMesh(opts)
%MAKEPERIODICVORONOIMESH Create seed points for a periodic Voronoi mesh.

    if ~ismember(opts.diffusionMode, {'hop', 'mixed'})
        mesh = struct('isActive', false, 'seedPositionsUm', zeros(0, 2), ...
            'nCompartments', 0, 'nominalCompartmentSizeUm', NaN, ...
            'meanAreaUm2', NaN);
        return;
    end

    areaUm2 = prod(opts.boxSizeUm);
    nSeeds = max(4, round(areaUm2 / opts.compartmentSizeUm^2));
    seeds = (rand(nSeeds, 2) - 0.5) .* opts.boxSizeUm;

    mesh = struct();
    mesh.isActive = true;
    mesh.seedPositionsUm = seeds;
    mesh.nCompartments = nSeeds;
    mesh.nominalCompartmentSizeUm = opts.compartmentSizeUm;
    mesh.meanAreaUm2 = areaUm2 / nSeeds;
    mesh.meanSqrtAreaUm = sqrt(mesh.meanAreaUm2);
end
