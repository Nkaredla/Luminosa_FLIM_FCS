function stack = normalizedStackAirInterfaceZPlanes( ...
        sim,coeffs,planeZ,x0,y0,z0)
%NORMALIZEDSTACKAIRINTERFACEZPLANES ISM bead scan at a glass/air interface.

    if nargin < 3 || isempty(planeZ), planeZ = [-0.5 0 0.5]; end
    if nargin < 4 || isempty(x0), x0 = 0; end
    if nargin < 5 || isempty(y0), y0 = 0; end
    if nargin < 6 || isempty(z0), z0 = 0; end

    stageZ = double(planeZ(:).')+z0;
    [objectSlices,heights] = airSurfaceBeadSlices(sim);
    nPlane = numel(stageZ);
    nHeight = numel(heights);
    stageConditions = repelem(stageZ(:),nHeight);
    heightConditions = repmat(heights(:),nPlane,1);

    simExc = sim;
    simExc.airInterfaceDirection = 'excitation';
    simExc = setVectorialPolarizationForPath(simExc, 'excitation');
    hExc = psfBesselAirInterface( ...
        simExc,coeffs,sim.lamExc,stageConditions,heightConditions);
    simDet = sim;
    simDet.airInterfaceDirection = 'collection';
    simDet = setVectorialPolarizationForPath(simDet, 'collection');
    hDet = psfBesselAirInterface( ...
        simDet,coeffs,sim.lamEm,stageConditions,heightConditions);
    eta = detectorCollectionEfficiencyExplicit(sim,hDet);

    nChannel = size(sim.detXY,1);
    stack = zeros(numel(sim.y),numel(sim.x),nChannel,nPlane);
    for ip = 1:nPlane
        for ih = 1:nHeight
            condition = (ip-1)*nHeight+ih;
            objectSlice = objectSlices(:,:,ih);
            for ch = 1:nChannel
                effectivePSF = hExc(:,:,condition).*eta(:,:,condition,ch);
                stack(:,:,ch,ip) = stack(:,:,ch,ip)+ ...
                    conv2(objectSlice,rot90(effectivePSF,2),'same');
            end
        end
    end

    if x0 ~= 0 || y0 ~= 0
        for ip = 1:nPlane
            for ch = 1:nChannel
                stack(:,:,ch,ip) = max(fourierShift2D( ...
                    stack(:,:,ch,ip),x0/sim.dx,y0/sim.dx),0);
            end
        end
    end

    total = sum(stack(:));
    if total > 0
        stack = stack/total;
    end
end
