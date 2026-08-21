function stack = normalizedStack(sim, coeffs, x0, y0)

    if nargin < 3, x0 = 0; end
    if nargin < 4, y0 = 0; end

    coeffs = coeffStruct(sim, coeffs);

    hExc = psfBessel(sim, coeffs, sim.lamExc);
    hDet = psfBessel(sim, coeffs, sim.lamEm);

    detPSF = detectorBlurStack(hDet, sim);

    chans = zeros(sim.nx, sim.nx, size(sim.detXY,1));

    for k = 1:size(sim.detXY,1)
        sx = sim.detXY(k,1);
        sy = sim.detXY(k,2);

        hEff = zeros(size(hExc));

        for iz = 1:size(hExc,3)
            shifted = max(fourierShift2D(detPSF(:,:,iz), sx/sim.dx, sy/sim.dx), 0);
            hEff(:,:,iz) = hExc(:,:,iz) .* shifted;
        end

        hEff = hEff / sum(hEff(:));

        img = zeros(sim.nx, sim.nx);
        for iz = 1:size(sim.obj,3)
            % Standard convolution (NOT correlation)
            img = img + conv2(sim.obj(:,:,iz), hEff(:,:,iz), 'same');
        end

        chans(:,:,k) = img;
    end

    stack = chans / sum(chans(:));

    if abs(x0) > 0 || abs(y0) > 0
        shiftedStack = zeros(size(stack));
        for k = 1:size(stack,3)
            shiftedStack(:,:,k) = max(fourierShift2D(stack(:,:,k), x0/sim.dx, y0/sim.dx), 0);
        end
        stack = shiftedStack / sum(shiftedStack(:));
    end
end
