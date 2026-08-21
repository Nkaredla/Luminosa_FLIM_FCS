function out = detectorBlurStack(stack, sim)
%DETECTORBLURSTACK Approximate finite detector-pixel integration by blur.
%
%   Uses sim.detectorPixelShape when available. The explicit detector model
%   is preferred for quantitative work; this helper keeps the older
%   normalizedStack approximation consistent with the configured aperture.

    pixelShape = 'square';
    if isfield(sim, 'detectorPixelShape') && ~isempty(sim.detectorPixelShape)
        pixelShape = lower(char(sim.detectorPixelShape));
    end

    switch pixelShape
        case {'hex','hexagon','hexagonal','honeycomb'}
            if isfield(sim, 'detectorHexRadius') && ~isempty(sim.detectorHexRadius)
                radius = sim.detectorHexRadius;
            else
                radius = sim.detSize / sqrt(3);
            end

            xHalf = sqrt(3)/2 * radius;
            yHalf = radius;
            nx = max(1, ceil(xHalf / sim.dx));
            ny = max(1, ceil(yHalf / sim.dx));
            xv = (-nx:nx) * sim.dx;
            yv = (-ny:ny) * sim.dx;
            [X, Y] = meshgrid(xv, yv);
            ker = double(detectorHexMask(X, Y, radius));

        case {'square','rect','rectangle'}
            detPix = max(1, round(sim.detSize / sim.dx));
            if mod(detPix,2)==0
                detPix = detPix + 1;
            end
            ker = ones(detPix, detPix);

        otherwise
            error('detectorBlurStack:UnknownPixelShape', ...
                'Unknown detectorPixelShape "%s". Use ''hex'' or ''square''.', pixelShape);
    end

    if sum(ker(:)) <= 0
        ker = 1;
    else
        ker = ker / sum(ker(:));
    end

    out = zeros(size(stack));
    for iz = 1:size(stack,3)
        out(:,:,iz) = conv2(stack(:,:,iz), ker, 'same');
    end
end
