function [DU, DV] = detectorPixelOffsets(sim, nSub)
%DETECTORPIXELOFFSETS Quadrature offsets inside one detector pixel.
%
%   The default for the Luminosa honeycomb model is a point-up hexagonal
%   detector cell. For legacy square-detector calculations, set
%   sim.detectorPixelShape = 'square'.

    if nargin < 2 || isempty(nSub)
        nSub = 5;
    end

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
            du = linspace(-xHalf, xHalf, nSub);
            dv = linspace(-yHalf, yHalf, nSub);
            [U, V] = meshgrid(du, dv);
            keep = detectorHexMask(U, V, radius);

            DU = U(keep);
            DV = V(keep);

        case {'square','rect','rectangle'}
            du = linspace(-sim.detSize/2, sim.detSize/2, nSub);
            dv = linspace(-sim.detSize/2, sim.detSize/2, nSub);
            [U, V] = meshgrid(du, dv);
            DU = U(:);
            DV = V(:);

        otherwise
            error('detectorPixelOffsets:UnknownPixelShape', ...
                'Unknown detectorPixelShape "%s". Use ''hex'' or ''square''.', pixelShape);
    end
end
