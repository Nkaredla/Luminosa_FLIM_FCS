function cube = photonCube(photonBins, pixelIndices, nx, ny, gateLength)
%PHOTONCUBE Build a memory-efficient uint16 B12 TCSPC cube.
% Bins are accumulated one at a time, avoiding a full temporary double cube.

    pixelCount = nx * ny;
    cube2d = zeros(pixelCount, gateLength, 'uint16');
    if isempty(photonBins)
        cube = reshape(cube2d, [nx ny gateLength]);
        return;
    end

    timeBin = round(double(photonBins(:)));
    pixel = round(double(pixelIndices(:)));
    valid = timeBin >= 1 & timeBin <= gateLength & pixel >= 1 & pixel <= pixelCount;
    timeBin = timeBin(valid);
    pixel = pixel(valid);
    [timeBin, order] = sort(timeBin);
    pixel = pixel(order);
    edges = [1; find(diff(timeBin) ~= 0) + 1; numel(timeBin) + 1];
    for group = 1:numel(edges)-1
        range = edges(group):(edges(group + 1) - 1);
        bin = timeBin(range(1));
        counts = accumarray(pixel(range), 1, [pixelCount 1], @sum, 0);
        cube2d(:, bin) = uint16(min(counts, double(intmax('uint16'))));
    end
    cube = reshape(cube2d, [nx ny gateLength]);
end
