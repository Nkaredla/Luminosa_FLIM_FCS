function image = photonImage(pixelIndices, nx, ny)
%PHOTONIMAGE Accumulate a photon list into an [x,y] intensity image.

    counts = accumarray(round(pixelIndices(:)), 1, [nx * ny 1], @sum, 0);
    image = reshape(counts, [nx ny]);
end
