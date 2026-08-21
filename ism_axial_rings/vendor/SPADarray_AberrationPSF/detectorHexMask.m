function mask = detectorHexMask(x, y, radius)
%DETECTORHEXMASK True inside a point-up regular hexagon centered at zero.
%
%   radius is the center-to-vertex distance. The horizontal flat-to-flat
%   width is sqrt(3)*radius, and neighboring rows touch at 1.5*radius.

    mask = abs(x) <= sqrt(3)/2 * radius & ...
           abs(y) <= radius - abs(x)/sqrt(3);
end
