function map = perceptualColormap(name, n)
%PERCEPTUALCOLORMAP Perceptually uniform colormaps for quantitative maps.
%
% map = perceptualColormap(name)        % 256 levels
% map = perceptualColormap(name, n)
%
% name: 'viridis' | 'magma' | 'cividis' | 'coolwarm' | 'gray'
%
% WHY NOT JET
%
% Jet is not perceptually uniform: its lightness rises, falls and rises again,
% so equal steps in value produce very unequal steps in apparent brightness. It
% invents sharp yellow/cyan bands where the data is smooth and flattens contrast
% in the green, which is exactly the failure mode for lifetime maps - a boundary
% that looks crisp in jet is often just the colormap's own cyan-to-yellow
% transition. It also collapses to near-uniform grey when printed and is not
% usable by red-green colourblind readers, which rules it out for publication.
%
% viridis and magma are monotonic in lightness, so a printed greyscale version
% still orders the data correctly; cividis is additionally optimised for
% colour-vision deficiency.
%
% coolwarm is diverging and is for SIGNED quantities only - the prior pull -
% where white marks zero and the eye should find the sign change. Using a
% diverging map on unsigned data would put a meaningless feature at its midpoint.
%
% Anchor values are the standard control points of each map, interpolated to n
% levels.

    if nargin < 2 || isempty(n); n = 256; end
    switch lower(char(name))
        case 'viridis'
            anchors = [
                0.267004 0.004874 0.329415
                0.282623 0.140926 0.457517
                0.253935 0.265254 0.529983
                0.206756 0.371758 0.553117
                0.163625 0.471133 0.558148
                0.127568 0.566949 0.550556
                0.134692 0.658636 0.517649
                0.266941 0.748751 0.440573
                0.477504 0.821444 0.318195
                0.741388 0.873449 0.149561
                0.993248 0.906157 0.143936];
        case 'magma'
            anchors = [
                0.001462 0.000466 0.013866
                0.078815 0.054184 0.211667
                0.232077 0.059889 0.437695
                0.390384 0.100379 0.501864
                0.550287 0.161158 0.505719
                0.716387 0.214982 0.475290
                0.868793 0.287728 0.409303
                0.967671 0.439703 0.359665
                0.994738 0.624350 0.427397
                0.996580 0.793550 0.545341
                0.987053 0.991438 0.749504];
        case 'cividis'
            anchors = [
                0.000000 0.135112 0.304751
                0.000000 0.212030 0.401800
                0.152000 0.285700 0.402400
                0.271000 0.356100 0.408800
                0.372000 0.427600 0.424700
                0.470000 0.500600 0.446300
                0.573000 0.576000 0.457000
                0.684000 0.655000 0.437000
                0.798000 0.738000 0.393000
                0.914000 0.825000 0.325000
                0.995737 0.909344 0.217772];
        case 'coolwarm'
            anchors = [
                0.2298 0.2987 0.7537
                0.4064 0.5375 0.9345
                0.6332 0.7304 0.9946
                0.8618 0.8618 0.8618
                0.9683 0.7176 0.6120
                0.8830 0.4730 0.3703
                0.7057 0.0156 0.1502];
        case 'gray'
            anchors = [0 0 0; 1 1 1];
        otherwise
            error('perceptualColormap:Unknown', ...
                ['Unknown colormap ''%s''. Use viridis, magma, cividis, ' ...
                 'coolwarm or gray.'], name);
    end
    x = linspace(0, 1, size(anchors, 1));
    xi = linspace(0, 1, n);
    map = min(max(interp1(x, anchors, xi, 'pchip'), 0), 1);
end
