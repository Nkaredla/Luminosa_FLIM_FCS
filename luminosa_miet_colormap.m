function [cmap, shadowFraction] = luminosa_miet_colormap(name, n, variant)
%LUMINOSA_MIET_COLORMAP Build MIET display and "nice image" colormaps.

    if nargin < 1 || isempty(name)
        name = 'MIET Bronze';
    end
    if nargin < 2 || isempty(n)
        n = 256;
    end
    if nargin < 3 || isempty(variant)
        variant = 'display';
    end

    n = max(2, round(double(n)));
    variant = lower(char(string(variant)));
    nameKey = lower(char(string(name)));
    nameKey = strrep(nameKey, ' ', '');
    nameKey = strrep(nameKey, '-', '');

    if strcmp(nameKey, 'mietbronze') || strcmp(nameKey, 'bronze') || strcmp(nameKey, 'miet')
        baseMap = mietBronzeBase();
        if strcmp(variant, 'nice')
            cmap = resizeColormap(baseMap, n);
            shadowFraction = 27 / size(baseMap, 1);
        else
            cmap = resizeColormap(baseMap(28:end, :), n);
            shadowFraction = 0;
        end
        return;
    end

    switch nameKey
        case 'parula'
            displayMap = parula(max(32, n));
        case 'hot'
            displayMap = hot(max(32, n));
        case 'copper'
            displayMap = copper(max(32, n));
        case 'gray'
            displayMap = gray(max(32, n));
        case 'bone'
            displayMap = bone(max(32, n));
        otherwise
            displayMap = jet(max(64, n + 32));
            displayMap = displayMap(17:end-16, :);
    end

    if strcmp(variant, 'nice')
        shadowLen = max(16, round(0.28 * n));
        displayLen = max(16, n - shadowLen);
        cmap = [buildShadowPrefix(shadowLen); resizeColormap(displayMap, displayLen)];
        shadowFraction = shadowLen / size(cmap, 1);
    else
        cmap = resizeColormap(displayMap, n);
        shadowFraction = 0;
    end
end

function cmap = resizeColormap(cmapIn, n)
    cmapIn = double(cmapIn);
    if size(cmapIn, 1) == n
        cmap = cmapIn;
        return;
    end

    xIn = linspace(0, 1, size(cmapIn, 1));
    xOut = linspace(0, 1, n);
    cmap = zeros(n, size(cmapIn, 2));
    for colIdx = 1:size(cmapIn, 2)
        cmap(:, colIdx) = interp1(xIn, cmapIn(:, colIdx), xOut, 'linear');
    end
    cmap = min(max(cmap, 0), 1);
end

function cmap = buildShadowPrefix(n)
    n = max(2, round(double(n)));
    grayVals = linspace(0.8, 0.3137, n - 1).';
    cmap = [grayVals grayVals grayVals; 0 0 0];
end

function cmap = mietBronzeBase()
    cmap = [ ...
        0.8 0.8 0.8
        0.7805 0.7805 0.7805
        0.7611 0.7611 0.7611
        0.7416 0.7416 0.7416
        0.7222 0.7222 0.7222
        0.7027 0.7027 0.7027
        0.6833 0.6833 0.6833
        0.6638 0.6638 0.6638
        0.6444 0.6444 0.6444
        0.6249 0.6249 0.6249
        0.6055 0.6055 0.6055
        0.5860 0.5860 0.5860
        0.5666 0.5666 0.5666
        0.5471 0.5471 0.5471
        0.5277 0.5277 0.5277
        0.5082 0.5082 0.5082
        0.4888 0.4888 0.4888
        0.4693 0.4693 0.4693
        0.4499 0.4499 0.4499
        0.4304 0.4304 0.4304
        0.4110 0.4110 0.4110
        0.3915 0.3915 0.3915
        0.3721 0.3721 0.3721
        0.3526 0.3526 0.3526
        0.3332 0.3332 0.3332
        0.3137 0.3137 0.3137
        0.0000 0.0000 0.0000
        0.07922 0.02275 0.01647
        0.1584 0.04549 0.03294
        0.2376 0.06824 0.04941
        0.3169 0.09098 0.06588
        0.3961 0.1137 0.08235
        0.4484 0.1366 0.1039
        0.5007 0.1595 0.1255
        0.5529 0.1824 0.1471
        0.6052 0.2052 0.1686
        0.6575 0.2281 0.1902
        0.7098 0.2510 0.2118
        0.7284 0.3103 0.2260
        0.7471 0.3696 0.2402
        0.7657 0.4289 0.2544
        0.7843 0.4882 0.2686
        0.8029 0.5475 0.2828
        0.8216 0.6069 0.2971
        0.8402 0.6662 0.3113
        0.8588 0.7255 0.3255
        0.8729 0.7518 0.3890
        0.8871 0.7780 0.4525
        0.9012 0.8043 0.5161
        0.9153 0.8306 0.5796
        0.9294 0.8569 0.6431
        0.9435 0.8831 0.7067
        0.9576 0.9094 0.7702
        0.9718 0.9357 0.8337
        0.9859 0.9620 0.8973
        1.0000 0.9882 0.9608
        1.0000 0.9897 0.9657
        1.0000 0.9912 0.9706
        1.0000 0.9926 0.9755
        1.0000 0.9941 0.9804
        1.0000 0.9956 0.9853
        1.0000 0.9971 0.9902
        1.0000 0.9985 0.9951
        1.0000 1.0000 1.0000];
end
