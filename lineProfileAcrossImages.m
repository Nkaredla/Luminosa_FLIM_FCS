function [profiles, distPix, pos] = lineProfileAcrossImages(imStack, refIdx, nSamples, imgTitles)
% lineProfileAcrossImages
%
% Draw one line on a reference image and extract intensity profiles
% along the SAME line from all images in the stack.
%
% INPUTS
% ------
% imStack    : H x W x N array
% refIdx     : reference image index for drawing the line (1..N), default = 1
% nSamples   : number of sampled points along the line, default = 500
% imgTitles  : cell array or string array of length N for image titles
%
% OUTPUTS
% -------
% profiles   : nSamples x N matrix, one column per image
% distPix    : nSamples x 1 distance axis in pixels
% pos        : 2 x 2 line endpoints [x1 y1; x2 y2]
%
% EXAMPLE
% -------
% titles = {'Raw','APR','APR deconv','ACO'};
% [profiles, distPix, pos] = lineProfileAcrossImages(ims, 2, 600, titles);

    if nargin < 2 || isempty(refIdx)
        refIdx = 1;
    end
    if nargin < 3 || isempty(nSamples)
        nSamples = 500;
    end

    if ndims(imStack) ~= 3
        error('imStack must be an H x W x N array.');
    end

    [~, ~, N] = size(imStack);

    if refIdx < 1 || refIdx > N
        error('refIdx must be between 1 and N.');
    end

    if nargin < 4 || isempty(imgTitles)
        imgTitles = arrayfun(@(k) sprintf('Image %d', k), 1:N, 'UniformOutput', false);
    else
        if isstring(imgTitles)
            imgTitles = cellstr(imgTitles);
        end
        if ~iscell(imgTitles) || numel(imgTitles) ~= N
            error('imgTitles must be a cell array or string array with one title per image.');
        end
    end

    imStack = double(imStack);

    % Draw line on chosen reference image
    figure('Color','w','Name','Draw line on reference image');
    imagesc(imStack(:,:,refIdx));
    axis image;
    colormap gray;
    colorbar;
    title(sprintf('Draw line on "%s", then double-click the line', imgTitles{refIdx}), ...
        'Interpreter', 'none');

    h = drawline('Color','r','LineWidth',1.5);
    wait(h);

    pos = h.Position;   % [x1 y1; x2 y2]
    x = pos(:,1);
    y = pos(:,2);

    % Distance axis in pixels
    lineLength = hypot(diff(x), diff(y));
    distPix = linspace(0, lineLength, nSamples).';

    % Extract profiles from all images
    profiles = zeros(nSamples, N);
    for k = 1:N
        profiles(:,k) = improfile(imStack(:,:,k), x, y, nSamples, 'bilinear');
    end

    % Show line on all images
    nCols = ceil(sqrt(N));
    nRows = ceil(N / nCols);

    figure('Color','w','Name','Line across all images');
    tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

    for k = 1:N
        nexttile;
        imagesc(imStack(:,:,k));
        axis image;
        colormap gray;
        hold on;
        plot(x, y, 'r-', 'LineWidth', 1.5);
        plot(x, y, 'ro', 'MarkerSize', 5, 'LineWidth', 1.2);
        hold off;
        title(imgTitles{k}, 'Interpreter', 'none');
    end

    % Plot all profiles
    figure('Color','w','Name','Intensity profiles');
    plot(distPix, profiles, 'LineWidth', 1.5);
    xlabel('Distance along line (pixels)');
    ylabel('Intensity');
    title('Intensity profiles along the same line');
    grid on;
    legend(imgTitles, 'Location', 'bestoutside', 'Interpreter', 'none');
end