function [profiles, distPix, pos] = lineProfileAcross4(imStack, refIdx, nSamples)
% lineProfileAcross4
%
% Draw one line on a reference image and extract intensity profiles
% along the SAME line from all 4 images.
%
% INPUTS
% ------
% imStack   : H x W x 4 array
% refIdx    : reference image index for drawing the line (1..4), default = 1
% nSamples  : number of sampled points along the line, default = 500
%
% OUTPUTS
% -------
% profiles  : nSamples x 4 matrix, one column per image
% distPix   : nSamples x 1 distance axis in pixels
% pos       : 2 x 2 line endpoints [x1 y1; x2 y2]
%
% EXAMPLE
% -------
% ims = cat(3, img1, img2, img3, img4);
% [profiles, distPix, pos] = lineProfileAcross4(ims, 1, 600);

    if nargin < 2 || isempty(refIdx)
        refIdx = 1;
    end
    if nargin < 3 || isempty(nSamples)
        nSamples = 500;
    end

    if ndims(imStack) ~= 3 || size(imStack,3) ~= 4
        error('imStack must be an H x W x 4 array.');
    end

    % Ensure all images are double for profiling/plotting
    imStack = double(imStack);

    % --- Show reference image and let user draw a line
    f1 = figure('Color','w','Name','Draw line on reference image');
    imagesc(imStack(:,:,refIdx));
    axis image;
    colormap gray;
    colorbar;
    title(sprintf('Draw line on image %d, then double-click the line', refIdx));

    h = drawline('Color','r','LineWidth',1.5);
    wait(h);  % user double-clicks when done

    pos = h.Position;   % [x1 y1; x2 y2]
    x = pos(:,1);
    y = pos(:,2);

    % Distance axis in pixels
    lineLength = hypot(diff(x), diff(y));
    distPix = linspace(0, lineLength, nSamples).';

    % --- Extract profiles from all 4 images
    profiles = zeros(nSamples, 4);
    for k = 1:4
        profiles(:,k) = improfile(imStack(:,:,k), x, y, nSamples, 'bilinear');
    end

    % --- Show the same line over all 4 images
    f2 = figure('Color','w','Name','Line across all 4 images');
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    for k = 1:4
        nexttile;
        imagesc(imStack(:,:,k));
        axis image;
        colormap gray;
        hold on;
        plot(x, y, 'r-', 'LineWidth', 1.5);
        plot(x, y, 'ro', 'MarkerSize', 5, 'LineWidth', 1.2);
        hold off;
        title(sprintf('Image %d', k));
    end

    % --- Plot the intensity profiles
    f3 = figure('Color','w','Name','Intensity profiles');
    plot(distPix, profiles(:,1), 'r', 'LineWidth', 1.5); hold on;
    plot(distPix, profiles(:,2), 'g', 'LineWidth', 1.5);
    plot(distPix, profiles(:,3), 'b', 'LineWidth', 1.5);
    plot(distPix, profiles(:,4), 'm', 'LineWidth', 1.5);
    hold off;
    xlabel('Distance along line (pixels)');
    ylabel('Intensity');
    legend({'Image 1','Image 2','Image 3','Image 4'}, 'Location','best');
    title('Intensity profiles along the same line');
    grid on;

end