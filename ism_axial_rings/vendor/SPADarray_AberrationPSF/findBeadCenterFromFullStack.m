function out = findBeadCenterFromFullStack(stackInput, varargin)
%FINDBEADCENTERFROMFULLSTACK Find bead center and focal plane from a full ISM stack.
%
%   out = findBeadCenterFromFullStack(stackInput)
%
%   Loads a detector-resolved ISM stack with loadFullStackISMData, sums the
%   detector channels into one 2-D image per z plane, and scores each plane by
%   both brightness and compactness:
%
%       focusScore = normalizedTotalSignal * normalizedCompactness
%
%   Compactness is the inverse of the weighted second spatial moment
%   varX + varY. The chosen focal plane is therefore bright and spatially
%   narrow. The reported bead center is the weighted centroid in that plane.
%
%   Output fields:
%       out.focusIndex
%       out.focusZUm
%       out.centerXY              [x y], subpixel centroid in image pixels
%       out.centerYX              rounded [row col] pixel
%       out.focusImage            detector-summed 2-D image at focus
%       out.metrics               per-plane table
%       out.data                  loader output, including rawCounts/stageZUm
%
%   Example:
%       addpath('SPADarray_AberrationPSF');
%       out = findBeadCenterFromFullStack( ...
%           'xz_yz_plots\x0_4uW_0_19collar_80mmlens_20260515_155744_frame_alignment.csv');
%       fprintf('focus plane %d, center x=%.2f y=%.2f\n', ...
%           out.focusIndex, out.centerXY(1), out.centerXY(2));

    if nargin < 1 || isempty(stackInput)
        stackInput = ['D:\Luminosa\Data\PSF_batch_outputs\ISM_Aberation2_73\' ...
            'xz_yz_plots\x0_4uW_0_19collar_80mmlens_20260515_155744_frame_alignment.csv'];
    end

    opts = parseInputs(varargin{:});

    data = loadFullStackISMData(stackInput, ...
        'channelIDs', opts.channelIDs, ...
        'ptuReaderFolder', opts.ptuReaderFolder, ...
        'darkFile', opts.darkFile, ...
        'darkScale', opts.darkScale, ...
        'backgroundMode', opts.backgroundMode, ...
        'xyPixelSizeUm', opts.xyPixelSizeUm, ...
        'reuseCache', opts.reuseCache, ...
        'cacheFile', opts.cacheFile, ...
        'verbose', opts.verbose);

    [ny, nx, ~, nPlane] = size(data.rawCounts);
    planeIndex = (1:nPlane).';
    stageZUm = data.stageZUm(:);
    centerX = nan(nPlane,1);
    centerY = nan(nPlane,1);
    totalSignal = nan(nPlane,1);
    peakSignal = nan(nPlane,1);
    varXPixel2 = nan(nPlane,1);
    varYPixel2 = nan(nPlane,1);
    secondMomentPixel2 = nan(nPlane,1);
    rmsRadiusPixel = nan(nPlane,1);
    rmsRadiusUm = nan(nPlane,1);

    for iz = 1:nPlane
        img = detectorSummedPlane(data, iz, opts.subtractBackground);
        m = weightedImageMoments(img, opts.thresholdFraction);

        centerX(iz) = m.centerXY(1);
        centerY(iz) = m.centerXY(2);
        totalSignal(iz) = m.totalSignal;
        peakSignal(iz) = m.peakSignal;
        varXPixel2(iz) = m.varX;
        varYPixel2(iz) = m.varY;
        secondMomentPixel2(iz) = m.secondMoment;
        rmsRadiusPixel(iz) = sqrt(max(m.secondMoment, 0));
        if ~isempty(data.xyPixelSizeUm) && isscalar(data.xyPixelSizeUm) && ...
                isfinite(data.xyPixelSizeUm)
            rmsRadiusUm(iz) = rmsRadiusPixel(iz) * data.xyPixelSizeUm;
        end
    end

    signalNorm = safeNormalizeToMax(totalSignal);
    compactness = safeInverseMoment(secondMomentPixel2);
    focusScore = signalNorm.^opts.signalPower .* compactness.^opts.compactnessPower;

    if isempty(opts.focusIndex)
        [~, focusIdx] = max(focusScore);
    else
        focusIdx = round(opts.focusIndex);
        if focusIdx < 1 || focusIdx > nPlane
            error('findBeadCenterFromFullStack:BadFocusIndex', ...
                'focusIndex must be between 1 and %d.', nPlane);
        end
    end

    focusImage = detectorSummedPlane(data, focusIdx, opts.subtractBackground);
    focusMoments = weightedImageMoments(focusImage, opts.thresholdFraction);
    centerXY = focusMoments.centerXY;
    centerYX = [ ...
        min(max(round(centerXY(2)), 1), ny), ...
        min(max(round(centerXY(1)), 1), nx)];

    metrics = table(planeIndex, stageZUm, totalSignal, peakSignal, ...
        signalNorm, compactness, focusScore, centerX, centerY, ...
        varXPixel2, varYPixel2, secondMomentPixel2, rmsRadiusPixel, ...
        rmsRadiusUm, ...
        'VariableNames', {'planeIndex','stageZUm','totalSignal', ...
        'peakSignal','normalizedTotalSignal','normalizedCompactness', ...
        'focusScore','centerXPixel','centerYPixel','varXPixel2', ...
        'varYPixel2','secondMomentPixel2','rmsRadiusPixel','rmsRadiusUm'});

    out = struct();
    out.focusIndex = focusIdx;
    out.focusZUm = stageZUm(focusIdx);
    out.centerXY = centerXY;
    out.centerYX = centerYX;
    out.focusImage = focusImage;
    out.metrics = metrics;
    out.data = data;
    out.options = opts;

    if opts.makeFigure || ~isempty(opts.outputDir)
        out.figure = plotDiagnostics(out, opts);
    end

    if ~isempty(opts.outputDir)
        if exist(opts.outputDir, 'dir') ~= 7
            mkdir(opts.outputDir);
        end
        writetable(metrics, fullfile(opts.outputDir, 'bead_focus_center_metrics.csv'));
        if isfield(out, 'figure') && ishghandle(out.figure)
            exportgraphics(out.figure, ...
                fullfile(opts.outputDir, 'bead_focus_center_diagnostic.png'), ...
                'Resolution', 180);
        end
    end

    if opts.verbose
        fprintf(['[findBeadCenterFromFullStack] focus plane %d/%d, z=%.4f um, ' ...
            'center x=%.2f y=%.2f, rms radius=%.2f px.\n'], ...
            focusIdx, nPlane, out.focusZUm, centerXY(1), centerXY(2), ...
            metrics.rmsRadiusPixel(focusIdx));
    end
end

function opts = parseInputs(varargin)
    p = inputParser;
    p.FunctionName = 'findBeadCenterFromFullStack';

    addParameter(p, 'channelIDs', 9:31);
    addParameter(p, 'ptuReaderFolder', defaultPtuReaderFolder());
    addParameter(p, 'darkFile', 'D:\Luminosa\Data\ISMdark_counts.ptu');
    addParameter(p, 'darkScale', 1);
    addParameter(p, 'backgroundMode', 'auto');
    addParameter(p, 'xyPixelSizeUm', []);
    addParameter(p, 'reuseCache', false);
    addParameter(p, 'cacheFile', '');
    addParameter(p, 'subtractBackground', true);
    addParameter(p, 'thresholdFraction', 0.10);
    addParameter(p, 'signalPower', 1);
    addParameter(p, 'compactnessPower', 1);
    addParameter(p, 'focusIndex', []);
    addParameter(p, 'makeFigure', true);
    addParameter(p, 'figureVisible', 'on');
    addParameter(p, 'outputDir', '');
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;
    opts.channelIDs = double(opts.channelIDs(:)).';
    opts.ptuReaderFolder = char(opts.ptuReaderFolder);
    opts.backgroundMode = lower(char(opts.backgroundMode));
    opts.cacheFile = char(opts.cacheFile);
    opts.darkFile = char(opts.darkFile);
    opts.figureVisible = char(opts.figureVisible);
    opts.outputDir = char(opts.outputDir);

    if opts.thresholdFraction < 0 || opts.thresholdFraction >= 1
        error('findBeadCenterFromFullStack:BadThreshold', ...
            'thresholdFraction must be in [0, 1).');
    end
    if opts.signalPower < 0 || opts.compactnessPower < 0
        error('findBeadCenterFromFullStack:BadScorePowers', ...
            'signalPower and compactnessPower must be nonnegative.');
    end
end

function folder = defaultPtuReaderFolder()
    thisDir = fileparts(mfilename('fullpath'));
    luminosaRoot = fileparts(fileparts(thisDir));
    folder = fullfile(luminosaRoot, 'Luminosa_FLIM_FCS');
end

function img = detectorSummedPlane(data, iz, subtractBackground)
    plane = double(data.rawCounts(:,:,:,iz));
    if subtractBackground
        bg = double(data.backgroundPerPixel(:,:,:,iz));
        plane = max(plane - bg, 0);
    end
    img = sum(plane, 3);
    img(~isfinite(img)) = 0;
    img = max(img, 0);
end

function m = weightedImageMoments(img, thresholdFraction)
    img = double(img);
    img(~isfinite(img)) = 0;
    img = max(img, 0);

    totalSignal = sum(img(:));
    peakSignal = max(img(:));

    weights = img;
    if peakSignal > 0 && thresholdFraction > 0
        weights(weights < thresholdFraction * peakSignal) = 0;
    end

    weightSum = sum(weights(:));
    if weightSum <= 0
        [~, idx] = max(img(:));
        [cy, cx] = ind2sub(size(img), idx);
        centerXY = [cx, cy];
        varX = NaN;
        varY = NaN;
        secondMoment = NaN;
    else
        [yy, xx] = ndgrid(1:size(img,1), 1:size(img,2));
        cx = sum(xx(:) .* weights(:)) / weightSum;
        cy = sum(yy(:) .* weights(:)) / weightSum;
        dx = xx - cx;
        dy = yy - cy;
        varX = sum((dx(:).^2) .* weights(:)) / weightSum;
        varY = sum((dy(:).^2) .* weights(:)) / weightSum;
        secondMoment = varX + varY;
        centerXY = [cx, cy];
    end

    m = struct();
    m.centerXY = centerXY;
    m.totalSignal = totalSignal;
    m.peakSignal = peakSignal;
    m.varX = varX;
    m.varY = varY;
    m.secondMoment = secondMoment;
end

function y = safeNormalizeToMax(x)
    x = double(x(:));
    y = zeros(size(x));
    valid = isfinite(x) & x >= 0;
    if any(valid)
        denom = max(x(valid));
        if denom > 0
            y(valid) = x(valid) / denom;
        end
    end
end

function compactness = safeInverseMoment(moment)
    moment = double(moment(:));
    compactness = zeros(size(moment));
    valid = isfinite(moment) & moment > 0;
    if any(valid)
        best = min(moment(valid));
        compactness(valid) = best ./ moment(valid);
    end
end

function fig = plotDiagnostics(out, opts)
    fig = figure('Color', 'w', 'Visible', opts.figureVisible, ...
        'Position', [80 80 1050 430]);
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    T = out.metrics;
    ax1 = nexttile(tl, 1);
    plot(ax1, T.stageZUm, T.normalizedTotalSignal, 'k-o', 'LineWidth', 1);
    hold(ax1, 'on');
    plot(ax1, T.stageZUm, T.normalizedCompactness, 'b-s', 'LineWidth', 1);
    plot(ax1, T.stageZUm, T.focusScore, 'r-^', 'LineWidth', 1.2);
    xline(ax1, out.focusZUm, 'r--', 'LineWidth', 1);
    hold(ax1, 'off');
    grid(ax1, 'on');
    xlabel(ax1, 'stage z [\mum]');
    ylabel(ax1, 'normalized metric');
    legend(ax1, {'total signal','compactness','focus score','chosen focus'}, ...
        'Location', 'best');
    title(ax1, 'Focus-plane score');

    ax2 = nexttile(tl, 2);
    imagesc(ax2, out.focusImage);
    axis(ax2, 'image');
    colormap(ax2, 'hot');
    colorbar(ax2);
    hold(ax2, 'on');
    plot(ax2, out.centerXY(1), out.centerXY(2), 'co', ...
        'MarkerSize', 9, 'LineWidth', 1.4);
    hold(ax2, 'off');
    xlabel(ax2, 'x pixel');
    ylabel(ax2, 'y pixel');
    title(ax2, sprintf('Chosen plane %d, z = %.3f \\mum', ...
        out.focusIndex, out.focusZUm));
end
