folderName = 'D:\Luminosa\Data\260323\t1_20260323-103604\';
pattern = 'Series_*.ptu';

files = dir(fullfile(folderName, pattern));
if isempty(files)
    error('No files matching %s in %s', pattern, folderName);
end

% Sort by numeric index in Series_xx.ptu
idx = nan(numel(files), 1);
for k = 1:numel(files)
    tok = regexp(files(k).name, 'Series_(\d+)\.ptu', 'tokens', 'once');
    if ~isempty(tok)
        idx(k) = str2double(tok{1});
    end
end
if all(isfinite(idx))
    [~, order] = sort(idx);
    files = files(order);
else
    [~, order] = sort({files.name});
    files = files(order);
end

% Use line-scan reader if only a single series file is present
useLineScan = (numel(files) == 1);

% Video settings
videoFile = fullfile(folderName, 'FLIM_Intensity_series.mp4');
fps = 5;
tauRangeNs = [3.5 6.5];

try
    vw = VideoWriter(videoFile, 'MPEG-4');
catch
    vw = VideoWriter(strrep(videoFile, '.mp4', '.avi'), 'Motion JPEG AVI');
end
vw.FrameRate = fps;
open(vw);

fig = figure('Color','w', 'Position', [100 100 1200 500]);
ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);

cmap = jet(256);
cmap = cmap(30:end-30,:);

for k = 1:numel(files)
    name = fullfile(folderName, files(k).name);
    if useLineScan
        [int2, tau2, head] = readLineScanToMaps(name);
    else
        head = PTU_Read_Head(name);
        out = PTU_MultiFrameScanReadFast(name); % read the channels

        int = out.tags;
        tauc = out.taus;

        % Reduce channel dimension if present
        if ndims(int) > 2
            int2 = sum(int, 3);
        else
            int2 = int;
        end
        if ndims(tauc) > 2
            tau2 = mean(tauc, 3);
        else
            tau2 = tauc;
        end
    end

    % intensity clipping and normalization
    lb = prctile(int2(:), 30); % lower bound, 30 for cells, 10 for tissue
    ub = prctile(int2(:), 99.5);
    int2 = min(max(int2, lb), ub);
    int2 = int2 - lb;
    int2 = int2 ./ max(int2(:));

    % Intensity panel
    imagesc(ax1, int2);
    axis(ax1, 'image');
    axis(ax1, 'off');
    colormap(ax1, gray);
    title(ax1, sprintf('Intensity: %s', files(k).name), 'Interpreter', 'none');

    % FLIM panel (tau color, intensity brightness)
    rgb = flim_rgb(tau2, int2, tauRangeNs, cmap);
    image(ax2, rgb);
    axis(ax2, 'image');
    axis(ax2, 'off');
    title(ax2, sprintf('FLIM: %s', files(k).name), 'Interpreter', 'none');

    % Optional scale bar on FLIM
    addPTUScaleBar(ax2.Children(1), head, size(tau2), 'bottomright', ...
        'BarLengthUm', 5, ...
        'Label', '5 \mum', ...
        'Color', 'w', ...
        'LineWidth', 2, ...
        'FontSize', 16);

    drawnow;
    frame = getframe(fig);
    writeVideo(vw, frame);
end

close(vw);

function rgb = flim_rgb(tau, inten, tauRange, cmap)
% Build FLIM RGB image from tau map and intensity brightness.
    tau = double(tau);
    inten = double(inten);
    tmin = tauRange(1);
    tmax = tauRange(2);
    if tmax <= tmin
        tmax = tmin + eps;
    end
    tau = min(max(tau, tmin), tmax);
    tnorm = (tau - tmin) / (tmax - tmin);
    idx = 1 + tnorm * (size(cmap,1) - 1);

    rgb = zeros([size(tau) 3]);
    for c = 1:3
        rgb(:,:,c) = interp1(1:size(cmap,1), cmap(:,c), idx, 'linear', 'extrap') .* inten;
    end
    rgb = max(min(rgb, 1), 0);
end

function [int2, tau2, head] = readLineScanToMaps(name)
% Read a single series file with line-scan reader and build intensity/tau maps.
    head = PTU_Read_Head(name);
    nRec = head.TTResult_NumberOfRecords;
    nx = double(head.ImgHdr_PixX);
    if isempty(nRec) || nRec <= 0
        [head, im_sync, im_tcspc, im_chan, im_line, im_col] = LPTU_LineScanRead(name, 1e6, nx); %#ok<ASGLU>
    else
        [head, im_sync, im_tcspc, im_chan, im_line, im_col] = LPTU_LineScanRead(name, [1 nRec], nx); %#ok<ASGLU>
    end

    if isempty(im_line) || isempty(im_col)
        int2 = zeros(1, nx);
        tau2 = zeros(1, nx);
        return;
    end

    nLines = double(max(im_line));
    dt_ns = head.MeasDesc_Resolution * 1e9;

    pixIdx = sub2ind([nLines, nx], double(im_line), double(im_col));
    nPix = nLines * nx;

    cnt  = accumarray(pixIdx, 1, [nPix, 1], @sum, 0);
    sum1 = accumarray(pixIdx, double(im_tcspc), [nPix, 1], @sum, 0);
    sum2 = accumarray(pixIdx, double(im_tcspc).^2, [nPix, 1], @sum, 0);

    denom = max(cnt, 1);
    m1 = sum1 ./ denom;
    m2 = sum2 ./ denom;
    tau_rms = sqrt(max(m2 - m1.^2, 0)) * dt_ns;

    int2 = reshape(double(cnt), [nLines, nx]);
    tau2 = reshape(double(tau_rms), [nLines, nx]);
end
