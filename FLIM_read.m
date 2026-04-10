folderName = 'D:\Luminosa\Data\260323\t1_20260323-104842\test\';

exportMode = 'frames';          % 'sum' or 'frames'
movieFormat = 'gif';         % 'gif', 'video', or 'both' (used only for exportMode='frames')
resultFolderName = 'FLIM_results';
fps = 5;
tauRangeNs = [3.5 6.5];
writeMat = true;

intensityLowerPct = 30;      % 30 for cells, 10 for tissue
intensityUpperPct = 99.5;
useSlidingTcspcBin2x2 = true;
slidingTcspcWindow = [2 2];  % [x y] pixels, applied with stride 1

scaleBarUm = 5;
scaleBarLabel = '5 \mum';
scaleBarColor = 'w';
scaleBarLineWidth = 3;
scaleBarFontSize = 20;
timeTagFontSize = 16;

validExportModes = {'sum', 'frames'};
validMovieFormats = {'gif', 'video', 'both'};
if ~any(strcmpi(exportMode, validExportModes))
    error('exportMode must be one of: %s', strjoin(validExportModes, ', '));
end
if strcmpi(exportMode, 'frames') && ~any(strcmpi(movieFormat, validMovieFormats))
    error('movieFormat must be one of: %s', strjoin(validMovieFormats, ', '));
end

ptuTargets = collectPTUTargets(folderName, resultFolderName);
if isempty(ptuTargets)
    error('No PTU inputs were found under %s', folderName);
end

opts = struct();
opts.photonsPerChunk = 5e6;
opts.computePerFrame = strcmpi(exportMode, 'frames');
opts.storeTcspcPix = false;
opts.storePhotonLists = useSlidingTcspcBin2x2;
opts.storeTimeCell = false;
opts.showWaitbar = false;

for k = 1:numel(ptuTargets)
    target = ptuTargets(k);
    name = target.ptuPath;

    fprintf('Processing %d/%d: %s\n', k, numel(ptuTargets), name);

    if ~exist(target.resultDir, 'dir')
        mkdir(target.resultDir);
    end

    out = PTU_FLIM_GPU(name, opts);
    head = out.head;
    headerTimeTag = getHeaderTimeTag(head);

    if writeMat
        save(fullfile(target.resultDir, [target.baseName '_FLIM.mat']), 'out', 'head', '-v7.3');
    end

    switch lower(exportMode)
        case 'sum'
            [intensityMap, tauMap] = buildDisplayMaps(out, [], useSlidingTcspcBin2x2, slidingTcspcWindow);
            overlayTag = composeOverlayTag(headerTimeTag, head, [], out.nFrames, 'sum');
            outPng = fullfile(target.resultDir, [target.baseName '_FLIM_sum.png']);
            exportFlimPanel(intensityMap, tauMap, head, outPng, tauRangeNs, overlayTag, ...
                scaleBarUm, scaleBarLabel, scaleBarColor, scaleBarLineWidth, scaleBarFontSize, timeTagFontSize, ...
                intensityLowerPct, intensityUpperPct);

        case 'frames'
            nFrames = getFrameCount(out);
            frameDir = fullfile(target.resultDir, 'frames');
            if ~exist(frameDir, 'dir')
                mkdir(frameDir);
            end

            movieState = initMovieWriters(target.resultDir, target.baseName, movieFormat, fps, nFrames);
            cleanupMovie = onCleanup(@() closeMovieWriters(movieState)); %#ok<NASGU>

            for frameIdx = 1:nFrames
                [intensityMap, tauMap] = buildDisplayMaps(out, frameIdx, useSlidingTcspcBin2x2, slidingTcspcWindow);
                overlayTag = composeOverlayTag(headerTimeTag, head, frameIdx, nFrames, 'frames');

                outPng = fullfile(frameDir, sprintf('%s_frame_%04d.png', target.baseName, frameIdx));
                exportFlimPanel(intensityMap, tauMap, head, outPng, tauRangeNs, overlayTag, ...
                    scaleBarUm, scaleBarLabel, scaleBarColor, scaleBarLineWidth, scaleBarFontSize, timeTagFontSize, ...
                    intensityLowerPct, intensityUpperPct);
                appendMovieFrame(movieState, outPng, frameIdx, fps);
            end

        otherwise
            error('Unknown exportMode: %s', exportMode);
    end
end


function targets = collectPTUTargets(folderName, resultFolderName)
targets = struct('ptuPath', {}, 'resultDir', {}, 'baseName', {});

dirEntries = dir(folderName);
for k = 1:numel(dirEntries)
    if ~dirEntries(k).isdir
        continue;
    end
    entryName = dirEntries(k).name;
    if strcmp(entryName, '.') || strcmp(entryName, '..')
        continue;
    end

    ptuPath = fullfile(folderName, entryName, 'RawImage.ptu');
    if exist(ptuPath, 'file') == 2
        targets(end+1).ptuPath = ptuPath; %#ok<AGROW>
        targets(end).resultDir = fullfile(folderName, entryName, resultFolderName);
        targets(end).baseName = entryName;
    end
end

ptuFiles = dir(fullfile(folderName, '*.ptu'));
for k = 1:numel(ptuFiles)
    ptuPath = fullfile(folderName, ptuFiles(k).name);
    [~, stem] = fileparts(ptuFiles(k).name);
    targets(end+1).ptuPath = ptuPath; %#ok<AGROW>
    targets(end).resultDir = fullfile(folderName, [stem '_' resultFolderName]);
    targets(end).baseName = stem;
end
end


function nFrames = getFrameCount(out)
nFrames = 1;

if isfield(out, 'nFrames') && ~isempty(out.nFrames) && isfinite(out.nFrames) && out.nFrames > 0
    nFrames = double(out.nFrames);
    return;
end

if isfield(out, 'tag') && ~isempty(out.tag) && ndims(out.tag) >= 4
    nFrames = size(out.tag, 4);
end
end


function [intensityMap, tauMap] = buildDisplayMaps(out, frameIdx, useSlidingTcspcBin2x2, slidingTcspcWindow)
if ~useSlidingTcspcBin2x2
    if isempty(frameIdx)
        intensityMap = collapseIntensity(out.tags);
        tauMap = collapseTau(out.taus);
    else
        intensityMap = collapseIntensity(out.tag(:,:,:,frameIdx));
        tauMap = collapseTau(out.tau(:,:,:,frameIdx));
    end
    return;
end

[cntCube, sum1Cube, sum2Cube] = photonMomentsFromOut(out, frameIdx);
if isempty(cntCube)
    intensityMap = [];
    tauMap = [];
    return;
end

[cntCube, sum1Cube, sum2Cube] = applySlidingWindowToMoments(cntCube, sum1Cube, sum2Cube, slidingTcspcWindow);
tauCube = tauFromMoments(sum1Cube, sum2Cube, cntCube, getResolutionNs(out));

intensityMap = sum(cntCube, 3);
tauMap = collapseTau(tauCube);
end


function intensityMap = collapseIntensity(tagMap)
if isempty(tagMap)
    intensityMap = [];
    return;
end

tagMap = double(tagMap);
if ndims(tagMap) > 2
    intensityMap = sum(tagMap, 3);
else
    intensityMap = tagMap;
end
end


function tauMap = collapseTau(tauCube)
if isempty(tauCube)
    tauMap = [];
    return;
end

tauCube = double(tauCube);
if ndims(tauCube) > 2
    tauMap = mean(tauCube, 3, 'omitnan');
else
    tauMap = tauCube;
end
end


function [cntCube, sum1Cube, sum2Cube] = photonMomentsFromOut(out, frameIdx)
cntCube = [];
sum1Cube = [];
sum2Cube = [];

requiredFields = {'im_col', 'im_line', 'im_chan', 'im_tcspc'};
for k = 1:numel(requiredFields)
    if ~isfield(out, requiredFields{k}) || isempty(out.(requiredFields{k}))
        return;
    end
end

nx = double(out.head.ImgHdr_PixX);
ny = double(out.head.ImgHdr_PixY);

if isfield(out, 'dind') && ~isempty(out.dind)
    dind = double(out.dind(:));
else
    dind = unique(double(out.im_chan(:)));
end
nCh = numel(dind);

x = double(out.im_col(:));
y = double(out.im_line(:));
ch = double(out.im_chan(:));
t = double(out.im_tcspc(:));

if ~isempty(frameIdx) && isfield(out, 'im_frame') && ~isempty(out.im_frame)
    keep = double(out.im_frame(:)) == double(frameIdx);
    x = x(keep);
    y = y(keep);
    ch = ch(keep);
    t = t(keep);
end

cntCube = zeros(nx, ny, nCh, 'double');
sum1Cube = zeros(nx, ny, nCh, 'double');
sum2Cube = zeros(nx, ny, nCh, 'double');

validBase = isfinite(x) & isfinite(y) & isfinite(ch) & isfinite(t) & ...
    x >= 1 & x <= nx & y >= 1 & y <= ny;

x = x(validBase);
y = y(validBase);
ch = ch(validBase);
t = t(validBase);

for kc = 1:nCh
    ind = ch == dind(kc);
    if ~any(ind)
        continue;
    end

    pixIdx = x(ind) + (y(ind) - 1) * nx;
    cnt = accumarray(pixIdx, 1, [nx * ny, 1], @sum, 0);
    sum1 = accumarray(pixIdx, t(ind), [nx * ny, 1], @sum, 0);
    sum2 = accumarray(pixIdx, t(ind).^2, [nx * ny, 1], @sum, 0);

    cntCube(:,:,kc) = reshape(cnt, [nx, ny]);
    sum1Cube(:,:,kc) = reshape(sum1, [nx, ny]);
    sum2Cube(:,:,kc) = reshape(sum2, [nx, ny]);
end
end


function [cntOut, sum1Out, sum2Out] = applySlidingWindowToMoments(cntCube, sum1Cube, sum2Cube, slidingTcspcWindow)
winX = max(1, round(slidingTcspcWindow(1)));
if numel(slidingTcspcWindow) >= 2
    winY = max(1, round(slidingTcspcWindow(2)));
else
    winY = winX;
end

kernel = ones(winX, winY);

cntOut = zeros(size(cntCube), 'double');
sum1Out = zeros(size(sum1Cube), 'double');
sum2Out = zeros(size(sum2Cube), 'double');

for kc = 1:size(cntCube, 3)
    cntOut(:,:,kc) = conv2(cntCube(:,:,kc), kernel, 'same');
    sum1Out(:,:,kc) = conv2(sum1Cube(:,:,kc), kernel, 'same');
    sum2Out(:,:,kc) = conv2(sum2Cube(:,:,kc), kernel, 'same');
end
end


function tauCube = tauFromMoments(sum1Cube, sum2Cube, cntCube, resolutionNs)
tauCube = zeros(size(cntCube), 'double');
valid = cntCube > 0;

mean1 = zeros(size(cntCube), 'double');
mean2 = zeros(size(cntCube), 'double');

mean1(valid) = (sum1Cube(valid) ./ cntCube(valid)) * resolutionNs;
mean2(valid) = (sum2Cube(valid) ./ cntCube(valid)) * (resolutionNs^2);

varT = max(mean2 - mean1.^2, 0);
tauCube(valid) = sqrt(varT(valid));
end


function resolutionNs = getResolutionNs(out)
resolutionNs = [];

if isfield(out, 'Resolution_ns') && ~isempty(out.Resolution_ns)
    resolutionNs = double(out.Resolution_ns);
elseif isfield(out, 'head') && isfield(out.head, 'MeasDesc_Resolution') && ~isempty(out.head.MeasDesc_Resolution)
    resolutionNs = double(out.head.MeasDesc_Resolution) * 1e9;
end

if isempty(resolutionNs) || ~isfinite(resolutionNs) || resolutionNs <= 0
    resolutionNs = 1;
end
end


function exportFlimPanel(intensityMap, tauMap, head, outPng, tauRangeNs, overlayTag, ...
    scaleBarUm, scaleBarLabel, scaleBarColor, scaleBarLineWidth, scaleBarFontSize, timeTagFontSize, ...
    intensityLowerPct, intensityUpperPct)

if isempty(intensityMap) || isempty(tauMap)
    warning('Skipping %s because no displayable FLIM maps were available.', outPng);
    return;
end

brightnessMap = normalizeIntensity(intensityMap, intensityLowerPct, intensityUpperPct);
brightnessMap = log10(0.001 + brightnessMap.^2);

cmap = jet(256);
cmap = cmap(30:end-30, :);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 760]);
cleanupFig = onCleanup(@() close(fig)); %#ok<NASGU>

hIm = cim(double(tauMap), brightnessMap, tauRangeNs, 'v', cmap);
ax = ancestor(hIm, 'axes');

set(fig, 'InvertHardcopy', 'off');
axes(ax);

addPTUScaleBar(hIm, head, size(tauMap), 'bottomright', ...
    'BarLengthUm', scaleBarUm, ...
    'Label', scaleBarLabel, ...
    'Color', scaleBarColor, ...
    'LineWidth', scaleBarLineWidth, ...
    'FontSize', scaleBarFontSize);

if ~isempty(overlayTag)
    text(ax, 0.98, 0.98, overlayTag, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'Color', 'w', ...
        'FontSize', timeTagFontSize, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none', ...
        'Clipping', 'off');
end

drawnow;
print(fig, outPng, '-dpng', '-r300');
end


function intensityNorm = normalizeIntensity(intensityMap, lowerPct, upperPct)
intensityMap = double(intensityMap);
intensityNorm = zeros(size(intensityMap));

validVals = intensityMap(isfinite(intensityMap));
if isempty(validVals)
    return;
end

lb = prctile(validVals, lowerPct);
ub = prctile(validVals, upperPct);

if ~isfinite(lb)
    lb = min(validVals);
end
if ~isfinite(ub)
    ub = max(validVals);
end
if ub <= lb
    ub = lb + eps;
end

intensityMap = min(max(intensityMap, lb), ub);
intensityMap = intensityMap - lb;
scale = max(intensityMap(:));
if scale > 0
    intensityNorm = intensityMap ./ scale;
end
end


function movieState = initMovieWriters(resultDir, baseName, movieFormat, fps, nFrames)
movieState = struct();
movieState.mode = lower(movieFormat);
movieState.gifFile = '';
movieState.videoWriter = [];
movieState.writeGif = strcmpi(movieFormat, 'gif') || strcmpi(movieFormat, 'both');
movieState.writeVideo = strcmpi(movieFormat, 'video') || strcmpi(movieFormat, 'both');

if nFrames <= 1
    movieState.writeGif = false;
    movieState.writeVideo = false;
    return;
end

if movieState.writeGif
    movieState.gifFile = fullfile(resultDir, [baseName '_FLIM_frames.gif']);
    if exist(movieState.gifFile, 'file') == 2
        delete(movieState.gifFile);
    end
end

if movieState.writeVideo
    videoFile = fullfile(resultDir, [baseName '_FLIM_frames.mp4']);
    try
        vw = VideoWriter(videoFile, 'MPEG-4');
    catch
        videoFile = fullfile(resultDir, [baseName '_FLIM_frames.avi']);
        vw = VideoWriter(videoFile, 'Motion JPEG AVI');
    end
    vw.FrameRate = fps;
    open(vw);
    movieState.videoWriter = vw;
end
end


function appendMovieFrame(movieState, pngFile, frameIdx, fps)
if ~movieState.writeGif && ~movieState.writeVideo
    return;
end

rgb = imread(pngFile);
if ndims(rgb) == 2
    rgb = repmat(rgb, [1 1 3]);
end

if movieState.writeGif
    [imind, cm] = rgb2ind(rgb, 256);
    delayTime = 1 / max(fps, eps);
    if frameIdx == 1
        imwrite(imind, cm, movieState.gifFile, 'gif', 'LoopCount', inf, 'DelayTime', delayTime);
    else
        imwrite(imind, cm, movieState.gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end

if movieState.writeVideo && ~isempty(movieState.videoWriter)
    writeVideo(movieState.videoWriter, rgb);
end
end


function closeMovieWriters(movieState)
if isfield(movieState, 'videoWriter') && ~isempty(movieState.videoWriter)
    try
        close(movieState.videoWriter);
    catch
    end
end
end


function overlayTag = composeOverlayTag(headerTimeTag, head, frameIdx, nFrames, exportMode)
lines = {};

if ~isempty(headerTimeTag)
    lines{end+1} = headerTimeTag; %#ok<AGROW>
end

frameTimeSec = estimateFrameTimeSeconds(head, nFrames);

switch lower(exportMode)
    case 'sum'
        if nargin >= 4 && ~isempty(nFrames) && nFrames > 1
            if isfinite(frameTimeSec)
                totalTimeSec = frameTimeSec * nFrames;
                lines{end+1} = sprintf('summed %d frames | %.2f s total', nFrames, totalTimeSec); %#ok<AGROW>
            else
                lines{end+1} = sprintf('summed %d frames', nFrames); %#ok<AGROW>
            end
        end

    case 'frames'
        if nargin >= 3 && ~isempty(frameIdx)
            if isfinite(frameTimeSec)
                relTimeSec = frameTimeSec * (frameIdx - 1);
                lines{end+1} = sprintf('frame %d/%d | t = %.2f s', frameIdx, nFrames, relTimeSec); %#ok<AGROW>
            else
                lines{end+1} = sprintf('frame %d/%d', frameIdx, nFrames); %#ok<AGROW>
            end
        end
end

overlayTag = strjoin(lines, newline);
end


function frameTimeSec = estimateFrameTimeSeconds(head, nFrames)
frameTimeSec = NaN;

if isfield(head, 'ImgHdr_TimePerPixel') && isfield(head, 'ImgHdr_PixX') && isfield(head, 'ImgHdr_PixY') ...
        && ~isempty(head.ImgHdr_TimePerPixel) && ~isempty(head.ImgHdr_PixX) && ~isempty(head.ImgHdr_PixY)
    frameTimeSec = double(head.ImgHdr_TimePerPixel) * double(head.ImgHdr_PixX) * double(head.ImgHdr_PixY) / 1e3;
end

if (~isfinite(frameTimeSec) || frameTimeSec <= 0) && isfield(head, 'TTResult_StopAfter') ...
        && ~isempty(head.TTResult_StopAfter) && nargin >= 2 && ~isempty(nFrames) && nFrames > 0
    frameTimeSec = double(head.TTResult_StopAfter) / double(nFrames) / 1e3;
end

if ~isfinite(frameTimeSec) || frameTimeSec <= 0
    frameTimeSec = NaN;
end
end


function timeTag = getHeaderTimeTag(head)
timeTag = '';

preferredFields = {'File_CreatingTime', 'Measurement_StartTime', 'Measurement_StopTime'};
for k = 1:numel(preferredFields)
    if isfield(head, preferredFields{k})
        timeTag = headerTimeOnly(head.(preferredFields{k}));
        if ~isempty(timeTag)
            return;
        end
    end
end

allFields = fieldnames(head);
for k = 1:numel(allFields)
    fieldName = allFields{k};
    lowerName = lower(fieldName);
    if contains(lowerName, 'time') || contains(lowerName, 'date')
        timeTag = headerTimeOnly(head.(fieldName));
        if ~isempty(timeTag)
            return;
        end
    end
end
end


function txt = headerTimeOnly(value)
txt = '';

if isstring(value)
    if isscalar(value)
        value = char(value);
    else
        return;
    end
end

if ischar(value)
    value = strtrim(value);
    if isempty(value)
        return;
    end

    tok = regexp(value, '\d{1,2}:\d{2}(:\d{2})?(\s*[APMapm]{2})?', 'match', 'once');
    if ~isempty(tok)
        txt = strtrim(tok);
        return;
    end

    try
        dt = datetime(value, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss', 'Locale', 'en_US');
        txt = char(string(dt, 'HH:mm:ss'));
        return;
    catch
    end

    try
        dt = datetime(value);
        txt = char(string(dt, 'HH:mm:ss'));
        return;
    catch
    end
end
end
