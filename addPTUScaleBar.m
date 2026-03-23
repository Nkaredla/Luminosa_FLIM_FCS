function h = addPTUScaleBar(axOrIm, head, imSize, location, varargin)
% addPTUScaleBar
% Draw a scale bar on the image axes using PTU header metadata.
%
% axOrIm can be:
%   - an image handle returned by cim(...)
%   - an axes handle
%
% If a colorbar axes is passed accidentally, the function will try to find
% the main image axes in the same figure.

if nargin < 4 || isempty(location)
    location = 'topright';
end

% ------------------------------------------------------------
% Resolve the correct image axes
% ------------------------------------------------------------
ax = resolveImageAxes(axOrIm);

p = inputParser;
p.addParameter('BarLengthUm', []);
p.addParameter('Color', 'w');
p.addParameter('LineWidth', 3);
p.addParameter('FontSize', 11);
p.addParameter('YOffsetFrac', 0.06);
p.addParameter('XOffsetFrac', 0.06);
p.addParameter('TextGapFrac', 0.02);
p.addParameter('Label', '');
p.parse(varargin{:});
opt = p.Results;

% --- PTU pixel size ---
if isfield(head, 'ImgHdr_PixResol') && ~isempty(head.ImgHdr_PixResol)
    pixSize = double(head.ImgHdr_PixResol);
else
    error('PTU header does not contain head.ImgHdr_PixResol.');
end

% Convert to microns if stored in meters
if pixSize < 1e-3
    pixSizeUm = pixSize * 1e6;
else
    pixSizeUm = pixSize;
end

% Native PTU dimensions if available
if isfield(head, 'ImgHdr_PixX') && ~isempty(head.ImgHdr_PixX)
    nativeNx = double(head.ImgHdr_PixX);
else
    nativeNx = imSize(2);
end
if isfield(head, 'ImgHdr_PixY') && ~isempty(head.ImgHdr_PixY)
    nativeNy = double(head.ImgHdr_PixY);
else
    nativeNy = imSize(1);
end

imgNy = imSize(1);
imgNx = imSize(2);

pixSizeXUm = pixSizeUm * nativeNx / imgNx;

% --- choose a nice bar length if not supplied ---
if isempty(opt.BarLengthUm)
    targetFrac = 0.18;
    targetUm = targetFrac * imgNx * pixSizeXUm;
    niceVals = [0.2 0.5 1 2 5 10 20 25 50 100 200 500 1000];
    [~, idx] = min(abs(log10(niceVals ./ targetUm)));
    barLengthUm = niceVals(idx);
else
    barLengthUm = opt.BarLengthUm;
end

barLengthPx = barLengthUm / pixSizeXUm;

xMargin = opt.XOffsetFrac * imgNx;
yMargin = opt.YOffsetFrac * imgNy;
textGap = opt.TextGapFrac * imgNy;

switch lower(location)
    case 'topright'
        x2 = imgNx - xMargin;
        x1 = x2 - barLengthPx;
        y  = 1 + yMargin;
        txtY = y - 0.03*imgNy;   % text above bar in image display
        txtVA = 'bottom';

    case 'topleft'
        x1 = 1 + xMargin;
        x2 = x1 + barLengthPx;
        y  = 1 + yMargin;
        txtY = y - 0.03*imgNy;
        txtVA = 'bottom';

    case 'bottomright'
        x2 = imgNx - xMargin;
        x1 = x2 - barLengthPx;
        y  = imgNy - yMargin;
        txtY = y - 0.03*imgNy;   % text above bar
        txtVA = 'bottom';

    case 'bottomleft'
        x1 = 1 + xMargin;
        x2 = x1 + barLengthPx;
        y  = imgNy - yMargin;
        txtY = y - 0.03*imgNy;
        txtVA = 'bottom';

    otherwise
        error('Unknown location: %s', location);
end

if isempty(opt.Label)
    if barLengthUm >= 1000
        lbl = sprintf('%g mm', barLengthUm/1000);
    else
        lbl = sprintf('%g \\mum', barLengthUm);
    end
else
    lbl = opt.Label;
end

hold(ax, 'on');
h.line = plot(ax, [x1 x2], [y y], '-', ...
    'Color', opt.Color, 'LineWidth', opt.LineWidth, 'Clipping', 'off');

h.text = text(ax, (x1+x2)/2, txtY, lbl, ...
    'Color', opt.Color, ...
    'FontSize', opt.FontSize, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', txtVA, ...
    'Interpreter', 'tex', ...
    'Clipping', 'off');
hold(ax, 'off');

h.axes = ax;
h.pixSizeUm = pixSizeUm;
h.pixSizeXUmDisplayed = pixSizeXUm;
h.barLengthUm = barLengthUm;
h.barLengthPx = barLengthPx;
end


function ax = resolveImageAxes(axOrIm)
% Accept image handle or axes handle; if the supplied axes is not the main
% image axes, try to find the image axes in the same figure.

if ~ishandle(axOrIm)
    error('First input must be an image handle or axes handle.');
end

typ = get(axOrIm, 'Type');

if strcmpi(typ, 'image')
    ax = ancestor(axOrIm, 'axes');
    return
elseif strcmpi(typ, 'axes')
    ax = axOrIm;
else
    error('First input must be an image handle or axes handle.');
end

% If this axes already contains an image, use it
if ~isempty(findobj(ax, 'Type', 'image'))
    return
end

% Otherwise search figure for the most likely image axes
fig = ancestor(ax, 'figure');
axs = findall(fig, 'Type', 'axes');

candidates = [];
areas = [];

for k = 1:numel(axs)
    if ~isempty(findobj(axs(k), 'Type', 'image'))
        candidates(end+1) = axs(k); %#ok<AGROW>
        pos = get(axs(k), 'Position');
        areas(end+1) = pos(3) * pos(4); %#ok<AGROW>
    end
end

if isempty(candidates)
    error('Could not find an image axes in this figure.');
end

[~, idx] = max(areas);
ax = candidates(idx);
end