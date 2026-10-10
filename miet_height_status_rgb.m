function rgb = miet_height_status_rgb(heightNm, statusCode, limits, cmap, varargin)
%MIET_HEIGHT_STATUS_RGB Colour a MIET height map without hiding out-of-range pixels.
%
% rgb = miet_height_status_rgb(heightNm, statusCode, limits, cmap)
% rgb = miet_height_status_rgb(..., 'Intensity', I, 'Background', [r g b])
%
% The ONE place that decides how a height map is coloured, shared by
% Luminosa_GUI (single file) and the batch figures (height maps, profiles,
% topography), so both paths always show the same thing.
%
% statusCode (same codes everywhere):
%   0  no lifetime / not evaluated   -> Background
%   1  inverted                      -> colormap at the pixel's height
%   2  below the calibration floor   -> grey  (more quenched than the curve)
%   3  above the calibration ceiling -> white (lifetime at or beyond the last
%                                        invertible value: FAR from the surface,
%                                        beyond MIET range)
%
% Why this exists: a pixel without a height used to be drawn with the
% colormap's LOWEST colour (min/max replace NaN by the bound), so the regions
% with the longest lifetime - the farthest from the metal - showed up as the
% closest contact. Never let a missing height borrow a colour from the scale.
%
% statusCode may be empty (results made before status was stored): finite
% heights are then treated as 1 and the rest as 0.
%
% 'Intensity' (same size, 0..1) scales every colour, as in an intensity-
% weighted FLIM image. 'Background' defaults to black.

    opts = struct('Intensity', [], 'Background', [0 0 0]);
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k + 1};
    end
    style = miet_out_of_range_style();

    heightNm = double(heightNm);
    if isempty(statusCode)
        statusCode = uint8(isfinite(heightNm));
    end
    statusCode = double(statusCode);
    inverted = statusCode == 1 & isfinite(heightNm);

    lo = limits(1);
    hi = max(limits(2), lo + eps);
    nMap = size(cmap, 1);
    scaled = (heightNm - lo) / (hi - lo);
    scaled(~inverted) = 0;
    idx = 1 + round(min(max(scaled, 0), 1) * (nMap - 1));

    rgb = zeros([size(heightNm) 3]);
    for c = 1:3
        channel = reshape(cmap(idx, c), size(heightNm));
        channel(~inverted) = opts.Background(c);
        channel(statusCode == 3) = style.aboveColour(c);
        channel(statusCode == 2) = style.belowColour(c);
        rgb(:, :, c) = channel;
    end

    if ~isempty(opts.Intensity)
        inten = double(opts.Intensity);
        inten(~isfinite(inten)) = 0;
        rgb = rgb .* inten;
    end
    rgb = max(min(rgb, 1), 0);
end
