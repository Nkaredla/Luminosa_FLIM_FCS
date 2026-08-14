function normalised = robustNormalise(image)
%ROBUSTNORMALISE Scale an image using its 1st and 99.5th percentiles.
% The explicit sorted-index implementation avoids a Statistics Toolbox
% dependency and limits the influence of isolated bright pixels.

    image = double(image);
    values = sort(image(isfinite(image)));
    if isempty(values)
        normalised = zeros(size(image));
        return;
    end

    lowValue = values(max(1, round(0.01 * numel(values))));
    highValue = values(max(1, round(0.995 * numel(values))));
    if highValue <= lowValue
        highValue = max(values);
    end
    normalised = min(max((image - lowValue) ./ max(highValue - lowValue, eps), 0), 1);
    normalised(~isfinite(normalised)) = 0;
end
