function image = maskedConstant(mask, value)
%MASKEDCONSTANT Create a constant single map with NaN outside a region.

    image = nan(size(mask), 'single');
    image(mask) = single(value);
end
