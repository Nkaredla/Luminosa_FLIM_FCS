function target = overwriteConstant(target, mask, value)
%OVERWRITECONSTANT Assign one scalar value at selected pixels.

    target(mask) = single(value);
end
