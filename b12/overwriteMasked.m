function target = overwriteMasked(target, source, mask)
%OVERWRITEMASKED Copy selected pixels without changing fallback pixels.

    source = single(source);
    target(mask) = source(mask);
end
