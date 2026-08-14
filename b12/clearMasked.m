function target = clearMasked(target, mask)
%CLEARMASKED Mark an unused component as NaN at selected pixels.

    target(mask) = NaN;
end
