function unwrapped = bnpUnwrapTrack(observed, opts)
    import membrane_tracking.focused_ism.*

    unwrapped = observed;
    for row = 2:size(observed, 1)
        increment = observed(row,:) - observed(row-1,:);
        increment = applyMinimumImage(increment, opts);
        unwrapped(row,:) = unwrapped(row-1,:) + increment;
    end
end
