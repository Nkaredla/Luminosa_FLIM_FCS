function globalIndices = mapLocalToGlobal(localIndices, globalLookup)
    import membrane_tracking.focused_ism.*

    globalIndices = zeros(size(localIndices));
    matched = localIndices > 0;
    globalIndices(matched) = globalLookup(localIndices(matched));
end
