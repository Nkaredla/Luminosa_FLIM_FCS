function ensureDir(d)
    import membrane_tracking.focused_ism.*

    if ~isempty(d) && exist(d, 'dir') ~= 7
        mkdir(d);
    end
end
