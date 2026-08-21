function ensureDir(directory)
    import membrane_tracking.curved_miet.*

    if exist(directory, 'dir') ~= 7
        mkdir(directory);
    end
end
