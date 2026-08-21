function ensureDir(directory)
    import membrane_tracking.fluctuating_miet.*

    if exist(directory,'dir') ~= 7
        mkdir(directory);
    end
end
