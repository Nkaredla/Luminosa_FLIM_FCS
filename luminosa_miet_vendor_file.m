function filePath = luminosa_miet_vendor_file(fileName)
% Resolve a support file shipped inside the vendored MIET-GUI tree.

    vendorRoot = luminosa_miet_vendor_root();
    filePath = '';
    if isempty(vendorRoot)
        return;
    end

    directCandidates = { ...
        fullfile(vendorRoot, fileName), ...
        fullfile(vendorRoot, 'subroutines', fileName), ...
        fullfile(vendorRoot, 'Documentation', fileName)};

    for idx = 1:numel(directCandidates)
        if isfile(directCandidates{idx})
            filePath = directCandidates{idx};
            return;
        end
    end

    matches = dir(fullfile(vendorRoot, '**', fileName));
    if ~isempty(matches)
        filePath = fullfile(matches(1).folder, matches(1).name);
    end
end
