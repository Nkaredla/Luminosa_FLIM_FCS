function analysisMat = immune_cell_MIET_biexp_resolve(source)
%IMMUNE_CELL_MIET_BIEXP_RESOLVE Accept a MAT path, or a folder to search.
    source = char(source);
    if isfolder(source)
        candidate = fullfile(source, 'immune_cell_MIET_640nm_red_analysis.mat');
        if isfile(candidate); analysisMat = candidate; return; end
        found = dir(fullfile(source, '**', ...
            'immune_cell_MIET_640nm_red_analysis.mat'));
        if isempty(found)
            error('immune_cell_MIET_biexp_resolve:NoMat', ...
                'No analysis MAT below %s', source);
        end
        analysisMat = fullfile(found(1).folder, found(1).name);
    elseif isfile(source)
        analysisMat = source;
    else
        error('immune_cell_MIET_biexp_resolve:NoSource', ...
            'Not a file or folder: %s', source);
    end
end
