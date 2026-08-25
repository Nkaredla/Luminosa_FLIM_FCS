function mask = immune_cell_MIET_biexp_mask(result, request, nRow, nCol)
%IMMUNE_CELL_MIET_BIEXP_MASK Which pixels to fit.
%
% 'cellFootprint' is the default because it is the set the three-model pipeline
% fitted, so the new maps are directly comparable with the old ones rather than
% differing partly because a different set of pixels was used.
    if islogical(request) || isnumeric(request)
        mask = logical(request);
        if ~isequal(size(mask), [nRow nCol])
            error('immune_cell_MIET_biexp_mask:Size', ...
                'A supplied mask must be %d x %d.', nRow, nCol);
        end
        return;
    end
    request = char(request);
    if strcmpi(request, 'all')
        mask = true(nRow, nCol);
        return;
    end
    value = immune_cell_MIET_explorer_field(result, ...
        ['segmentation.masks.' request]);
    if isempty(value)
        error('immune_cell_MIET_biexp_mask:Unknown', ...
            ['segmentation.masks.%s is not present. Available options are ' ...
             'the mask names in the result, or ''all''.'], request);
    end
    mask = logical(value);
end
