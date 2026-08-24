function value = immune_cell_MIET_explorer_field(root, path)
%IMMUNE_CELL_MIET_EXPLORER_FIELD Fetch a dotted field path, or [] if absent.
%
% value = immune_cell_MIET_explorer_field(result, 'bayesian.maps.tauMeanArithmetic')
% value = immune_cell_MIET_explorer_field(result, 'bayesian.maps.componentLifetimeNs(2)')
%
% Returns [] rather than erroring when any level is missing. That is the point:
% which stages and maps a result contains depends on how the pipeline was
% configured for that run, so the map catalogue has to probe rather than assume.
% An error here would mean a result missing one optional map could not be opened
% at all.
%
% A trailing (k) selects slice k along the THIRD dimension of a 3-D array,
% which is how MATLAB holds the pipeline's per-component maps: [y, x, k]. An
% h5py listing of the same variable shows (k, y, x) because HDF5 reports
% dimensions in reverse order - reading that as the MATLAB layout is what made
% the map catalogue enumerate 602 components instead of 3.

    value = [];
    if ~isstruct(root) || isempty(path); return; end

    index = [];
    token = regexp(path, '^(.*)\((\d+)\)$', 'tokens', 'once');
    if ~isempty(token)
        path = token{1};
        index = str2double(token{2});
    end

    parts = strsplit(path, '.');
    node = root;
    for k = 1:numel(parts)
        name = parts{k};
        if ~isstruct(node) || ~isscalar(node) || ~isfield(node, name)
            return;
        end
        node = node.(name);
    end

    if ~isempty(index)
        if ndims(node) ~= 3 || index < 1 || index > size(node, 3)
            return;
        end
        node = node(:, :, index);
    end
    value = node;
end
