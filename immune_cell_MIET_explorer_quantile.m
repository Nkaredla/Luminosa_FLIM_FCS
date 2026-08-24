function value = immune_cell_MIET_explorer_quantile(values, level)
%IMMUNE_CELL_MIET_EXPLORER_QUANTILE Quantile without the Statistics Toolbox.
%
% The rest of this project deliberately avoids assuming the Statistics Toolbox
% is licensed, so prctile is not used for display limits.

    values = sort(double(values(:)));
    if isempty(values); value = NaN; return; end
    if isscalar(values); value = values; return; end
    position = min(max(level * numel(values) + 0.5, 1), numel(values));
    low = floor(position); high = ceil(position);
    if low == high
        value = values(low);
    else
        value = values(low) + (position - low) * (values(high) - values(low));
    end
end
