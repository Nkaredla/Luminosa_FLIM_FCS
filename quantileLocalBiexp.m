function value = quantileLocalBiexp(values, level)
%QUANTILELOCALBIEXP Quantile without the Statistics Toolbox.
% The rest of this project deliberately avoids assuming that toolbox is
% licensed, so prctile is not used.
    values = sort(double(values(:)));
    values = values(isfinite(values));
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
