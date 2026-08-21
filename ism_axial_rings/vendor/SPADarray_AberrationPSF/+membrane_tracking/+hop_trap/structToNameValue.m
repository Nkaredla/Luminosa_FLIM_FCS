function args = structToNameValue(s)
%STRUCTTONAMEVALUE Convert a scalar structure to name-value arguments.

    names = fieldnames(s);
    args = cell(1, 2 * numel(names));
    for k = 1:numel(names)
        args{2*k-1} = names{k};
        args{2*k} = s.(names{k});
    end
end
