function text = maskLabel(request)
%MASKLABEL Human-readable name for a pixel-mask request.
    if islogical(request) || isnumeric(request)
        text = sprintf('supplied mask, %d px', nnz(request));
    else
        text = char(request);
    end
end
