function data = applyDisplayMask(data, out)
%APPLYDISPLAYMASK Blank map pixels outside the display region.
%
% data = applyDisplayMask(data, out)
%
% The fit now covers the WHOLE file, because including the bare SLB makes the
% anchor self-checking - tau1 there must come back at the measured value. But the
% bare SLB is not what the figures are about, and showing it stretches every
% colour scale over a region whose answer is already known, flattening the
% contrast inside the cell where the structure is.
%
% So display and fit are separated: everything is fitted, only the cell is drawn.
% Pixels outside out.maps.displayMask are set to NaN, which drawMap renders
% transparent, so the colour limits are set by the cell alone.

    if ~isstruct(out) || ~isfield(out, 'maps') || ...
            ~isfield(out.maps, 'displayMask') || isempty(out.maps.displayMask)
        return;
    end
    keep = logical(out.maps.displayMask);
    if ~isequal(size(keep), size(data))
        return;
    end
    data(~keep) = NaN;
end
