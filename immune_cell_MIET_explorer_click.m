function immune_cell_MIET_explorer_click(fig, event)
%IMMUNE_CELL_MIET_EXPLORER_CLICK Turn a click on the map into a selection.
%
% immune_cell_MIET_explorer_click(fig, event)
%
% Kept as its own file so the click handler and the coordinate boxes reach the
% same selection code rather than each having their own rounding.

    if ~isprop(event, 'IntersectionPoint') || numel(event.IntersectionPoint) < 2
        return;
    end
    state = fig.UserData;
    col = round(event.IntersectionPoint(1));
    row = round(event.IntersectionPoint(2));
    col = min(max(col, 1), state.imageSize(2));
    row = min(max(row, 1), state.imageSize(1));
    immune_cell_MIET_explorer_select(fig, row, col);
end
