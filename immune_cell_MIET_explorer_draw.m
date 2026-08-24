function immune_cell_MIET_explorer_draw(fig)
%IMMUNE_CELL_MIET_EXPLORER_DRAW Redraw the map for the current dropdown choice.
%
% immune_cell_MIET_explorer_draw(fig)
%
% Separate from the selection callback so changing the displayed map never
% re-reads the TCSPC cube or refits anything: the decay panel keeps showing the
% region you already chose while you flip between maps looking at it. That is
% the point of having a dropdown rather than one figure per map.

    state = fig.UserData;
    index = find(strcmp({state.catalogue.label}, state.mapMenu.Value), 1);
    if isempty(index); index = 1; end
    entry = state.catalogue(index);

    ax = state.mapAxes;
    delete(findobj(ax, 'Type', 'image'));
    data = entry.data;
    image = imagesc(ax, data);
    uistack(image, 'bottom');
    axis(ax, 'image');
    set(ax, 'YDir', 'reverse');

    if entry.isTau
        colormap(ax, parula);
    else
        colormap(ax, gray);
    end

    finite = data(isfinite(data) & data ~= 0);
    if ~isempty(entry.limits)
        caxis(ax, entry.limits);
    elseif numel(finite) > 10
        % Robust limits: a handful of hot pixels otherwise flatten the whole
        % map to one colour and make it useless for navigating.
        lo = immune_cell_MIET_explorer_quantile(finite, 0.02);
        hi = immune_cell_MIET_explorer_quantile(finite, 0.98);
        if hi > lo; caxis(ax, [lo hi]); end
    end
    set(image, 'AlphaData', isfinite(data));

    % Clicking the image is how a region is chosen, so the callback has to be
    % re-attached every time the image object is replaced.
    image.ButtonDownFcn = @(src, evt) immune_cell_MIET_explorer_click(fig, evt);
    ax.ButtonDownFcn = @(src, evt) immune_cell_MIET_explorer_click(fig, evt);

    title(ax, sprintf('%s  [%s]', entry.label, entry.units));
    xlabel(ax, 'column (pixels)');
    ylabel(ax, 'row (pixels)');

    % The marker must be recreated above the new image.
    if isfield(state, 'marker') && isgraphics(state.marker)
        delete(state.marker);
    end
    hold(ax, 'on');
    state.marker = plot(ax, state.colBox.Value, state.rowBox.Value, ...
        '+', 'Color', [1 0.2 0.2], 'MarkerSize', 14, 'LineWidth', 1.5, ...
        'HitTest', 'off', 'PickableParts', 'none');
    hold(ax, 'off');
    fig.UserData = state;
end
