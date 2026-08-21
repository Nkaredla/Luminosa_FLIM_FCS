function [xPlot, yPlot] = breakPeriodicSegments(x, y, opts)
    import membrane_tracking.focused_ism.*

    x = x(:);
    y = y(:);
    if ~strcmp(opts.boundaryCondition, 'periodic') || numel(x) < 2
        xPlot = x;
        yPlot = y;
        return;
    end

    startsNewSegment = [false; ...
        abs(diff(x)) > opts.boxSizeUm(1) / 2 | ...
        abs(diff(y)) > opts.boxSizeUm(2) / 2];
    outputLength = numel(x) + sum(startsNewSegment);
    xPlot = nan(outputLength, 1);
    yPlot = nan(outputLength, 1);
    cursor = 0;
    for k = 1:numel(x)
        if startsNewSegment(k)
            cursor = cursor + 1;
        end
        cursor = cursor + 1;
        xPlot(cursor) = x(k);
        yPlot(cursor) = y(k);
    end
end
