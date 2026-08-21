function [prob, dPdx, dPdy] = gaussianPixelModel(x0, y0, sigma, xEdges, yEdges)
    import membrane_tracking.focused_ism.*

    xLeft = (xEdges(1:end-1) - x0) / sigma;
    xRight = (xEdges(2:end) - x0) / sigma;
    yBottom = (yEdges(1:end-1) - y0) / sigma;
    yTop = (yEdges(2:end) - y0) / sigma;

    px = normalCdf(xRight) - normalCdf(xLeft);
    py = normalCdf(yTop) - normalCdf(yBottom);
    dpx = (normalPdf(xLeft) - normalPdf(xRight)) / sigma;
    dpy = (normalPdf(yBottom) - normalPdf(yTop)) / sigma;

    px = max(px, 0);
    py = max(py, 0);
    prob = py(:) * px(:).';

    if nargout > 1
        dPdx = py(:) * dpx(:).';
        dPdy = dpy(:) * px(:).';
    end
end
