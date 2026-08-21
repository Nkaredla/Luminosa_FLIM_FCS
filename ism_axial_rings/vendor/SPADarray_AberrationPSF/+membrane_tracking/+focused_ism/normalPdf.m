function y = normalPdf(x)
    import membrane_tracking.focused_ism.*

    y = exp(-0.5 * x.^2) / sqrt(2*pi);
end
