function relIllum = illuminationRelative(pos, opts)
    import membrane_tracking.focused_ism.*

    switch opts.illuminationMode
        case 'uniform'
            relIllum = ones(size(pos, 1), 1);

        case 'gaussian'
            r2 = sum(pos.^2, 2);
            relIllum = exp(-2 * r2 / opts.laserWaistUm^2);
    end
end
