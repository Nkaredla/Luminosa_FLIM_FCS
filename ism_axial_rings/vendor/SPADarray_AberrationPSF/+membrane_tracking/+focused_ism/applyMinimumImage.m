function delta = applyMinimumImage(delta, opts)
    import membrane_tracking.focused_ism.*

    if strcmp(opts.boundaryCondition, 'periodic')
        delta(:,1) = delta(:,1) - opts.boxSizeUm(1) * ...
            round(delta(:,1) / opts.boxSizeUm(1));
        delta(:,2) = delta(:,2) - opts.boxSizeUm(2) * ...
            round(delta(:,2) / opts.boxSizeUm(2));
    end
end
