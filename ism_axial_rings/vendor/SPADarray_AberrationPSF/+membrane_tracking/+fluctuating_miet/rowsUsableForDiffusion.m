function usable = rowsUsableForDiffusion(tbl, opts)
    import membrane_tracking.fluctuating_miet.*

    usable = true(height(tbl),1);
    if opts.minMahalanobisAcceptance > 0
        kernelValue = exp(-0.5*tbl.mahalanobisSquared);
        usable = usable & (isnan(tbl.mahalanobisSquared) | ...
            kernelValue >= opts.minMahalanobisAcceptance);
    end
end
