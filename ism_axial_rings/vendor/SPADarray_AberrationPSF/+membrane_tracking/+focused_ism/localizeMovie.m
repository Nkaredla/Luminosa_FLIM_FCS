function localizationTable = localizeMovie(frames, detector, opts, positionsUm, ...
        photonCountsInFov)
    import membrane_tracking.focused_ism.*

    switch opts.localizationMethod
        case 'jointISMPoissonBayes'
            localizationTable = localizeFocusedISMMovie( ...
                frames, detector, opts, positionsUm, photonCountsInFov);
    end
end
