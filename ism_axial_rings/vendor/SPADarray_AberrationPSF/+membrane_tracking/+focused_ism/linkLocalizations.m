function trackTable = linkLocalizations(localizationTable, opts)
    import membrane_tracking.focused_ism.*

    switch opts.trackingMethod
        case 'jpda'
            trackTable = linkLocalizationsJpda(localizationTable, opts);
        case 'nearestNeighbor'
            trackTable = linkLocalizationsNearestNeighbor( ...
                localizationTable, opts);
    end
end
