function safeExportFigure(fig, outputFile, resolution)
    import membrane_tracking.curved_miet.*

    try
        exportgraphics(fig, outputFile, 'Resolution', resolution);
    catch
        print(fig, outputFile, '-dpng', sprintf('-r%d', resolution));
    end
end
