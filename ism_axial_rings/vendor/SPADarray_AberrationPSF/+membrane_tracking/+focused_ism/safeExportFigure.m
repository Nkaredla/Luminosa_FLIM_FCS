function safeExportFigure(fig, outFile, resolution)
    import membrane_tracking.focused_ism.*

    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, outFile, 'Resolution', resolution);
    else
        print(fig, outFile, '-dpng', sprintf('-r%d', resolution));
    end
end
