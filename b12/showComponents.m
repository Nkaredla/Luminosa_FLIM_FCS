function figureHandle = showComponents(componentLifetime, componentAmplitude, ...
        componentFraction, backgroundAmplitude, visibility)
%SHOWCOMPONENTS Display both lifetime components and their photon amplitudes.

    figureHandle = figure('Name', 'B12 lifetime components and amplitudes', ...
        'Color', 'w', 'Visible', visibility);
    tiledlayout(figureHandle, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    showMap(componentLifetime(:, :, 1), 'Component 1 lifetime', 'ns', turbo(256));
    nexttile;
    showMap(componentAmplitude(:, :, 1), 'Component 1 amplitude', 'photons', parula(256));
    nexttile;
    showMap(componentFraction(:, :, 1), 'Component 1 fraction', '', parula(256));
    nexttile;
    showMap(componentLifetime(:, :, 2), 'Component 2 lifetime', 'ns', turbo(256));
    nexttile;
    showMap(componentAmplitude(:, :, 2), 'Component 2 amplitude', 'photons', parula(256));
    nexttile;
    showMap(backgroundAmplitude, 'Background amplitude', 'photons', gray(256));
end
