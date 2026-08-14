function figureHandle = showLifetimeHistograms(mitoTauNs, mitoModelOrder, ...
        lysoTauNs, lysoModelOrder, lifetimeBoundsNs, binCount, visibility)
%SHOWLIFETIMEHISTOGRAMS Compare fitted lifetime components by compartment.
% Each connected region contributes at most one value to each component.
% Component 1 is the shorter fitted lifetime because fit outputs are sorted.

    lowerBound = lifetimeBoundsNs(1);
    upperBound = lifetimeBoundsNs(2);
    binEdges = linspace(lowerBound, upperBound, max(2, round(binCount)) + 1);
    component1Colour = [0.05 0.55 0.85];
    component2Colour = [0.95 0.35 0.10];

    validMito1 = mitoModelOrder >= 1 & isfinite(mitoTauNs(:, 1));
    validMito2 = mitoModelOrder >= 2 & isfinite(mitoTauNs(:, 2));
    validLyso1 = lysoModelOrder >= 1 & isfinite(lysoTauNs(:, 1));
    validLyso2 = lysoModelOrder >= 2 & isfinite(lysoTauNs(:, 2));

    figureHandle = figure('Name', 'Organelle lifetime histograms', ...
        'Color', 'w', 'Visible', visibility);
    tiledlayout(figureHandle, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    histogram(mitoTauNs(validMito1, 1), binEdges, 'FaceColor', component1Colour, ...
        'EdgeColor', component1Colour, 'FaceAlpha', 0.55);
    hold on;
    histogram(mitoTauNs(validMito2, 2), binEdges, 'FaceColor', component2Colour, ...
        'EdgeColor', component2Colour, 'FaceAlpha', 0.55);
    xlabel('Lifetime (ns)');
    ylabel('Connected-region count');
    title(sprintf('Mitochondria: C1 n=%d, C2 n=%d', nnz(validMito1), nnz(validMito2)));
    legend('Component 1 (short)', 'Component 2 (long)', 'Location', 'best');
    grid on;
    xlim([lowerBound upperBound]);

    nexttile;
    histogram(lysoTauNs(validLyso1, 1), binEdges, 'FaceColor', component1Colour, ...
        'EdgeColor', component1Colour, 'FaceAlpha', 0.55);
    hold on;
    histogram(lysoTauNs(validLyso2, 2), binEdges, 'FaceColor', component2Colour, ...
        'EdgeColor', component2Colour, 'FaceAlpha', 0.55);
    xlabel('Lifetime (ns)');
    ylabel('Connected-region count');
    title(sprintf('Lysosomes: C1 n=%d, C2 n=%d', nnz(validLyso1), nnz(validLyso2)));
    legend('Component 1 (short)', 'Component 2 (long)', 'Location', 'best');
    grid on;
    xlim([lowerBound upperBound]);
end
