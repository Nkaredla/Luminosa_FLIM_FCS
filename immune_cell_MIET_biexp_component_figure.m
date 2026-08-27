function name = immune_cell_MIET_biexp_component_figure(out, outputDir)
%IMMUNE_CELL_MIET_BIEXP_COMPONENT_FIGURE Publication figure of both components.
%
% name = immune_cell_MIET_biexp_component_figure(out, outputDir)
%
% Four panels: the two component lifetimes, the photon share of the long one,
% and the photon-weighted mean lifetime. Written at 300 dpi with no axis
% furniture, so it can go into a figure as it is.
%
% COLOUR CHOICES, WHICH ARE PART OF THE RESULT
%
% Sequential quantities use viridis (lifetimes) and magma (fractions). Both are
% monotonic in lightness, so a reader printing in greyscale still sees the
% ordering, and neither invents features.
%
% Jet does invent them. Its lightness rises, falls and rises again, so wherever
% the data crosses the middle of the scale it produces a sharp cyan-to-yellow
% edge - on a lifetime map that reads as a membrane boundary which is not in the
% data. It also compresses contrast through the greens, and is unusable by
% red-green colourblind readers. It is kept available through opts only for
% matching older figures.
%
% Each panel is scaled to the 2nd-98th percentile of the DISPLAYED region, so a
% few extreme pixels cannot flatten the contrast everywhere else.
%
% WHY tau1 IS SHOWN EVEN WHEN IT IS PINNED
%
% Because a pinned parameter that has drifted is the first sign the anchor is
% wrong, and that has to be visible rather than hidden by the fact that it was
% not free.

    maps = out.maps;
    tau1Map = applyDisplayMask(maps.tau1Ns, out);
    tau2Map = applyDisplayMask(maps.tau2Ns, out);
    shareMap = applyDisplayMask(maps.photonFraction2, out);
    meanMap = applyDisplayMask(maps.tauMeanNs, out);

    h = figure('Color', 'w', 'Visible', 'off', 'Position', [40 40 1450 1180]);
    layout = tiledlayout(h, 2, 2, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    fs = 12;

    fixed = isfield(out, 'fixSlbTau') && ~isempty(out.fixSlbTau) && out.fixSlbTau;
    v1 = tau1Map(isfinite(tau1Map));
    if fixed || (~isempty(v1) && (max(v1) - min(v1)) < 1e-9)
        t1text = sprintf('\\tau_1  SLB, fixed at %.4f ns', out.opts.slbTauNs);
    else
        t1text = sprintf('\\tau_1  SLB  (median %.4f ns)', median(v1));
    end

    publicationMapPanel(nexttile(layout), tau1Map, 'viridis', t1text, 'ns', fs);
    publicationMapPanel(nexttile(layout), tau2Map, 'viridis', ...
        sprintf('\\tau_2  long component  (median %.3f ns)', ...
        median(tau2Map(isfinite(tau2Map)))), 'ns', fs);
    publicationMapPanel(nexttile(layout), shareMap, 'magma', ...
        sprintf('photon share of \\tau_2  (median %.3f)', ...
        median(shareMap(isfinite(shareMap)))), 'fraction', fs);
    publicationMapPanel(nexttile(layout), meanMap, 'viridis', ...
        sprintf('photon-weighted \\langle\\tau\\rangle  (median %.3f ns)', ...
        median(meanMap(isfinite(meanMap)))), 'ns', fs);

    binText = 'no binning';
    if out.opts.binSize > 1
        binText = sprintf('%dx%d sliding', out.opts.binSize, out.opts.binSize);
    end
    photons = median(out.amplitude1 + out.amplitude2);
    [~, acq] = fileparts(fileparts(fileparts(out.analysisMat)));
    title(layout, sprintf('%s   -   %s   -   %d pixels, median %.0f photons', ...
        acq, binText, numel(out.pixelIndex), photons), ...
        'FontSize', fs + 4, 'FontWeight', 'bold', 'Interpreter', 'none');

    name = fullfile(outputDir, 'biexp_slb_component_lifetimes.png');
    drawnow;
    exportgraphics(h, name, 'Resolution', 300);
    close(h);
end
