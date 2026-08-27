function outFile = export_biexp_for_explorer(mapsMat, analysisMat)
%EXPORT_BIEXP_FOR_EXPLORER Make a biexp result loadable by the existing GUI.
%
% outFile = export_biexp_for_explorer(biexpMapsMat)
% outFile = export_biexp_for_explorer(biexpMapsMat, analysisMat)
%
% immune_cell_MIET_explorer reads an ANALYSIS mat: it takes its map list from
% fields of `result` and pools per-pixel decays from the `tcspc_pix` cube in the
% same file. The biexponential results live in a separate file with a different
% layout, so the GUI cannot see them.
%
% Rather than fork the GUI, this writes a small companion file in the GUI's own
% shape: a `result` struct carrying the biexp maps under the paths the map
% catalogue already understands, and a reference to the original cube so the
% decay panel still has photons to pool. The GUI is then pointed at this file
% and needs no change.
%
% The cube is REFERENCED, not copied. It is 113 MB per acquisition and this
% machine runs with C: chronically full, so duplicating it per binning would
% cost gigabytes for nothing; the explorer opens it through matfile anyway.

    if nargin < 2 || isempty(analysisMat)
        L = load(mapsMat, 'out');
        analysisMat = L.out.analysisMat;
    else
        L = load(mapsMat, 'out');
    end
    out = L.out;
    if ~isfile(analysisMat)
        error('export_biexp_for_explorer:NoAnalysis', ...
            'The source analysis MAT is missing: %s', analysisMat);
    end

    m = out.maps;
    result = struct();
    % Paths chosen to match what immune_cell_MIET_explorer_maps already looks
    % for, so the dropdown populates without touching the catalogue.
    result.bayesian.maps.tauMeanArithmetic = m.tauMeanNs;
    result.bayesian.maps.membraneFraction = m.photonFraction2;
    result.biexp.tau1Ns = m.tau1Ns;
    result.biexp.tau2Ns = m.tau2Ns;
    result.biexp.photonFraction2 = m.photonFraction2;
    result.biexp.speciesFraction2 = m.speciesFraction2;
    result.biexp.tauMeanNs = m.tauMeanNs;
    result.biexp.reducedDeviance = m.reducedDeviance;
    if isfield(m, 'slbPhotons'); result.biexp.slbPhotons = m.slbPhotons; end
    if isfield(m, 'longPhotons'); result.biexp.longPhotons = m.longPhotons; end
    if isfield(m, 'residualAcf1')
        result.biexp.residualAcf1 = m.residualAcf1;
    end
    result.redMeanFlim.intensity = m.intensity;

    % Timing and IRF the decay panel needs.
    result.bayesian.compact.dtNs = out.dtNs;
    result.bayesian.compact.pulsePeriodNs = out.periodNs;
    result.irf.curve = out.irf;
    result.bayesian.compact.fixedSlbLifetimeNs = out.opts.slbTauNs;

    result.explorerSource = struct('analysisMat', analysisMat, ...
        'biexpMaps', mapsMat, 'binSize', out.opts.binSize, ...
        'method', out.method);

    outFile = strrep(mapsMat, '.mat', '_explorer.mat');
    save(outFile, 'result', '-v7.3');
    fprintf('  wrote %s\n', outFile);
    fprintf(['    open with: immune_cell_MIET_explorer(''%s'')\n' ...
        '    (the TCSPC cube is read from %s)\n'], outFile, analysisMat);
end
