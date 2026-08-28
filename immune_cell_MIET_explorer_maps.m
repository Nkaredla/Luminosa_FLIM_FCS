function catalogue = immune_cell_MIET_explorer_maps(result)
%IMMUNE_CELL_MIET_EXPLORER_MAPS Catalogue of displayable maps in an analysis.
%
% catalogue = immune_cell_MIET_explorer_maps(result)
%
% Returns a struct array, one entry per map that is actually present in this
% result, with fields:
%   label    text for the dropdown, grouped by stage
%   data     602x602 (or whatever the image size is) numeric map
%   units    axis/colourbar label
%   limits   [lo hi] display range, or [] to autoscale
%   isTau    true for lifetime-valued maps, which get a different colormap
%
% Maps are collected from every stage the pipeline may have written - native,
% 2x2 sliding and 4x4 sliding - and entries whose source is absent are skipped
% rather than erroring, because which stages exist depends on the run.
%
% The pooled TCSPC shown alongside always comes from the NATIVE per-pixel cube,
% whichever map is on screen. That is deliberate: the binned maps are anchored
% sliding sums of the native cube, so pooling the native cube over the region
% the user selects reproduces them exactly and there is only one code path that
% can be wrong. The binned maps are for navigating to a feature, not a separate
% source of photons.

    if ~isstruct(result)
        error('immune_cell_MIET_explorer_maps:Result', ...
            'result must be the struct saved by immune_cell_MIET.');
    end

    entries = {
        % label,                                    path,                                                  units,            isTau
        'red mean FLIM - arrival time',             'redMeanFlim.meanArrivalNs',                            'ns',              true
        'red mean FLIM - intensity',                'redMeanFlim.intensity',                                'photons',         false
        'blue mean FLIM - arrival time',            'blueMeanFlim.meanArrivalNs',                           'ns',              true
        'blue mean FLIM - intensity',               'blueMeanFlim.intensity',                               'photons',         false
        'reassigned - mean arrival',                'reassigned.total.meanArrival',                         'ns',              true
        'reassigned - tau mean',                    'reassigned.total.tauMean',                             'ns',              true
        'ISM - ACO average',                        'ism.acoAverage',                                       'photons',         false
        'ISM - deconvolved',                        'ism.deconvolvedImage',                                 'photons',         false
        'native fit - tau mean',                    'bayesian.maps.tauMeanArithmetic',                      'ns',              true
        'native fit - tau posterior sd',            'bayesian.maps.tauPosteriorStd',                        'ns',              true
        'native fit - P(biexponential)',            'bayesian.maps.probabilityBiexponential',               'probability',     false
        'native fit - P(triexponential)',           'bayesian.maps.probabilityTriexponential',              'probability',     false
        'native fit - P(SLB only)',                 'bayesian.maps.probabilityFixedSlbOnly',                'probability',     false
        'native fit - membrane fraction',           'bayesian.maps.membraneFraction',                       'fraction',        false
        'native fit - signal fraction',             'bayesian.maps.signalFraction',                         'fraction',        false
        'native fit - log10 BF M3 vs M2',           'bayesian.maps.logBayesFactorM3VsM2',                   'log Bayes factor', false
        '2x2 mean FLIM - arrival time',             'spatialBinning.meanFlim.meanArrivalNs',                'ns',              true
        '2x2 mean FLIM - intensity',                'spatialBinning.meanFlim.intensity',                    'photons',         false
        '2x2 fit - tau mean',                       'spatialBinning.bayesian.maps.tauMeanArithmetic',       'ns',              true
        '2x2 fit - P(triexponential)',              'spatialBinning.bayesian.maps.probabilityTriexponential', 'probability',   false
        '4x4 mean FLIM - arrival time',             'spatialBinning4x4.meanFlim.meanArrivalNs',             'ns',              true
        '4x4 mean FLIM - intensity',                'spatialBinning4x4.meanFlim.intensity',                 'photons',         false
        '4x4 fit - tau mean',                       'spatialBinning4x4.bayesian.maps.tauMeanArithmetic',    'ns',              true
        '4x4 fit - P(triexponential)',              'spatialBinning4x4.bayesian.maps.probabilityTriexponential', 'probability', false
        'segmentation - label map',                 'segmentation.masks.labelMap',                          'label',           false
        % ---- soft-SLB biexponential results -----------------------------
        % Present only in the companion file an anchored biexp run writes
        % (biexp_slb_maps_explorer.mat). Entries whose source is absent are
        % skipped, so listing them here costs nothing on an ordinary
        % analysis MAT and means the biexp maps are browsable with the same
        % tool instead of only as PNGs.
        'biexp - tau1 (SLB)',                       'biexp.tau1Ns',                                          'ns',              true
        'biexp - tau2 (long component)',            'biexp.tau2Ns',                                          'ns',              true
        'biexp - photon-weighted mean tau',         'biexp.tauMeanNs',                                       'ns',              true
        'biexp - photon share of tau2',             'biexp.photonFraction2',                                 'fraction',        false
        'biexp - species share of tau2',            'biexp.speciesFraction2',                                'fraction',        false
        'biexp - SLB photons per pixel',            'biexp.slbPhotons',                                      'photons',         false
        'biexp - long-component photons',           'biexp.longPhotons',                                     'photons',         false
        'biexp - reduced deviance',                 'biexp.reducedDeviance',                                 'deviance',        false
        'biexp - residual autocorrelation',         'biexp.residualAcf1',                                    'correlation',     false
        };

    catalogue = struct('label', {}, 'data', {}, 'units', {}, ...
        'limits', {}, 'isTau', {}, 'path', {});
    for k = 1:size(entries, 1)
        value = immune_cell_MIET_explorer_field(result, entries{k, 2});
        if isempty(value) || ~isnumeric(value) || ~ismatrix(value) || ...
                any(size(value) < 2)
            continue;
        end
        catalogue(end + 1) = struct( ... %#ok<AGROW>
            'label', entries{k, 1}, 'data', double(value), ...
            'units', entries{k, 3}, 'limits', [], ...
            'isTau', logical(entries{k, 4}), 'path', entries{k, 2});
    end

    % Per-component lifetime and photon maps are 3-D (component, y, x); expose
    % each component separately, since the sorted components are the whole
    % point of the three-model fit.
    stages = {'', 'spatialBinning.', 'spatialBinning4x4.'};
    stageNames = {'native', '2x2', '4x4'};
    for s = 1:numel(stages)
        stack = immune_cell_MIET_explorer_field(result, ...
            [stages{s} 'bayesian.maps.componentLifetimeNs']);
        shares = immune_cell_MIET_explorer_field(result, ...
            [stages{s} 'bayesian.maps.componentPhotonFraction']);
        if isempty(stack) || ndims(stack) ~= 3; continue; end
        % Components are the THIRD dimension in MATLAB: [y, x, component].
        % An h5py listing of the same array shows (3, 602, 602) because HDF5
        % reports dimensions in the reverse order, and taking that at face
        % value made this loop run 602 times per stage.
        for c = 1:size(stack, 3)
            catalogue(end + 1) = struct( ... %#ok<AGROW>
                'label', sprintf('%s fit - component %d lifetime', ...
                stageNames{s}, c), ...
                'data', double(stack(:, :, c)), 'units', 'ns', ...
                'limits', [], 'isTau', true, ...
                'path', sprintf('%sbayesian.maps.componentLifetimeNs(%d)', ...
                stages{s}, c));
            if ~isempty(shares) && ndims(shares) == 3 && size(shares, 3) >= c
                catalogue(end + 1) = struct( ... %#ok<AGROW>
                    'label', sprintf('%s fit - component %d photon share', ...
                    stageNames{s}, c), ...
                    'data', double(shares(:, :, c)), ...
                    'units', 'fraction', 'limits', [0 1], 'isTau', false, ...
                    'path', sprintf(['%sbayesian.maps.' ...
                    'componentPhotonFraction(%d)'], stages{s}, c));
            end
        end
    end

    if isempty(catalogue)
        error('immune_cell_MIET_explorer_maps:NoMaps', ...
            'No displayable 2-D maps were found in this result.');
    end
end
