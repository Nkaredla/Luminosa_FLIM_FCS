function out = immune_cell_MIET_biexp_backfill(mapsMat)
%IMMUNE_CELL_MIET_BIEXP_BACKFILL Recompute the derived maps in a saved fit.
%
% out = immune_cell_MIET_biexp_backfill(pathToBiexpSlbMapsMat)
%
% The photon and species fractions and the photon-weighted mean lifetime were
% originally computed with an extra factor of tau. That was wrong because
% biexp_slb_pattern.m normalises every basis column to UNIT SUM, which makes a
% fitted amplitude equal to that component's PHOTON COUNT already - confirmed
% numerically: a decay synthesised from 200 and 400 component photons returns
% amplitudes 200.0000 and 400.0000. Weighting those by tau again reported a
% photon share of 0.921 where the truth was 0.667, and a mean lifetime of
% 1.94 ns where the truth was 1.50 ns.
%
% Nothing needs refitting. tau1, tau2, the amplitudes, the background and the
% deviance are all correct and already stored, so this rewrites only the
% quantities derived from them. The original file is left in place as
% biexp_slb_maps_preFractionFix.mat the first time this runs, so the old numbers
% remain auditable rather than being silently replaced.
%
% Also fills in tau1AtGridEdge, which the fit did not record before the grid-edge
% diagnostic existed. It is reconstructed by comparing the fitted tau1 against
% the ends of the saved tau1Grid, which is exact: tau1 can only ever take a grid
% value.

    if ~isfile(mapsMat)
        error('immune_cell_MIET_biexp_backfill:NoFile', ...
            'No such file: %s', mapsMat);
    end
    loaded = load(mapsMat, 'out');
    out = loaded.out;

    backup = strrep(mapsMat, '.mat', '_preFractionFix.mat');
    if ~isfile(backup)
        copyfile(mapsMat, backup);
    end

    idx = out.pixelIndex;
    a1 = out.amplitude1(:);
    a2 = out.amplitude2(:);
    t1 = out.tau1Ns(:);
    t2 = out.tau2Ns(:);
    imageSize = out.imageSize;

    % ---- photon fractions: the amplitudes ARE photon counts ------------
    photonTotal = a1 + a2;
    out.maps.photonFraction1 = nan(imageSize);
    out.maps.photonFraction1(idx) = a1 ./ max(photonTotal, eps);
    out.maps.photonFraction2 = nan(imageSize);
    out.maps.photonFraction2(idx) = a2 ./ max(photonTotal, eps);

    % ---- species fractions: divide by tau for a unit-sum basis ---------
    s1 = a1 ./ max(t1, eps);
    s2 = a2 ./ max(t2, eps);
    speciesTotal = s1 + s2;
    out.maps.speciesFraction1 = nan(imageSize);
    out.maps.speciesFraction1(idx) = s1 ./ max(speciesTotal, eps);
    out.maps.speciesFraction2 = nan(imageSize);
    out.maps.speciesFraction2(idx) = s2 ./ max(speciesTotal, eps);

    % ---- photon-weighted mean lifetime ---------------------------------
    out.maps.tauMeanNs = nan(imageSize);
    out.maps.tauMeanNs(idx) = ...
        (a1 .* t1 + a2 .* t2) ./ max(photonTotal, eps);

    % ---- grid-edge flag, reconstructed exactly from the saved grid -----
    grid1 = out.tau1Grid(:);
    tol = 1e-9 * max(1, max(grid1));
    out.maps.tau1AtGridEdge = false(imageSize);
    out.maps.tau1AtGridEdge(idx) = ...
        abs(t1 - min(grid1)) < tol | abs(t1 - max(grid1)) < tol;

    out.fractionFixApplied = true;
    save(mapsMat, 'out', '-v7.3');

    valid = isfinite(out.maps.tau1Ns);
    fprintf(['  %s\n    photon share tau2 %.4f | species share tau2 %.4f | ' ...
        'mean tau %.4f ns | at edge %.1f%%\n'], mapsMat, ...
        median(out.maps.photonFraction2(valid)), ...
        median(out.maps.speciesFraction2(valid)), ...
        median(out.maps.tauMeanNs(valid)), ...
        100 * mean(out.maps.tau1AtGridEdge(valid)));
end
