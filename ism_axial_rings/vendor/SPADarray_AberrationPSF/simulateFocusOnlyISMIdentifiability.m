function results = simulateFocusOnlyISMIdentifiability(varargin)
%SIMULATEFOCUSONLYISMIDENTIFIABILITY Test focal ISM magnitude/sign information.
%
%   results = simulateFocusOnlyISMIdentifiability()
%
%   The test has two parts:
%     1. Compare the 23-pixel detector microimages for +a and -a of each
%        fitted Zernike mode at focus and at a +500 nm diversity plane.
%     2. Retrieve a mixed aberration using the focus-anchored estimator and
%        compare focus-only magnitude error with final signed error.
%
%   Radially symmetric modes can have indistinguishable +/- focal
%   microimages. The diversity-plane distance measures how strongly the
%   second plane breaks that ambiguity.

    opts = parseOptions(varargin{:});
    addRequiredPaths();

    if ~isempty(opts.randomSeed)
        rng(opts.randomSeed);
    end

    sim = configureSim(opts);
    planeZ = [0 opts.defocusUm];
    sim = prepareSim(sim, planeZ, opts);

    signTable = makeSignContrastTable(sim, planeZ, opts);
    [retrievalTable, trials] = runMixedRetrievalTrials(sim, planeZ, opts);

    results = struct();
    results.options = opts;
    results.sim = sim;
    results.planeZ = planeZ;
    results.signContrastTable = signTable;
    results.retrievalTable = retrievalTable;
    results.trials = trials;
    results.outputDir = resolveOutputDir(opts);

    if opts.writeOutputs
        writeOutputs(results);
    end
    if opts.verbose
        printSummary(results);
    end
end

function opts = parseOptions(varargin)
    p = inputParser;
    p.FunctionName = 'simulateFocusOnlyISMIdentifiability';

    addParameter(p, 'defocusUm', 0.5);
    addParameter(p, 'signAmplitudeWaves', 0.04);
    addParameter(p, 'fitModes', ...
        {'defocus','astig_x','astig_y','coma_x','coma_y','spherical'});
    addParameter(p, 'truthCoeffs', defaultMixedTruth());
    addParameter(p, 'fitXY', true);
    addParameter(p, 'trueXY', [0 0]);
    addParameter(p, 'trueZOffset', 0);

    addParameter(p, 'photonsPerPlane', 5e4);
    addParameter(p, 'addNoise', true);
    addParameter(p, 'nNoiseReplicates', 3);
    addParameter(p, 'randomSeed', 1);

    addParameter(p, 'centerNormalization', 'perPlane');
    addParameter(p, 'centerSampleMode', 'subpixel');
    addParameter(p, 'modelSampleXY', [0 0]);
    addParameter(p, 'focusWeight', 1);
    addParameter(p, 'diversityWeight', 0.20);
    addParameter(p, 'focusSeedAmplitude', 0.05);
    addParameter(p, 'signSelectionModes', ...
        {'defocus','astig_x','astig_y','spherical'});
    addParameter(p, 'maxIter', 8);

    addParameter(p, 'sim', []);
    addParameter(p, 'modelDz', 0.05);
    addParameter(p, 'modelZPadding', 0.50);
    addParameter(p, 'outputDir', '');
    addParameter(p, 'writeOutputs', true);
    addParameter(p, 'verbose', true);

    parse(p, varargin{:});
    opts = p.Results;
    if ischar(opts.fitModes) || isstring(opts.fitModes)
        opts.fitModes = cellstr(opts.fitModes);
    end
    opts.fitModes = opts.fitModes(:).';
    if ischar(opts.signSelectionModes) || isstring(opts.signSelectionModes)
        opts.signSelectionModes = cellstr(opts.signSelectionModes);
    end
    opts.signSelectionModes = opts.signSelectionModes(:).';
    opts.centerNormalization = lower(char(opts.centerNormalization));
    opts.centerSampleMode = lower(char(opts.centerSampleMode));
    opts.trueXY = double(opts.trueXY(:)).';
    opts.modelSampleXY = double(opts.modelSampleXY(:)).';
    opts.nNoiseReplicates = max(1, round(opts.nNoiseReplicates));

    if ~isscalar(opts.defocusUm) || ~isfinite(opts.defocusUm) || opts.defocusUm <= 0
        error('simulateFocusOnlyISMIdentifiability:BadDefocus', ...
            'defocusUm must be a positive finite scalar.');
    end
    if numel(opts.trueXY) ~= 2 || numel(opts.modelSampleXY) ~= 2
        error('simulateFocusOnlyISMIdentifiability:BadXY', ...
            'trueXY and modelSampleXY must be [x y].');
    end
end

function truth = defaultMixedTruth()
    truth = struct();
    truth.defocus = 0.030;
    truth.astig_x = 0.025;
    truth.astig_y = -0.020;
    truth.coma_x = 0.018;
    truth.coma_y = -0.015;
    truth.spherical = 0.028;
end

function addRequiredPaths()
    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    parentDir = fileparts(thisDir);
    if exist(parentDir, 'dir') == 7
        addpath(parentDir);
    end
end

function sim = configureSim(opts)
    if isempty(opts.sim)
        sim = defaultParams();
    else
        sim = opts.sim;
    end
end

function sim = prepareSim(sim, planeZ, opts)
    zMin = min([planeZ(:); 0]) - opts.modelZPadding;
    zMax = max([planeZ(:); 0]) + opts.modelZPadding;
    nZ = ceil((zMax - zMin) / opts.modelDz) + 1;
    if mod(nZ, 2) == 0
        nZ = nZ + 1;
    end
    sim.z = linspace(zMin, zMax, nZ);
    sim.nz = numel(sim.z);
    sim.nzRange = zMax - zMin;
    sim.obj = beadObject3D(sim);
end

function T = makeSignContrastTable(sim, planeZ, opts)
    n = numel(opts.fitModes);
    mode = opts.fitModes(:);
    amplitudeWaves = repmat(opts.signAmplitudeWaves, n, 1);
    focusPlusMinusDistance = zeros(n,1);
    diversityPlusMinusDistance = zeros(n,1);
    twoPlanePlusMinusDistance = zeros(n,1);
    focusToDiversityRatio = zeros(n,1);

    for k = 1:n
        plus = struct();
        minus = struct();
        plus.(opts.fitModes{k}) = opts.signAmplitudeWaves;
        minus.(opts.fitModes{k}) = -opts.signAmplitudeWaves;
        plusValues = simulateCenterValues(sim, plus, planeZ, opts);
        minusValues = simulateCenterValues(sim, minus, planeZ, opts);

        focusPlusMinusDistance(k) = norm(plusValues(:,1) - minusValues(:,1));
        diversityPlusMinusDistance(k) = norm(plusValues(:,2) - minusValues(:,2));
        twoPlanePlusMinusDistance(k) = norm(plusValues(:) - minusValues(:));
        focusToDiversityRatio(k) = focusPlusMinusDistance(k) / ...
            max(diversityPlusMinusDistance(k), eps);
    end

    focusSignAmbiguous = focusPlusMinusDistance <= ...
        max(1e-10, 1e-3 * diversityPlusMinusDistance);
    T = table(mode, amplitudeWaves, focusPlusMinusDistance, ...
        diversityPlusMinusDistance, twoPlanePlusMinusDistance, ...
        focusToDiversityRatio, focusSignAmbiguous);
end

function [T, trials] = runMixedRetrievalTrials(sim, planeZ, opts)
    truthVec = coeffVector(sim, opts.truthCoeffs);
    fitMask = ismember(sim.modeOrder, opts.fitModes);
    symmetricMask = ismember(sim.modeOrder, opts.signSelectionModes) & ...
        abs(truthVec) > 1e-12;
    clean = simulateCenterValues(sim, opts.truthCoeffs, planeZ, opts);

    n = opts.nNoiseReplicates;
    trials = repmat(struct(), n, 1);
    replicate = (1:n).';
    focusMagnitudeRmseNm = nan(n,1);
    focusSignedRmseNm = nan(n,1);
    finalSignedRmseNm = nan(n,1);
    focusRank = nan(n,1);
    focusConditionNumber = nan(n,1);
    finalRank = nan(n,1);
    finalConditionNumber = nan(n,1);
    symmetricSignsCorrect = false(n,1);
    selectedSignLabel = cell(n,1);
    signResidualRatio = nan(n,1);

    for k = 1:n
        observed = addNoise(clean, opts);
        args = retrievalArgs(opts, sim, planeZ);
        fitResult = estimateCenterPointISMWavefront(observed, [], args{:});
        focusVec = fitResult.focusAnchored.focusFit.estCoeffVector(:).';
        finalVec = fitResult.fit.estCoeffVector(:).';

        focusMagnitudeRmseNm(k) = rmsError(abs(focusVec(fitMask)), ...
            abs(truthVec(fitMask))) * sim.lamRef * 1000;
        focusSignedRmseNm(k) = rmsError(focusVec(fitMask), ...
            truthVec(fitMask)) * sim.lamRef * 1000;
        finalSignedRmseNm(k) = rmsError(finalVec(fitMask), ...
            truthVec(fitMask)) * sim.lamRef * 1000;
        focusRank(k) = fitResult.focusAnchored.focusSufficiency.rank;
        focusConditionNumber(k) = ...
            fitResult.focusAnchored.focusSufficiency.conditionNumber;
        finalRank(k) = fitResult.sufficiency.rank;
        finalConditionNumber(k) = fitResult.sufficiency.conditionNumber;
        symmetricSignsCorrect(k) = all(sign(finalVec(symmetricMask)) == ...
            sign(truthVec(symmetricMask)));
        selectedSignLabel{k} = fitResult.focusAnchored.selectedSignLabel;
        signResidualRatio(k) = fitResult.focusAnchored.signResidualRatio;

        trials(k).observed = observed;
        trials(k).fitResult = fitResult;
    end

    T = table(replicate, focusMagnitudeRmseNm, focusSignedRmseNm, ...
        finalSignedRmseNm, symmetricSignsCorrect, selectedSignLabel, ...
        signResidualRatio, focusRank, focusConditionNumber, ...
        finalRank, finalConditionNumber);
end

function values = simulateCenterValues(sim, coeffs, planeZ, opts)
    stack = normalizedStackExplicitDetectorZPlanes(sim, coeffs, planeZ, ...
        opts.trueXY(1), opts.trueXY(2), opts.trueZOffset);
    values = sampleModelStack(stack, sim, opts.modelSampleXY, ...
        opts.centerSampleMode, numel(planeZ));
    values = normalizeValues(values, opts.centerNormalization);
end

function values = sampleModelStack(stack, sim, sampleXY, mode, nPlane)
    x = min(max(sampleXY(1), min(sim.x(:))), max(sim.x(:)));
    y = min(max(sampleXY(2), min(sim.y(:))), max(sim.y(:)));
    switch mode
        case {'subpixel','linear','bilinear','interp'}
            values = zeros(sim.nDet, nPlane);
            for ip = 1:nPlane
                for ch = 1:sim.nDet
                    values(ch,ip) = interp2(sim.x, sim.y, ...
                        double(stack(:,:,ch,ip)), x, y, 'linear', 0);
                end
            end
        case {'nearest','round','pixel'}
            [~, ix] = min(abs(sim.x - x));
            [~, iy] = min(abs(sim.y - y));
            values = squeeze(stack(iy,ix,:,:));
            values = reshape(double(values), sim.nDet, nPlane);
        otherwise
            error('simulateFocusOnlyISMIdentifiability:BadSampleMode', ...
                'centerSampleMode must be subpixel or nearest.');
    end
end

function values = normalizeValues(values, mode)
    values = max(double(values), 0);
    switch mode
        case {'perplane','plane','eachplane'}
            scale = sum(values, 1);
            scale(scale <= 0) = 1;
            values = values ./ reshape(scale, 1, []);
        case {'global','all'}
            scale = sum(values(:));
            if scale <= 0
                scale = 1;
            end
            values = values / scale;
        case {'none','raw'}
        otherwise
            error('simulateFocusOnlyISMIdentifiability:BadNormalization', ...
                'Unknown centerNormalization mode.');
    end
end

function observed = addNoise(clean, opts)
    expected = clean .* opts.photonsPerPlane;
    if opts.addNoise
        observed = poissonSample(expected);
    else
        observed = expected;
    end
end

function args = retrievalArgs(opts, sim, planeZ)
    args = { ...
        'planeZ', planeZ, ...
        'sim', sim, ...
        'fitModes', opts.fitModes, ...
        'fitStrategy', 'focusAnchored', ...
        'fitXY', opts.fitXY, ...
        'fitZ', false, ...
        'focusWeight', opts.focusWeight, ...
        'diversityWeight', opts.diversityWeight, ...
        'focusSeedAmplitude', opts.focusSeedAmplitude, ...
        'signSelectionModes', opts.signSelectionModes, ...
        'maxIter', opts.maxIter, ...
        'centerNormalization', opts.centerNormalization, ...
        'centerSampleMode', opts.centerSampleMode, ...
        'modelSampleXY', opts.modelSampleXY, ...
        'darkMode', 'none', ...
        'subtractBoundary', false, ...
        'writeOutputs', false, ...
        'verbose', false};
end

function vec = coeffVector(sim, coeffs)
    vec = zeros(1, numel(sim.modeOrder));
    for k = 1:numel(sim.modeOrder)
        if isfield(coeffs, sim.modeOrder{k})
            vec(k) = coeffs.(sim.modeOrder{k});
        end
    end
end

function value = rmsError(a, b)
    d = double(a) - double(b);
    value = sqrt(mean(d.^2));
end

function outDir = resolveOutputDir(opts)
    if ~isempty(opts.outputDir)
        outDir = char(opts.outputDir);
        return;
    end
    thisDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(thisDir);
    outDir = fullfile(rootDir, 'output_matlab', ...
        'focus_only_ism_identifiability');
end

function writeOutputs(results)
    outDir = results.outputDir;
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    writetable(results.signContrastTable, ...
        fullfile(outDir, 'focus_vs_500nm_sign_contrast.csv'));
    writetable(results.retrievalTable, ...
        fullfile(outDir, 'focus_anchored_retrieval_trials.csv'));
    save(fullfile(outDir, 'focus_only_ism_identifiability.mat'), ...
        'results', '-v7.3');
    writeFigure(results, fullfile(outDir, ...
        'focus_only_ism_identifiability_summary.png'));
end

function writeFigure(results, outFile)
    S = results.signContrastTable;
    R = results.retrievalTable;
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [80 80 1150 620]);
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', ...
        'TileSpacing', 'compact');

    ax = nexttile(tl);
    x = 1:height(S);
    semilogy(ax, x, max(S.focusPlusMinusDistance, eps), '-o', ...
        x, max(S.diversityPlusMinusDistance, eps), '-s', ...
        'LineWidth', 1.2);
    set(ax, 'XTick', x, 'XTickLabel', S.mode, ...
        'XTickLabelRotation', 35);
    ylabel(ax, 'normalized +/- signature distance');
    title(ax, 'Sign information by plane');
    legend(ax, {'focus','+500 nm'}, 'Location', 'best');
    grid(ax, 'on');

    ax = nexttile(tl);
    bar(ax, [R.focusMagnitudeRmseNm R.finalSignedRmseNm]);
    xlabel(ax, 'noise replicate');
    ylabel(ax, 'coefficient RMSE (nm RMS)');
    title(ax, 'Magnitude from focus, signed result from two planes');
    legend(ax, {'focus magnitude','final signed'}, 'Location', 'best');
    grid(ax, 'on');

    try
        exportgraphics(fig, outFile, 'Resolution', 180);
    catch
        print(fig, outFile, '-dpng', '-r180');
    end
    close(fig);
end

function printSummary(results)
    S = results.signContrastTable;
    R = results.retrievalTable;
    fprintf('\nFocus-only ISM identifiability simulation\n');
    fprintf('  diversity plane: %.0f nm\n', 1000*results.planeZ(2));
    for k = 1:height(S)
        fprintf('  %-10s +/- distance: focus %.3g, diversity %.3g%s\n', ...
            S.mode{k}, S.focusPlusMinusDistance(k), ...
            S.diversityPlusMinusDistance(k), ...
            ternary(S.focusSignAmbiguous(k), ' [focus sign ambiguous]', ''));
    end
    fprintf('  mean focus magnitude RMSE: %.3g nm\n', ...
        mean(R.focusMagnitudeRmseNm, 'omitnan'));
    fprintf('  mean final signed RMSE: %.3g nm\n', ...
        mean(R.finalSignedRmseNm, 'omitnan'));
    fprintf('  symmetric signs correct: %d/%d trials\n', ...
        nnz(R.symmetricSignsCorrect), height(R));
    fprintf('  outputs: %s\n\n', results.outputDir);
end

function out = ternary(condition, yesValue, noValue)
    if condition
        out = yesValue;
    else
        out = noValue;
    end
end
