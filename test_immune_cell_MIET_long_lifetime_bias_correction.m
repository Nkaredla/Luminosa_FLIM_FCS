function test_immune_cell_MIET_long_lifetime_bias_correction
%TEST_IMMUNE_CELL_MIET_LONG_LIFETIME_BIAS_CORRECTION Regression test.
    tauAxis = [0.7 0.9 1.2 2.0];
    slbAxis = [250 4000];
    longAxis = [50 500];
    rowCount = numel(tauAxis) * numel(slbAxis) * numel(longAxis);
    trueTau = zeros(rowCount, 1);
    slbPhotons = zeros(rowCount, 1);
    longPhotons = zeros(rowCount, 1);
    forwardMean = zeros(rowCount, 1);
    selectionRate = ones(rowCount, 1);
    empiricalStd = 0.03 * ones(rowCount, 1);
    row = 0;
    for iSlb = 1:numel(slbAxis)
        for iLong = 1:numel(longAxis)
            for iTau = 1:numel(tauAxis)
                row = row + 1;
                trueTau(row) = tauAxis(iTau);
                slbPhotons(row) = slbAxis(iSlb);
                longPhotons(row) = longAxis(iLong);
                forwardMean(row) = tauAxis(iTau) + 0.02;
            end
        end
    end
    bad = slbPhotons == 4000 & longPhotons == 50;
    forwardMean(bad) = [1.70; 1.43; 1.31; 2.05];
    selectionRate(bad) = [0.2; 0.5; 0.9; 1.0];
    empiricalStd(bad) = 0.3;
    meanPosteriorExpectedLongPhotons = longPhotons;
    summary = table(trueTau, slbPhotons, longPhotons, ...
        meanPosteriorExpectedLongPhotons, ...
        ones(rowCount, 1), forwardMean, selectionRate, empiricalStd, ...
        0.304 * ones(rowCount, 1), ...
        'VariableNames', {'trueLongLifetimeNs','slbPhotons','longPhotons', ...
        'meanPosteriorExpectedLongPhotons', ...
        'slbPriorScale','posteriorMeanEstimateNs','displaySelectionRate', ...
        'posteriorEmpiricalStdNs','fixedSlbLifetimeNs'});
    analysis = struct('summary', summary, ...
        'instrument', struct('tauSlbNs', 0.304, 'dtNs', 0.16, ...
            'irf', [0.7; 0.2; 0.1]), ...
        'estimator', struct('lifetimeGridNs', [0.624 1 2 5]), ...
        'config', struct('slbPriorRelativeStd', 0.2, ...
            'slbCountPriorNodes', 5));

    good = immune_cell_MIET_correct_long_lifetime_bias(1.22, 250, 500, analysis);
    assert(good.reasonCode == 1);
    assert(abs(double(good.correctedLifetimeNs) - 1.2) < 1e-6);
    badCorrection = immune_cell_MIET_correct_long_lifetime_bias( ...
        1.70, 4000, 50, analysis);
    assert(badCorrection.reasonCode == 5);
    assert(isnan(badCorrection.correctedLifetimeNs));

    testDir = tempname;
    mkdir(testDir);
    calibrationFile = fullfile(testDir, 'calibration.mat');
    save(calibrationFile, 'analysis');
    result = mockResult(analysis);
    resultFile = fullfile(testDir, 'result.mat');
    save(resultFile, 'result');
    outputFile = fullfile(testDir, 'correction.mat');
    cfg = struct('writeFigures', false, 'outputDir', testDir, ...
        'outputMatFile', outputFile);
    correction = immune_cell_MIET_apply_long_lifetime_bias_correction( ...
        resultFile, calibrationFile, cfg);
    assert(isfile(outputFile));
    assert(correction.stages.native.available);
    assert(correction.stages.native.inverse.summary.inputDisplayPixelCount == 4);
    assert(correction.stages.native.inverse.summary.reliableDisplayPixelCount == 3);
    fprintf('test_immune_cell_MIET_long_lifetime_bias_correction: PASS\n');
end

function result = mockResult(analysis)
    imageSize = [2 2];
    raw = single([1.22 1.70; 2.02 0.92]);
    slb = single([250 4000; 250 250]);
    long = single([500 50; 500 500]);
    layers = struct('imageSize', imageSize, ...
        'ungated', struct('secondLifetimeNs', raw), ...
        'expectedPhotonCount', struct('secondComponent', long), ...
        'masks', struct('secondDisplay', true(imageSize)));
    prior = struct('countMarginalised', true, ...
        'targetPhotonCount', slb(:), 'photonCountStd', 0.2 * slb(:), ...
        'priorNodes', 5);
    compact = struct('fixedSlbCountPrior', prior, ...
        'membraneTauGridNs', analysis.estimator.lifetimeGridNs);
    result = struct();
    result.config = struct('maxMembraneStates', 1);
    result.slbReference = struct('fixedLifetimeNs', 0.304, ...
        'amplitudeConstraint', struct( ...
            'priorMeanPhotonCountMapNative', slb));
    result.channel = struct('dtNs', 0.16);
    result.irf = struct('curve', analysis.instrument.irf);
    result.bayesian = struct('compact', compact, ...
        'orderedComponents', layers, 'membranePixelIndices', (1:4).');
    result.spatialBinning = struct('status', 'disabled');
    result.spatialBinning4x4 = struct('status', 'disabled');
end
