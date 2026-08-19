function result = immune_cell_MIET(ptuFile, cfg)
%IMMUNE_CELL_MIET Two-colour PIE ISM-FLIM analysis for immune-cell MIET.
%
% result = immune_cell_MIET(ptuFile, cfg)
%
% The pipeline keeps every ISM detector in the selected
% 640-nm PIE excitation window, performs ACO-ISM/pixel reassignment, and
% displays a red-channel mean-photon-arrival-time FLIM image. It segments
% the cell from outside SLB, fixes the IRF-fitted SLB lifetime, and compares
% one-, two-, and three-exponential Bayesian models in membrane pixels.
%
% Important cfg fields (all optional):
%   excitationPulseIndex  PIE pulse number for 640 nm (default 2)
%   excitationNm          excitation wavelength label (default 640)
%   ismWavelengthNm       emission wavelength used by ACO-ISM (default 690)
%   tcspcBinNs            coarse reader resolution in ns (default 0.16)
%   minPhotonsPerPixel    display mask threshold (default 5)
%   lifetimeRangeNs       fixed display limits; [] chooses percentiles
%   showFigure            show FLIM and Bayesian figures (default true)
%   saveOutputs           save PNG and MAT beside the PTU (default true)
%   outputDir             output directory (default <PTU>/immune_cell_MIET)
%   useGPU                use a GPU in ACO-ISM when available (default false)
%   ism                   overrides fields in the ACO-ISM parameter struct
%   segmentation          overrides immune_cell_MIET_segment options
%   bayes                 overrides flim_bayes_fixed_slb options
%   fixSlbAmplitude       fix expected SLB photons/pixel from outside fit (true)
%   slbAmplitudeScale     multiplier for that outside calibration (default 1)
%   maxSlbClippedFraction warning threshold when observed N is below the
%                         fixed SLB count (default 0.10)
%   componentMaps         posterior/photon gates for sorted tau2/tau3 layers
%   spatialBinning        optional overlapping TCSPC binning; the supported
%                         mode is 2x2 windows with one-pixel steps
%   spatialBinning4x4     the corresponding 4x4, one-pixel-step analysis

    if nargin < 2 || isempty(cfg)
        cfg = struct();
    end
    cfg = fillDefaults(cfg);
    if ~(ischar(ptuFile) || (isstring(ptuFile) && isscalar(ptuFile)))
        error('immune_cell_MIET:PTUInput', 'ptuFile must be one file name.');
    end
    ptuFile = char(ptuFile);
    if ~isfile(ptuFile)
        error('immune_cell_MIET:MissingPTU', 'PTU file was not found: %s', ptuFile);
    end
    outputDir = cfg.outputDir;
    if isempty(outputDir)
        outputDir = fullfile(fileparts(ptuFile), 'immune_cell_MIET');
    end
    if cfg.saveOutputs && ~isfolder(outputDir)
        mkdir(outputDir);
    end
    preliminaryMatFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_preliminary.mat');
    segmentationMatFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_segmentation.mat');
    analysisMatFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_analysis.mat');

    required = {'PTU_Read_Head', 'PTU_FLIM_GPU', ...
        'run_ism_reconstruction_from_ptu', ...
        'immune_cell_MIET_reassigned_moments', 'immune_cell_MIET_segment', ...
        'immune_cell_MIET_reassigned_tcspc', ...
        'immune_cell_MIET_reassigned_sliding_tcspc', 'flim_bayes_fixed_slb', ...
        'immune_cell_MIET_sorted_components'};
    for k = 1:numel(required)
        if exist(required{k}, 'file') ~= 2
            error('immune_cell_MIET:MissingDependency', ...
                'Required function %s.m is not on the MATLAB path.', required{k});
        end
    end

    helperDir = fullfile(fileparts(mfilename('fullpath')), 'b12');
    addedHelperPath = ~contains([pathsep path pathsep], [pathsep helperDir pathsep]);
    if addedHelperPath
        addpath(helperDir, '-begin');
        helperCleanup = onCleanup(@() rmpath(helperDir));
    end

    pipelineSourceFile = [mfilename('fullpath') '.m'];
    pipelineSourceInfo = dir(pipelineSourceFile);
    if isempty(pipelineSourceInfo)
        pipelineSourceStamp = 'unknown';
    else
        pipelineSourceStamp = pipelineSourceInfo(1).date;
    end
    fprintf(['immune_cell_MIET: stages enabled - native 1; 2x2 sliding %d; ' ...
        '4x4 sliding %d (source modified %s).\n'], ...
        cfg.spatialBinning.enabled, cfg.spatialBinning4x4.enabled, ...
        pipelineSourceStamp);

    fprintf('immune_cell_MIET: reading %s\n', ptuFile);
    previewHead = PTU_Read_Head(ptuFile);
    [readerGateOk, readerGateStartNative, readerGateLengthNative] = ...
        nativeReaderGate(previewHead, cfg);
    readOpts = struct('photonsPerChunk', cfg.photonsPerChunk, ...
        'computePerFrame', false, 'storeTcspcPix', false, ...
        'storePhotonLists', true, 'storeTimeCell', false, ...
        'showWaitbar', cfg.showWaitbar, 'useGPU', cfg.useGPU, ...
        'minLifetimeBin_ns', cfg.tcspcBinNs, 'maxNgate', cfg.maxNgate, ...
        'computeGlobalMaps', false, 'storePhotonSync', false, ...
        'storePhotonFrame', false, 'lowMemoryPhotonLists', true, ...
        'photonBlockSize', cfg.photonBlockSize);
    if readerGateOk
        readOpts.photonGateStartNative = readerGateStartNative;
        readOpts.photonGateLengthNative = readerGateLengthNative;
        fprintf(['immune_cell_MIET: prefiltering the PTU reader to the ' ...
            '%g-nm PIE window.\n'], cfg.excitationNm);
    end
    ptuRaw = PTU_FLIM_GPU(ptuFile, readOpts);
    assertPhotonLists(ptuRaw);

    detectorIds = activeDetectorIds(ptuRaw.im_chan, ...
        getFieldOr(ptuRaw, 'dind', []), cfg.photonBlockSize);
    detectorHistogramCounts = detectorHistogram(ptuRaw.im_tcspc, ...
        ptuRaw.im_chan, detectorIds, double(ptuRaw.Ngate), ...
        cfg.photonBlockSize);
    dtNs = resolveDtNs(ptuRaw);
    preGateBins = max(1, round(cfg.gatePreNs / dtNs));
    [headerGateOk, gateStarts, gateStops, gateLength, riseBins] = ...
        gatesFromPtuHeader(ptuRaw.head, double(ptuRaw.Ngate), ...
        cfg.piePulseCount, preGateBins);
    gateProfile = max(mean(double(detectorHistogramCounts), 2), 0);
    smoothBins = max(3, min(21, 2 * floor(numel(gateProfile) / 200) + 1));
    smoothGateProfile = conv(gateProfile, ones(smoothBins, 1) / smoothBins, 'same');
    sortedProfile = sort(smoothGateProfile);
    gateBaseline = mean(sortedProfile(1:max(1, round(0.2 * numel(sortedProfile)))));
    gateThreshold = gateBaseline + cfg.gateThresholdFrac * ...
        (max(smoothGateProfile) - gateBaseline);
    if headerGateOk
        [peakBins, peakHeights] = peaksWithinHeaderGates( ...
            smoothGateProfile, riseBins, gateLength);
        gateSource = 'PTU header PIETimeGate';
        fprintf('immune_cell_MIET: using %d PIE windows recorded in the PTU header.\n', ...
            numel(gateStarts));
    else
        try
            [gateStarts, gateStops, gateLength, gateProfile, smoothGateProfile, ...
                gateThreshold, gateBaseline, riseBins, peakBins, peakHeights] = ...
                detectGates(detectorHistogramCounts, cfg.piePulseCount, ...
                cfg.gateThresholdFrac, preGateBins, cfg.minGateSeparationBins, ...
                cfg.weakPeakFraction);
        catch gateError
            if contains(gateError.message, 'too weak for reliable gating')
                warning('immune_cell_MIET:WeakPiePulse', ...
                    ['Automatic PIE detection found a weak pulse. Retrying without ' ...
                     'the relative peak-height rejection because no header gates exist.']);
                [gateStarts, gateStops, gateLength, gateProfile, smoothGateProfile, ...
                    gateThreshold, gateBaseline, riseBins, peakBins, peakHeights] = ...
                    detectGates(detectorHistogramCounts, cfg.piePulseCount, ...
                    cfg.gateThresholdFrac, preGateBins, cfg.minGateSeparationBins, 0);
            else
                rethrow(gateError);
            end
        end
        gateSource = 'TCSPC histogram detection';
    end

    pulseIndex = cfg.excitationPulseIndex;
    if pulseIndex > numel(gateStarts)
        error('immune_cell_MIET:PulseIndex', ...
            'Requested PIE pulse %d, but only %d pulses were detected.', ...
            pulseIndex, numel(gateStarts));
    end
    recordedWavelengths = resolveRecordedWavelengths(ptuRaw.head);
    if numel(recordedWavelengths) >= pulseIndex && isfinite(recordedWavelengths(pulseIndex))
        recordedExcitationNm = recordedWavelengths(pulseIndex);
        if abs(recordedExcitationNm - cfg.excitationNm) > cfg.wavelengthToleranceNm
            error('immune_cell_MIET:WavelengthMismatch', ...
                'PIE pulse %d is recorded as %g nm, not the requested %g nm.', ...
                pulseIndex, recordedExcitationNm, cfg.excitationNm);
        end
    else
        recordedExcitationNm = NaN;
        warning('immune_cell_MIET:NoWavelengthMetadata', ...
            'No laser-wavelength metadata was available; pulse %d is assumed to be %g nm.', ...
            pulseIndex, cfg.excitationNm);
    end
    fprintf('immune_cell_MIET: using pulse %d (%g nm), bins %d:%d; retaining %d detectors.\n', ...
        pulseIndex, cfg.excitationNm, gateStarts(pulseIndex), ...
        gateStops(pulseIndex), numel(detectorIds));

    ptu640 = isolateExcitationGate(ptuRaw, gateStarts(pulseIndex), ...
        gateLength, detectorIds, cfg.photonBlockSize);
    clear ptuRaw
    ismParams = defaultIsmParams(ptu640, cfg);
    ismParams = mergeStruct(ismParams, cfg.ism);

    fprintf('immune_cell_MIET: running ACO-ISM reconstruction.\n');
    ismResult = run_ism_reconstruction_from_ptu(ptu640, ismParams);
    fprintf('immune_cell_MIET: applying detector and TCSPC reassignment.\n');
    reassignmentParams = struct('minCounts', cfg.minPhotonsPerPixel, ...
        'photonBlockSize', cfg.photonBlockSize, 'oversampleXY', 2);
    reassigned = immune_cell_MIET_reassigned_moments( ...
        ptu640, ismResult, reassignmentParams);
    [meanArrivalNs, intensity, tcspcSource] = resolveReassignedMean(reassigned, dtNs);
    meanArrivalNs(double(intensity) < cfg.minPhotonsPerPixel) = NaN;
    figureHandle = [];
    preliminaryPngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_preliminary_mean_FLIM.png');
    if cfg.showFigure || cfg.saveOutputs
        visibility = figureVisibility(cfg.showFigure);
        figureHandle = showRedMeanFlim(meanArrivalNs, intensity, ...
            cfg.lifetimeRangeNs, cfg.excitationNm, visibility);
        drawnow;
        if cfg.saveOutputs
            exportgraphics(figureHandle, preliminaryPngFile, 'Resolution', 300);
            fprintf('immune_cell_MIET: saved preliminary mean FLIM: %s\n', ...
                preliminaryPngFile);
        end
    end
    if cfg.saveOutputs
        preliminary = struct();
        preliminary.pipeline = 'immune_cell_MIET';
        preliminary.stage = 'preliminary reassigned mean FLIM';
        preliminary.sourceFile = ptuFile;
        preliminary.config = cfg;
        preliminary.channel = struct('excitationNm', cfg.excitationNm, ...
            'piePulseIndex', pulseIndex, 'gateStartBin', gateStarts(pulseIndex), ...
            'gateStopBin', gateStops(pulseIndex), 'gateLengthBins', gateLength, ...
            'detectorIds', detectorIds(:).', 'dtNs', dtNs);
        preliminary.meanArrivalNs = single(meanArrivalNs);
        preliminary.intensity = uint32(min(max(round(double(intensity)), 0), ...
            double(intmax('uint32'))));
        preliminary.ism = ismResult;
        preliminary.reassignedMoments = reassigned;
        preliminary.containsPhotonStream = false;
        save(preliminaryMatFile, 'preliminary', '-v7.3');
        clear preliminary
        fprintf('immune_cell_MIET: saved preliminary MAT without photon stream: %s\n', ...
            preliminaryMatFile);
    end

    fprintf('immune_cell_MIET: segmenting cell membrane from surrounding SLB.\n');
    segmentationCfg = cfg.segmentation;
    if ~isfield(segmentationCfg, 'minPhotonsPerPixel') || ...
            isempty(segmentationCfg.minPhotonsPerPixel)
        segmentationCfg.minPhotonsPerPixel = cfg.minPhotonsPerPixel;
    end
    [masks, segmentationDiagnostics] = immune_cell_MIET_segment( ...
        intensity, meanArrivalNs, segmentationCfg);
    segmentationPngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_segmentation_check.png');
    if cfg.showFigure || cfg.saveOutputs
        segmentationFigureHandle = showSegmentationCheck(meanArrivalNs, ...
            intensity, masks, segmentationDiagnostics, ...
            figureVisibility(cfg.showFigure));
        drawnow;
        if cfg.saveOutputs
            exportgraphics(segmentationFigureHandle, segmentationPngFile, ...
                'Resolution', 300);
            fprintf('immune_cell_MIET: saved segmentation check: %s\n', ...
                segmentationPngFile);
        end
    end
    if cfg.saveOutputs
        segmentation = struct('pipeline', 'immune_cell_MIET', ...
            'stage', 'cell-SLB segmentation', 'sourceFile', ptuFile, ...
            'masks', masks, 'diagnostics', segmentationDiagnostics, ...
            'meanArrivalNs', single(meanArrivalNs), ...
            'intensity', uint32(min(max(round(double(intensity)), 0), ...
                double(intmax('uint32')))), ...
            'containsPhotonStream', false);
        save(segmentationMatFile, 'segmentation', '-v7.3');
        clear segmentation
        fprintf('immune_cell_MIET: saved segmentation MAT without photon stream: %s\n', ...
            segmentationMatFile);
    end
    if ~strcmp(segmentationDiagnostics.status, 'ok')
        warning('immune_cell_MIET:SegmentationStatus', ...
            'Segmentation status is %s. Inspect result.segmentation before interpretation.', ...
            segmentationDiagnostics.status);
        if cfg.requireValidSegmentation
            error('immune_cell_MIET:InvalidSegmentation', ...
                ['FLIM-only segmentation status is %s. Bayesian fitting was ' ...
                 'stopped so an ' ...
                 'implausible footprint cannot be interpreted as a cell.'], ...
                segmentationDiagnostics.status);
        end
    end
    if ~any(masks.cellMembrane(:))
        error('immune_cell_MIET:NoMembrane', ...
            'No cell-membrane pixels were segmented for Bayesian analysis.');
    end
    if ~any(masks.slbReference(:))
        error('immune_cell_MIET:NoSLBReference', ...
            'No outside-SLB reference pixels were segmented.');
    end

    fprintf('immune_cell_MIET: collecting compact reassigned membrane TCSPC data.\n');
    storeFullTcspc = cfg.saveOutputs && cfg.saveTcspcPix;
    [membraneCube, membranePixelIndices, slbCurve, compactMeta, tcspc_pix] = ...
        immune_cell_MIET_reassigned_tcspc(ptu640, ismResult, ...
        masks.cellMembrane, masks.slbReference, storeFullTcspc, ...
        cfg.photonBlockSize);
    tcspcPixSaved = false;
    if storeFullTcspc
        if ~isa(tcspc_pix, 'uint16')
            error('immune_cell_MIET:TcspcStorageClass', ...
                'The reassigned tcspc_pix cube must be uint16 before saving.');
        end
        save(analysisMatFile, 'tcspc_pix', '-v7.3');
        tcspcPixSaved = true;
        fprintf(['immune_cell_MIET: saved reassigned detector-summed ' ...
            'uint16 tcspc_pix: %s\n'], analysisMatFile);
        clear tcspc_pix
    end
    if sum(slbCurve) < cfg.minSlbPhotons
        error('immune_cell_MIET:TooFewSLBPhotons', ...
            'The outside-SLB reference has %d photons; at least %d are required.', ...
            round(sum(slbCurve)), cfg.minSlbPhotons);
    end
    if abs(compactMeta.dtNs - dtNs) > max(1e-9, 1e-6 * dtNs)
        error('immune_cell_MIET:TimingMetadata', ...
            'Reassigned TCSPC bin width (%g ns) differs from reader metadata (%g ns).', ...
            compactMeta.dtNs, dtNs);
    end

    pulsePeriodNs = resolvePulsePeriodNs(ptu640);
    [irf, irfSource, irfParameters, irfTau, irfFit, irfChi] = ...
        estimateIrf(compactMeta.globalCurve, dtNs, cfg.irf, ...
        cfg.minPhotonsForIRF, cfg.useGPU, cfg.irfNCasc, cfg.irfNSub, ...
        cfg.irfTauSeedsNs);
    [slbTauFit, slbTauMean, slbCoefficients, slbFraction, slbFitCounts, ...
        slbDeviance, slbBic, slbReducedChiSquare] = fitFixedOrder( ...
        slbCurve, irf, dtNs, 1, cfg.slbTauSeedNs, ...
        cfg.membraneTauSeedsNs, cfg.lifetimeBoundsNs, ...
        cfg.includeBackground, cfg.fitMaxIter);
    tauSlbNs = slbTauFit(1);
    fprintf('immune_cell_MIET: fixed outside-SLB lifetime = %.4g ns.\n', tauSlbNs);
    slbAmplitude = estimateSlbAmplitudeConstraint(slbCurve, slbFitCounts, ...
        slbCoefficients, intensity, masks.slbReference, ...
        masks.outsideSlbGeometry, ...
        cfg.slbAmplitudeScale);
    slbAmplitude.constraintApplied = cfg.fixSlbAmplitude;
    if cfg.fixSlbAmplitude
        fprintf(['immune_cell_MIET: fixing the expected SLB amplitude to ' ...
            '%.3g fitted photons/pixel (outside-population dispersion %.3g; ' ...
            '%d outside pixels).\n'], ...
            slbAmplitude.signalPhotonsPerPixel, ...
            slbAmplitude.dispersionPhotonsPerPixel, ...
            slbAmplitude.amplitudePopulationPixelCount);
    end

    bayesOpts = cfg.bayes;
    bayesOpts.analysisMask = true(size(membraneCube, 1), 1);
    bayesOpts.minPhotons = cfg.bayesMinPhotons;
    bayesOpts.useGPU = cfg.useGPU;
    if cfg.fixSlbAmplitude
        bayesOpts.fixedSlbPhotonCount = slbAmplitude.signalPhotonsPerPixel;
        bayesOpts.fixedSlbPhotonCountStd = ...
            slbAmplitude.dispersionPhotonsPerPixel;
    else
        % The top-level switch is authoritative even if stale count fields
        % remain in a reused cfg.bayes struct.
        bayesOpts.fixedSlbPhotonCount = [];
        bayesOpts.fixedSlbPhotonCountStd = [];
    end
    fprintf(['immune_cell_MIET: running fixed-SLB Bayesian comparison ' ...
        '(SLB, biexponential, triexponential).\n']);
    bayesCompact = flim_bayes_fixed_slb(membraneCube, irf, ...
        pulsePeriodNs, dtNs, tauSlbNs, bayesOpts);
    clear membraneCube
    slbAmplitude.bayesianClippedPixelCount = 0;
    slbAmplitude.bayesianClippedPixelFraction = 0;
    bayesianStatus = 'complete';
    bayesianValidForInterpretation = true;
    if cfg.fixSlbAmplitude && ...
            isfield(bayesCompact, 'fixedSlbPhotonConstraint')
        slbAmplitude.bayesianClippedPixelCount = ...
            bayesCompact.fixedSlbPhotonConstraint.clippedPixelCount;
        slbAmplitude.bayesianClippedPixelFraction = ...
            bayesCompact.fixedSlbPhotonConstraint.clippedPixelFraction;
        if slbAmplitude.bayesianClippedPixelFraction > ...
                cfg.maxSlbClippedFraction
            bayesianStatus = 'invalid_fixed_slb_count_clipping';
            bayesianValidForInterpretation = false;
            warning('immune_cell_MIET:FixedSlbCountClipping', ...
                ['The outside-derived SLB count exceeded the observed total ' ...
                 'in %.1f%% of fitted cell pixels. Those pixels contain no ' ...
                 'residual photons for membrane-state inference; Bayesian maps ' ...
                 'are saved for diagnosis but are marked invalid.'], ...
                100 * slbAmplitude.bayesianClippedPixelFraction);
        end
    end
    bayesMaps = scatterBayesMaps(bayesCompact, membranePixelIndices, size(intensity));
    orderedComponents = immune_cell_MIET_sorted_components(bayesCompact, ...
        membranePixelIndices, size(intensity), masks.cellMembrane, intensity, ...
        cfg.componentMaps);
    modelSelection = bayesianModelSelectionMetadata(bayesCompact, ...
        cfg.componentMaps);

    spatialBinningResult = struct('status', 'disabled', ...
        'config', cfg.spatialBinning, 'containsPhotonStream', false);
    binnedMeanFigureHandle = [];
    binnedBayesFigureHandle = [];
    binnedSecondLifetimeFigureHandle = [];
    binnedThirdLifetimeFigureHandle = [];
    binnedLayeredLifetimeFigureHandle = [];
    binnedMeanPngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_2x2_step1_mean_FLIM.png');
    binnedBayesPngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_2x2_step1_Bayes_maps.png');
    binnedSecondLifetimePngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_2x2_step1_second_sorted_lifetime.png');
    binnedThirdLifetimePngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_2x2_step1_third_sorted_lifetime.png');
    binnedLayeredLifetimePngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_2x2_step1_sorted_lifetime_layers.png');
    if cfg.spatialBinning.enabled
        fprintf(['immune_cell_MIET: collecting overlapping 2x2 reassigned ' ...
            'TCSPC windows (one-pixel step).\n']);
        [binnedCube, binnedPixelIndices, binnedIntensity, binnedTcspcMeta] = ...
            immune_cell_MIET_reassigned_sliding_tcspc(ptu640, ismResult, ...
            masks.cellMembrane, cfg.spatialBinning.binSize, ...
            cfg.photonBlockSize);
        if isempty(binnedPixelIndices)
            error('immune_cell_MIET:NoBinnedMembrane', ...
                'No complete 2x2 membrane-anchored TCSPC windows were available.');
        end
        if abs(binnedTcspcMeta.dtNs - dtNs) > max(1e-9, 1e-6 * dtNs)
            error('immune_cell_MIET:BinnedTimingMetadata', ...
                ['Sliding-window TCSPC bin width (%g ns) differs from reader ' ...
                 'metadata (%g ns).'], binnedTcspcMeta.dtNs, dtNs);
        end
        binnedMeanArrivalNs = compactMeanArrivalMap(binnedCube, ...
            binnedPixelIndices, size(intensity), compactMeta.globalCurve, ...
            dtNs, cfg.minPhotonsPerPixel);
        binnedSlbAmplitude = slidingSlbAmplitude(slbAmplitude, ...
            cfg.spatialBinning.binSize, cfg.spatialBinning.step);

        binnedBayesOpts = cfg.bayes;
        binnedBayesOpts.analysisMask = true(size(binnedCube, 1), 1);
        binnedBayesOpts.minPhotons = cfg.bayesMinPhotons;
        binnedBayesOpts.useGPU = cfg.useGPU;
        if cfg.fixSlbAmplitude
            binnedBayesOpts.fixedSlbPhotonCount = ...
                binnedSlbAmplitude.signalPhotonsPerPixel;
            binnedBayesOpts.fixedSlbPhotonCountStd = ...
                binnedSlbAmplitude.dispersionPhotonsPerPixel;
        else
            binnedBayesOpts.fixedSlbPhotonCount = [];
            binnedBayesOpts.fixedSlbPhotonCountStd = [];
        end
        fprintf(['immune_cell_MIET: running the 2x2/step-1 fixed-SLB ' ...
            'Bayesian comparison.\n']);
        binnedBayesCompact = flim_bayes_fixed_slb(binnedCube, irf, ...
            pulsePeriodNs, dtNs, tauSlbNs, binnedBayesOpts);
        [binnedBayesianStatus, binnedBayesianValidForInterpretation, ...
            binnedSlbAmplitude] = assessSlbCountClipping( ...
            binnedBayesCompact, binnedSlbAmplitude, cfg, ...
            'immune_cell_MIET:BinnedFixedSlbCountClipping', ...
            '2x2 sliding-window');
        binnedBayesMaps = scatterBayesMaps(binnedBayesCompact, ...
            binnedPixelIndices, size(intensity));
        binnedOrderedComponents = immune_cell_MIET_sorted_components( ...
            binnedBayesCompact, binnedPixelIndices, size(intensity), ...
            masks.cellMembrane, binnedIntensity, cfg.componentMaps);
        binnedModelSelection = bayesianModelSelectionMetadata( ...
            binnedBayesCompact, cfg.componentMaps);
        validAnchorMask = false(size(intensity));
        validAnchorMask(binnedPixelIndices) = true;
        spatialBinningResult = struct();
        spatialBinningResult.status = 'complete';
        spatialBinningResult.config = cfg.spatialBinning;
        spatialBinningResult.method = ...
            'overlapping reassigned 2x2 TCSPC sums with one-pixel stride';
        spatialBinningResult.meanFlim = struct( ...
            'meanArrivalNs', single(binnedMeanArrivalNs), ...
            'intensity', binnedIntensity, ...
            'source', 'summed reassigned TCSPC window histograms');
        spatialBinningResult.validAnchorMask = validAnchorMask;
        spatialBinningResult.bayesian = struct( ...
            'status', binnedBayesianStatus, ...
            'validForInterpretation', binnedBayesianValidForInterpretation, ...
            'compact', binnedBayesCompact, 'maps', binnedBayesMaps, ...
            'orderedComponents', binnedOrderedComponents, ...
            'modelSelection', binnedModelSelection, ...
            'pixelLinearIndices', binnedPixelIndices, ...
            'tcspcMeta', binnedTcspcMeta, 'candidateOrders', [1 2 3], ...
            'fixedSlbAmplitude', binnedSlbAmplitude);
        spatialBinningResult.containsPhotonStream = false;
        spatialBinningResult.compactTcspcCubeSaved = false;
        clear binnedCube
    end

    spatialBinning4x4Result = struct('status', 'disabled', ...
        'config', cfg.spatialBinning4x4, 'containsPhotonStream', false);
    binned4x4MeanFigureHandle = [];
    binned4x4BayesFigureHandle = [];
    binned4x4SecondLifetimeFigureHandle = [];
    binned4x4ThirdLifetimeFigureHandle = [];
    binned4x4LayeredLifetimeFigureHandle = [];
    binned4x4MeanPngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_4x4_step1_mean_FLIM.png');
    binned4x4BayesPngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_4x4_step1_Bayes_maps.png');
    binned4x4SecondLifetimePngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_4x4_step1_second_sorted_lifetime.png');
    binned4x4ThirdLifetimePngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_4x4_step1_third_sorted_lifetime.png');
    binned4x4LayeredLifetimePngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_4x4_step1_sorted_lifetime_layers.png');
    if cfg.spatialBinning4x4.enabled
        fprintf(['immune_cell_MIET: collecting overlapping 4x4 reassigned ' ...
            'TCSPC windows (one-pixel step).\n']);
        [binned4x4Cube, binned4x4PixelIndices, binned4x4Intensity, ...
            binned4x4TcspcMeta] = ...
            immune_cell_MIET_reassigned_sliding_tcspc(ptu640, ismResult, ...
            masks.cellMembrane, cfg.spatialBinning4x4.binSize, ...
            cfg.photonBlockSize);
        if isempty(binned4x4PixelIndices)
            error('immune_cell_MIET:NoBinned4x4Membrane', ...
                'No complete 4x4 membrane-anchored TCSPC windows were available.');
        end
        if abs(binned4x4TcspcMeta.dtNs - dtNs) > max(1e-9, 1e-6 * dtNs)
            error('immune_cell_MIET:Binned4x4TimingMetadata', ...
                ['Sliding-window TCSPC bin width (%g ns) differs from reader ' ...
                 'metadata (%g ns).'], binned4x4TcspcMeta.dtNs, dtNs);
        end
        binned4x4MeanArrivalNs = compactMeanArrivalMap(binned4x4Cube, ...
            binned4x4PixelIndices, size(intensity), compactMeta.globalCurve, ...
            dtNs, cfg.minPhotonsPerPixel);
        binned4x4SlbAmplitude = slidingSlbAmplitude(slbAmplitude, ...
            cfg.spatialBinning4x4.binSize, cfg.spatialBinning4x4.step);

        binned4x4BayesOpts = cfg.bayes;
        binned4x4BayesOpts.analysisMask = true(size(binned4x4Cube, 1), 1);
        binned4x4BayesOpts.minPhotons = cfg.bayesMinPhotons;
        binned4x4BayesOpts.useGPU = cfg.useGPU;
        if cfg.fixSlbAmplitude
            binned4x4BayesOpts.fixedSlbPhotonCount = ...
                binned4x4SlbAmplitude.signalPhotonsPerPixel;
            binned4x4BayesOpts.fixedSlbPhotonCountStd = ...
                binned4x4SlbAmplitude.dispersionPhotonsPerPixel;
        else
            binned4x4BayesOpts.fixedSlbPhotonCount = [];
            binned4x4BayesOpts.fixedSlbPhotonCountStd = [];
        end
        fprintf(['immune_cell_MIET: running the 4x4/step-1 fixed-SLB ' ...
            'Bayesian comparison.\n']);
        binned4x4BayesCompact = flim_bayes_fixed_slb(binned4x4Cube, irf, ...
            pulsePeriodNs, dtNs, tauSlbNs, binned4x4BayesOpts);
        [binned4x4BayesianStatus, ...
            binned4x4BayesianValidForInterpretation, ...
            binned4x4SlbAmplitude] = assessSlbCountClipping( ...
            binned4x4BayesCompact, binned4x4SlbAmplitude, cfg, ...
            'immune_cell_MIET:Binned4x4FixedSlbCountClipping', ...
            '4x4 sliding-window');
        binned4x4BayesMaps = scatterBayesMaps(binned4x4BayesCompact, ...
            binned4x4PixelIndices, size(intensity));
        binned4x4OrderedComponents = immune_cell_MIET_sorted_components( ...
            binned4x4BayesCompact, binned4x4PixelIndices, size(intensity), ...
            masks.cellMembrane, binned4x4Intensity, cfg.componentMaps);
        binned4x4ModelSelection = bayesianModelSelectionMetadata( ...
            binned4x4BayesCompact, cfg.componentMaps);
        valid4x4AnchorMask = false(size(intensity));
        valid4x4AnchorMask(binned4x4PixelIndices) = true;
        spatialBinning4x4Result = struct();
        spatialBinning4x4Result.status = 'complete';
        spatialBinning4x4Result.config = cfg.spatialBinning4x4;
        spatialBinning4x4Result.method = ...
            'overlapping reassigned 4x4 TCSPC sums with one-pixel stride';
        spatialBinning4x4Result.meanFlim = struct( ...
            'meanArrivalNs', single(binned4x4MeanArrivalNs), ...
            'intensity', binned4x4Intensity, ...
            'source', 'summed reassigned TCSPC window histograms');
        spatialBinning4x4Result.validAnchorMask = valid4x4AnchorMask;
        spatialBinning4x4Result.bayesian = struct( ...
            'status', binned4x4BayesianStatus, ...
            'validForInterpretation', ...
                binned4x4BayesianValidForInterpretation, ...
            'compact', binned4x4BayesCompact, 'maps', binned4x4BayesMaps, ...
            'orderedComponents', binned4x4OrderedComponents, ...
            'modelSelection', binned4x4ModelSelection, ...
            'pixelLinearIndices', binned4x4PixelIndices, ...
            'tcspcMeta', binned4x4TcspcMeta, 'candidateOrders', [1 2 3], ...
            'fixedSlbAmplitude', binned4x4SlbAmplitude);
        spatialBinning4x4Result.containsPhotonStream = false;
        spatialBinning4x4Result.compactTcspcCubeSaved = false;
        fprintf(['immune_cell_MIET: 4x4/step-1 stage complete; %d complete ' ...
            'windows; second-lifetime pixels %d; third-lifetime pixels %d.\n'], ...
            numel(binned4x4PixelIndices), ...
            binned4x4OrderedComponents.summary.secondDisplayPixelCount, ...
            binned4x4OrderedComponents.summary.thirdDisplayPixelCount);
        clear binned4x4Cube
    else
        fprintf(['immune_cell_MIET: 4x4/step-1 stage SKIPPED ' ...
            '(cfg.spatialBinning4x4.enabled is false).\n']);
    end
    if cfg.spatialBinning4x4.enabled && ...
            ~strcmp(spatialBinning4x4Result.status, 'complete')
        error('immune_cell_MIET:Incomplete4x4Stage', ...
            ['The 4x4 sliding-window stage is enabled but its status is ' ...
             '%s. Refusing to continue with an incomplete analysis.'], ...
            spatialBinning4x4Result.status);
    end
    clear ptu640

    bayesFigureHandle = [];
    secondLifetimeFigureHandle = [];
    thirdLifetimeFigureHandle = [];
    layeredLifetimeFigureHandle = [];
    secondLifetimePngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_second_sorted_lifetime.png');
    thirdLifetimePngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_third_sorted_lifetime.png');
    layeredLifetimePngFile = fullfile(outputDir, ...
        'immune_cell_MIET_640nm_red_sorted_lifetime_layers.png');
    if cfg.showFigure || cfg.saveOutputs
        visibility = figureVisibility(cfg.showFigure);
        bayesFigureHandle = showBayesSummary(meanArrivalNs, intensity, ...
            masks, bayesMaps, tauSlbNs, visibility);
        if cfg.componentMaps.enabled
            lifetimeLimits = orderedLifetimeLimits(orderedComponents, ...
                cfg.componentMaps.lifetimeRangeNs);
            secondLifetimeFigureHandle = showOrderedLifetimeMap( ...
                orderedComponents, masks, 2, lifetimeLimits, visibility, ...
                tauSlbNs, slbAmplitude);
            thirdLifetimeFigureHandle = showOrderedLifetimeMap( ...
                orderedComponents, masks, 3, lifetimeLimits, visibility, ...
                tauSlbNs, slbAmplitude);
            layeredLifetimeFigureHandle = showOrderedLifetimeLayers( ...
                orderedComponents, masks, lifetimeLimits, visibility, ...
                tauSlbNs, slbAmplitude);
        end
        if cfg.spatialBinning.enabled
            binnedMeanFigureHandle = showRedMeanFlim( ...
                double(spatialBinningResult.meanFlim.meanArrivalNs), ...
                spatialBinningResult.meanFlim.intensity, cfg.lifetimeRangeNs, ...
                cfg.excitationNm, visibility, ...
                '2x2 sliding TCSPC, one-pixel step');
            binnedBayesFigureHandle = showBayesSummary( ...
                double(spatialBinningResult.meanFlim.meanArrivalNs), ...
                spatialBinningResult.meanFlim.intensity, masks, ...
                spatialBinningResult.bayesian.maps, tauSlbNs, visibility, ...
                '2x2 sliding TCSPC, one-pixel step; native FLIM segmentation reused');
            if cfg.componentMaps.enabled
                binnedLifetimeLimits = orderedLifetimeLimits( ...
                    spatialBinningResult.bayesian.orderedComponents, ...
                    cfg.componentMaps.lifetimeRangeNs);
                binnedSecondLifetimeFigureHandle = showOrderedLifetimeMap( ...
                    spatialBinningResult.bayesian.orderedComponents, masks, 2, ...
                    binnedLifetimeLimits, visibility, tauSlbNs, ...
                    spatialBinningResult.bayesian.fixedSlbAmplitude, ...
                    '2x2 sliding TCSPC, one-pixel step');
                binnedThirdLifetimeFigureHandle = showOrderedLifetimeMap( ...
                    spatialBinningResult.bayesian.orderedComponents, masks, 3, ...
                    binnedLifetimeLimits, visibility, tauSlbNs, ...
                    spatialBinningResult.bayesian.fixedSlbAmplitude, ...
                    '2x2 sliding TCSPC, one-pixel step');
                binnedLayeredLifetimeFigureHandle = showOrderedLifetimeLayers( ...
                    spatialBinningResult.bayesian.orderedComponents, masks, ...
                    binnedLifetimeLimits, visibility, tauSlbNs, ...
                    spatialBinningResult.bayesian.fixedSlbAmplitude, ...
                    '2x2 sliding TCSPC, one-pixel step');
            end
        end
        if cfg.spatialBinning4x4.enabled
            binned4x4MeanFigureHandle = showRedMeanFlim( ...
                double(spatialBinning4x4Result.meanFlim.meanArrivalNs), ...
                spatialBinning4x4Result.meanFlim.intensity, ...
                cfg.lifetimeRangeNs, cfg.excitationNm, visibility, ...
                '4x4 sliding TCSPC, one-pixel step');
            binned4x4BayesFigureHandle = showBayesSummary( ...
                double(spatialBinning4x4Result.meanFlim.meanArrivalNs), ...
                spatialBinning4x4Result.meanFlim.intensity, masks, ...
                spatialBinning4x4Result.bayesian.maps, tauSlbNs, visibility, ...
                '4x4 sliding TCSPC, one-pixel step; native FLIM segmentation reused');
            if cfg.componentMaps.enabled
                binned4x4LifetimeLimits = orderedLifetimeLimits( ...
                    spatialBinning4x4Result.bayesian.orderedComponents, ...
                    cfg.componentMaps.lifetimeRangeNs);
                binned4x4SecondLifetimeFigureHandle = showOrderedLifetimeMap( ...
                    spatialBinning4x4Result.bayesian.orderedComponents, ...
                    masks, 2, binned4x4LifetimeLimits, visibility, tauSlbNs, ...
                    spatialBinning4x4Result.bayesian.fixedSlbAmplitude, ...
                    '4x4 sliding TCSPC, one-pixel step');
                binned4x4ThirdLifetimeFigureHandle = showOrderedLifetimeMap( ...
                    spatialBinning4x4Result.bayesian.orderedComponents, ...
                    masks, 3, binned4x4LifetimeLimits, visibility, tauSlbNs, ...
                    spatialBinning4x4Result.bayesian.fixedSlbAmplitude, ...
                    '4x4 sliding TCSPC, one-pixel step');
                binned4x4LayeredLifetimeFigureHandle = ...
                    showOrderedLifetimeLayers( ...
                    spatialBinning4x4Result.bayesian.orderedComponents, ...
                    masks, binned4x4LifetimeLimits, visibility, tauSlbNs, ...
                    spatialBinning4x4Result.bayesian.fixedSlbAmplitude, ...
                    '4x4 sliding TCSPC, one-pixel step');
            end
        end
    end

    result = struct();
    result.pipeline = 'immune_cell_MIET';
    result.analysisSchemaVersion = 6;
    result.stage = '640 nm red-channel fixed-SLB Bayesian ISM-FLIM';
    result.sourceFile = ptuFile;
    result.config = cfg;
    result.channel = struct('excitationNm', cfg.excitationNm, ...
        'piePulseIndex', pulseIndex, 'gateStartBin', gateStarts(pulseIndex), ...
        'gateStopBin', gateStops(pulseIndex), 'gateLengthBins', gateLength, ...
        'recordedExcitationNm', recordedExcitationNm, ...
        'ismWavelengthNm', cfg.ismWavelengthNm, ...
        'detectorIds', detectorIds(:).', 'dtNs', dtNs);
    result.gateDetection = struct('source', gateSource, ...
        'starts', gateStarts, 'stops', gateStops, ...
        'profile', gateProfile, 'smoothProfile', smoothGateProfile, ...
        'threshold', gateThreshold, 'baseline', gateBaseline, ...
        'riseBins', riseBins, 'peakBins', peakBins, 'peakHeights', peakHeights);
    result.ism = ismResult;
    result.reassigned = reassigned;
    result.redMeanFlim = struct('meanArrivalNs', meanArrivalNs, ...
        'intensity', intensity, 'tcspcSource', tcspcSource);
    result.blueMeanFlim = result.redMeanFlim; % Backward-compatible first-stage field.
    result.segmentation = struct('masks', masks, ...
        'diagnostics', segmentationDiagnostics);
    result.irf = struct('curve', irf, 'source', irfSource, ...
        'parameters', irfParameters, 'estimatedTauNs', irfTau, ...
        'fitCounts', irfFit, 'chi', irfChi);
    result.slbReference = struct('fixedLifetimeNs', tauSlbNs, ...
        'tauMeanNs', slbTauMean, 'curve', slbCurve, ...
        'fitCounts', slbFitCounts, 'coefficients', slbCoefficients, ...
        'speciesFraction', slbFraction, 'deviance', slbDeviance, ...
        'bic', slbBic, 'reducedChiSquare', slbReducedChiSquare, ...
        'photonCount', sum(slbCurve), 'amplitudeConstraint', slbAmplitude);
    result.bayesian = struct('status', bayesianStatus, ...
        'validForInterpretation', bayesianValidForInterpretation, ...
        'compact', bayesCompact, 'maps', bayesMaps, ...
        'orderedComponents', orderedComponents, ...
        'modelSelection', modelSelection, ...
        'membranePixelIndices', membranePixelIndices, ...
        'tcspcMeta', compactMeta, 'candidateOrders', [1 2 3]);
    result.spatialBinning = spatialBinningResult;
    result.spatialBinning4x4 = spatialBinning4x4Result;
    result.tcspcPix = struct('saved', tcspcPixSaved, ...
        'matVariable', 'tcspc_pix', 'class', 'uint16', ...
        'size', compactMeta.fullCubeSize, ...
        'detectorsSummed', true, 'spatiallyReassigned', true, ...
        'containsPhotonStream', false);
    result.containsPhotonStream = false;
    result.outputFiles = struct('mat', '', ...
        'preliminaryMat', preliminaryMatFile, ...
        'segmentationMat', segmentationMatFile, ...
        'preliminaryMeanFlimPng', preliminaryPngFile, ...
        'segmentationPng', segmentationPngFile, ...
        'meanFlimPng', '', 'bayesPng', '', ...
        'secondSortedLifetimePng', secondLifetimePngFile, ...
        'thirdSortedLifetimePng', thirdLifetimePngFile, ...
        'sortedLifetimeLayersPng', layeredLifetimePngFile, ...
        'sliding2x2Step1MeanFlimPng', binnedMeanPngFile, ...
        'sliding2x2Step1BayesPng', binnedBayesPngFile, ...
        'sliding2x2Step1SecondSortedLifetimePng', ...
            binnedSecondLifetimePngFile, ...
        'sliding2x2Step1ThirdSortedLifetimePng', ...
            binnedThirdLifetimePngFile, ...
        'sliding2x2Step1SortedLifetimeLayersPng', ...
            binnedLayeredLifetimePngFile, ...
        'sliding4x4Step1MeanFlimPng', binned4x4MeanPngFile, ...
        'sliding4x4Step1BayesPng', binned4x4BayesPngFile, ...
        'sliding4x4Step1SecondSortedLifetimePng', ...
            binned4x4SecondLifetimePngFile, ...
        'sliding4x4Step1ThirdSortedLifetimePng', ...
            binned4x4ThirdLifetimePngFile, ...
        'sliding4x4Step1SortedLifetimeLayersPng', ...
            binned4x4LayeredLifetimePngFile);

    if cfg.saveOutputs
        matFile = analysisMatFile;
        pngFile = fullfile(outputDir, 'immune_cell_MIET_640nm_red_mean_FLIM.png');
        result.outputFiles.mat = matFile;
        bayesPngFile = fullfile(outputDir, 'immune_cell_MIET_640nm_red_Bayes_maps.png');
        result.outputFiles.meanFlimPng = pngFile;
        result.outputFiles.bayesPng = bayesPngFile;
        if ~isempty(figureHandle)
            exportgraphics(figureHandle, pngFile, 'Resolution', 300);
        end
        if ~isempty(bayesFigureHandle)
            exportgraphics(bayesFigureHandle, bayesPngFile, 'Resolution', 300);
        end
        if ~isempty(secondLifetimeFigureHandle)
            exportgraphics(secondLifetimeFigureHandle, secondLifetimePngFile, ...
                'Resolution', 300);
        end
        if ~isempty(thirdLifetimeFigureHandle)
            exportgraphics(thirdLifetimeFigureHandle, thirdLifetimePngFile, ...
                'Resolution', 300);
        end
        if ~isempty(layeredLifetimeFigureHandle)
            exportgraphics(layeredLifetimeFigureHandle, layeredLifetimePngFile, ...
                'Resolution', 300);
        end
        if ~isempty(binnedMeanFigureHandle)
            exportgraphics(binnedMeanFigureHandle, binnedMeanPngFile, ...
                'Resolution', 300);
        end
        if ~isempty(binnedBayesFigureHandle)
            exportgraphics(binnedBayesFigureHandle, binnedBayesPngFile, ...
                'Resolution', 300);
        end
        if ~isempty(binnedSecondLifetimeFigureHandle)
            exportgraphics(binnedSecondLifetimeFigureHandle, ...
                binnedSecondLifetimePngFile, 'Resolution', 300);
        end
        if ~isempty(binnedThirdLifetimeFigureHandle)
            exportgraphics(binnedThirdLifetimeFigureHandle, ...
                binnedThirdLifetimePngFile, 'Resolution', 300);
        end
        if ~isempty(binnedLayeredLifetimeFigureHandle)
            exportgraphics(binnedLayeredLifetimeFigureHandle, ...
                binnedLayeredLifetimePngFile, 'Resolution', 300);
        end
        if ~isempty(binned4x4MeanFigureHandle)
            exportgraphics(binned4x4MeanFigureHandle, binned4x4MeanPngFile, ...
                'Resolution', 300);
        end
        if ~isempty(binned4x4BayesFigureHandle)
            exportgraphics(binned4x4BayesFigureHandle, ...
                binned4x4BayesPngFile, 'Resolution', 300);
        end
        if ~isempty(binned4x4SecondLifetimeFigureHandle)
            exportgraphics(binned4x4SecondLifetimeFigureHandle, ...
                binned4x4SecondLifetimePngFile, 'Resolution', 300);
        end
        if ~isempty(binned4x4ThirdLifetimeFigureHandle)
            exportgraphics(binned4x4ThirdLifetimeFigureHandle, ...
                binned4x4ThirdLifetimePngFile, 'Resolution', 300);
        end
        if ~isempty(binned4x4LayeredLifetimeFigureHandle)
            exportgraphics(binned4x4LayeredLifetimeFigureHandle, ...
                binned4x4LayeredLifetimePngFile, 'Resolution', 300);
        end
        if cfg.spatialBinning4x4.enabled
            expected4x4Figures = {binned4x4MeanPngFile, binned4x4BayesPngFile};
            if cfg.componentMaps.enabled
                expected4x4Figures = [expected4x4Figures, ...
                    {binned4x4SecondLifetimePngFile, ...
                     binned4x4ThirdLifetimePngFile, ...
                     binned4x4LayeredLifetimePngFile}];
            end
            for figureIndex = 1:numel(expected4x4Figures)
                if ~isfile(expected4x4Figures{figureIndex})
                    error('immune_cell_MIET:Missing4x4Figure', ...
                        ['The 4x4 stage is enabled but its figure %s was ' ...
                         'not written.'], expected4x4Figures{figureIndex});
                end
                fprintf('immune_cell_MIET: saved 4x4 figure %s\n', ...
                    expected4x4Figures{figureIndex});
            end
        end
        if tcspcPixSaved
            save(matFile, 'result', '-append');
        else
            save(matFile, 'result', '-v7.3');
        end
        fprintf('immune_cell_MIET: saved %s\n', outputDir);
    end
end

function visibility = figureVisibility(showFigure)
    visibility = 'off';
    if showFigure
        visibility = 'on';
    end
end

function cfg = fillDefaults(cfg)
    componentMapDefaults = struct('enabled', true, ...
        'posteriorThreshold', [0.8 0.95], ...
        'minExpectedPhotons', [10 10], 'lifetimeRangeNs', [], ...
        'probabilityContourLevels', [0.5 0.7 0.9]);
    spatialBinningDefaults = struct('enabled', false, ...
        'binSize', [2 2], 'step', [1 1]);
    spatialBinning4x4Defaults = struct('enabled', false, ...
        'binSize', [4 4], 'step', [1 1]);
    defaults = struct('excitationPulseIndex', 2, 'excitationNm', 640, ...
        'ismWavelengthNm', 690, 'wavelengthToleranceNm', 10, ...
        'piePulseCount', 2, 'tcspcBinNs', 0.16, 'maxNgate', 1024, ...
        'photonsPerChunk', 5e6, 'photonBlockSize', 1e6, ...
        'showWaitbar', false, 'useGPU', false, ...
        'gateThresholdFrac', 0.15, 'gatePreNs', 1.5, ...
        'minGateSeparationBins', 50, 'weakPeakFraction', 0.05, ...
        'minPhotonsPerPixel', 5, 'lifetimeRangeNs', [], ...
        'showFigure', true, 'saveOutputs', true, 'saveTcspcPix', true, ...
        'requireValidSegmentation', true, ...
        'outputDir', '', ...
        'minSlbPhotons', 500, 'minPhotonsForIRF', 500, ...
        'irf', [], 'irfNCasc', 4, 'irfNSub', 6, ...
        'irfTauSeedsNs', [0.4 2.5], 'slbTauSeedNs', 0.5, ...
        'membraneTauSeedsNs', [1.2 3.0], 'lifetimeBoundsNs', [0.03 10], ...
        'includeBackground', true, 'fitMaxIter', 500, ...
        'fixSlbAmplitude', true, 'slbAmplitudeScale', 1, ...
        'maxSlbClippedFraction', 0.10, ...
        'bayesMinPhotons', 10, 'ism', struct(), 'segmentation', struct(), ...
        'componentMaps', componentMapDefaults, ...
        'spatialBinning', spatialBinningDefaults, ...
        'spatialBinning4x4', spatialBinning4x4Defaults, ...
        'bayes', struct('batchSize', 256, 'includeBackground', true, ...
            'signalGrid', [0.25 0.5 0.75 1], 'membraneTauCount', 10, ...
            'fractionStep', 0.2, 'minimumMembraneFraction', 0.1));
    cfg = mergeStruct(defaults, cfg);
    cfg.componentMaps = mergeStruct(componentMapDefaults, cfg.componentMaps);
    cfg.spatialBinning = mergeStruct(spatialBinningDefaults, cfg.spatialBinning);
    cfg.spatialBinning4x4 = mergeStruct( ...
        spatialBinning4x4Defaults, cfg.spatialBinning4x4);
    validateattributes(cfg.excitationPulseIndex, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.piePulseCount, {'numeric'}, {'scalar','integer','>=',2});
    validateattributes(cfg.fixSlbAmplitude, {'numeric','logical'}, {'scalar'});
    cfg.fixSlbAmplitude = logical(cfg.fixSlbAmplitude);
    validateattributes(cfg.requireValidSegmentation, ...
        {'numeric','logical'}, {'scalar'});
    cfg.requireValidSegmentation = logical(cfg.requireValidSegmentation);
    validateattributes(cfg.slbAmplitudeScale, {'numeric'}, ...
        {'real','finite','scalar','positive'});
    validateattributes(cfg.maxSlbClippedFraction, {'numeric'}, ...
        {'real','finite','scalar','>=',0,'<=',1});
    validateattributes(cfg.componentMaps.enabled, {'numeric','logical'}, {'scalar'});
    cfg.componentMaps.enabled = logical(cfg.componentMaps.enabled);
    validateattributes(cfg.spatialBinning.enabled, {'numeric','logical'}, {'scalar'});
    cfg.spatialBinning.enabled = logical(cfg.spatialBinning.enabled);
    validateattributes(cfg.spatialBinning.binSize, {'numeric'}, ...
        {'real','finite','integer','positive','vector','numel',2});
    validateattributes(cfg.spatialBinning.step, {'numeric'}, ...
        {'real','finite','integer','positive','vector','numel',2});
    cfg.spatialBinning.binSize = double(cfg.spatialBinning.binSize(:)).';
    cfg.spatialBinning.step = double(cfg.spatialBinning.step(:)).';
    if ~isequal(cfg.spatialBinning.binSize, [2 2]) || ...
            ~isequal(cfg.spatialBinning.step, [1 1])
        error('immune_cell_MIET:SpatialBinning', ...
            ['The full-size sliding analysis currently supports only 2x2 ' ...
             'TCSPC windows with a one-pixel step.']);
    end
    validateattributes(cfg.spatialBinning4x4.enabled, ...
        {'numeric','logical'}, {'scalar'});
    cfg.spatialBinning4x4.enabled = logical(cfg.spatialBinning4x4.enabled);
    validateattributes(cfg.spatialBinning4x4.binSize, {'numeric'}, ...
        {'real','finite','integer','positive','vector','numel',2});
    validateattributes(cfg.spatialBinning4x4.step, {'numeric'}, ...
        {'real','finite','integer','positive','vector','numel',2});
    cfg.spatialBinning4x4.binSize = ...
        double(cfg.spatialBinning4x4.binSize(:)).';
    cfg.spatialBinning4x4.step = double(cfg.spatialBinning4x4.step(:)).';
    if ~isequal(cfg.spatialBinning4x4.binSize, [4 4]) || ...
            ~isequal(cfg.spatialBinning4x4.step, [1 1])
        error('immune_cell_MIET:SpatialBinning4x4', ...
            ['The second sliding analysis supports only 4x4 TCSPC windows ' ...
             'with a one-pixel step.']);
    end
end

function amplitude = estimateSlbAmplitudeConstraint(slbCurve, slbFitCounts, ...
        coefficients, intensity, referenceMask, amplitudeMask, scale)
    fitReferencePixelCount = nnz(referenceMask);
    if fitReferencePixelCount < 1
        error('immune_cell_MIET:NoSLBAmplitudeReference', ...
            'At least one outside-SLB pixel is required to calibrate amplitude.');
    end
    amplitudePopulationPixelCount = nnz(amplitudeMask);
    usedReferenceFallback = amplitudePopulationPixelCount < 1;
    if usedReferenceFallback
        amplitudeMask = referenceMask;
        amplitudePopulationPixelCount = fitReferencePixelCount;
    end
    signalPhotonTotal = max(double(coefficients(2)), 0);
    fittedPhotonTotal = sum(max(double(slbFitCounts(:)), 0));
    observedPhotonTotal = sum(max(double(slbCurve(:)), 0));
    signalFraction = signalPhotonTotal / max(fittedPhotonTotal, eps);
    stableReferenceSignalPhotonsPerPixel = ...
        scale * signalPhotonTotal / fitReferencePixelCount;

    % The pooled, lifetime-stable reference determines the signal fraction.
    % Apply that fraction to the broader exterior-SLB population to obtain the
    % constant expected SLB signal count used in every membrane-pixel model.
    outsideCounts = max(double(intensity(amplitudeMask)), 0);
    proxySignalCounts = scale * signalFraction * outsideCounts;
    signalPhotonsPerPixel = mean(proxySignalCounts);
    if ~(isfinite(signalPhotonsPerPixel) && signalPhotonsPerPixel > 0)
        error('immune_cell_MIET:InvalidSLBAmplitude', ...
            'The fitted outside-SLB signal amplitude is not positive.');
    end
    proxyMedian = median(proxySignalCounts);
    robustSigma = 1.4826 * median(abs(proxySignalCounts - proxyMedian));
    poissonSigma = sqrt(max(signalPhotonsPerPixel, 1));
    dispersion = max([robustSigma, poissonSigma, 1]);

    amplitude = struct();
    amplitude.constraintApplied = true;
    amplitude.method = ['pooled stable outside-SLB IRF fit for signal fraction; ' ...
        'mean signal proxy over the geometric outside-SLB population, ' ...
        'including low-count and zero-count pixels'];
    amplitude.referencePixelCount = fitReferencePixelCount;
    amplitude.fitReferencePixelCount = fitReferencePixelCount;
    amplitude.amplitudePopulationPixelCount = amplitudePopulationPixelCount;
    amplitude.usedReferenceFallbackForAmplitudePopulation = usedReferenceFallback;
    amplitude.scale = scale;
    amplitude.signalPhotonsPerPixel = signalPhotonsPerPixel;
    amplitude.stableReferenceSignalPhotonsPerPixel = ...
        stableReferenceSignalPhotonsPerPixel;
    amplitude.dispersionPhotonsPerPixel = dispersion;
    amplitude.poissonSigmaPhotonsPerPixel = poissonSigma;
    amplitude.robustSigmaPhotonsPerPixel = robustSigma;
    amplitude.relativeDispersion = dispersion / signalPhotonsPerPixel;
    amplitude.proxyMedianPhotonsPerPixel = proxyMedian;
    amplitude.meanObservedPhotonsPerPixel = ...
        observedPhotonTotal / fitReferencePixelCount;
    amplitude.meanAmplitudePopulationObservedPhotonsPerPixel = mean(outsideCounts);
    amplitude.fittedSignalPhotonFraction = signalFraction;
    amplitude.fittedBackgroundPhotonsPerPixel = ...
        max(fittedPhotonTotal - signalPhotonTotal, 0) / fitReferencePixelCount;
end

function amplitude = slidingSlbAmplitude(nativeAmplitude, windowSize, step)
    area = prod(double(windowSize));
    amplitude = nativeAmplitude;
    amplitude.nativeSignalPhotonsPerPixel = ...
        nativeAmplitude.signalPhotonsPerPixel;
    amplitude.nativeDispersionPhotonsPerPixel = ...
        nativeAmplitude.dispersionPhotonsPerPixel;
    amplitude.signalPhotonsPerPixel = area * ...
        nativeAmplitude.signalPhotonsPerPixel;
    % This width is retained for calibration metadata. The direct fixed-count
    % likelihood uses the expected count itself, not this uncertainty value.
    amplitude.dispersionPhotonsPerPixel = sqrt(area) * ...
        nativeAmplitude.dispersionPhotonsPerPixel;
    amplitude.poissonSigmaPhotonsPerPixel = ...
        sqrt(max(amplitude.signalPhotonsPerPixel, 1));
    amplitude.robustSigmaPhotonsPerPixel = sqrt(area) * ...
        nativeAmplitude.robustSigmaPhotonsPerPixel;
    amplitude.relativeDispersion = amplitude.dispersionPhotonsPerPixel / ...
        amplitude.signalPhotonsPerPixel;
    fieldsScaledByArea = {'stableReferenceSignalPhotonsPerPixel', ...
        'proxyMedianPhotonsPerPixel', 'meanObservedPhotonsPerPixel', ...
        'meanAmplitudePopulationObservedPhotonsPerPixel', ...
        'fittedBackgroundPhotonsPerPixel'};
    for k = 1:numel(fieldsScaledByArea)
        name = fieldsScaledByArea{k};
        if isfield(nativeAmplitude, name)
            amplitude.(name) = area * nativeAmplitude.(name);
        end
    end
    amplitude.method = [nativeAmplitude.method '; expected window amplitude ' ...
        'equals native photons/pixel multiplied by the complete window area'];
    amplitude.windowSize = double(windowSize(:)).';
    amplitude.windowArea = area;
    amplitude.step = double(step(:)).';
    amplitude.windowCentreOffsetPixels = (double(windowSize(:)).' - 1) / 2;
    amplitude.dispersionUsedForInference = false;
    amplitude.bayesianClippedPixelCount = 0;
    amplitude.bayesianClippedPixelFraction = 0;
end

function meanArrival = compactMeanArrivalMap(pixelCube, pixelIndices, ...
        imageSize, globalCurve, dtNs, minimumPhotons)
    pixelCount = size(pixelCube, 1);
    gateLength = size(pixelCube, 3);
    counts = zeros(pixelCount, 1, 'double');
    firstMoment = zeros(pixelCount, 1, 'double');
    [~, t0Bin] = max(double(globalCurve(:)));
    for timeIndex = 1:gateLength
        timeCounts = double(pixelCube(:, 1, timeIndex));
        counts = counts + timeCounts;
        firstMoment = firstMoment + ...
            timeCounts * ((timeIndex - t0Bin) * dtNs);
    end
    compactMean = nan(pixelCount, 1);
    valid = counts >= minimumPhotons;
    compactMean(valid) = max(firstMoment(valid) ./ counts(valid), 0);
    meanArrival = nan(imageSize);
    meanArrival(pixelIndices) = compactMean;
end

function [status, validForInterpretation, amplitude] = ...
        assessSlbCountClipping(bayesCompact, amplitude, cfg, ...
        warningIdentifier, analysisLabel)
    amplitude.bayesianClippedPixelCount = 0;
    amplitude.bayesianClippedPixelFraction = 0;
    status = 'complete';
    validForInterpretation = true;
    if ~cfg.fixSlbAmplitude || ...
            ~isfield(bayesCompact, 'fixedSlbPhotonConstraint')
        return;
    end
    constraint = bayesCompact.fixedSlbPhotonConstraint;
    amplitude.bayesianClippedPixelCount = constraint.clippedPixelCount;
    amplitude.bayesianClippedPixelFraction = constraint.clippedPixelFraction;
    if amplitude.bayesianClippedPixelFraction <= cfg.maxSlbClippedFraction
        return;
    end
    status = 'invalid_fixed_slb_count_clipping';
    validForInterpretation = false;
    warning(warningIdentifier, ...
        ['The outside-derived SLB count exceeded the observed total in %.1f%% ' ...
         'of fitted %s pixels. Those pixels contain no residual photons for ' ...
         'membrane-state inference; maps are saved for diagnosis but are ' ...
         'marked invalid.'], 100 * amplitude.bayesianClippedPixelFraction, ...
        analysisLabel);
end

function metadata = bayesianModelSelectionMetadata(bayesCompact, componentCfg)
    posteriorThreshold = double(componentCfg.posteriorThreshold(:)).';
    minimumExpectedPhotons = double(componentCfg.minExpectedPhotons(:)).';
    if isscalar(posteriorThreshold)
        posteriorThreshold = [posteriorThreshold posteriorThreshold];
    end
    if isscalar(minimumExpectedPhotons)
        minimumExpectedPhotons = ...
            [minimumExpectedPhotons minimumExpectedPhotons];
    end
    metadata = struct();
    metadata.criterion = ...
        'posterior model probability from discrete marginal likelihood';
    metadata.bicUsedForPixelModelOrder = false;
    metadata.likelihood = ...
        'multinomial TCSPC shape likelihood conditional on observed photon total';
    metadata.evidence = ['uniform-prior average of likelihood over each ' ...
        'model''s configured lifetime/fraction/background grid'];
    metadata.modelNames = bayesCompact.modelNames;
    metadata.modelPrior = bayesCompact.modelPrior;
    if isfield(bayesCompact, 'gridInfo') && ...
            isfield(bayesCompact.gridInfo, 'modelStateCount')
        metadata.modelStateCount = bayesCompact.gridInfo.modelStateCount;
    else
        metadata.modelStateCount = [];
    end
    metadata.secondDisplayCriterion = sprintf( ...
        ['P(M2)+P(M3) >= %.4g and posterior-expected component-2 ' ...
         'photons >= %.4g'], posteriorThreshold(1), ...
        minimumExpectedPhotons(1));
    metadata.thirdDisplayCriterion = sprintf( ...
        ['P(M3 | M2 or M3) >= %.4g and posterior-expected component-3 ' ...
         'photons >= %.4g; also inside component-2 display support'], ...
        posteriorThreshold(2), ...
        minimumExpectedPhotons(2));
    metadata.incrementalM3VsM2 = struct( ...
        'logBayesFactorMapField', 'maps.logBayesFactorM3VsM2', ...
        'conditionalProbabilityMapField', ...
            'maps.probabilityM3GivenM2orM3', ...
        'usedForCurrentDisplayGate', true);
    metadata.note = ['Model probabilities are conditional on the configured ' ...
        'grids, priors, fixed SLB calibration, and IRF; they are not BIC ' ...
        'scores or prior-independent probabilities of physical species.'];
end

function periodNs = resolvePulsePeriodNs(ptu)
    periodNs = double(ptu.head.MeasDesc_GlobalResolution) * 1e9;
    if ~(isfinite(periodNs) && periodNs > 0)
        periodNs = double(ptu.Ngate) * resolveDtNs(ptu);
    end
end

function out = mergeStruct(out, overrides)
    if isempty(overrides)
        return;
    end
    names = fieldnames(overrides);
    for k = 1:numel(names)
        out.(names{k}) = overrides.(names{k});
    end
end

function assertPhotonLists(ptu)
    fields = {'head', 'Ngate', 'im_tcspc', 'im_chan', 'im_col', 'im_line'};
    for k = 1:numel(fields)
        if ~isfield(ptu, fields{k}) || isempty(ptu.(fields{k}))
            error('immune_cell_MIET:PTUContent', ...
                'PTU_FLIM_GPU did not return required field %s.', fields{k});
        end
    end
end

function value = getFieldOr(s, name, fallback)
    value = fallback;
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    end
end

function dtNs = resolveDtNs(ptu)
    dtNs = getFieldOr(ptu, 'Resolution_ns', []);
    if isempty(dtNs)
        dtNs = double(ptu.head.MeasDesc_Resolution) * 1e9;
    end
    dtNs = double(dtNs);
end

function wavelengths = resolveRecordedWavelengths(head)
    wavelengths = [];
    candidateFields = {'LaserWL', 'LaserWavelength', 'LaserWavelengths'};
    for k = 1:numel(candidateFields)
        name = candidateFields{k};
        if isfield(head, name) && ~isempty(head.(name))
            value = head.(name);
            if iscell(value)
                value = cellfun(@double, value);
            end
            if isnumeric(value)
                wavelengths = double(value(:)).';
                return;
            end
        end
    end
end

function [ok, gateStartNative, gateLengthNative] = nativeReaderGate(head, cfg)
    ok = false;
    gateStartNative = [];
    gateLengthNative = [];
    requiredFields = {'MeasDesc_Resolution', ...
        'MeasDesc_GlobalResolution', 'PIETimeGate'};
    if ~all(isfield(head, requiredFields))
        return;
    end
    rawResolutionSeconds = double(head.MeasDesc_Resolution);
    periodSeconds = double(head.MeasDesc_GlobalResolution);
    if ~(isfinite(rawResolutionSeconds) && rawResolutionSeconds > 0 && ...
            isfinite(periodSeconds) && periodSeconds > 0)
        return;
    end
    requestedSeconds = max(rawResolutionSeconds, cfg.tcspcBinNs * 1e-9);
    channelDivision = max(1, round(requestedSeconds / rawResolutionSeconds));
    actualResolutionSeconds = channelDivision * rawResolutionSeconds;
    periodBins = min(cfg.maxNgate, ...
        ceil(periodSeconds / actualResolutionSeconds) + 1);
    gateHead = head;
    gateHead.TCSPC_chDiv = channelDivision;
    preGateBins = max(1, round(cfg.gatePreNs / ...
        (actualResolutionSeconds * 1e9)));
    [headerOk, starts, ~, gateLength] = gatesFromPtuHeader( ...
        gateHead, periodBins, cfg.piePulseCount, preGateBins);
    if ~headerOk || cfg.excitationPulseIndex > numel(starts)
        return;
    end
    nativePeriodBins = max(1, round(periodSeconds / rawResolutionSeconds));
    gateStartNative = mod((starts(cfg.excitationPulseIndex) - 1) * ...
        channelDivision, nativePeriodBins);
    gateLengthNative = gateLength * channelDivision;
    if gateLengthNative < 1 || gateLengthNative > nativePeriodBins
        gateStartNative = [];
        gateLengthNative = [];
        return;
    end
    ok = true;
end

function [ok, gateStarts, gateStops, gateLength, nominalStarts] = ...
        gatesFromPtuHeader(head, periodBins, requestedPulseCount, preGateBins)
    ok = false;
    gateStarts = [];
    gateStops = [];
    gateLength = [];
    nominalStarts = [];
    if ~isfield(head, 'PIETimeGate') || isempty(head.PIETimeGate)
        return;
    end
    timeGatesNative = double(head.PIETimeGate(:));
    if isfield(head, 'PIENumPIEWindows') && ~isempty(head.PIENumPIEWindows)
        recordedCount = round(double(head.PIENumPIEWindows(1)));
    else
        recordedCount = numel(timeGatesNative);
    end
    pulseCount = min([requestedPulseCount, recordedCount, numel(timeGatesNative)]);
    if pulseCount < requestedPulseCount
        return;
    end
    timeGatesNative = timeGatesNative(1:pulseCount);
    if isfield(head, 'TCSPC_chDiv') && ~isempty(head.TCSPC_chDiv)
        channelDivision = double(head.TCSPC_chDiv(1));
    elseif isfield(head, 'MeasDesc_Resolution_Original') && ...
            ~isempty(head.MeasDesc_Resolution_Original)
        channelDivision = round(double(head.MeasDesc_Resolution) / ...
            double(head.MeasDesc_Resolution_Original));
    else
        channelDivision = 1;
    end
    channelDivision = max(1, round(channelDivision));
    nominalStarts = floor(timeGatesNative / channelDivision) + 1;
    if any(~isfinite(nominalStarts)) || any(diff(nominalStarts) <= 0) || ...
            nominalStarts(1) < 1 || nominalStarts(end) > periodBins
        nominalStarts = [];
        return;
    end
    circularIntervals = diff([nominalStarts; nominalStarts(1) + periodBins]);
    if any(circularIntervals < 1)
        nominalStarts = [];
        return;
    end
    gateLength = floor(min(circularIntervals));
    gateStarts = mod(nominalStarts - preGateBins - 1, periodBins) + 1;
    gateStops = mod(gateStarts + gateLength - 2, periodBins) + 1;
    gateStarts = round(gateStarts(:));
    gateStops = round(gateStops(:));
    nominalStarts = round(nominalStarts(:));
    ok = true;
end

function [peakBins, peakHeights] = peaksWithinHeaderGates(profile, nominalStarts, gateLength)
    periodBins = numel(profile);
    peakBins = zeros(numel(nominalStarts), 1);
    peakHeights = zeros(numel(nominalStarts), 1);
    for pulse = 1:numel(nominalStarts)
        indices = mod((nominalStarts(pulse) - 1) + (0:gateLength-1), ...
            periodBins) + 1;
        [peakHeights(pulse), relativeIndex] = max(profile(indices));
        peakBins(pulse) = indices(relativeIndex);
    end
end

function gated = isolateExcitationGate(ptu, gateStart, gateLength, ...
        detectorIds, photonBlockSize)
    gated = ptu;
    periodBins = double(ptu.Ngate);
    n = min([numel(ptu.im_tcspc), numel(ptu.im_chan), ...
        numel(ptu.im_col), numel(ptu.im_line)]);
    relativeBin = zeros(n, 1, 'uint16');
    keep = true(n, 1);
    for first = 1:photonBlockSize:n
        last = min(first + photonBlockSize - 1, n);
        indices = first:last;
        blockTime = round(double(ptu.im_tcspc(indices)));
        blockRelative = mod(blockTime - gateStart, periodBins) + 1;
        relativeBin(indices) = uint16(blockRelative);
        keep(indices) = blockRelative <= gateLength;
    end
    alignedFields = {'im_chan', 'im_col', 'im_line', 'im_sync', ...
        'im_frame', 'im_tcspc_native'};
    if all(keep)
        gated.im_tcspc = relativeBin;
    else
        gated.im_tcspc = relativeBin(keep);
        for k = 1:numel(alignedFields)
            name = alignedFields{k};
            if isfield(ptu, name) && numel(ptu.(name)) >= n
                values = ptu.(name);
                values = values(:);
                gated.(name) = values(1:n);
                gated.(name) = gated.(name)(keep);
            end
        end
    end
    gated.Ngate = gateLength;
    gated.dind = detectorIds(:).';
    gated.Resolution_ns = resolveDtNs(ptu);
    gated.head.MeasDesc_Resolution = gated.Resolution_ns * 1e-9;

    nx = double(ptu.head.ImgHdr_PixX);
    ny = double(ptu.head.ImgHdr_PixY);
    nDetectors = numel(detectorIds);
    tags = zeros(nx, ny, nDetectors, 'uint32');
    detectorLut = zeros(256, 1, 'uint16');
    detectorLut(double(detectorIds(:)) + 1) = uint16(1:nDetectors);
    photonCount = min([numel(gated.im_chan), numel(gated.im_col), ...
        numel(gated.im_line)]);
    for first = 1:photonBlockSize:photonCount
        last = min(first + photonBlockSize - 1, photonCount);
        indices = first:last;
        column = double(gated.im_col(indices));
        line = double(gated.im_line(indices));
        detector = double(detectorLut(double(gated.im_chan(indices)) + 1));
        valid = detector >= 1 & column >= 1 & column <= nx & ...
            line >= 1 & line <= ny;
        linear = column(valid) + (line(valid) - 1) * nx + ...
            (detector(valid) - 1) * nx * ny;
        tags = addCountsUint32(tags, linear);
    end
    gated.tags = tags;
    gated.tag = tags;
    if isfield(gated, 'tcspc_pix')
        gated.tcspc_pix = [];
    end
    if isfield(gated, 'tcspc_pix_mt')
        gated.tcspc_pix_mt = [];
    end
end

function array = addCountsUint32(array, linearIndex)
    if isempty(linearIndex)
        return;
    end
    [uniqueIndex, ~, group] = unique(double(linearIndex(:)));
    additions = accumarray(group, 1);
    updated = double(array(uniqueIndex)) + additions;
    if any(updated > double(intmax('uint32')))
        error('immune_cell_MIET:DetectorCountOverflow', ...
            'Detector image counts exceed uint32 range.');
    end
    array(uniqueIndex) = uint32(updated);
end

function params = defaultIsmParams(ptu, cfg)
    params = struct('imageSource', 'tags', 'smoothSigma', 1, ...
        'useWindow', true, 'normalizeImages', true, 'upsampleReg', 20, ...
        'nIter', 500, 'checkEvery', 25, 'stopTol', 1e-7, ...
        'minIter', 50, 'pixelSize', double(ptu.head.ImgHdr_PixResol) * 1e3, ...
        'lambda', cfg.ismWavelengthNm, 'NA', 1.45, 'showPlots', false, ...
        'useGPU', cfg.useGPU);
end

function [tau, intensity, source] = resolveReassignedMean(flim, dtNs)
    if isfield(flim, 'reassigned') && isfield(flim.reassigned, 'total')
        total = flim.reassigned.total;
    elseif isfield(flim, 'total')
        total = flim.total;
    else
        error('immune_cell_MIET:ReassignedOutput', ...
            'reassigned_flim returned no total reassigned FLIM data.');
    end
    if isfield(total, 'xyT') && ~isempty(total.xyT)
        cube = double(total.xyT);
        intensity = sum(cube, 3);
        globalDecay = squeeze(sum(sum(cube, 1), 2));
        [~, t0] = max(globalDecay);
        if isfield(total, 'tAxisNs') && numel(total.tAxisNs) == size(cube, 3)
            timeNs = double(total.tAxisNs(:));
        else
            timeNs = (0:size(cube, 3)-1).' * dtNs;
        end
        timeNs = timeNs - timeNs(t0);
        tau = sum(cube .* reshape(timeNs, 1, 1, []), 3) ./ max(intensity, 1);
        source = 'reassigned.total.xyT';
    elseif isfield(total, 'tauMean') && isfield(total, 'tag')
        tau = double(total.tauMean);
        intensity = double(total.tag);
        source = 'reassigned.total.tauMean';
    else
        error('immune_cell_MIET:ReassignedOutput', ...
            'No reassigned TCSPC cube or mean-arrival map was returned.');
    end
end

function h = showRedMeanFlim(tau, intensity, requestedRange, excitationNm, ...
        visibility, analysisLabel)
    if nargin < 6
        analysisLabel = '';
    end
    h = figure('Name', 'immune_cell_MIET: red-channel mean FLIM', ...
        'Color', 'w', 'Visible', visibility);
    ax = axes(h);
    imagesc(ax, tau.');
    axis(ax, 'image', 'off');
    set(ax, 'YDir', 'normal', 'Color', 'k');
    % The acquisition channel is red; the lifetime itself uses a blue-to-red
    % scale so short outside-SLB and long cell-membrane arrivals separate.
    colormap(ax, turbo(256));
    alpha = isfinite(tau) & double(intensity) > 0;
    set(ax.Children(1), 'AlphaData', alpha.');
    lifetimeRange = requestedRange;
    finiteTau = tau(isfinite(tau));
    if isempty(lifetimeRange) && ~isempty(finiteTau)
        lifetimeRange = finitePercentiles(finiteTau, [0.02 0.98]);
    end
    if numel(lifetimeRange) == 2 && lifetimeRange(2) > lifetimeRange(1)
        clim(ax, lifetimeRange);
    end
    cb = colorbar(ax);
    cb.Label.String = 'Mean photon arrival time (ns)';
    titleText = sprintf('%g nm excitation: ACO-ISM reassigned mean FLIM', ...
        excitationNm);
    if ~isempty(analysisLabel)
        titleText = sprintf('%s | %s', titleText, analysisLabel);
    end
    title(ax, titleText, 'Interpreter', 'none');
end

function h = showSegmentationCheck(meanArrival, ~, masks, diagnostics, visibility)
    h = figure('Name', 'immune_cell_MIET: segmentation check', ...
        'Color', 'w', 'Visible', visibility, 'Position', [100 80 1100 850]);
    layout = tiledlayout(h, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(layout);
    showImage(ax, meanArrival, turbo(256), 'Preliminary reassigned mean FLIM (ns)');
    overlayBoundary(ax, masks.cellBoundary);

    ax = nexttile(layout);
    showImage(ax, diagnostics.smoothedMeanArrivalNs, turbo(256), ...
        'Photon-weighted smoothed mean FLIM (ns)');
    overlayBoundary(ax, masks.cellBoundary);

    ax = nexttile(layout);
    showImage(ax, diagnostics.lifetimeExcessZ, parula(256), ...
        'Lifetime excess above outside SLB (robust \sigma)');
    overlayBoundary(ax, masks.cellBoundary);

    ax = nexttile(layout);
    imagesc(ax, masks.labelMap.');
    axis(ax, 'image', 'off');
    set(ax, 'YDir', 'normal');
    colormap(ax, [0 0 0; 0.1 0.35 1; 1 0.75 0; 1 0.1 0.1]);
    colorbar(ax, 'Ticks', 0:3, ...
        'TickLabels', {'unclassified','outside SLB','transition','membrane'});
    title(ax, 'Derived cell footprint and SLB reference');

    threshold = diagnostics.thresholds;
    title(layout, sprintf(['Segmentation: %s | SLB %.3g ns | core %.3g ns | ' ...
        'cell %d px | SLB ref %d px'], diagnostics.status, ...
        threshold.slbMeanArrivalNs, threshold.cellCoreMeanArrivalNs, ...
        diagnostics.pixelCounts.cellMembrane, diagnostics.pixelCounts.slbReference), ...
        'Interpreter', 'none');
    if isfield(diagnostics, 'fallback') && diagnostics.fallback.used
        subtitle(layout, diagnostics.fallback.method, 'Interpreter', 'none');
    end
end

function overlayBoundary(ax, boundary)
    if any(boundary(:))
        hold(ax, 'on');
        contour(ax, boundary.', [0.5 0.5], 'w', 'LineWidth', 1.2);
        hold(ax, 'off');
    end
end

function maps = scatterBayesMaps(compact, indices, imageSize)
    scalarFields = {'probabilityFixedSlbOnly', 'probabilityBiexponential', ...
        'probabilityTriexponential', 'tauMeanArithmetic', 'tauPosteriorStd', ...
        'fixedSlbFraction', 'membraneFraction', 'backgroundFraction', ...
        'signalFraction', 'fixedSlbPhotonFraction', ...
        'fixedSlbExpectedPhotonCount', ...
        'membrane1PhotonFraction', 'membrane2PhotonFraction'};
    maps = struct();
    for k = 1:numel(scalarFields)
        name = scalarFields{k};
        maps.(name) = scatterScalarMap(compact.(name), indices, imageSize, NaN);
    end
    maps.modelProbability = scatterStack(compact.modelProbability, ...
        indices, imageSize, NaN);
    maps.logEvidence = scatterStack(compact.logEvidence, indices, imageSize, NaN);
    maps.logBayesFactorM3VsM2 = ...
        maps.logEvidence(:, :, 3) - maps.logEvidence(:, :, 2);
    membraneModelProbability = maps.probabilityBiexponential + ...
        maps.probabilityTriexponential;
    maps.probabilityM3GivenM2orM3 = nan(imageSize, 'single');
    validMembraneModel = isfinite(membraneModelProbability) & ...
        membraneModelProbability > 0;
    maps.probabilityM3GivenM2orM3(validMembraneModel) = ...
        maps.probabilityTriexponential(validMembraneModel) ./ ...
        membraneModelProbability(validMembraneModel);
    maps.componentLifetimeNs = scatterStack(compact.componentLifetimeNs, ...
        indices, imageSize, NaN);
    maps.componentSpeciesFraction = scatterStack( ...
        compact.componentSpeciesFraction, indices, imageSize, NaN);
    maps.componentPhotonFraction = scatterStack( ...
        compact.componentPhotonFraction, indices, imageSize, NaN);
    maps.completeExponentialCountMAP = scatterScalarMap( ...
        compact.completeExponentialCountMAP, indices, imageSize, 0);
    maps.membraneComponentCountMAP = scatterScalarMap( ...
        compact.membraneComponentCountMAP, indices, imageSize, 0);
    maps.fixedSlbLifetimeNs = compact.fixedSlbLifetimeNs;

    selectedFields = {'tauMeanArithmetic', 'fixedSlbFraction', ...
        'membraneLifetime1Ns', 'membraneLifetime2Ns', ...
        'membraneFraction1', 'membraneFraction2', 'backgroundFraction', ...
        'signalFraction', 'fixedSlbPhotonFraction', ...
        'membrane1PhotonFraction', 'membrane2PhotonFraction'};
    maps.selected = struct('fixedSlbLifetimeNs', compact.fixedSlbLifetimeNs);
    for k = 1:numel(selectedFields)
        name = selectedFields{k};
        maps.selected.(name) = scatterScalarMap( ...
            compact.selected.(name), indices, imageSize, NaN);
    end
end

function fullMap = scatterScalarMap(compactMap, indices, imageSize, fillValue)
    if isa(compactMap, 'uint8')
        fullMap = repmat(uint8(fillValue), imageSize);
    else
        fullMap = repmat(single(fillValue), imageSize);
    end
    fullMap(indices) = reshape(compactMap, [], 1);
end

function fullStack = scatterStack(compactStack, indices, imageSize, fillValue)
    planeCount = size(compactStack, 3);
    fullStack = repmat(single(fillValue), [imageSize planeCount]);
    for plane = 1:planeCount
        map = fullStack(:, :, plane);
        compactPlane = compactStack(:, :, plane);
        map(indices) = compactPlane(:);
        fullStack(:, :, plane) = map;
    end
end

function h = showBayesSummary(meanArrival, intensity, masks, bayes, tauSlb, ...
        visibility, analysisLabel)
    if nargin < 7
        analysisLabel = '';
    end
    h = figure('Name', 'immune_cell_MIET: fixed-SLB Bayesian FLIM', ...
        'Color', 'w', 'Visible', visibility, 'Position', [80 80 1250 720]);
    layout = tiledlayout(h, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(layout);
    showImage(ax, meanArrival, parula(256), 'Mean-arrival segmentation');
    hold(ax, 'on');
    contour(ax, masks.cellBoundary.', [0.5 0.5], 'w', 'LineWidth', 0.8);
    hold(ax, 'off');

    ax = nexttile(layout);
    imagesc(ax, masks.labelMap.');
    axis(ax, 'image', 'off'); set(ax, 'YDir', 'normal');
    colormap(ax, [0 0 0; 0.1 0.35 1; 1 0.75 0; 1 0.1 0.1]);
    colorbar(ax, 'Ticks', 0:3, ...
        'TickLabels', {'unclassified','outside SLB','transition','membrane'});
    title(ax, sprintf('Segmentation; fixed \tau_{SLB}=%.3g ns', tauSlb));

    ax = nexttile(layout);
    showImage(ax, bayes.tauMeanArithmetic, turbo(256), ...
        'Bayesian model-averaged lifetime (ns)');

    ax = nexttile(layout);
    showImage(ax, bayes.fixedSlbFraction, parula(256), ...
        'Fixed-SLB species fraction', [0 1]);

    ax = nexttile(layout);
    showImage(ax, bayes.probabilityBiexponential, parula(256), ...
        'P(SLB + one membrane state)', [0 1]);

    ax = nexttile(layout);
    showImage(ax, bayes.probabilityTriexponential, hot(256), ...
        'P(SLB + two membrane states)', [0 1]);

    titleText = sprintf('Fixed-SLB Bayesian FLIM | membrane photons: %.3g', ...
        sum(double(intensity(masks.cellMembrane)), 'all'));
    if ~isempty(analysisLabel)
        titleText = sprintf('%s | %s', titleText, analysisLabel);
    end
    title(layout, titleText, 'Interpreter', 'none');
end

function limits = orderedLifetimeLimits(layers, requested)
    if ~isempty(requested)
        validateattributes(requested, {'numeric'}, ...
            {'real','finite','vector','numel',2}, ...
            mfilename, 'cfg.componentMaps.lifetimeRangeNs');
        limits = double(requested(:)).';
        if limits(2) <= limits(1)
            error('immune_cell_MIET:ComponentLifetimeRange', ...
                'cfg.componentMaps.lifetimeRangeNs must be increasing.');
        end
        return;
    end
    values = [double(layers.display.secondLifetimeNs(:)); ...
        double(layers.display.thirdLifetimeNs(:))];
    values = values(isfinite(values));
    if isempty(values)
        values = [double(layers.ungated.secondLifetimeNs(:)); ...
            double(layers.ungated.thirdLifetimeNs(:))];
        values = values(isfinite(values));
    end
    if isempty(values)
        limits = [0 1];
        return;
    end
    limits = finitePercentiles(values, [0.02 0.98]);
    if ~(limits(2) > limits(1))
        padding = max(0.05 * abs(limits(1)), 0.05);
        limits = limits(1) + [-padding padding];
    end
end

function h = showOrderedLifetimeMap(layers, ~, component, limits, ...
        visibility, tauSlbNs, slbAmplitude, analysisLabel)
    if nargin < 8
        analysisLabel = '';
    end
    if component == 2
        data = layers.display.secondLifetimeNs;
        displayMask = layers.masks.secondDisplay;
        threshold = layers.config.posteriorThreshold(1);
        expectedThreshold = layers.config.minExpectedPhotons(1);
        label = 'Second sorted lifetime';
        probabilityLabel = 'P(M2 or M3)';
        count = layers.summary.secondDisplayPixelCount;
    elseif component == 3
        data = layers.display.thirdLifetimeNs;
        displayMask = layers.masks.thirdDisplay;
        threshold = layers.config.posteriorThreshold(2);
        expectedThreshold = layers.config.minExpectedPhotons(2);
        label = 'Third sorted (longest) lifetime';
        probabilityLabel = 'P(M3 | M2 or M3)';
        count = layers.summary.thirdDisplayPixelCount;
    else
        error('immune_cell_MIET:OrderedComponent', ...
            'Ordered component must be 2 or 3.');
    end

    h = figure('Name', ['immune_cell_MIET: ' label], 'Color', 'w', ...
        'Visible', visibility, 'Position', [120 80 850 760]);
    ax = axes(h);
    object = imagesc(ax, data.');
    axis(ax, 'image', 'off');
    set(ax, 'YDir', 'normal', 'Color', [0.015 0.015 0.025]);
    % The separate component figures are deliberately plain lifetime maps.
    % layers.display and displayMask already apply the Bayesian presence and
    % expected-photon gates, so no support/probability contours are needed.
    set(object, 'AlphaData', displayMask.');
    colormap(ax, turbo(256));
    clim(ax, limits);
    colourBar = colorbar(ax);
    colourBar.Label.String = 'Lifetime [ns]';
    title(ax, sprintf(['%s (component %d) | %s >= %.2f, ' ...
        'expected photons >= %.1f | %d pixels'], label, component, ...
        probabilityLabel, threshold, expectedThreshold, count));
    subtitleText = sprintf(['Fixed SLB: tau %.3g ns, expected amplitude ' ...
        '%.3g photons/pixel; %.1f%% clipped'], tauSlbNs, ...
        slbAmplitude.signalPhotonsPerPixel, ...
        100 * slbAmplitude.bayesianClippedPixelFraction);
    if ~isempty(analysisLabel)
        subtitleText = sprintf('%s | %s', analysisLabel, subtitleText);
    end
    subtitle(ax, subtitleText, 'Interpreter', 'none');
end

function h = showOrderedLifetimeLayers(layers, masks, limits, visibility, ...
        tauSlbNs, slbAmplitude, analysisLabel)
    if nargin < 7
        analysisLabel = '';
    end
    h = figure('Name', 'immune_cell_MIET: ordered lifetime contour layers', ...
        'Color', 'w', 'Visible', visibility, 'Position', [150 90 900 760]);
    ax = axes(h);
    base = layers.display.secondLifetimeNs;
    longest = layers.display.thirdLifetimeNs;
    baseObject = imagesc(ax, base.');
    axis(ax, 'image', 'off');
    set(ax, 'YDir', 'normal', 'Color', [0.015 0.015 0.025]);
    baseAlpha = orderedLifetimeAlpha(base, ...
        layers.probability.secondComponentPresent, ...
        layers.config.posteriorThreshold(1));
    set(baseObject, 'AlphaData', (0.72 * baseAlpha).');
    colormap(ax, turbo(256));
    clim(ax, limits);
    hold(ax, 'on');
    longestObject = imagesc(ax, longest.');
    longestAlpha = orderedLifetimeAlpha(longest, ...
        layers.probability.thirdComponentPresent, ...
        layers.config.posteriorThreshold(2));
    set(longestObject, 'AlphaData', longestAlpha.');
    contour(ax, masks.cellBoundary.', [0.5 0.5], 'w', ...
        'LineWidth', 1.0, 'LineStyle', '--');
    if any(layers.masks.secondDisplay(:))
        contour(ax, layers.masks.secondDisplay.', [0.5 0.5], ...
            'LineColor', [0.1 0.95 1], 'LineWidth', 1.5);
    end
    if any(layers.masks.thirdDisplay(:))
        contour(ax, layers.masks.thirdDisplay.', [0.5 0.5], ...
            'LineColor', [1 0.15 0.8], 'LineWidth', 1.5);
    end
    contourProbability = layers.probability.thirdComponentPresent;
    contourSupport = layers.expectedPhotonCount.thirdComponent >= ...
        layers.config.minExpectedPhotons(2) & ...
        layers.expectedPhotonCount.secondComponent >= ...
        layers.config.minExpectedPhotons(1);
    contourProbability(~contourSupport) = NaN;
    drawProbabilityContours(ax, contourProbability, layers.config, [1 1 1]);
    hold(ax, 'off');
    colourBar = colorbar(ax);
    colourBar.Label.String = 'Lifetime [ns]';
    title(ax, sprintf(['Ordered lifetime layers: cyan = component 2 support; ' ...
        'magenta/white = component 3 contours | %d / %d pixels'], ...
        layers.summary.secondDisplayPixelCount, ...
        layers.summary.thirdDisplayPixelCount));
    subtitleText = sprintf(['Fixed SLB: tau %.3g ns, expected amplitude ' ...
        '%.3g photons/pixel; %.1f%% clipped'], tauSlbNs, ...
        slbAmplitude.signalPhotonsPerPixel, ...
        100 * slbAmplitude.bayesianClippedPixelFraction);
    if ~isempty(analysisLabel)
        subtitleText = sprintf('%s | %s', analysisLabel, subtitleText);
    end
    subtitle(ax, subtitleText, 'Interpreter', 'none');
end

function alpha = orderedLifetimeAlpha(data, probability, threshold)
    alpha = zeros(size(data));
    valid = isfinite(data) & isfinite(probability) & probability >= threshold;
    if threshold < 1
        confidence = (double(probability(valid)) - threshold) / (1 - threshold);
    else
        confidence = ones(nnz(valid), 1);
    end
    alpha(valid) = 0.4 + 0.6 * min(max(confidence, 0), 1);
end

function drawProbabilityContours(ax, probability, cfg, colour)
    if ~isfield(cfg, 'probabilityContourLevels') || ...
            isempty(cfg.probabilityContourLevels)
        return;
    end
    finiteProbability = double(probability(isfinite(probability)));
    if isempty(finiteProbability)
        return;
    end
    levels = unique(double(cfg.probabilityContourLevels(:)).');
    levels = levels(isfinite(levels) & levels > min(finiteProbability) & ...
        levels < max(finiteProbability));
    if isempty(levels)
        return;
    end
    if isscalar(levels)
        levels = [levels levels];
    end
    contour(ax, probability.', levels, 'LineColor', colour, 'LineWidth', 0.8);
end

function showImage(ax, imageData, colourMap, titleText, limits)
    object = imagesc(ax, imageData.');
    axis(ax, 'image', 'off');
    set(ax, 'YDir', 'normal', 'Color', 'k');
    set(object, 'AlphaData', isfinite(imageData).');
    colormap(ax, colourMap);
    if nargin >= 5 && numel(limits) == 2 && limits(2) > limits(1)
        clim(ax, limits);
    else
        finiteValues = imageData(isfinite(imageData));
        if ~isempty(finiteValues)
            robustLimits = finitePercentiles(finiteValues, [0.02 0.98]);
            if robustLimits(2) > robustLimits(1)
                clim(ax, robustLimits);
            end
        end
    end
    colorbar(ax);
    title(ax, titleText);
end

function values = finitePercentiles(data, probabilities)
    data = sort(double(data(isfinite(data))));
    values = nan(size(probabilities));
    if isempty(data)
        return;
    end
    for k = 1:numel(probabilities)
        position = 1 + (numel(data) - 1) * min(max(probabilities(k), 0), 1);
        low = floor(position);
        high = ceil(position);
        fraction = position - low;
        values(k) = data(low) * (1 - fraction) + data(high) * fraction;
    end
end
