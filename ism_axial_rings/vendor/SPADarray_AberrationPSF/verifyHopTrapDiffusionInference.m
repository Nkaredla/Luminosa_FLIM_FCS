function report = verifyHopTrapDiffusionInference()
%VERIFYHOPTRAPDIFFUSIONINFERENCE Numerical checks without photon simulation.
%
%   report = verifyHopTrapDiffusionInference()
%
%   This verifies periodic mesh indexing, free-diffusion recovery, trap-HMM
%   evidence, and hop-MSD evidence using synthetic inputs. It does not replace
%   an end-to-end focused-ISM simulation.

    projectDirectory = fileparts(mfilename('fullpath'));
    addpath(projectDirectory);
    import membrane_tracking.hop_trap.*

    rng(91);
    opts = parseInputs('makeFigure', false, 'verbose', false, ...
        'runStickyBnp', false, 'minimumInferenceSteps', 30);

    mesh = makePeriodicVoronoiMesh(opts);
    probe = (rand(80, 2) - 0.5) .* opts.boxSizeUm;
    originalCell = periodicVoronoiCell(probe, mesh, opts.boxSizeUm);
    shiftedCell = periodicVoronoiCell(probe + opts.boxSizeUm, ...
        mesh, opts.boxSizeUm);
    periodicMeshPassed = isequal(originalCell, shiftedCell);

    nStep = 1600;
    dt = opts.dtS * ones(nStep, 1);
    localizationSigma = 0.008;
    noiseVariance = 2 * localizationSigma^2 * ones(nStep, 1);
    trueFreeD = 0.50;
    stepSigma = sqrt(2 * trueFreeD * dt + noiseVariance);
    dx = stepSigma .* randn(nStep, 1);
    dy = stepSigma .* randn(nStep, 1);
    sequenceId = ones(nStep, 1);
    frameStart = (1:2:2*nStep).';
    frameEnd = frameStart + 1;
    timeCenterS = 0.5 * (frameStart + frameEnd - 2) * opts.dtS;
    trackId = ones(nStep, 1);
    noiseCovXYUm2 = zeros(nStep, 1);
    trueTrapState = ones(nStep, 1);
    freeSteps = table(sequenceId, trackId, frameStart, frameEnd, ...
        timeCenterS, dt, dx, dy, noiseVariance, noiseVariance, ...
        noiseCovXYUm2, trueTrapState, 'VariableNames', ...
        {'sequenceId', 'trackId', 'frameStart', 'frameEnd', ...
        'timeCenterS', 'stepDtS', 'dxUm', 'dyUm', 'noiseVarXUm2', ...
        'noiseVarYUm2', 'noiseCovXYUm2', 'trueTrapState'});
    freeFit = fitFreeStepModel(freeSteps, opts);
    freeRecoveryError = abs(freeFit.diffusionUm2PerS / trueFreeD - 1);

    transition = [0.96 0.04; 0.10 0.90];
    state = ones(nStep, 1);
    for t = 2:nStep
        if rand > transition(state(t-1), state(t-1))
            state(t) = 3 - state(t-1);
        else
            state(t) = state(t-1);
        end
    end
    stateDiffusion = [0.015 0.50];
    trapSigma = sqrt(2 * stateDiffusion(state).' .* dt + noiseVariance);
    trapSteps = freeSteps;
    trapSteps.dxUm = trapSigma .* randn(nStep, 1);
    trapSteps.dyUm = trapSigma .* randn(nStep, 1);
    trapSteps.trueTrapState = state;
    trapFreeFit = fitFreeStepModel(trapSteps, opts);
    trapFit = fitTrapHmmBic(trapSteps, trapFreeFit, opts);
    deltaTrapBic = trapFreeFit.bic - trapFit.bic;

    lagFrame = (1:60).';
    lagTimeS = lagFrame * opts.dtS;
    hopTruth = [0.035, 0.10^2/3, 0.012];
    hopMean = hopMsdModel(hopTruth, lagTimeS);
    standardErrorUm2 = 0.025 * max(hopMean, max(hopMean) / 20);
    correctedMsdUm2 = hopMean + standardErrorUm2 .* ...
        randn(size(hopMean));
    rawMsdUm2 = correctedMsdUm2;
    nPairs = 500 * ones(size(lagFrame));
    msdTable = table(lagFrame, lagTimeS, rawMsdUm2, ...
        correctedMsdUm2, standardErrorUm2, nPairs);
    hopFit = fitHopMsdQuasiBic(msdTable, opts);

    report = struct();
    report.periodicMeshPassed = periodicMeshPassed;
    report.trueFreeDiffusionUm2PerS = trueFreeD;
    report.fittedFreeDiffusionUm2PerS = freeFit.diffusionUm2PerS;
    report.freeRelativeError = freeRecoveryError;
    report.deltaTrapBic = deltaTrapBic;
    report.fittedTrapDiffusionUm2PerS = trapFit.diffusionUm2PerS;
    report.deltaHopQuasiBic = hopFit.deltaHopQuasiBic;
    report.fittedCompartmentSizeUm = hopFit.effectiveCompartmentSizeUm;
    report.passed = periodicMeshPassed && freeRecoveryError < 0.12 && ...
        deltaTrapBic > opts.bicEvidenceThreshold && ...
        hopFit.deltaHopQuasiBic > opts.bicEvidenceThreshold;

    fprintf(['[hop/trap verification] periodic=%d, free error=%.2f%%, ' ...
        'DeltaBIC_trap=%.2f, DeltaqBIC_hop=%.2f, passed=%d.\n'], ...
        periodicMeshPassed, 100 * freeRecoveryError, deltaTrapBic, ...
        hopFit.deltaHopQuasiBic, report.passed);
    if ~report.passed
        error('verifyHopTrapDiffusionInference:Failed', ...
            'At least one hop/trap inference verification check failed.');
    end
end
