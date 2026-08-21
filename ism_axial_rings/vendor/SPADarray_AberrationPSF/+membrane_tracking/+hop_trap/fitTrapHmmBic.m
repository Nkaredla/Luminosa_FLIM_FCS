function model = fitTrapHmmBic(stepTable, freeModel, opts)
%FITTRAPHMMBIC Fit a two-state diffusion HMM by multiple-start EM.

    import membrane_tracking.hop_trap.*

    model = emptyTrapHmmResult(stepTable);
    if ~freeModel.estimateAvailable || ...
            height(stepTable) < opts.minimumInferenceSteps
        return;
    end

    nStart = opts.trapHmmStarts;
    bestLogLikelihood = -Inf;
    best = struct();
    freeD = freeModel.diffusionUm2PerS;
    for start = 1:nStart
        fraction = (start - 1) / max(nStart - 1, 1);
        diffusion = sort([max(1e-8, freeD * 10^(-2 + fraction)), ...
            max(1e-7, freeD * (1.2 + 1.8 * fraction))]);
        persistence = 0.85 + 0.14 * fraction;
        transition = [persistence, 1-persistence; ...
            1-persistence, persistence];
        initial = [0.5 0.5];
        previous = -Inf;

        for iteration = 1:opts.trapHmmMaximumIterations
            logEmission = [logGaussianStepDensity(stepTable, diffusion(1)), ...
                logGaussianStepDensity(stepTable, diffusion(2))];
            [logLikelihood, gamma, transitionCounts] = ...
                forwardBackwardHmm(logEmission, stepTable.sequenceId, ...
                initial, transition);

            sequenceStarts = [true; diff(stepTable.sequenceId) ~= 0];
            initial = sum(gamma(sequenceStarts,:), 1) + 1e-3;
            initial = initial / sum(initial);
            transition = transitionCounts + 1e-3;
            transition = transition ./ sum(transition, 2);
            for state = 1:2
                weights = gamma(:,state);
                objective = @(logD) -sum(weights .* ...
                    logGaussianStepDensity(stepTable, exp(logD)));
                diffusion(state) = exp(fminbnd(objective, log(1e-8), ...
                    log(max(100, 20 * opts.freeDiffusionUm2PerS)), ...
                    optimset('Display', 'off', 'TolX', 1e-6)));
            end

            [diffusion, order] = sort(diffusion);
            transition = transition(order, order);
            initial = initial(order);
            if iteration > 2 && abs(logLikelihood - previous) < ...
                    1e-6 * (1 + abs(logLikelihood))
                break;
            end
            previous = logLikelihood;
        end

        logEmission = [logGaussianStepDensity(stepTable, diffusion(1)), ...
            logGaussianStepDensity(stepTable, diffusion(2))];
        [logLikelihood, gamma, transitionCounts] = ...
            forwardBackwardHmm(logEmission, stepTable.sequenceId, ...
            initial, transition);
        if logLikelihood > bestLogLikelihood
            bestLogLikelihood = logLikelihood;
            best.diffusion = diffusion;
            best.transition = transition;
            best.initial = initial;
            best.gamma = gamma;
            best.transitionCounts = transitionCounts;
            best.iterations = iteration;
        end
    end

    observationSpacing = estimateStepObservationSpacing(stepTable);
    dwell = observationSpacing ./ max(1 - diag(best.transition), eps);
    model.status = 'ok';
    model.estimateAvailable = true;
    model.diffusionUm2PerS = best.diffusion;
    model.diffusionRatioSlowToFast = ...
        best.diffusion(1) / max(best.diffusion(2), realmin);
    model.transitionMatrix = best.transition;
    model.initialProbability = best.initial;
    model.posteriorSlowProbability = best.gamma(:,1);
    model.posteriorState = 1 + (best.gamma(:,2) > best.gamma(:,1));
    model.stateOccupancy = mean(best.gamma, 1);
    model.meanDwellTimeS = dwell(:).';
    model.observationSpacingS = observationSpacing;
    model.logLikelihood = bestLogLikelihood;
    model.bic = model.nParameters * log(height(stepTable)) - ...
        2 * model.logLikelihood;
    model.iterations = best.iterations;
end
