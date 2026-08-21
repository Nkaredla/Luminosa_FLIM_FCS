function model = fitFreeStepModel(stepTable, opts)
%FITFREESTEPMODEL Fit one Brownian diffusion coefficient to independent steps.

    import membrane_tracking.hop_trap.*

    model = struct('status', 'insufficient steps', 'estimateAvailable', false, ...
        'diffusionUm2PerS', NaN, 'logLikelihood', NaN, 'bic', NaN, ...
        'nParameters', 1, 'nSteps', height(stepTable));
    if height(stepTable) < opts.minimumInferenceSteps
        return;
    end

    lower = log(1e-8);
    upper = log(max(100, 20 * opts.freeDiffusionUm2PerS));
    objective = @(logD) -sum(logGaussianStepDensity( ...
        stepTable, exp(logD)));
    settings = optimset('Display', 'off', 'TolX', 1e-7);
    [logD, negativeLogLikelihood] = fminbnd( ...
        objective, lower, upper, settings);

    model.status = 'ok';
    model.estimateAvailable = true;
    model.diffusionUm2PerS = exp(logD);
    model.logLikelihood = -negativeLogLikelihood;
    model.bic = model.nParameters * log(height(stepTable)) - ...
        2 * model.logLikelihood;
end
