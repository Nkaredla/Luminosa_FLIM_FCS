function model = emptyTrapHmmResult(stepTable)
%EMPTYTRAPHMMRESULT Construct an unavailable two-state HMM result.

    model = struct();
    model.status = 'insufficient steps';
    model.estimateAvailable = false;
    model.diffusionUm2PerS = [NaN NaN];
    model.diffusionRatioSlowToFast = NaN;
    model.transitionMatrix = nan(2);
    model.initialProbability = [NaN NaN];
    model.posteriorSlowProbability = nan(height(stepTable), 1);
    model.posteriorState = nan(height(stepTable), 1);
    model.stateOccupancy = [NaN NaN];
    model.meanDwellTimeS = [NaN NaN];
    model.observationSpacingS = NaN;
    model.logLikelihood = NaN;
    model.bic = NaN;
    model.nParameters = 5;
    model.nSteps = height(stepTable);
    model.iterations = 0;
end
