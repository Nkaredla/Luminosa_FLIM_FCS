function modelComparison = fitDiffusionNatureModels( ...
        stepTable, msdTable, opts)
%FITDIFFUSIONNATUREMODELS Fit free, trapping, and hop-diffusion evidence.

    import membrane_tracking.hop_trap.*

    freeModel = fitFreeStepModel(stepTable, opts);
    trapModel = fitTrapHmmBic(stepTable, freeModel, opts);
    hopComparison = fitHopMsdQuasiBic(msdTable, opts);

    modelComparison = struct();
    modelComparison.freeStep = freeModel;
    modelComparison.trapHmm = trapModel;
    modelComparison.hopMsd = hopComparison;
    if freeModel.estimateAvailable && trapModel.estimateAvailable
        modelComparison.deltaTrapBic = freeModel.bic - trapModel.bic;
    else
        modelComparison.deltaTrapBic = NaN;
    end
    modelComparison.deltaHopQuasiBic = ...
        hopComparison.deltaHopQuasiBic;
end
