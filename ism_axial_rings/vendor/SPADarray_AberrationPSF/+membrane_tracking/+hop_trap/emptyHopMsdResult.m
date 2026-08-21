function comparison = emptyHopMsdResult()
%EMPTYHOPMSDRESULT Construct an unavailable hop-MSD fit result.

    comparison = struct();
    comparison.status = 'not run';
    comparison.estimateAvailable = false;
    comparison.nLagPoints = 0;
    comparison.validRows = zeros(0, 1);
    comparison.lagTimeS = zeros(0, 1);
    comparison.observedMsdUm2 = zeros(0, 1);
    comparison.standardErrorUm2 = zeros(0, 1);
    comparison.freeDiffusionUm2PerS = NaN;
    comparison.freePredictionUm2 = zeros(0, 1);
    comparison.freeLogLikelihood = NaN;
    comparison.freeQuasiBic = NaN;
    comparison.hopParameters = [NaN NaN NaN];
    comparison.hopPredictionUm2 = zeros(0, 1);
    comparison.hopLogLikelihood = NaN;
    comparison.hopQuasiBic = NaN;
    comparison.deltaHopQuasiBic = NaN;
    comparison.longRangeDiffusionUm2PerS = NaN;
    comparison.shortRangeDiffusionUm2PerS = NaN;
    comparison.shortToLongDiffusionRatio = NaN;
    comparison.effectiveCompartmentSizeUm = NaN;
    comparison.confinementTimeS = NaN;
end
