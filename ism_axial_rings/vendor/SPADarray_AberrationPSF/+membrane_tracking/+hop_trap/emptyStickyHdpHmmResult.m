function result = emptyStickyHdpHmmResult(stepTable)
%EMPTYSTICKYHDPHMMRESULT Construct an unavailable sticky BNP result.

    result = struct();
    result.status = 'insufficient steps';
    result.estimateAvailable = false;
    result.method = 'finite sticky HDP-HMM approximation';
    result.maximumStates = NaN;
    result.nPosteriorSamples = 0;
    result.stateCountTrace = zeros(0, 1);
    result.posteriorModeStateCount = NaN;
    result.posteriorMeanStateCount = NaN;
    result.diffusionTraceUm2PerS = zeros(0, 0);
    result.occupancyTrace = zeros(0, 0);
    result.posteriorRankProbability = nan(height(stepTable), 0);
    result.posteriorSlowStateProbability = nan(height(stepTable), 1);
    result.posteriorMedianDiffusionUm2PerS = zeros(1, 0);
    result.posteriorMedianOccupancy = zeros(1, 0);
    result.stateCountProbability = zeros(1, 0);
end
