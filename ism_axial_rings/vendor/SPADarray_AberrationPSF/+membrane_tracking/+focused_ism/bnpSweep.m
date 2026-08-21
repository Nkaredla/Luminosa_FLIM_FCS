function [state, diagnostics] = bnpSweep(state, detector, opts, freezeDiffusion)
    import membrane_tracking.focused_ism.*

    [state, acceptedPositions, proposedPositions] = ...
        bnpUpdatePositions(state, detector, opts);
    [state, indicatorChanges] = bnpUpdateIndicators(state, opts);
    [state, informativeStepCount] = ...
        bnpUpdateDiffusion(state, opts, freezeDiffusion);
    [state, backgroundAccepted] = bnpUpdateBackground(state, opts);
    [state, brightnessAccepted] = bnpUpdateBrightness(state, opts);
    % Rebuilding suppresses roundoff drift after many add/remove operations.
    state.meanCounts = bnpRebuildMeanCounts(state);

    diagnostics = struct();
    if proposedPositions > 0
        diagnostics.positionAcceptance = ...
            acceptedPositions / proposedPositions;
    else
        diagnostics.positionAcceptance = NaN;
    end
    diagnostics.indicatorChanges = indicatorChanges;
    diagnostics.informativeStepCount = informativeStepCount;
    diagnostics.backgroundAccepted = backgroundAccepted;
    diagnostics.brightnessAccepted = brightnessAccepted;
end
