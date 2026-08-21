function [state, nInformative] = bnpUpdateDiffusion( ...
        state, opts, freezeDiffusion)
    import membrane_tracking.focused_ism.*

    squaredDisplacementSum = 0;
    nInformative = 0;
    for candidate = find(state.active).'
        path = reshape(state.positionsUm(candidate,:,:), state.nFrames, 2);
        contribution = reshape(state.contributions(candidate,:,:), ...
            state.nFrames, state.nChannels);
        expectedSignal = sum(contribution, 2);
        insideModelDisc = hypot(path(:,1), path(:,2)) <= ...
            opts.maxLocalizationRadiusUm;
        informative = expectedSignal(1:end-1) >= ...
            opts.minPhotonsPerLocalization & ...
            expectedSignal(2:end) >= opts.minPhotonsPerLocalization & ...
            insideModelDisc(1:end-1) & insideModelDisc(2:end);
        displacement = path(2:end,:) - path(1:end-1,:);
        displacement = applyMinimumImage(displacement, opts);
        squaredDisplacement = sum(displacement.^2, 2);
        squaredDisplacementSum = squaredDisplacementSum + ...
            sum(squaredDisplacement(informative));
        nInformative = nInformative + sum(informative);
    end

    if freezeDiffusion || nInformative < opts.bnpMinimumInformativeSteps
        return;
    end
    posteriorShape = opts.bnpDiffusionPriorShape + nInformative;
    posteriorScale = opts.bnpDiffusionPriorScale + ...
        squaredDisplacementSum / (4 * opts.dtS);
    state.diffusionUm2PerS = posteriorScale / ...
        max(randg(posteriorShape), realmin);
end
