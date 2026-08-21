function [state, nAccepted, nProposed] = ...
        bnpUpdatePositions(state, detector, opts)
    import membrane_tracking.focused_ism.*

    nAccepted = 0;
    nProposed = 0;
    for parity = 1:2
        for candidate = 1:state.nCandidates
            isActive = state.active(candidate);
            for frame = parity:2:state.nFrames
                currentPath = reshape(state.positionsUm(candidate,:,:), ...
                    state.nFrames, 2);
                proposal = bnpBrownianBridgeProposal(currentPath, frame, ...
                    state.diffusionUm2PerS, opts);
                proposedContribution = bnpEmitterContribution( ...
                    proposal, detector, opts, state.brightnessScale);

                accept = true;
                if isActive
                    nProposed = nProposed + 1;
                    oldContribution = reshape( ...
                        state.contributions(candidate,frame,:), ...
                        1, state.nChannels);
                    proposedMean = state.meanCounts(frame,:) - ...
                        oldContribution + proposedContribution;
                    logRatio = bnpPoissonLogLikelihood( ...
                        state.data(frame,:), proposedMean, opts) - ...
                        bnpPoissonLogLikelihood(state.data(frame,:), ...
                        state.meanCounts(frame,:), opts);
                    accept = log(rand) < min(0, logRatio);
                end

                if accept
                    state.positionsUm(candidate,frame,:) = ...
                        reshape(proposal, 1, 1, 2);
                    state.contributions(candidate,frame,:) = ...
                        reshape(proposedContribution, 1, 1, state.nChannels);
                    if isActive
                        state.meanCounts(frame,:) = proposedMean;
                        nAccepted = nAccepted + 1;
                    end
                end
            end
        end
    end
end
