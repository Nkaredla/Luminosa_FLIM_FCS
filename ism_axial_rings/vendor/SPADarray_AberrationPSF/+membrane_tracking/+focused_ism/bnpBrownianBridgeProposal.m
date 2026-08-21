function proposal = bnpBrownianBridgeProposal(path, frame, ...
        diffusionUm2PerS, opts)
    import membrane_tracking.focused_ism.*

    nFrames = size(path, 1);
    transitionVariance = 2 * diffusionUm2PerS * opts.dtS;
    if nFrames == 1
        proposal = (rand(1, 2) - 0.5) .* opts.boxSizeUm;
        return;
    elseif frame == 1
        bridgeMean = path(2,:);
        proposalVariance = transitionVariance;
    elseif frame == nFrames
        bridgeMean = path(nFrames-1,:);
        proposalVariance = transitionVariance;
    else
        neighborDisplacement = path(frame+1,:) - path(frame-1,:);
        neighborDisplacement = applyMinimumImage(neighborDisplacement, opts);
        bridgeMean = path(frame-1,:) + 0.5 * neighborDisplacement;
        proposalVariance = transitionVariance / 2;
    end
    proposal = bridgeMean + sqrt(max(proposalVariance, realmin)) * ...
        randn(1, 2);
    proposal = applyBoundary(proposal, opts.boxSizeUm, ...
        opts.boundaryCondition);
end
