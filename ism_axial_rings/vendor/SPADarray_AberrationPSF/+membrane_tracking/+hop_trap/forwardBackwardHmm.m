function [logLikelihood, gamma, transitionCounts] = ...
        forwardBackwardHmm(logEmission, sequenceId, initialProbability, ...
        transitionMatrix)
%FORWARDBACKWARDHMM Scaled HMM forward-backward recursion for many tracks.

    nObservation = size(logEmission, 1);
    nState = size(logEmission, 2);
    gamma = zeros(nObservation, nState);
    transitionCounts = zeros(nState);
    logLikelihood = 0;
    sequences = unique(sequenceId, 'stable').';

    for sequence = sequences
        rows = find(sequenceId == sequence);
        n = numel(rows);
        localEmission = logEmission(rows,:);
        emissionMaximum = max(localEmission, [], 2);
        emission = exp(localEmission - emissionMaximum);
        alpha = zeros(n, nState);
        scale = zeros(n, 1);

        alpha(1,:) = initialProbability(:).' .* emission(1,:);
        scale(1) = max(sum(alpha(1,:)), realmin);
        alpha(1,:) = alpha(1,:) / scale(1);
        for t = 2:n
            alpha(t,:) = (alpha(t-1,:) * transitionMatrix) .* ...
                emission(t,:);
            scale(t) = max(sum(alpha(t,:)), realmin);
            alpha(t,:) = alpha(t,:) / scale(t);
        end
        logLikelihood = logLikelihood + ...
            sum(log(scale) + emissionMaximum);

        beta = ones(n, nState);
        for t = n-1:-1:1
            beta(t,:) = (transitionMatrix * ...
                (emission(t+1,:) .* beta(t+1,:)).').';
            beta(t,:) = beta(t,:) / scale(t+1);
        end
        localGamma = alpha .* beta;
        localGamma = localGamma ./ max(sum(localGamma, 2), realmin);
        gamma(rows,:) = localGamma;

        for t = 1:n-1
            xi = (alpha(t,:).' * ...
                (emission(t+1,:) .* beta(t+1,:))) .* transitionMatrix;
            xi = xi / max(sum(xi, 'all'), realmin);
            transitionCounts = transitionCounts + xi;
        end
    end
end
