function state = sampleHmmStatePaths(logEmission, sequenceId, ...
        initialProbability, transitionMatrix)
%SAMPLEHMMSTATEPATHS Forward-filter backward-sample all HMM sequences.

    import membrane_tracking.hop_trap.*

    nObservation = size(logEmission, 1);
    nState = size(logEmission, 2);
    state = ones(nObservation, 1);
    sequences = unique(sequenceId, 'stable').';

    for sequence = sequences
        rows = find(sequenceId == sequence);
        n = numel(rows);
        filtered = zeros(n, nState);
        logFiltered = log(max(initialProbability, realmin)) + ...
            logEmission(rows(1),:);
        logFiltered = logFiltered - logSumExp(logFiltered, 2);
        filtered(1,:) = exp(logFiltered);
        for t = 2:n
            predicted = filtered(t-1,:) * transitionMatrix;
            logFiltered = log(max(predicted, realmin)) + ...
                logEmission(rows(t),:);
            logFiltered = logFiltered - logSumExp(logFiltered, 2);
            filtered(t,:) = exp(logFiltered);
        end

        localState = ones(n, 1);
        localState(n) = sampleDiscreteLogWeights( ...
            log(max(filtered(n,:), realmin)));
        for t = n-1:-1:1
            weights = filtered(t,:) .* ...
                transitionMatrix(:,localState(t+1)).';
            localState(t) = sampleDiscreteLogWeights( ...
                log(max(weights, realmin)));
        end
        state(rows) = localState;
    end
end
