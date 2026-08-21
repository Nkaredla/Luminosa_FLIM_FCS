function tableCounts = sampleHdpTableCounts( ...
        transitionCounts, beta, alpha, kappa)
%SAMPLEHDPTABLECOUNTS Sample sticky-HDP Chinese-restaurant table counts.
%
%   Self-transition tables are thinned into global-beta and sticky-mass
%   contributions before the global weights are updated.

    nState = size(transitionCounts, 1);
    tableCounts = zeros(nState);
    for row = 1:nState
        for column = 1:nState
            customerCount = round(transitionCounts(row, column));
            concentration = alpha * beta(column);
            if row == column
                concentration = concentration + kappa;
            end
            if customerCount <= 0 || concentration <= 0
                continue;
            end
            probability = concentration ./ ...
                (concentration + (0:customerCount-1));
            nTables = sum(rand(1, customerCount) < probability);
            if row == column && nTables > 0
                globalFraction = alpha * beta(column) / concentration;
                nTables = sum(rand(1, nTables) < globalFraction);
            end
            tableCounts(row, column) = nTables;
        end
    end
end
